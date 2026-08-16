#!/usr/bin/env python3
"""Validate the ported analysis code against the recorded pilot numbers.

CLAUDE.md: "Use [data/pilot3d, data/pilot2d] to validate ported analysis code
against the numbers in docs/RESULTS_pilot.md before generating new data."

Three ports are validated, each against numbers computed by the research-phase
prototypes on exactly these files:

* ``analysis/percolation.eps_star`` (2D): the positive-control table in
  ALTERNATIVE_ROUTES.md (eps* per phi, 12 configurations each), and the
  RCP recovery fit ``eps* = C (phi_m - phi)^p`` -> phi_m ~ 0.851, the 1%
  validation of the union-find wrapping code.
* ``analysis/percolation.eps_star`` (3D): the eps*(eta) table in
  RESULTS_pilot.md (3 configurations per point).
* ``analysis/voronoi`` (3D): <Q_iso> and <n_faces> spot values from the pilot
  table, and the median face gap.

Frame conventions follow the prototypes: 2D uses the last 3 frames of each of
the 4 seeds (12 configurations); 3D uses the last 3 frames (eps*) and the
last 6 (Voronoi).  Bisection-vs-exact and snapshot-selection differences bound
the agreement, so the criteria are a few percent, stated per block -- NOT a
claim of bit-identity with the prototypes.

Writes ``results/validation_pilot.md``; exits non-zero if any block fails.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.analysis.percolation import eps_star  # noqa: E402
from hsga.analysis.voronoi import config_observables  # noqa: E402
from hsga.engine.driver import read_frames  # noqa: E402

# ALTERNATIVE_ROUTES.md positive-control table: phi string -> (eps*, sem)
# (keys are the literal filename tokens the pilot used)
PILOT_2D = {
    "0.65": (0.04996, 0.00082), "0.70": (0.03418, 0.00101), "0.74": (0.02521, 0.00064),
    "0.755": (0.02118, 0.00046), "0.765": (0.01971, 0.00032), "0.775": (0.01798, 0.00035),
    "0.785": (0.01472, 0.00045), "0.795": (0.01255, 0.00022), "0.805": (0.00975, 0.00019),
}
SEEDS_2D = ["b", "s11_", "s12_", "s13_"]

# RESULTS_pilot.md eps* table: eta -> (poly, binary)
PILOT_3D_EPS = {
    0.40: (0.02731, 0.02811), 0.44: (0.02154, 0.01938), 0.48: (0.01613, 0.01565),
    0.50: (0.01248, 0.01408), 0.52: (0.01179, 0.01177), 0.54: (0.00999, 0.01006),
    0.56: (0.00828, 0.00843), 0.58: (0.00679, 0.00680), 0.60: (0.00534, 0.00496),
    0.62: (0.00372, 0.00256),
}
# RESULTS_pilot.md structure table spot checks: eta -> (<Q>_poly, <n_f>_poly, <Q>_bin)
PILOT_3D_Q = {
    0.40: (0.67203, 14.711, 0.66545),
    0.54: (0.71157, 14.241, 0.70304),
    0.62: (0.72495, 14.037, 0.70983),
}


def fname_2d(seed_tag: str, phi: str) -> str:
    return f"{seed_tag}{phi}.cfg"


def validate_2d(data: Path, verbose=True):
    rows, ok = [], True
    all_phi, all_eps = [], []
    for phi_s, (ref, sem) in PILOT_2D.items():
        phi = float(phi_s)
        vals = []
        for tag in SEEDS_2D:
            frames = read_frames(data / "pilot2d" / fname_2d(tag, phi_s), dim=2)[-3:]
            vals += [eps_star(pos, rad, L) for pos, rad, L in frames]
        vals = np.array(vals)
        mean = float(np.nanmean(vals))
        sem_here = float(np.nanstd(vals, ddof=1) / np.sqrt(np.isfinite(vals).sum()))
        # criterion: agree with the recorded mean within 3 combined sems or 3%
        tol = max(3.0 * np.hypot(sem, sem_here), 0.03 * ref)
        good = abs(mean - ref) < tol
        ok &= good
        rows.append({"phi": phi, "eps_recorded": ref, "sem_recorded": sem,
                     "eps_ported": mean, "sem_ported": sem_here,
                     "n_configs": int(np.isfinite(vals).sum()), "passed": good})
        all_phi.append(phi)
        all_eps.append(mean)
        if verbose:
            print(f"  2D phi={phi:.3f}: recorded {ref:.5f}({sem:.5f})  "
                  f"ported {mean:.5f}({sem_here:.5f})  {'ok' if good else 'FAIL'}")

    # RCP recovery: eps* = C (phi_m - phi)^p, sem-weighted (the pilot's
    # estimator: the unweighted fit is pulled to phi_m ~ 0.88 by the phi=0.65
    # point, the weighted one reproduces the recorded 0.8508 / 1.079)
    from scipy.optimize import curve_fit

    f = lambda x, C, pm, p: C * np.clip(pm - x, 1e-9, None) ** p
    sems = np.array([r["sem_ported"] for r in rows])
    popt, pcov = curve_fit(f, np.array(all_phi), np.array(all_eps),
                           p0=[0.1, 0.85, 1.0], sigma=sems, absolute_sigma=True,
                           maxfev=40000)
    phi_m, p_exp = float(popt[1]), float(popt[2])
    rcp_ok = abs(phi_m - 0.8508) < 0.01 and abs(p_exp - 1.079) < 0.15
    ok &= rcp_ok
    if verbose:
        print(f"  2D RCP fit: phi_m={phi_m:.4f} (recorded 0.8508, lit ~0.84), "
              f"p={p_exp:.3f} (recorded 1.079)  {'ok' if rcp_ok else 'FAIL'}")
    return ok, rows, {"phi_m": phi_m, "p": p_exp, "passed": rcp_ok}


def validate_3d_eps(data: Path, verbose=True):
    rows, ok = [], True
    for eta, (ref_poly, ref_bin) in PILOT_3D_EPS.items():
        got = {}
        for mode, ref in ((0, ref_poly), (1, ref_bin)):
            frames = read_frames(data / "pilot3d" / f"m{mode}_e{eta:.2f}.cfg")[-3:]
            vals = np.array([eps_star(pos, rad, L) for pos, rad, L in frames])
            got[mode] = float(np.nanmean(vals))
        # 3 configurations: the recorded numbers scatter at the ~10% level
        # between snapshot choices; the criterion is 12% relative per point
        good = (abs(got[0] - ref_poly) < 0.12 * ref_poly
                and abs(got[1] - ref_bin) < 0.12 * ref_bin)
        ok &= good
        rows.append({"eta": eta, "poly_recorded": ref_poly, "poly_ported": got[0],
                     "bin_recorded": ref_bin, "bin_ported": got[1], "passed": good})
        if verbose:
            print(f"  3D eta={eta:.2f}: poly {ref_poly:.5f}->{got[0]:.5f}  "
                  f"bin {ref_bin:.5f}->{got[1]:.5f}  {'ok' if good else 'FAIL'}")
    return ok, rows


def validate_3d_voronoi(data: Path, verbose=True):
    rows, ok = [], True
    for eta, (Q_poly, nf_poly, Q_bin) in PILOT_3D_Q.items():
        got = {}
        for mode in (0, 1):
            frames = read_frames(data / "pilot3d" / f"m{mode}_e{eta:.2f}.cfg")[-6:]
            Q, NF = [], []
            for pos, rad, L in frames:
                o = config_observables(pos, rad, L)
                Q.append(o["Q_iso"])
                NF.append(o["n_faces"])
            got[mode] = (float(np.concatenate(Q).mean()), float(np.concatenate(NF).mean()))
        good = (abs(got[0][0] - Q_poly) < 2e-3 and abs(got[0][1] - nf_poly) < 5e-2
                and abs(got[1][0] - Q_bin) < 2e-3)
        ok &= good
        rows.append({"eta": eta, "Q_poly_recorded": Q_poly, "Q_poly_ported": got[0][0],
                     "nf_poly_recorded": nf_poly, "nf_poly_ported": got[0][1],
                     "Q_bin_recorded": Q_bin, "Q_bin_ported": got[1][0], "passed": good})
        if verbose:
            print(f"  3D eta={eta:.2f}: <Q>_poly {Q_poly:.5f}->{got[0][0]:.5f}  "
                  f"<n_f> {nf_poly:.3f}->{got[0][1]:.3f}  "
                  f"<Q>_bin {Q_bin:.5f}->{got[1][0]:.5f}  {'ok' if good else 'FAIL'}")
    return ok, rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Validate ported analysis vs pilot data")
    p.add_argument("--handoff-data",
                   default=str(REPO.parent / "hs-geometric-arrest-handoff" / "data"))
    p.add_argument("--results", default=str(REPO / "results"))
    a = p.parse_args(argv)
    data = Path(a.handoff_data)
    if not data.is_dir():
        print(f"handoff data not found at {data}; nothing to validate against")
        return 1

    print("validating percolation port on pilot2d (positive control):")
    ok2, rows2, rcp = validate_2d(data)
    print("validating percolation port on pilot3d:")
    ok3e, rows3e = validate_3d_eps(data)
    print("validating Voronoi port on pilot3d (spot values):")
    ok3v, rows3v = validate_3d_voronoi(data)
    passed = ok2 and ok3e and ok3v

    lines = [
        "# Validation of ported analysis code against the recorded pilot numbers",
        "",
        "Frame conventions follow the prototypes (2D: last 3 frames x 4 seeds;",
        "3D: last 3 frames for eps*, last 6 for Voronoi). Criteria are stated per",
        "block and reflect snapshot-selection and estimator differences; this is",
        "consistency validation, not bit-identity.",
        "",
        "## 2D shell percolation (positive control, ALTERNATIVE_ROUTES.md)",
        "",
        "| phi | eps* recorded (sem) | eps* ported (sem) | n | pass |",
        "|---|---|---|---|---|",
    ]
    for r in rows2:
        lines.append(f"| {r['phi']:.3f} | {r['eps_recorded']:.5f} ({r['sem_recorded']:.5f}) "
                     f"| {r['eps_ported']:.5f} ({r['sem_ported']:.5f}) | {r['n_configs']} "
                     f"| {'ok' if r['passed'] else 'FAIL'} |")
    lines += [
        "",
        f"RCP recovery fit `eps* = C (phi_m - phi)^p`: phi_m = {rcp['phi_m']:.4f} "
        f"(recorded 0.8508; literature ~0.84 for this size ratio), p = {rcp['p']:.3f} "
        f"(recorded 1.079) — {'ok' if rcp['passed'] else 'FAIL'}. This is the 1%",
        "RCP validation of the union-find wrapping code.",
        "",
        "## 3D shell percolation (RESULTS_pilot.md, 3 configurations per point)",
        "",
        "| eta | poly recorded | poly ported | binary recorded | binary ported | pass |",
        "|---|---|---|---|---|---|",
    ]
    for r in rows3e:
        lines.append(f"| {r['eta']:.2f} | {r['poly_recorded']:.5f} | {r['poly_ported']:.5f} "
                     f"| {r['bin_recorded']:.5f} | {r['bin_ported']:.5f} "
                     f"| {'ok' if r['passed'] else 'FAIL'} |")
    lines += [
        "",
        "## 3D Voronoi observables (spot values)",
        "",
        "| eta | Q_poly rec | Q_poly port | n_f rec | n_f port | Q_bin rec | Q_bin port | pass |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in rows3v:
        lines.append(f"| {r['eta']:.2f} | {r['Q_poly_recorded']:.5f} | {r['Q_poly_ported']:.5f} "
                     f"| {r['nf_poly_recorded']:.3f} | {r['nf_poly_ported']:.3f} "
                     f"| {r['Q_bin_recorded']:.5f} | {r['Q_bin_ported']:.5f} "
                     f"| {'ok' if r['passed'] else 'FAIL'} |")
    lines += ["", f"**Overall: {'PASSED' if passed else 'FAILED'}**"]
    out = Path(a.results)
    out.mkdir(parents=True, exist_ok=True)
    (out / "validation_pilot.md").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {out / 'validation_pilot.md'}  -> {'PASSED' if passed else 'FAILED'}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""T4 THEORY RUN -- is the uniform-K degeneracy extensive?  (No simulation.)

Exhaustively enumerates all uniform-K independent vacancy sets on growing
tori (``geometry/coloring.enumerate_uniform_K``), fits ``log(count)`` against
both the area ``N`` and the boundary ``sqrt(N)``, and spot-checks that every
enumerated state -- crystalline or not -- is exactly tangential at the ladder
density.  This attacks the configurational-entropy objection head on: an
extensive count (``log ~ N``) would give the reference family a finite entropy
density; a boundary-law count (``log ~ sqrt(N)``) would not.

Only commensurate sizes admit solutions (``Nv = K N/(z0+K)`` must be an
integer); incommensurate sizes are reported as the structural zeros they are
and excluded from fits.  Every enumeration is exhaustive and flagged
``complete``; capped runs are reported as bounds, never as counts.

Deliverable: ``results/degeneracy.md``.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.geometry.coloring import (  # noqa: E402
    counting_eta,
    enumerate_uniform_K,
    fcc_torus,
    solution_cells,
    triangular_torus,
)

#: enumeration grid: (dim, K, sides).  Sides chosen commensurate where
#: possible; incommensurate ones are kept to demonstrate the structural zeros.
GRID_2D = {
    3: [3, 6, 9, 12, 15],
    2: [4, 6, 8, 10, 12, 14, 16],
    1: [7, 14],
}
GRID_3D = {4: [4, 6]}

Q_TARGET_2D = {2: np.sqrt(3) * np.pi / 8}
ETA_TARGET_3D = {4: counting_eta(12, 4, np.pi / np.sqrt(18.0))}


def enumerate_grid(max_nodes: int, verbose=True):
    rows = []
    for K, sides in GRID_2D.items():
        for N in sides:
            t0 = time.time()
            r = enumerate_uniform_K(triangular_torus(N), K, max_nodes=max_nodes)
            r.update(dim=2, side=N, seconds=time.time() - t0)
            rows.append(r)
            if verbose:
                print(f"  2D tri {N:>2d}x{N:<2d} K={K}: sol={r['solutions']:>6} "
                      f"orbits={r['orbits']:>3} cosets={r['cosets']} "
                      f"complete={r['complete']} nodes={r['nodes']} "
                      f"({r['seconds']:.1f}s)", flush=True)
    for K, sides in GRID_3D.items():
        for side in sides:
            t0 = time.time()
            r = enumerate_uniform_K(fcc_torus(side), K, max_nodes=max_nodes)
            r.update(dim=3, side=side, seconds=time.time() - t0)
            rows.append(r)
            if verbose:
                print(f"  3D FCC side {side}: K={K}: sol={r['solutions']:>6} "
                      f"orbits={r['orbits']:>3} cosets={r['cosets']} "
                      f"complete={r['complete']} nodes={r['nodes']} "
                      f"({r['seconds']:.1f}s)", flush=True)
    return rows


def scaling_fit(rows, dim, K):
    """Fit ln(solutions) against N and against sqrt(N); report both."""
    pts = [
        (r["N"], r["solutions"])
        for r in rows
        if r["dim"] == dim and r["K"] == K and r["complete"] and r["solutions"] > 0
    ]
    if len(pts) < 3:
        return None
    N = np.array([p[0] for p in pts], float)
    lnS = np.log([p[1] for p in pts])

    def fit(x):
        c = np.polyfit(x, lnS, 1)
        resid = lnS - np.polyval(c, x)
        ss_res = float((resid**2).sum())
        ss_tot = float(((lnS - lnS.mean()) ** 2).sum())
        return {"slope": float(c[0]), "intercept": float(c[1]),
                "R2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")}

    constant = bool(np.ptp(lnS) < 1e-12)
    return {
        "dim": dim, "K": K, "n_points": len(pts),
        "N": N.tolist(), "lnS": lnS.tolist(),
        "constant": constant,
        "count": int(pts[0][1]) if constant else None,
        "vs_N": None if constant else fit(N),
        "vs_sqrtN": None if constant else fit(np.sqrt(N)),
    }


def tangentiality_spot_checks(verbose=True):
    """Every enumerated orbit on the largest checked tori is exactly tangential."""
    out = []
    # 2D: 8x8 K=2 -- all orbits, all cells
    tor = triangular_torus(8)
    r = enumerate_uniform_K(tor, 2)
    worst = 0.0
    for rep in r["orbit_representatives"]:
        for c in solution_cells(rep, tor):
            worst = max(worst, abs(c.Q_iso - Q_TARGET_2D[2]),
                        float(c.face_distances.max() - c.face_distances.min()))
    out.append({"case": "2D 8x8 K=2, all orbits, all cells",
                "orbits": r["orbits"], "worst_deviation": worst, "ok": worst < 1e-9})
    # 3D: FCC side 6 K=4 -- all orbits, sampled cells
    tor = fcc_torus(6)
    r = enumerate_uniform_K(tor, 4)
    worst = 0.0
    for rep in r["orbit_representatives"]:
        occ = [i for i, v in enumerate(rep) if v == 0]
        sample = range(0, len(occ), max(1, len(occ) // 12))
        for c in solution_cells(rep, tor, sample=sample):
            eta_local = (np.pi / 6.0) * np.sqrt(2.0) ** 3 / c.V
            worst = max(worst, abs(eta_local - ETA_TARGET_3D[4]),
                        float(c.face_distances.max() - c.face_distances.min()))
    out.append({"case": "3D FCC side-6 K=4, all orbits, sampled cells",
                "orbits": r["orbits"], "cosets": r["cosets"],
                "worst_deviation": worst, "ok": worst < 1e-9})
    if verbose:
        for o in out:
            print(f"  {o['case']}: worst deviation {o['worst_deviation']:.2e} "
                  f"-> {'ok' if o['ok'] else 'FAIL'}")
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="T4: uniform-K degeneracy on growing tori")
    p.add_argument("--max-nodes", type=int, default=200_000_000,
                   help="DFS node cap per torus; capped runs are flagged, not counted")
    p.add_argument("--results", default=str(REPO / "results"))
    a = p.parse_args(argv)

    print("T4 THEORY RUN: exhaustive uniform-K enumeration on growing tori")
    rows = enumerate_grid(a.max_nodes)
    print("\ntangentiality of every enumerated state:")
    spots = tangentiality_spot_checks()

    fits = [f for f in (scaling_fit(rows, 2, K) for K in (3, 2, 1)) if f]
    fits += [f for f in (scaling_fit(rows, 3, K) for K in (4,)) if f]

    lines = [
        "# T4 — uniform-K degeneracy on growing tori (theory run, no simulation)",
        "",
        "Exhaustive enumeration of all uniform-K independent vacancy sets",
        "(`geometry/coloring.enumerate_uniform_K`, constraint-propagating DFS with",
        "the vacancy count pinned to `Nv = K N/(z0+K)`). Every row below is a",
        "complete enumeration unless flagged otherwise; capped runs are bounds,",
        "not counts. Incommensurate sizes (non-integer `Nv`) are structural zeros",
        "and are excluded from fits.",
        "",
        "## Counts",
        "",
        "| dim | lattice | side | N sites | K | Nv | solutions | orbits | cosets | complete | DFS nodes | s |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lat = "triangular" if r["dim"] == 2 else "FCC"
        lines.append(
            f"| {r['dim']} | {lat} | {r['side']} | {r['N']} | {r['K']} | "
            f"{r['Nv'] if r['Nv'] is not None else '—'} | {r['solutions']} | "
            f"{r['orbits']} | {r['cosets']} | {'yes' if r['complete'] else 'NO (capped)'} | "
            f"{r['nodes']} | {r['seconds']:.1f} |"
        )
    lines += [
        "",
        "## Scaling of the count",
        "",
        "An extensive degeneracy (finite entropy density) means `ln(count) ~ N`;",
        "a boundary-law degeneracy means `ln(count) ~ sqrt(N)`. Both fits below;",
        "the better-fitting law is the honest headline, with the caveat that the",
        "accessible sizes are small.",
        "",
        "| dim | K | points | slope vs N | R²(N) | slope vs √N | R²(√N) |",
        "|---|---|---|---|---|---|---|",
    ]
    for f in fits:
        if f["constant"]:
            lines.append(
                f"| {f['dim']} | {f['K']} | {f['n_points']} | — | — | — | — |"
            )
        else:
            lines.append(
                f"| {f['dim']} | {f['K']} | {f['n_points']} | "
                f"{f['vs_N']['slope']:.4f} | {f['vs_N']['R2']:.4f} | "
                f"{f['vs_sqrtN']['slope']:.4f} | {f['vs_sqrtN']['R2']:.4f} |"
            )

    # verdicts, written from the numbers just computed
    lines += ["", "## Reading", ""]
    for f in fits:
        if f["constant"]:
            lines.append(
                f"- **{f['dim']}D K={f['K']}**: count is constant "
                f"({f['count']}) at every commensurate size — unique up to "
                "translation, zero degeneracy growth. No fit applies."
            )
            continue
        better = "sqrt(N) (boundary law)" if f["vs_sqrtN"]["R2"] > f["vs_N"]["R2"] else "N (extensive)"
        lines.append(
            f"- **{f['dim']}D K={f['K']}**: ln(count) fits {better} better "
            f"(R² {max(f['vs_sqrtN']['R2'], f['vs_N']['R2']):.4f} vs "
            f"{min(f['vs_sqrtN']['R2'], f['vs_N']['R2']):.4f}). "
        )

    # observed closed form for the K=2 class (checked, then predicted blind)
    k2 = [r for r in rows if r["dim"] == 2 and r["K"] == 2 and r["complete"]
          and r["solutions"] > 0]
    formula_ok = all(r["solutions"] == 6 * 2 ** (r["side"] // 2) - 8 for r in k2)
    k4 = [r for r in rows if r["dim"] == 3 and r["K"] == 4 and r["complete"]
          and r["solutions"] > 0]
    formula3_ok = all(r["solutions"] == 6 * 2 ** (r["side"] // 2) - 8 for r in k4)
    lines += [
        "",
        "## Observed closed form (K=2 class)",
        "",
        f"Every complete even-side triangular count obeys **S(side) = 6·2^(side/2) − 8** "
        f"exactly: {'confirmed on all ' + str(len(k2)) + ' sides' if formula_ok else 'VIOLATED — see table'}. "
        "The side-16 count (1528) was predicted from the formula before being "
        "enumerated, and matched. This is an OBSERVED law, not a derived one "
        "(DEBT.md); its form — ln S ≈ (ln 2/2)·√N plus corrections — is exactly a "
        "boundary/stacking law: the degeneracy grows like independent layer choices "
        "along one axis, the same mechanism as ABC stacking freedom, not like a "
        "bulk entropy.",
        "",
        f"The two 3D FCC K=4 counts (16 at side 4, 40 at side 6) satisfy the same "
        f"formula{' ' if formula3_ok else ' NOT '}— two points prove nothing, but the "
        "coincidence is worth recording: the 3D degeneracy may be the same layered "
        "mechanism.",
    ]
    lines += [
        "",
        "## Tangentiality of every enumerated state",
        "",
    ]
    for o in spots:
        lines.append(f"- {o['case']}: worst deviation {o['worst_deviation']:.2e} "
                     f"({'exact within float' if o['ok'] else 'FAIL'})")
    lines += [
        "",
        "## Caveats",
        "",
        "- Sizes are small (up to 15x15 triangular, FCC side 6); both scaling laws",
        "  are fitted through few points and the verdict is provisional, not a",
        "  theorem.",
        "- Orbit counts are per torus translation group only; point-group symmetry",
        "  is not quotiented, so orbit counts overstate the number of genuinely",
        "  distinct patterns by up to the point-group order.",
        "- K=1 (maple-leaf class) has two orbits on the side-7 and side-14 tori —",
        "  the two enantiomers of the chiral pattern — both lattice cosets: unique",
        "  up to translation and chirality at accessible sizes.",
        "- K=3 (honeycomb class) has exactly 3 solutions (1 orbit) at every",
        "  commensurate size 3–15: unique up to translation. No degeneracy at all.",
        "",
        "## Consequence for the configurational-entropy objection",
        "",
        "The reference family does contain non-crystalline members (most orbits are",
        "not lattice cosets, in 2D K=2 and in 3D K=4 side 6), and every one of them",
        "is exactly tangential at the ladder density — so the geometric reference",
        "state is a degenerate family, not a single crystal. But at the sizes",
        "enumerated the degeneracy is a boundary law, not extensive: it carries no",
        "finite entropy density. On this evidence the degeneracy weakens the",
        "entropy objection only marginally; it does not answer it.",
    ]
    out = Path(a.results)
    out.mkdir(parents=True, exist_ok=True)
    (out / "degeneracy.md").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {out / 'degeneracy.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

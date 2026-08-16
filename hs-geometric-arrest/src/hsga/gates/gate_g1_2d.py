"""Gate G1-2D (task T5) -- the 2D engine against known physics.

Pass criteria:

* monodisperse EOS sanity: ``Z = 1 + 2 phi g(sigma+)`` within 2% of Henderson
  at ``phi = 0.55`` and ``0.60`` (Henderson itself is good to ~1% there; the
  criterion honestly reflects both errors and is fixed here, not tuned);
* the reference-composition binary (R^-1 = 1.4) at ``phi = 0.80``: anneal
  reaches exactly zero overlap, the final audit is zero, and the shell
  percolation threshold is finite -- the full production path of the 2D
  campaign works at its densest state point.

Everything is small (N ~ 500, seconds per run); the gate is meant to run on a
laptop every time the engine changes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ..analysis.eos import compressibility, contact_value, radial_distribution
from ..analysis.percolation import eps_star, eps_star_bisect
from ..engine.driver import RunSpec, build_engine, read_frames, run_sweep, write_manifest

PHIS_EOS = (0.55, 0.60)
PHI_DENSE = 0.80
TOL = 0.02


def run(
    *,
    data_dir: str = "data/gate_g1_2d",
    ncell: int = 12,
    eq: int = 5000,
    prod: int = 60000,
    nsnap: int = 100,
    seeds: tuple[int, ...] = (1, 2, 3, 4),
    nproc: int | None = None,
    verbose: bool = True,
) -> dict:
    exe = build_engine(name="hsmc2d")
    Path(data_dir).mkdir(parents=True, exist_ok=True)

    specs = [
        RunSpec(eta=phi, prefix=f"{data_dir}/mono_p{phi:.2f}_s{s}", dim=2,
                ncell=ncell, mode=1, seed=s, eq=eq, prod=prod, nsnap=nsnap, melt=eq)
        for phi in PHIS_EOS for s in seeds
    ] + [
        RunSpec(eta=PHI_DENSE, prefix=f"{data_dir}/bin_p{PHI_DENSE:.2f}_s{s}", dim=2,
                ncell=ncell, mode=0, seed=s, eq=eq, prod=prod // 4, nsnap=6, melt=eq)
        for s in seeds[:2]
    ]
    runs = run_sweep(specs, exe, nproc=nproc)
    write_manifest(f"{data_dir}/manifest.json", runs, gate="G1-2D",
                   phis=list(PHIS_EOS), phi_dense=PHI_DENSE, tol=TOL)

    rows, ok = [], True
    for phi in PHIS_EOS:
        frames = []
        for s in seeds:
            frames += read_frames(f"{data_dir}/mono_p{phi:.2f}_s{s}.cfg", dim=2)
        gr = radial_distribution(frames, dim=2)
        cv = contact_value(gr["r"], gr["g"], gr["sigma"], counts=gr["counts"])
        comp = compressibility(phi, cv["g_contact"], cv["err"], dim=2)
        good = comp["rel_dev"] < TOL
        ok &= good
        rows.append({**comp, "g_contact": cv["g_contact"],
                     "g_err_window": cv["err_window"], "n_frames": len(frames),
                     "passed": good})
        if verbose:
            print(f"  phi={phi:.2f}  g(sigma+)={cv['g_contact']:.5f} "
                  f"(window sys {cv['err_window']:.5f})  Z={comp['Z']:.5f}  "
                  f"Z_Henderson={comp['Z_ref']:.5f}  "
                  f"dev={comp['rel_dev'] * 100:.3f}%  {'ok' if good else 'FAIL'}")

    # dense binary: audits + a working percolation threshold, both estimators
    dense_ok = True
    eps_vals = []
    for s in seeds[:2]:
        log = next(r["log"] for r in runs
                   if r["spec"]["prefix"].endswith(f"bin_p{PHI_DENSE:.2f}_s{s}"))
        dense_ok &= log.get("anneal_final_energy", 1.0) == 0.0
        dense_ok &= log.get("final_overlap_audit", 1) == 0
        frames = read_frames(f"{data_dir}/bin_p{PHI_DENSE:.2f}_s{s}.cfg", dim=2)
        pos, rad, L = frames[-1]
        e_exact = eps_star(pos, rad, L)
        e_bis = eps_star_bisect(pos, rad, L)
        eps_vals.append(e_exact)
        dense_ok &= np.isfinite(e_exact) and abs(e_exact - e_bis) / e_exact < 1e-3
    ok &= dense_ok
    if verbose:
        print(f"  binary R^-1=1.4 at phi={PHI_DENSE}: zero-overlap anneal + audit, "
              f"eps* = {', '.join(f'{e:.5f}' for e in eps_vals)} "
              f"(exact == bisect)  {'ok' if dense_ok else 'FAIL'}")
        print(f"  GATE G1-2D {'PASSED' if ok else 'FAILED'}")
    return {"gate": "G1-2D", "passed": bool(ok), "rows": rows,
            "eps_star_dense": eps_vals}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Gate G1-2D: 2D EOS sanity + phi=0.80 audits")
    p.add_argument("--data-dir", default="data/gate_g1_2d")
    p.add_argument("--nproc", type=int, default=None)
    a = p.parse_args(argv)
    print("GATE G1-2D: hsmc2d vs Henderson EOS; zero overlaps at phi=0.80")
    return 0 if run(data_dir=a.data_dir, nproc=a.nproc)["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

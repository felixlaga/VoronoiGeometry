"""Gate G4 -- finite size (IMPLEMENTED; NOT EXECUTED at gate scale).

Pass criterion (spec Sec. 3): ``N = 864 / 4000 / 10^4`` at three densities;
``eps*(eta)`` consistent within error (2 combined sems between consecutive
sizes).  ``N = 4 ncell^3`` gives ``ncell = 6, 10, 14`` -> ``864, 4000, 10976``.

Percolation thresholds drift with box size on general grounds, so this is a
real test: if the drift between N=4000 and N=10976 exceeds the error, any
feature in ``eps*(eta)`` must be established at fixed N and via the
finite-size crossing analysis (``analysis.pwrap.finite_size_crossings``),
never quoted from one size.

Cost: the N=10976 runs are hours each; NOT run in this repository (DEBT.md).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ..analysis.percolation import eps_star
from ..engine.driver import RunSpec, build_engine, read_frames, run_sweep, write_manifest

ETAS = (0.50, 0.56, 0.60)
NCELLS = (6, 10, 14)          # N = 864, 4000, 10976


def run(
    *,
    data_dir: str = "data/gate_s4",
    mode: int = 0,
    eq: int = 30000,
    prod: int = 40000,
    nsnap: int = 10,
    seeds: tuple[int, ...] = (1, 2, 3, 4),
    etas=ETAS,
    ncells=NCELLS,
    nproc: int | None = None,
    verbose: bool = True,
) -> dict:
    exe = build_engine(name="hsmc")
    Path(data_dir).mkdir(parents=True, exist_ok=True)
    specs = [
        RunSpec(eta=eta, prefix=f"{data_dir}/n{nc}_e{eta:.2f}_s{s}", dim=3,
                ncell=nc, mode=mode, seed=s, eq=eq, prod=prod, nsnap=nsnap,
                melt=20000)
        for eta in etas for nc in ncells for s in seeds
    ]
    runs = run_sweep(specs, exe, nproc=nproc)
    write_manifest(f"{data_dir}/manifest.json", runs, gate="G4",
                   etas=list(etas), ncells=list(ncells))

    rows, ok = [], True
    for eta in etas:
        summary = {}
        for nc in ncells:
            per_seed = []
            for s in seeds:
                frames = read_frames(f"{data_dir}/n{nc}_e{eta:.2f}_s{s}.cfg")
                es = [eps_star(p_, r_, L_) for p_, r_, L_ in frames]
                es = [e for e in es if np.isfinite(e)]
                if es:
                    per_seed.append(float(np.mean(es)))
            summary[nc] = {
                "N": 4 * nc**3,
                "mean": float(np.mean(per_seed)),
                "sem": float(np.std(per_seed, ddof=1) / np.sqrt(len(per_seed))),
            }
        for a, b in zip(ncells[:-1], ncells[1:]):
            d = abs(summary[a]["mean"] - summary[b]["mean"])
            se = float(np.hypot(summary[a]["sem"], summary[b]["sem"]))
            good = d <= 2.0 * se
            ok &= good
            rows.append({"eta": eta, "N_a": summary[a]["N"], "N_b": summary[b]["N"],
                         "eps_a": summary[a]["mean"], "eps_b": summary[b]["mean"],
                         "diff": d, "combined_sem": se, "passed": bool(good)})
            if verbose:
                print(f"  eta={eta:.2f}  N {summary[a]['N']:>6d}->{summary[b]['N']:<6d} "
                      f"eps* {summary[a]['mean']:.5f}->{summary[b]['mean']:.5f} "
                      f"|d|={d:.5f} ({d / se if se > 0 else float('inf'):.1f} sigma) "
                      f"{'ok' if good else 'FAIL'}")
    if verbose:
        print(f"  GATE G4 {'PASSED' if ok else 'FAILED'}")
    return {"gate": "G4", "passed": bool(ok), "rows": rows}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Gate G4: eps* vs system size")
    p.add_argument("--data-dir", default="data/gate_s4")
    p.add_argument("--mode", type=int, default=0)
    p.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4])
    p.add_argument("--ncells", type=int, nargs="+", default=list(NCELLS))
    p.add_argument("--nproc", type=int, default=None)
    a = p.parse_args(argv)
    print("GATE G4: eps*(eta) vs system size (N = 864 / 4000 / 10976)")
    res = run(data_dir=a.data_dir, mode=a.mode, seeds=tuple(a.seeds),
              ncells=tuple(a.ncells), nproc=a.nproc)
    return 0 if res["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

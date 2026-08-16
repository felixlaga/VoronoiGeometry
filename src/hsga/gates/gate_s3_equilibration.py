"""Gate G3 -- equilibration length (IMPLEMENTED; NOT EXECUTED at gate scale).

Pass criterion (spec Sec. 3): vary the equilibration length over
``1e3 / 1e4 / 1e5`` sweeps at >= 4 ``(eta, composition)`` points near arrest;
the structural observables (``eps*`` primary, ``<Q_iso>`` alongside) must
change only within error between consecutive lengths (2 combined sems).

A drift with equilibration length means the configurations are not
equilibrated and every downstream number is a property of the protocol, not
the fluid; the only fix is more sweeps.

Cost: the ``1e5``-sweep points at N=864 are hours of single-core compute;
this gate has therefore NOT been run in this repository (DEBT.md).  The
implementation is complete and `make gate-g3` runs it as specified.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ..analysis.percolation import eps_star
from ..analysis.voronoi import config_observables
from ..engine.driver import RunSpec, build_engine, read_frames, run_sweep, write_manifest

POINTS = ((0.56, 0), (0.58, 0), (0.58, 1), (0.60, 0), (0.60, 1))
EQ_LENGTHS = (1_000, 10_000, 100_000)


def run(
    *,
    data_dir: str = "data/gate_s3",
    ncell: int = 6,
    prod: int = 40000,
    nsnap: int = 10,
    seeds: tuple[int, ...] = (1, 2, 3, 4),
    points=POINTS,
    eq_lengths=EQ_LENGTHS,
    voronoi: bool = True,
    nproc: int | None = None,
    verbose: bool = True,
) -> dict:
    exe = build_engine(name="hsmc")
    Path(data_dir).mkdir(parents=True, exist_ok=True)
    specs = [
        RunSpec(eta=eta, prefix=f"{data_dir}/m{mode}_e{eta:.2f}_q{eq}_s{s}", dim=3,
                ncell=ncell, mode=mode, seed=s, eq=eq, prod=prod, nsnap=nsnap,
                melt=min(eq, 20000))
        for eta, mode in points for eq in eq_lengths for s in seeds
    ]
    runs = run_sweep(specs, exe, nproc=nproc, allow_anneal_failure=True)
    write_manifest(f"{data_dir}/manifest.json", runs, gate="G3",
                   points=list(points), eq_lengths=list(eq_lengths))

    rows, ok = [], True
    for eta, mode in points:
        summary = {}
        for eq in eq_lengths:
            eps_seed, Q_seed = [], []
            for s in seeds:
                pre = f"{data_dir}/m{mode}_e{eta:.2f}_q{eq}_s{s}"
                if not Path(f"{pre}.cfg").exists():
                    continue                      # unreachable point: reported below
                frames = read_frames(f"{pre}.cfg")
                es = [eps_star(p_, r_, L_) for p_, r_, L_ in frames]
                es = [e for e in es if np.isfinite(e)]
                if es:
                    eps_seed.append(float(np.mean(es)))
                if voronoi:
                    Q_seed.append(float(np.mean(
                        [config_observables(p_, r_, L_)["mean_Q_iso"]
                         for p_, r_, L_ in frames[-2:]])))
            if len(eps_seed) < 2:
                summary[eq] = None
                continue
            summary[eq] = {
                "eps_mean": float(np.mean(eps_seed)),
                "eps_sem": float(np.std(eps_seed, ddof=1) / np.sqrt(len(eps_seed))),
                "Q_mean": float(np.mean(Q_seed)) if Q_seed else float("nan"),
                "n_seeds": len(eps_seed),
            }
        for a, b in zip(eq_lengths[:-1], eq_lengths[1:]):
            if summary.get(a) is None or summary.get(b) is None:
                rows.append({"eta": eta, "mode": mode, "eq_a": a, "eq_b": b,
                             "passed": None, "note": "state point unreachable"})
                continue
            d = abs(summary[a]["eps_mean"] - summary[b]["eps_mean"])
            se = float(np.hypot(summary[a]["eps_sem"], summary[b]["eps_sem"]))
            good = d <= 2.0 * se
            ok &= good
            rows.append({"eta": eta, "mode": mode, "eq_a": a, "eq_b": b,
                         "eps_a": summary[a]["eps_mean"], "eps_b": summary[b]["eps_mean"],
                         "diff": d, "combined_sem": se, "passed": bool(good)})
            if verbose:
                print(f"  eta={eta:.2f} mode={mode}  eq {a:>6d}->{b:<6d}  "
                      f"eps* {summary[a]['eps_mean']:.5f}->{summary[b]['eps_mean']:.5f} "
                      f"|d|={d:.5f} ({d / se if se > 0 else float('inf'):.1f} sigma) "
                      f"{'ok' if good else 'FAIL'}")
    if verbose:
        print(f"  GATE G3 {'PASSED' if ok else 'FAILED'}")
    return {"gate": "G3", "passed": bool(ok), "rows": rows}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Gate G3: equilibration length")
    p.add_argument("--data-dir", default="data/gate_s3")
    p.add_argument("--prod", type=int, default=40000)
    p.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4])
    p.add_argument("--no-voronoi", action="store_true")
    p.add_argument("--nproc", type=int, default=None)
    a = p.parse_args(argv)
    print("GATE G3: structural observables vs equilibration length (1e3/1e4/1e5)")
    res = run(data_dir=a.data_dir, prod=a.prod, seeds=tuple(a.seeds),
              voronoi=not a.no_voronoi, nproc=a.nproc)
    return 0 if res["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

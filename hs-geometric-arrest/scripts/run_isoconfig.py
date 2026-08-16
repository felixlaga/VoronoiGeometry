#!/usr/bin/env python3
"""T10 driver -- isoconfigurational propensity and the held-out baseline test.

For each selected equilibrated configuration this launches ``M`` production
runs from the SAME configuration with different seeds (the engines' restart
mode; swap moves off throughout, since propensity is a dynamical quantity),
reads per-particle displacements at the snapshot lags, and evaluates whether
the structural fields {s_k, marked-cluster size, persistence} predict
propensity beyond the baseline {local eta, tetrahedrality/psi6, anisotropy,
V_cell} on held-out configurations (``analysis.isoconfig.held_out_r2``).

If the structural fields add nothing, the hypothesis fails de Graaf's own
standard -- that verdict is printed verbatim, not softened.

Cost: M x (configurations) full production runs; the real evaluation needs
the T6/T9 campaign data and has NOT been executed in this repository
(DEBT.md).  ``--smoke`` runs the identical path at toy scale.
"""

from __future__ import annotations

import argparse
import sys
import zlib
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.analysis.isoconfig import (  # noqa: E402
    baseline_features, held_out_r2, propensity_from_frames, structural_features,
)
from hsga.analysis.voronoi import config_cells  # noqa: E402
from hsga.engine.driver import (  # noqa: E402
    RunSpec, build_engine, read_frames, run_sweep, write_manifest,
)


def propensity_for_config(cfg_path: str, dim: int, *, M: int, prod: int,
                          nsnap: int, out_dir: Path, nproc, tag: str):
    """Launch M restarts from the last frame of ``cfg_path``; return propensity."""
    exe = build_engine(name="hsmc" if dim == 3 else "hsmc2d")
    frames = read_frames(cfg_path, dim=dim)
    start_pos, start_rad, L = frames[-1]
    specs = [
        RunSpec(eta=0.0, prefix=str(out_dir / f"{tag}_m{m}"), dim=dim,
                seed=zlib.crc32(f"iso|{tag}|{m}".encode()) % (2**31),
                eq=0, prod=prod, nsnap=nsnap, swap=0,
                infile=cfg_path, inframe=len(frames) - 1)
        for m in range(M)
    ]
    runs = run_sweep(specs, exe, nproc=nproc)
    ensembles = [read_frames(str(out_dir / f"{tag}_m{m}.cfg"), dim=dim)
                 for m in range(M)]
    prop = propensity_from_frames(start_pos, ensembles, L)
    return start_pos, start_rad, L, prop, runs


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Isoconfigurational propensity (T10)")
    p.add_argument("--configs", nargs="+",
                   help=".cfg files of equilibrated configurations (campaign output)")
    p.add_argument("--dim", type=int, choices=(2, 3), default=2)
    p.add_argument("--M", type=int, default=32, help="restarts per configuration")
    p.add_argument("--prod", type=int, default=40000)
    p.add_argument("--nsnap", type=int, default=6, help="lag times per restart")
    p.add_argument("--lag-index", type=int, default=-1,
                   help="which lag's propensity to evaluate (default: longest)")
    p.add_argument("--data-dir", default=str(REPO / "data" / "isoconfig"))
    p.add_argument("--results", default=str(REPO / "results"))
    p.add_argument("--nproc", type=int, default=None)
    p.add_argument("--smoke", action="store_true",
                   help="toy-scale end-to-end wiring check (generates its own configs)")
    a = p.parse_args(argv)

    out_dir = Path(a.data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if a.smoke:
        exe = build_engine(name="hsmc2d")
        cfgs = []
        for k in range(3):
            spec = RunSpec(eta=0.74, prefix=str(out_dir / f"smokebase_{k}"), dim=2,
                           ncell=6, mode=0, seed=100 + k, eq=800, prod=1500, nsnap=2,
                           melt=800)
            run_sweep([spec], exe, nproc=1)
            cfgs.append(f"{spec.prefix}.cfg")
        a.configs, a.dim, a.M, a.prod, a.nsnap = cfgs, 2, 6, 1500, 3
    if not a.configs:
        p.error("--configs is required (or --smoke)")

    blocks, props, all_runs = [], [], []
    for cfg in a.configs:
        tag = Path(cfg).stem
        pos, rad, L, prop, runs = propensity_for_config(
            cfg, a.dim, M=a.M, prod=a.prod, nsnap=a.nsnap,
            out_dir=out_dir, nproc=a.nproc, tag=tag)
        all_runs += runs
        cells = config_cells(pos, rad, L)
        Xb, base_names = baseline_features(cells, pos, rad, L)
        Xs, struct_names = structural_features(cells, pos, rad, L)
        blocks.append((Xb, Xs))
        props.append(prop[a.lag_index])
        print(f"  {tag}: N={len(pos)}  M={a.M}  "
              f"<propensity>={prop[a.lag_index].mean():.4f}")

    res = held_out_r2(blocks, props)
    print(f"\nheld-out delta R^2 = {res['delta_r2_mean']:+.4f} "
          f"+- {res['delta_r2_spread']:.4f}  "
          f"(baseline: {base_names}; structural: {struct_names})")
    print(res["verdict_note"])

    write_manifest(out_dir / "manifest.json", all_runs, task="T10",
                   configs=a.configs, M=a.M, lag_index=a.lag_index)
    out = Path(a.results)
    out.mkdir(parents=True, exist_ok=True)
    import json

    (out / ("isoconfig_smoke.json" if a.smoke else "isoconfig.json")).write_text(
        json.dumps({"delta_r2_mean": res["delta_r2_mean"],
                    "delta_r2_spread": res["delta_r2_spread"],
                    "folds": res["folds"], "baseline": base_names,
                    "structural": struct_names, "smoke": a.smoke,
                    "verdict_note": res["verdict_note"]}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

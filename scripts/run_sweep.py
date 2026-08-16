#!/usr/bin/env python3
"""Generate campaign configurations (2D replication first; 3D gated on it).

Spec Sec. 4.  The 2D replication (T6) is the decision node for every 3D
physics claim: N = 2016 (ny=24, nx=42 -- the closest ncell realisation of the
reference study's N=2048), >= 350 configurations per density, grid
``Delta_phi = 5e-4`` within ``+-0.006`` of each exact target

    0.604600 (honeycomb, K=3), 0.680175 (kagome, K=2 -- the novel rung),
    0.777343 (maple-leaf, K=1 -- de Graaf's value),

coarse elsewhere.  The 3D campaign (T9) runs ONLY if T6 produces a feature;
its grid is coarse over [0.40, 0.70] plus the same fine windows around every
admissible value.

Statistics per state point are ``replicas x nsnap``; only the replica axis is
independent and both counts are recorded.  State points at or beyond their
jamming density fail their anneal, are recorded as unreachable, produce no
data, and are never interpolated over.

HEAVY: the full presets are cluster-scale.  ``--dry-run`` prints the plan and
an honest cost estimate; the ``*-smoke`` presets exercise the identical code
path at toy scale.  Neither full preset has been executed in this repository
(DEBT.md, Execution status).
"""

from __future__ import annotations

import argparse
import sys
import time
import zlib
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.engine.driver import (  # noqa: E402
    MODES_2D, MODES_3D, RunSpec, build_engine, run_sweep, write_manifest,
)

#: exact fine-grid targets
TARGETS_2D = [0.6045997880780726, 0.6801747615878316, 0.7773425846718076]
TARGETS_3D = [0.5235987755982988, 0.5553603672697957, 0.5923843917544487,
              0.6045997880780726, 0.6346975625940521, 0.6835204520243637]


def grid_with_fine_windows(lo, hi, coarse, targets, fine=5e-4, half=0.006):
    """Coarse grid over [lo, hi] plus a fine grid within +-half of each target."""
    pts = set(np.round(np.arange(lo, hi + 1e-12, coarse), 6))
    for t in targets:
        if lo - half <= t <= hi + half:
            fine_pts = np.arange(t - half, t + half + 1e-12, fine)
            pts.update(np.round(fine_pts, 6))
    return sorted(p for p in pts if lo - half - 1e-9 <= p <= hi + half + 1e-9)


PRESETS = {
    # T6 -- the decision node. ny=24 -> nx=42 -> N=2016.
    "2d-replication": dict(
        dim=2, ncell=24, modes=[0], replicas=44, nsnap=8,
        eq=60_000, prod=120_000,
        grid=lambda: grid_with_fine_windows(0.58, 0.82, 0.01, TARGETS_2D),
    ),
    "2d-smoke": dict(
        dim=2, ncell=8, modes=[0], replicas=2, nsnap=2, eq=1000, prod=2000,
        grid=lambda: [0.70, 0.775, 0.80],
    ),
    # T9 -- gated on T6. N=4000 (ncell=10); the G4 sizes run separately.
    "3d-campaign": dict(
        dim=3, ncell=10, modes=[0, 1, 2], replicas=50, nsnap=10,
        eq=100_000, prod=100_000,
        grid=lambda: grid_with_fine_windows(0.40, 0.70, 0.01, TARGETS_3D),
    ),
    "3d-smoke": dict(
        dim=3, ncell=4, modes=[0], replicas=2, nsnap=2, eq=1000, prod=2000,
        grid=lambda: [0.45, 0.55],
    ),
}


def build_specs(a, preset) -> list[RunSpec]:
    grid = preset["grid"]()
    specs = []
    for mode in preset["modes"]:
        for eta in grid:
            for rep in range(a.replicas or preset["replicas"]):
                seed = zlib.crc32(f"{a.preset}|{mode}|{eta:.6f}|{rep}".encode()) % (2**31)
                specs.append(RunSpec(
                    eta=float(eta),
                    prefix=f"{a.data_dir}/m{mode}/e{eta:.6f}_r{rep}",
                    dim=preset["dim"], ncell=a.ncell or preset["ncell"],
                    mode=mode, seed=seed,
                    eq=a.eq or preset["eq"], prod=a.prod or preset["prod"],
                    nsnap=a.nsnap or preset["nsnap"],
                    melt=min(a.eq or preset["eq"], 20_000),
                ))
    return specs


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Campaign configuration generator")
    p.add_argument("--preset", choices=sorted(PRESETS), required=True)
    p.add_argument("--data-dir", default=None)
    p.add_argument("--replicas", type=int, default=None)
    p.add_argument("--nsnap", type=int, default=None)
    p.add_argument("--ncell", type=int, default=None)
    p.add_argument("--eq", type=int, default=None)
    p.add_argument("--prod", type=int, default=None)
    p.add_argument("--nproc", type=int, default=None)
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args(argv)

    preset = PRESETS[a.preset]
    if a.data_dir is None:
        a.data_dir = str(REPO / "data" / a.preset)
    specs = build_specs(a, preset)
    n_cfg = (a.replicas or preset["replicas"]) * (a.nsnap or preset["nsnap"])
    grid = preset["grid"]()
    sweeps = sum(s.N * (s.eq + s.prod + s.melt) for s in specs)

    print(f"preset       : {a.preset} (dim={preset['dim']})")
    print(f"modes        : {preset['modes']} "
          f"({', '.join((MODES_2D if preset['dim'] == 2 else MODES_3D)[m] for m in preset['modes'])})")
    print(f"densities    : {len(grid)} points in [{grid[0]:.4f}, {grid[-1]:.4f}] "
          f"(fine windows of 5e-4 around the exact targets)")
    print(f"N            : {specs[0].N}")
    print(f"per point    : {a.replicas or preset['replicas']} replicas x "
          f"{a.nsnap or preset['nsnap']} snapshots = {n_cfg} configurations "
          f"(only replicas are independent)")
    if n_cfg < 350 and "smoke" not in a.preset:
        print(f"  WARNING: {n_cfg} configurations/point is below the calibrated "
              "requirement of ~350-500 (ALTERNATIVE_ROUTES.md)")
    print(f"total runs   : {len(specs)}")
    print(f"particle-sweeps total: {sweeps:.3e} "
          f"(~{sweeps / 5e9:.0f} core-hours at the measured ~1.4e6 particle-sweeps/s)")
    if a.dry_run:
        print("\n(dry run: nothing executed)")
        return 0

    exe = build_engine(name="hsmc2d" if preset["dim"] == 2 else "hsmc")
    t0 = time.time()
    runs = run_sweep(specs, exe, nproc=a.nproc, allow_anneal_failure=True)
    dt = time.time() - t0
    write_manifest(f"{a.data_dir}/manifest.json", runs,
                   preset=a.preset, grid=[float(g) for g in grid],
                   configurations_per_point=n_cfg, wall_seconds=dt)
    ok = [r for r in runs if not r["unreachable"]]
    bad = [r for r in runs if r["unreachable"]]
    print(f"\n{len(ok)}/{len(runs)} runs completed in {dt / 60:.1f} min; "
          f"overlap audits all zero: "
          f"{all(r['log'].get('final_overlap_audit') == 0 for r in ok)}")
    if bad:
        pts = sorted({(r['spec']['mode'], r['spec']['eta']) for r in bad})
        print(f"{len(bad)} runs were at unreachable state points (at/beyond jamming); "
              "no data written for:")
        for mode, eta in pts:
            print(f"    mode {mode}  eta={eta:.6f}")
        print("The analysis reports these as missing; it never interpolates over them.")
    print(f"manifest: {a.data_dir}/manifest.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())

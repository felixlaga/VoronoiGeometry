#!/usr/bin/env python3
"""Kagome seeded-stability analysis (follow-up group ``stab``).

Question (DEBT: thermodynamic relevance at T > 0): started IN the exact
kagome packing at its own density, does the system (a) persist anomalously
-- evidence the state is a dynamical attractor, (b) convert to the
hexagonal crystal, or (c) melt into the fluid on the same timescale as any
other overpacked start?

Runs: 3 seeds x {kagome N=2340, hexagonal N=3120 at the SAME phi=0.6783}
restarted from constructed configurations (see run_followup.py, incl. the
rectangular-superlattice landmine), plus 3 fresh monodisperse fluid
baselines.  Metrics per snapshot: MSD (from the engine, relative to the
seed), and the K2-refscore mark fraction -- the fraction of cells the
frozen pre-registered classifier assigns to the kagome class -- plus the
same for the hexagonal reference on the hex runs.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.analysis.refscore import load_frozen, refscores  # noqa: E402
from hsga.analysis.voronoi import config_cells  # noqa: E402
from hsga.engine.driver import read_frames, read_msd  # noqa: E402

STAB = REPO / "data" / "followup" / "stab"
FRAME_PICKS = [0, 2, 5, 10, 19, 29]


def metrics(pos, rad, L):
    """K2-mark fraction (kagome memory) and mean |psi6| (hexagonal order)."""
    from hsga.analysis.isoconfig import psi6_2d

    cells = config_cells(pos, rad, L)
    sk, names = refscores(cells, rad, dim=2)
    thr = load_frozen()["mark_threshold"]
    k2 = float(np.mean(sk[:, list(names).index("K2_kagome")] < thr))
    psi6 = float(np.mean(np.abs(psi6_2d(pos, L, cells))))
    return k2, psi6


def series(prefix):
    frames = read_frames(f"{prefix}.cfg", dim=2)
    t, msd = read_msd(f"{prefix}.msd")
    picks = [k for k in FRAME_PICKS if k < len(frames)]
    vals = [metrics(*frames[k]) for k in picks]
    return {"frames": picks, "k2_mark": [v[0] for v in vals],
            "psi6": [v[1] for v in vals],
            "msd_final": float(msd[-1]), "msd": [float(m) for m in msd[::5]]}


def main() -> int:
    out = {}
    for name in ("kagome", "hex"):
        rs = []
        for c in sorted(STAB.glob(f"{name}_r*.cfg")):
            try:
                rs.append(series(str(c)[:-4]))
            except Exception as e:  # noqa: BLE001 -- report, never invent
                rs.append({"error": str(e)})
        out[name] = rs
        ok = [r for r in rs if "error" not in r]
        if ok:
            k0 = np.mean([r["k2_mark"][0] for r in ok])
            k1 = np.mean([r["k2_mark"][-1] for r in ok])
            p0 = np.mean([r["psi6"][0] for r in ok])
            p1 = np.mean([r["psi6"][-1] for r in ok])
            msd = np.mean([r["msd_final"] for r in ok])
            print(f"{name:7s}: K2-mark {k0:.3f} -> {k1:.3f}   "
                  f"|psi6| {p0:.3f} -> {p1:.3f}   "
                  f"final MSD {msd:.3f} sigma^2   ({len(ok)} runs)")
    fluid = []
    for c in sorted((STAB / "m1").glob("e*_r*.cfg")):
        pos, rad, L = read_frames(c, dim=2)[-1]
        k2, psi6 = metrics(pos, rad, L)
        fluid.append({"k2_mark": k2, "psi6": psi6})
    if fluid:
        out["fluid_baseline"] = fluid
        print(f"fluid  : K2-mark {np.mean([f['k2_mark'] for f in fluid]):.4f}"
              f"   |psi6| {np.mean([f['psi6'] for f in fluid]):.3f}"
              f"   (equilibrium at the same phi)")
    (REPO / "results/campaign_2d/stability.json").write_text(
        json.dumps(out, indent=2))
    print("-> results/campaign_2d/stability.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Follow-up analyses beyond the main re-test: composition independence,
finite size, equilibration defense, and the dilution ladder.

* comp : kagome window at binary R^-1 = 1.7 (mode 2).  A geometric feature
  is pinned at sqrt(3)pi/8 for EVERY composition; the calibrated kink test
  runs on this dataset exactly as on the main one.
* fs12/fs32 : kagome window at N = 504 and N = 3520.  Kink test per size
  (a size-suppressed feature grows with N) and the eps* level itself vs N.
* eq : sixteen 3x-equilibration replicas at five critical densities; each
  mean must agree with the campaign value (two-sample z), else the
  campaign was not equilibrated and every conclusion is suspect.
* stab2 : dilution ladder -- melted-or-rigid verdict per structure and
  shrink factor (MSD threshold 0.1 sigma^2 over the run), the direct
  comparison kagome vs hexagonal vs the NON-tangential z=4 square control.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.analysis.featuretest import (  # noqa: E402
    calibrate, empirical_p, injection_power, kink_at_target)
from hsga.engine.driver import parse_log, read_msd  # noqa: E402
from rerun_feature_test import reduce_group  # noqa: E402

T_KAG = 0.6801747615878316
FOLLOW = REPO / "data" / "followup"
OUT = REPO / "results" / "campaign_2d"


def kagome_kink(gdir, label):
    eta, eps, sem, reps = reduce_group(gdir)
    r = kink_at_target(eta, eps, sem, T_KAG, degree=1, band=0.012)
    if r.get("tested"):
        null = calibrate(eta, eps, sem, T_KAG, degree=1, band=0.012)
        r["p_emp"] = empirical_p(r["z_kink"], null["z_kink_null"])
        r["a95_rel"] = injection_power(eta, eps, sem, T_KAG, n_synth=400,
                                       degree=1, band=0.012).get("a95_rel")
    print(f"  {label:14s} kink z={r.get('z_kink', float('nan')):6.2f} "
          f"p={r.get('p_emp', float('nan')):.3f} "
          f"a95={100*r.get('a95_rel', float('nan')):.1f}%  "
          f"eps*(kagome)={np.interp(T_KAG, eta, eps):.5f}  "
          f"replicas {int(reps.min())}-{int(reps.max())}")
    return {"label": label, "kink": r,
            "eps_at_target": float(np.interp(T_KAG, eta, eps)),
            "eta": eta.tolist(), "eps": eps.tolist(), "sem": sem.tolist()}


def eq_defense():
    import csv
    camp = {float(r["eta"]): (float(r["eps_star"]), float(r["sem"]))
            for r in csv.DictReader(open(OUT / "m0_eps_star.csv"))}
    eta, eps, sem, reps = reduce_group(FOLLOW / "eq")
    rows = []
    print("  equilibration defense (eq x3 vs campaign):")
    for e, y, s in zip(eta, eps, sem):
        key = min(camp, key=lambda c: abs(c - e))
        if abs(key - e) > 1e-6:
            continue
        y0, s0 = camp[key]
        z = (y - y0) / np.hypot(s, s0)
        rows.append({"eta": float(e), "eps_eq3": float(y), "eps_campaign": y0,
                     "z": float(z)})
        print(f"    eta={e:.6f}: {y:.5f} vs {y0:.5f}  z={z:+.2f}")
    worst = max(abs(r["z"]) for r in rows)
    print(f"    worst |z| = {worst:.2f}  -> "
          f"{'PASS (<3)' if worst < 3 else 'FAIL: not equilibrated'}")
    return {"rows": rows, "worst_abs_z": worst, "pass": bool(worst < 3)}


def dilution():
    rows = []
    print("  dilution ladder (rigid if final MSD < 0.1 sigma^2):")
    for cfg in sorted((FOLLOW / "stab2").glob("*_r*.msd")):
        name, shrink, rep = cfg.stem.split("_")
        t, m = read_msd(cfg)
        log = parse_log(str(cfg)[:-4] + ".log")
        rows.append({"structure": name, "shrink": float(shrink[1:]),
                     "rep": int(rep[1:]), "msd_final": float(m[-1]),
                     "phi": log.get("phi"), "rigid": bool(m[-1] < 0.1)})
    out = {}
    for name in ("kagome", "hex", "square"):
        rs = [r for r in rows if r["structure"] == name]
        by = {}
        for r in rs:
            by.setdefault(r["shrink"], []).append(r)
        line = []
        for s in sorted(by, reverse=True):
            msd = np.mean([r["msd_final"] for r in by[s]])
            rigid = all(r["rigid"] for r in by[s])
            line.append(f"s={s}: {'RIGID' if rigid else f'melts (MSD {msd:.1f})'}")
        print(f"    {name:7s} " + " | ".join(line))
        out[name] = rows and [r for r in rows if r["structure"] == name]
    return {"rows": rows}


def main() -> int:
    res = {}
    print("composition & finite-size (kagome window, calibrated kink test):")
    res["comp_R1.7"] = kagome_kink(FOLLOW / "comp", "R^-1=1.7 N=2016")
    res["fs_N504"] = kagome_kink(FOLLOW / "fs12", "R^-1=1.4 N=504")
    res["fs_N3520"] = kagome_kink(FOLLOW / "fs32", "R^-1=1.4 N=3520")
    res["eq"] = eq_defense()
    res["dilution"] = dilution()
    (OUT / "followup.json").write_text(json.dumps(res, indent=2, default=float))
    print(f"-> {OUT}/followup.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())

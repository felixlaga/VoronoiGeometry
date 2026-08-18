#!/usr/bin/env python3
"""Redesigned T6 feature test (see analysis/featuretest.py for why).

Frozen design (committed before the calibration/extension data existed):

* offset test : sideband hold-out, polynomial degree 3, band 0.045
* kink test   : linear + hinge pinned at the target, band 0.012
                (requires the near-sideband extension grid; until that data
                exists the degree-3/0.045 variant is reported as fallback)
* calibration : 2000 synthetic smooth curves per window with the measured
                sem structure; empirical p-values only
* sensitivity : injection a95 -- the smallest slope change detected with
                95% power at the calibrated 5% threshold
* multiplicity: Bonferroni over the three theory windows; the two control
                windows (0.6300, 0.7250) and the honeycomb fluid window are
                negative controls -- the procedure is invalid if it fires
                there
* replicas    : eps* per density = mean over 44 replicas; sem over replicas

Inputs: results/campaign_2d/m0_eps_star.csv (T6) plus, when present,
data/followup/{sb,ctrl}/ raw runs (reduced here with the same estimator).
Output: results/campaign_2d/retest.json and a human-readable addendum.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.analysis.featuretest import (  # noqa: E402
    ALPHA, calibrate, empirical_p, injection_power, kink_at_target,
    sideband_offset,
)
from hsga.analysis.percolation import eps_star  # noqa: E402
from hsga.engine.driver import read_frames  # noqa: E402

THEORY = {"honeycomb_K3": 0.6045997880780726,
          "kagome_K2": 0.6801747615878316,
          "maple_leaf_K1": 0.7773425846718076}
CONTROLS = {"control_0.63": 0.63, "control_0.725": 0.725}
KINK_BAND, KINK_DEG = 0.012, 1
FALL_BAND, FALL_DEG = 0.045, 3


def load_t6():
    import csv
    rows = list(csv.DictReader(open(REPO / "results/campaign_2d/m0_eps_star.csv")))
    eta = np.array([float(r["eta"]) for r in rows])
    eps = np.array([float(r["eps_star"]) for r in rows])
    sem = np.array([float(r["sem"]) for r in rows])
    return eta, eps, sem


def reduce_group(gdir: Path, nproc=None):
    """eps* mean/sem per density from a followup group's raw runs."""
    from collections import defaultdict
    from concurrent.futures import ProcessPoolExecutor

    cfgs = sorted(gdir.glob("m*/e*_r*.cfg"))
    per = defaultdict(list)
    def one(c):
        eta = float(c.name.split("_")[0][1:])
        vals = [eps_star(p, r, L) for p, r, L in read_frames(c, dim=2)]
        vals = [v for v in vals if np.isfinite(v)]
        return eta, float(np.mean(vals)) if vals else np.nan
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        for eta, v in ex.map(one, cfgs, chunksize=8):
            if np.isfinite(v):
                per[eta].append(v)
    eta = np.array(sorted(per))
    eps = np.array([np.mean(per[e]) for e in eta])
    sem = np.array([np.std(per[e], ddof=1) / np.sqrt(len(per[e])) for e in eta])
    reps = np.array([len(per[e]) for e in eta])
    return eta, eps, sem, reps


def merge(base, extra):
    """Merge (eta, eps, sem) sets; extra densities not in base are appended."""
    e0, y0, s0 = base
    e1, y1, s1 = extra
    keep = ~np.isin(np.round(e1, 6), np.round(e0, 6))
    eta = np.concatenate([e0, e1[keep]])
    order = np.argsort(eta)
    return (eta[order], np.concatenate([y0, y1[keep]])[order],
            np.concatenate([s0, s1[keep]])[order])


def run_window(name, t, eta, eps, sem, *, is_control, n_theory):
    have_sb = np.sum((np.abs(eta - t) > 0.0061) & (np.abs(eta - t) <= KINK_BAND)) >= 4
    res = {"window": name, "target": t, "is_control": is_control,
           "sb_extension_present": bool(have_sb)}
    res["offset"] = sideband_offset(eta, eps, sem, t, degree=FALL_DEG, band=FALL_BAND)
    kink_kw = (dict(degree=KINK_DEG, band=KINK_BAND) if have_sb
               else dict(degree=FALL_DEG, band=FALL_BAND))
    res["kink"] = kink_at_target(eta, eps, sem, t, **kink_kw)
    res["kink"]["design"] = "primary (deg1/0.012)" if have_sb else "fallback (deg3/0.045)"
    null = calibrate(eta, eps, sem, t, degree=kink_kw["degree"], band=kink_kw["band"])
    if res["offset"].get("tested"):
        # offset null must use the offset design
        null_off = calibrate(eta, eps, sem, t, degree=FALL_DEG, band=FALL_BAND)
        res["offset"]["p_emp"] = empirical_p(res["offset"]["z_offset"],
                                             null_off["z_offset_null"])
    if res["kink"].get("tested"):
        res["kink"]["p_emp"] = empirical_p(res["kink"]["z_kink"], null["z_kink_null"])
    alpha_w = ALPHA / n_theory
    for key in ("offset", "kink"):
        r = res[key]
        if r.get("tested") and np.isfinite(r.get("p_emp", np.nan)):
            r["significant"] = bool(r["p_emp"] < (ALPHA if is_control else alpha_w))
    res["sensitivity"] = injection_power(eta, eps, sem, t, n_synth=400,
                                         degree=kink_kw["degree"], band=kink_kw["band"])
    return res


def main(argv=None) -> int:
    followup = REPO / "data" / "followup"
    base = load_t6()
    used = ["results/campaign_2d/m0_eps_star.csv"]
    for grp in ("sb", "ctrl"):
        g = followup / grp
        if g.exists() and list(g.glob("m*/e*_r*.cfg")):
            print(f"reducing followup group '{grp}' ...")
            e, y, s, reps = reduce_group(g)
            print(f"  {len(e)} densities, replicas {reps.min()}–{reps.max()}")
            base = merge(base, (e, y, s))
            used.append(f"data/followup/{grp}")
    eta, eps, sem = base

    n_theory = len(THEORY)
    out = {"design": {"kink": f"deg{KINK_DEG}/band{KINK_BAND}",
                      "offset": f"deg{FALL_DEG}/band{FALL_BAND}",
                      "alpha_family": ALPHA, "n_theory_windows": n_theory},
           "inputs": used, "windows": []}
    print(f"\n{'window':16s} {'z_off':>7s} {'p_off':>7s} {'z_kink':>7s} "
          f"{'p_kink':>7s} {'a95(rel)':>9s} verdict")
    for name, t in list(THEORY.items()) + list(CONTROLS.items()):
        is_ctrl = name in CONTROLS or name == "honeycomb_K3"
        if name in CONTROLS and not any("ctrl" in u for u in used):
            continue
        r = run_window(name, t, eta, eps, sem, is_control=is_ctrl,
                       n_theory=n_theory)
        out["windows"].append(r)
        zo = r["offset"].get("z_offset", float("nan"))
        po = r["offset"].get("p_emp", float("nan"))
        zk = r["kink"].get("z_kink", float("nan"))
        pk = r["kink"].get("p_emp", float("nan"))
        a95 = r["sensitivity"].get("a95_rel", float("nan"))
        sig = (r["offset"].get("significant") or r["kink"].get("significant"))
        verdict = ("CONTROL-FIRED (procedure invalid!)" if sig and is_ctrl else
                   "SIGNIFICANT" if sig else "null")
        print(f"{name:16s} {zo:7.2f} {po:7.3f} {zk:7.2f} {pk:7.3f} "
              f"{a95*100:8.1f}% {verdict}"
              + ("   [control]" if is_ctrl else ""))
    ctrl_fired = any(w["is_control"] and (w["offset"].get("significant")
                                          or w["kink"].get("significant"))
                     for w in out["windows"])
    out["controls_clean"] = not ctrl_fired
    (REPO / "results/campaign_2d/retest.json").write_text(
        json.dumps(out, indent=2, default=float))
    print(f"\ncontrols clean: {not ctrl_fired}   "
          f"-> results/campaign_2d/retest.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())

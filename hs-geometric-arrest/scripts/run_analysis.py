#!/usr/bin/env python3
"""Campaign analysis: does anything happen at the exact targets?  (T6 / T9.)

For each mode in a campaign data directory this computes, per density:

* ``eps*`` with errors over independent replicas only;
* the smooth null: the sem-weighted fit ``eps* = C (eta_m - eta)^p`` over the
  whole sweep -- the positive control showed a naive local exponent "jump" is
  entirely explained by this law, so a feature claim is ONLY an excursion
  above it.  Standardised residuals, plus de Graaf's own diagnostic
  ``d log(eta_m - eta)/d log eps`` (flat for a pure power law);
* the feature test at each exact target: the maximum standardised residual
  inside the fine window against a look-elsewhere-corrected threshold, with a
  bootstrap over replicas for significance and location, and a minimum-replica
  guard (a bootstrap over a handful of replicas underestimates its own spread
  and manufactures z-scores);
* ``P_wrap^(k)``, susceptibility and correlation length per reference
  (subsampled configurations -- cell construction is the expensive step);
* the large-particle q-distribution submode weight near the target q;
* the MSD arrest density as a consistency check only, window-scanned and
  labelled MC.

Verdict logic (T6): a resolvable, look-elsewhere-surviving feature at
0.777343 validates the pipeline; anything at 0.680175 (kagome) is a new
result either way; nothing anywhere means STOP -- the 3D campaign is not run
and the negative result is written up as such.  "Could not be tested" is
reported as absent, never as negative.

This script has NOT been executed on campaign data (none exists); it is
exercised end to end by the smoke presets.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.analysis.dynamics import (  # noqa: E402
    diffusion_coefficient, fit_eta_a, threshold_sweep,
)
from hsga.analysis.percolation import eps_star  # noqa: E402
from hsga.analysis.pwrap import p_wrap_all_references  # noqa: E402
from hsga.analysis.voronoi import config_cells, q_submode_weight  # noqa: E402
from hsga.engine.driver import read_frames, read_msd  # noqa: E402

TARGETS = {
    2: {"honeycomb_K3": 0.6045997880780726,
        "kagome_K2": 0.6801747615878316,
        "maple_leaf_K1": 0.7773425846718076},
    3: {"simple_cubic": 0.5235987755982988, "K4": 0.5553603672697957,
        "K3": 0.5923843917544487, "simple_hexagonal": 0.6045997880780726,
        "K2": 0.6346975625940521, "K1": 0.6835204520243637},
}
FINE_HALF = 0.006


# --------------------------------------------------------------------------- #
# per-run reduction (worker processes)
# --------------------------------------------------------------------------- #
def reduce_run(task):
    prefix, dim, structural, n_struct_frames = task
    frames = read_frames(f"{prefix}.cfg", dim=dim)
    eps = [eps_star(pos, rad, L) for pos, rad, L in frames]
    out = {"prefix": prefix, "eps": [e for e in eps if np.isfinite(e)],
           "n_frames": len(frames),
           "n_eps_failed": int(sum(not np.isfinite(e) for e in eps))}
    if structural and frames:
        pos, rad, L = frames[-1]
        cells = config_cells(pos, rad, L)
        pw = p_wrap_all_references(cells, pos, rad, L, dim=dim)
        out["pwrap"] = {k: {"wraps": v["wraps"], "marked": v["marked_fraction"],
                            "chi": v["chi"], "xi": v["xi"]}
                        for k, v in pw.items()}
        Q = np.array([c.Q_iso for c in cells])
        out["submode"] = {
            name: q_submode_weight(Q, rad, t)
            for name, t in TARGETS[dim].items()
        }
    return out


# --------------------------------------------------------------------------- #
# the smooth null and the feature test
# --------------------------------------------------------------------------- #
def smooth_null(eta, eps, sem):
    """Sem-weighted ``eps* = C (eta_m - eta)^p`` and standardised residuals."""
    from scipy.optimize import curve_fit

    eta = np.asarray(eta, float)
    eps = np.asarray(eps, float)
    sem = np.asarray(sem, float)
    ok = np.isfinite(eps) & np.isfinite(sem) & (sem > 0)
    f = lambda x, C, em, p: C * np.clip(em - x, 1e-9, None) ** p
    popt, pcov = curve_fit(f, eta[ok], eps[ok], p0=[0.1, eta[ok].max() + 0.05, 1.0],
                           sigma=sem[ok], absolute_sigma=True, maxfev=60000)
    resid = np.full(len(eta), np.nan)
    resid[ok] = (eps[ok] - f(eta[ok], *popt)) / sem[ok]
    # de Graaf's diagnostic: d log(eta_m - eta) / d log eps == 1/p, flat
    lx = np.log(np.clip(popt[1] - eta[ok], 1e-12, None))
    ly = np.log(eps[ok])
    diag = np.gradient(lx, ly)
    return {
        "C": float(popt[0]), "eta_m": float(popt[1]), "p": float(popt[2]),
        "eta_m_err": float(np.sqrt(pcov[1, 1])),
        "residuals": resid, "diagnostic_flat": diag,
        "note": "eta_m is a nuisance parameter of the null; the form is known "
                "to be wrong near contact percolation (DEBT.md) and eta_m is "
                "NEVER quoted as a physical vanishing point",
    }


def feature_test(eta, per_replica, targets, *, n_boot=400, min_replicas=8,
                 alpha=0.05, rng_seed=20260816):
    """Excess above the smooth null inside each fine target window.

    ``per_replica[j]`` is a replica's eps* curve over ``eta``.  The statistic
    per target is the maximum standardised residual inside the +-FINE_HALF
    window; its null distribution comes from bootstrapping replicas; the
    threshold is Bonferroni-corrected for the number of windows scanned.
    Below ``min_replicas`` no significance is claimed at all.
    """
    from scipy.stats import norm

    eta = np.asarray(eta, float)
    reps = np.asarray(per_replica, float)
    n_rep = len(reps)
    if n_rep < min_replicas:
        return {name: {"tested": False,
                       "note": f"only {n_rep} replicas; >= {min_replicas} needed"}
                for name in targets}

    def curve_resid(sample):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mean = np.nanmean(sample, axis=0)
            sem = np.nanstd(sample, axis=0, ddof=1) / np.sqrt(
                np.sum(np.isfinite(sample), axis=0).clip(1))
        try:
            return smooth_null(eta, mean, sem)["residuals"]
        except Exception:
            return None

    obs = curve_resid(reps)
    if obs is None:
        return {name: {"tested": False, "note": "null fit failed"} for name in targets}

    rng = np.random.default_rng(rng_seed)
    boot = []
    for _ in range(n_boot):
        pick = rng.integers(0, n_rep, n_rep)
        r = curve_resid(reps[pick])
        if r is not None:
            boot.append(r)
    boot = np.array(boot)
    n_windows = len(targets)
    z_thresh = float(norm.isf(0.5 * alpha / max(1, n_windows)))

    out = {}
    for name, t in targets.items():
        win = np.abs(eta - t) <= FINE_HALF
        if not win.any() or not np.isfinite(obs[win]).any():
            out[name] = {"tested": False, "note": "no data in the fine window"}
            continue
        k = np.nanargmax(np.where(win, obs, -np.inf))
        z_obs = float(obs[k])
        spread = float(np.nanstd(boot[:, k], ddof=1)) if len(boot) else float("nan")
        z_eff = z_obs / spread if spread > 0 else float("nan")
        out[name] = {
            "tested": True,
            "target": float(t),
            "eta_at_max": float(eta[k]),
            "z_raw": z_obs,
            "bootstrap_spread": spread,
            "z_effective": float(z_eff),
            "z_threshold": z_thresh,
            "resolvable": bool(np.isfinite(z_eff) and z_eff > z_thresh),
        }
    return out


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def analyse_mode(data_dir: Path, mode_dir: Path, dim: int, a):
    prefixes = sorted({str(p.with_suffix("")) for p in mode_dir.glob("*.cfg")})
    if not prefixes:
        return None
    tasks = []
    for i, pre in enumerate(prefixes):
        structural = a.structural and (i % max(1, a.structural_every) == 0)
        tasks.append((pre, dim, structural, 1))
    if a.nproc == 1:
        results = [reduce_run(t) for t in tasks]
    else:
        with ProcessPoolExecutor(max_workers=a.nproc) as pool:
            results = list(pool.map(reduce_run, tasks))

    by_eta: dict[float, list] = {}
    for r in results:
        eta = float(Path(r["prefix"]).name.split("_")[0][1:])
        by_eta.setdefault(eta, []).append(r)
    etas = np.array(sorted(by_eta))
    n_rep = max(len(v) for v in by_eta.values())
    per_replica = np.full((n_rep, len(etas)), np.nan)
    rows = []
    for i, e in enumerate(etas):
        means = [float(np.mean(r["eps"])) for r in by_eta[e] if r["eps"]]
        per_replica[: len(means), i] = means
        rows.append({
            "eta": float(e),
            "eps_star": float(np.mean(means)) if means else float("nan"),
            "sem": float(np.std(means, ddof=1) / np.sqrt(len(means)))
                   if len(means) > 1 else float("nan"),
            "n_replicas": len(means),
            "n_configurations": int(sum(len(r["eps"]) for r in by_eta[e])),
        })

    # attempted-but-unreachable densities: log exists, no cfg
    attempted = {float(p.name.split("_")[0][1:]) for p in mode_dir.glob("*.log")}
    missing = sorted(attempted - set(np.round(etas, 6)))

    null = None
    sems = np.array([r["sem"] for r in rows])
    means = np.array([r["eps_star"] for r in rows])
    if np.isfinite(sems).sum() >= 5:
        try:
            null = smooth_null(etas, means, sems)
        except Exception as exc:
            null = {"error": str(exc)}

    features = feature_test(etas, per_replica, TARGETS[dim],
                            n_boot=a.n_boot, min_replicas=a.min_replicas)

    # structural summaries (subsampled)
    pwrap_rows, submode_rows = [], []
    for e in etas:
        recs = [r for r in by_eta[e] if "pwrap" in r]
        if not recs:
            continue
        for k in recs[0]["pwrap"]:
            pwrap_rows.append({
                "eta": float(e), "reference": k,
                "P_wrap": float(np.mean([r["pwrap"][k]["wraps"] for r in recs])),
                "marked": float(np.mean([r["pwrap"][k]["marked"] for r in recs])),
                "chi": float(np.mean([r["pwrap"][k]["chi"] for r in recs])),
                "xi": float(np.mean([r["pwrap"][k]["xi"] for r in recs])),
                "n": len(recs),
            })
        for name in recs[0]["submode"]:
            fr = [r["submode"][name]["flank_ratio"] for r in recs]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                fr_mean = float(np.nanmean(fr)) if np.isfinite(fr).any() else float("nan")
            submode_rows.append({
                "eta": float(e), "target": name,
                "flank_ratio": fr_mean,
                "weight": float(np.mean([r["submode"][name]["weight"] for r in recs])),
                "n": len(recs),
            })

    # MSD consistency check (never primary)
    msd_fit = None
    Ds = {}
    for f in sorted(mode_dir.glob("*.msd")):
        eta = float(f.name.split("_")[0][1:])
        t, m = read_msd(f)
        d = diffusion_coefficient(t, m, dim=dim)
        if np.isfinite(d["D"]) and d["D"] > 0 and d["diffusive"]:
            Ds.setdefault(eta, []).append(d["D"])
    if len(Ds) >= 4:
        e_arr = np.array(sorted(Ds))
        msd_fit = fit_eta_a(e_arr, np.array([np.mean(Ds[e]) for e in e_arr]),
                            dim=dim, dynamics_label="MC (Brownian proxy)")

    # observation-time-threshold sweep (Kramers-vs-pinned discriminator):
    # does the apparent arrest density move with the tau_alpha criterion?
    # Error bars from a bootstrap over the replica axis -- the only
    # independent axis -- never over snapshots.
    tsweep = None
    if len(Ds) >= 4:
        e_arr = np.array(sorted(Ds))
        tsweep = threshold_sweep(e_arr,
                                 np.array([np.mean(Ds[e]) for e in e_arr]))
        rng = np.random.default_rng(20260818)
        drifts, cross_boot = [], {}
        for _ in range(400):
            Db = [np.mean(rng.choice(Ds[e], size=len(Ds[e]), replace=True))
                  for e in e_arr]
            rb = threshold_sweep(e_arr, np.array(Db))
            if np.isfinite(rb["drift_per_decade"]):
                drifts.append(rb["drift_per_decade"])
            for c in rb["crossings"]:
                cross_boot.setdefault(c["decade"], []).append(c["eta_x"])
        tsweep["drift_boot_sem"] = (float(np.std(drifts)) if drifts
                                    else float("nan"))
        for c in tsweep["crossings"]:
            bs = cross_boot.get(c["decade"], [])
            c["eta_x_sem"] = float(np.std(bs)) if len(bs) > 1 else float("nan")
        tsweep["dynamics_label"] = "MC (Brownian proxy)"

    return {"rows": rows, "missing": missing, "null": null, "features": features,
            "pwrap": pwrap_rows, "submode": submode_rows, "msd_fit": msd_fit,
            "tsweep": tsweep, "n_replicas": int(n_rep)}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Campaign analysis and verdict")
    p.add_argument("--data-dir", required=True)
    p.add_argument("--dim", type=int, choices=(2, 3), required=True)
    p.add_argument("--results", default=None)
    p.add_argument("--structural", action="store_true",
                   help="also compute P_wrap and submode (expensive cell builds)")
    p.add_argument("--structural-every", type=int, default=5,
                   help="structural analysis on every k-th run only")
    p.add_argument("--n-boot", type=int, default=400)
    p.add_argument("--min-replicas", type=int, default=8)
    p.add_argument("--nproc", type=int, default=None)
    a = p.parse_args(argv)

    data_dir = Path(a.data_dir)
    if a.results is None:
        a.results = str(REPO / "results" / f"campaign_{data_dir.name}")
    out_dir = Path(a.results)
    out_dir.mkdir(parents=True, exist_ok=True)

    verdicts, all_reports = [], {}
    for mode_dir in sorted(data_dir.glob("m*")):
        rep = analyse_mode(data_dir, mode_dir, a.dim, a)
        if rep is None:
            continue
        all_reports[mode_dir.name] = rep
        print(f"\n== {mode_dir.name}: {len(rep['rows'])} densities, "
              f"{rep['n_replicas']} replicas ==")
        if rep["missing"]:
            print(f"  unreachable densities (no data, not interpolated): "
                  f"{[f'{m:.4f}' for m in rep['missing']]}")
        if rep["null"] and "eta_m" in rep["null"]:
            print(f"  smooth null: eta_m={rep['null']['eta_m']:.4f} "
                  f"p={rep['null']['p']:.3f} (nuisance parameters, not physics)")
        for name, f in rep["features"].items():
            if not f.get("tested"):
                print(f"  {name}: NOT TESTED ({f.get('note')}) — absent, not negative")
            else:
                print(f"  {name} @ {f['target']:.4f}: z_eff={f['z_effective']:.2f} "
                      f"(threshold {f['z_threshold']:.2f}) "
                      f"-> {'RESOLVABLE' if f['resolvable'] else 'no feature'}")
                verdicts.append((mode_dir.name, name, f))
        if rep["msd_fit"] is not None and np.isfinite(rep["msd_fit"]["eta_a"]):
            m = rep["msd_fit"]
            print(f"  MSD consistency: eta_a={m['eta_a']:.4f} +- {m['err']:.4f} "
                  f"[{m['note'].split(';')[0]}]")
        ts = rep.get("tsweep")
        if ts and ts["crossings"]:
            xs = ", ".join(f"10^-{c['decade']}: {c['eta_x']:.4f}"
                           + (f"+-{c['eta_x_sem']:.4f}"
                              if np.isfinite(c.get("eta_x_sem", float("nan")))
                              else "")
                           for c in ts["crossings"])
            print(f"  threshold sweep (D/D0 criterion): {xs}")
            if np.isfinite(ts["drift_per_decade"]):
                print(f"    arrest-location drift {ts['drift_per_decade']:+.4f} "
                      f"+- {ts['drift_boot_sem']:.4f} per decade "
                      f"[{ts['dynamics_label']}] — a smooth Kramers barrier "
                      f"predicts steady drift; a pinned geometric feature "
                      f"predicts the structural observables stay put")

    if not all_reports:
        print("no campaign data found; run scripts/run_sweep.py first")
        return 1

    # ---- decision text (T6 semantics) ----------------------------------- #
    lines = [f"# Campaign analysis — {data_dir.name} (dim={a.dim})", ""]
    tested = [v for v in verdicts]
    resolvable = [v for v in tested if v[2].get("resolvable")]
    if a.dim == 2:
        at_maple = [v for v in resolvable if "maple" in v[1]]
        at_kagome = [v for v in resolvable if "kagome" in v[1]]
        if not tested:
            decision = ("INCONCLUSIVE: no target window could be tested. This is an "
                        "absent result, not a negative one; the T6 decision is NOT "
                        "reached and the 3D campaign stays blocked.")
        elif at_maple:
            decision = ("FEATURE AT 0.777343: the pipeline resolves de Graaf's "
                        "transition. T6 decision: proceed to the 3D campaign (T9). "
                        + ("Kagome rung ALSO shows a feature — a new result."
                           if at_kagome else "No kagome feature."))
        elif at_kagome:
            decision = ("No feature at 0.777343 but a feature at the kagome rung "
                        "0.680175: a new result either way; the maple-leaf null "
                        "must be reported alongside it. 3D remains blocked unless "
                        "the maple-leaf control is understood.")
        else:
            decision = ("NO FEATURE at any target at these statistics. T6 decision: "
                        "STOP — do not run the 3D campaign; write the negative "
                        "result (it is publishable and must not be dressed up).")
    else:
        decision = (f"{len(resolvable)} of {len(tested)} tested target windows show "
                    "a resolvable feature. The admissible 3D ladder is dense: any "
                    "match must be reported with the number of ladder values inside "
                    "the uncertainty, never as bare agreement.")
    print("\nDECISION: " + decision)
    lines += ["## Decision", "", decision, ""]

    # ---- persist -------------------------------------------------------- #
    for mode, rep in all_reports.items():
        lines += [f"## {mode}", "",
                  "| eta | eps* | sem | replicas | configs |", "|---|---|---|---|---|"]
        for r in rep["rows"]:
            lines.append(f"| {r['eta']:.6f} | {r['eps_star']:.6f} | {r['sem']:.6f} "
                         f"| {r['n_replicas']} | {r['n_configurations']} |")
        if rep["missing"]:
            lines += ["", "Unreachable (no data, never interpolated): "
                      + ", ".join(f"{m:.4f}" for m in rep["missing"])]
        lines += [""]
        with open(out_dir / f"{mode}_eps_star.csv", "w") as fh:
            fh.write("eta,eps_star,sem,n_replicas,n_configurations\n")
            for r in rep["rows"]:
                fh.write(f"{r['eta']},{r['eps_star']},{r['sem']},"
                         f"{r['n_replicas']},{r['n_configurations']}\n")
        if rep["pwrap"]:
            with open(out_dir / f"{mode}_pwrap.csv", "w") as fh:
                fh.write("eta,reference,P_wrap,marked,chi,xi,n\n")
                for r in rep["pwrap"]:
                    fh.write(f"{r['eta']},{r['reference']},{r['P_wrap']},"
                             f"{r['marked']},{r['chi']},{r['xi']},{r['n']}\n")
        if rep["submode"]:
            with open(out_dir / f"{mode}_submode.csv", "w") as fh:
                fh.write("eta,target,flank_ratio,weight,n\n")
                for r in rep["submode"]:
                    fh.write(f"{r['eta']},{r['target']},{r['flank_ratio']},"
                             f"{r['weight']},{r['n']}\n")
    (out_dir / "features.json").write_text(json.dumps(
        {m: {"features": rep["features"],
             "null": {k: v for k, v in (rep["null"] or {}).items()
                      if not isinstance(v, np.ndarray)},
             "msd_fit": ({k: v for k, v in rep["msd_fit"].items() if k != "windows"}
                         if rep["msd_fit"] else None),
             "tsweep": rep.get("tsweep")}
         for m, rep in all_reports.items()}, indent=2, default=str))
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {out_dir}/report.md, features.json and per-mode CSVs")
    return 0


if __name__ == "__main__":
    sys.exit(main())

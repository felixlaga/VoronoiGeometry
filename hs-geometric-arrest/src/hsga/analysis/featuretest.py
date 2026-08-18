"""Redesigned feature test for fine-window campaigns (T6 re-analysis).

Why this module exists: the pre-registered test compared each fine window
against a single GLOBAL three-parameter power law.  At campaign precision
(sem ~ 5e-5) that form is falsified by its own fit (chi2/dof = 137) and the
max-residual statistic tracks the misfit -- demonstrated by an 8.4-sigma
"feature" in the normal fluid (results/campaign_2d/DIAGNOSTICS.md).  The
tests here use NO global functional form:

* ``sideband_offset``   -- the null is measured from the data in the
  sidebands of the window alone (sem-weighted local polynomial); the
  window's precision-weighted mean deviation from the sideband prediction
  is standardised with the FULL covariance (measurement variance plus
  propagated prediction covariance, which is correlated across the window).
* ``kink_at_target``    -- one extra degree of freedom over a local smooth
  polynomial: a slope discontinuity at EXACTLY the predicted density.  The
  theories under test predict where, so the test is pinned there; a free-
  breakpoint search would pay an unnecessary look-elsewhere penalty.
* ``calibrate``         -- empirical null distributions for both statistics
  from synthetic smooth curves carrying the campaign's own noise structure;
  p-values are read from these, never from a Gaussian table.
* ``injection_power``   -- the minimum kink amplitude detectable at a given
  power, so a negative result states what it excludes.

Frozen analysis constants (chosen before the calibration data existed;
changing them after seeing calibration or re-test output would void the
design): HALF = 0.006 (window half-width, matches the campaign grid),
BAND = 0.045 (sideband reach), DEGREE = 3 (local polynomial), GEN_DEGREE
= 4 (synthetic-null generator), N_SYNTH = 2000, ALPHA = 0.05 family-wise
over the tested windows (Bonferroni).
"""

from __future__ import annotations

import numpy as np

HALF = 0.006
BAND = 0.045
DEGREE = 3
GEN_DEGREE = 4
N_SYNTH = 2000
ALPHA = 0.05


def _wls(X, y, w):
    """Weighted least squares with absolute sigmas: returns (beta, cov)."""
    Xw = X * np.sqrt(w)[:, None]
    yw = y * np.sqrt(w)
    XtX = Xw.T @ Xw
    cov = np.linalg.inv(XtX)
    beta = cov @ (Xw.T @ yw)
    return beta, cov


def _poly_design(x, degree):
    return np.vander(x, degree + 1, increasing=True)


def sideband_offset(eta, eps, sem, target, *, half=HALF, band=BAND,
                    degree=DEGREE) -> dict:
    """Standardised window offset against a sideband-only local polynomial."""
    eta = np.asarray(eta, float); eps = np.asarray(eps, float)
    sem = np.asarray(sem, float)
    d = eta - target
    win = np.abs(d) <= half + 1e-12
    side = (np.abs(d) > half + 1e-12) & (np.abs(d) <= band + 1e-12)
    n_side = int(side.sum())
    if n_side < degree + 2 or not win.any():
        return {"tested": False,
                "note": f"{n_side} sideband points for degree {degree}"}
    w = 1.0 / sem[side] ** 2
    Xs = _poly_design(d[side], degree)
    beta, covb = _wls(Xs, eps[side], w)
    chi2_side = float(np.sum(w * (eps[side] - Xs @ beta) ** 2))
    Xw = _poly_design(d[win], degree)
    pred = Xw @ beta
    dev = eps[win] - pred
    C = np.diag(sem[win] ** 2) + Xw @ covb @ Xw.T
    Ci = np.linalg.inv(C)
    one = np.ones(len(dev))
    denom = float(one @ Ci @ one)
    offset = float(one @ Ci @ dev) / denom
    z = offset * np.sqrt(denom)
    return {"tested": True, "target": float(target), "offset": offset,
            "offset_err": float(1.0 / np.sqrt(denom)), "z_offset": float(z),
            "n_side": n_side, "n_window": int(win.sum()),
            "chi2_side_dof": chi2_side / max(1, n_side - degree - 1)}


def kink_at_target(eta, eps, sem, target, *, band=BAND, degree=DEGREE) -> dict:
    """Slope discontinuity at exactly ``target`` over a local polynomial."""
    eta = np.asarray(eta, float); eps = np.asarray(eps, float)
    sem = np.asarray(sem, float)
    d = eta - target
    sel = np.abs(d) <= band + 1e-12
    x = d[sel]
    if (x < 0).sum() < 3 or (x > 0).sum() < 3 or len(x) < degree + 3:
        return {"tested": False, "note": "too few points on a side"}
    X = np.column_stack([_poly_design(x, degree), np.clip(x, 0.0, None)])
    w = 1.0 / sem[sel] ** 2
    beta, cov = _wls(X, eps[sel], w)
    resid = eps[sel] - X @ beta
    chi2 = float(np.sum(w * resid ** 2))
    dof = len(x) - X.shape[1]
    c, sc = float(beta[-1]), float(np.sqrt(cov[-1, -1]))
    slope_at = float(beta[1])            # d eps / d eta just below target
    return {"tested": True, "target": float(target), "kink": c,
            "kink_err": sc, "z_kink": c / sc, "chi2_dof": chi2 / max(1, dof),
            "slope_below": slope_at,
            "kink_rel": c / abs(slope_at) if slope_at else float("nan"),
            "n_points": int(len(x))}


def _generator(eta, eps, sem, target, *, band=BAND, gen_degree=GEN_DEGREE):
    """Smooth local generator for synthetic nulls.  Fitted through ALL local
    data (window included), so any true localised feature is partially
    absorbed into the generator -- which makes the calibration CONSERVATIVE
    for claiming positives, and is documented as such."""
    d = np.asarray(eta, float) - target
    sel = np.abs(d) <= band + 1e-12
    w = 1.0 / np.asarray(sem, float)[sel] ** 2
    X = _poly_design(d[sel], gen_degree)
    beta, _ = _wls(X, np.asarray(eps, float)[sel], w)
    return d[sel] + target, X @ beta, np.asarray(sem, float)[sel]


def calibrate(eta, eps, sem, target, *, n_synth=N_SYNTH, seed=20260818,
              injection=0.0, **kw) -> dict:
    """Empirical null (or signal, if ``injection``) distributions of both
    statistics on synthetic curves: generator + N(0, sem) per point."""
    ge, gy, gs = _generator(eta, eps, sem, target, band=kw.get("band", BAND))
    rng = np.random.default_rng(seed)
    y0 = gy + injection * np.clip(ge - target, 0.0, None)
    zo, zk = [], []
    for _ in range(n_synth):
        y = y0 + rng.normal(0.0, gs)
        r1 = sideband_offset(ge, y, gs, target, **kw)
        r2 = kink_at_target(ge, y, gs, target, **kw)
        if r1.get("tested"):
            zo.append(r1["z_offset"])
        if r2.get("tested"):
            zk.append(r2["z_kink"])
    return {"z_offset_null": np.array(zo), "z_kink_null": np.array(zk)}


def empirical_p(z_obs: float, z_null: np.ndarray) -> float:
    """Two-sided empirical p-value, CENTERED on the null median, +1 rule.

    Centering matters: smooth local curvature leaks deterministically into
    the hinge coefficient, so the calibrated null is shifted away from
    zero.  The question "is the observed kink unusual UNDER THE SMOOTH
    NULL" is answered relative to that shifted distribution; comparing
    |z| to |null| uncentered would both lose power and mistake curvature
    for signal."""
    if not np.isfinite(z_obs) or len(z_null) == 0:
        return float("nan")
    med = float(np.median(z_null))
    return float((np.sum(np.abs(z_null - med) >= abs(z_obs - med)) + 1)
                 / (len(z_null) + 1))


def injection_power(eta, eps, sem, target, *, power=0.95, z_crit=None,
                    n_synth=400, seed=1, **kw) -> dict:
    """Minimum kink amplitude (slope-change units of eps per unit eta)
    detected with ``power`` at threshold ``z_crit`` (from the calibrated
    null).  Bisection over amplitude; returns NaN if even the largest
    scanned amplitude is not detectable (never fabricates sensitivity)."""
    null = calibrate(eta, eps, sem, target, n_synth=n_synth, seed=seed, **kw)
    med = float(np.median(null["z_kink_null"]))
    if z_crit is None:
        z_crit = float(np.quantile(np.abs(null["z_kink_null"] - med), 1 - ALPHA))

    def detect_rate(a):
        r = calibrate(eta, eps, sem, target, n_synth=n_synth, seed=seed + 7,
                      injection=a, **kw)
        z = r["z_kink_null"]
        return float(np.mean(np.abs(z - med) >= z_crit)) if len(z) else 0.0

    slope = abs(kink_at_target(eta, eps, sem, target, **kw).get("slope_below",
                                                                np.nan))
    lo, hi = 0.0, 2.0 * slope if np.isfinite(slope) and slope > 0 else 1.0
    if detect_rate(hi) < power:
        return {"a95": float("nan"), "a95_rel": float("nan"),
                "z_crit": z_crit, "note": "not reachable within scan range"}
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        if detect_rate(mid) >= power:
            hi = mid
        else:
            lo = mid
    return {"a95": hi, "a95_rel": hi / slope if slope else float("nan"),
            "z_crit": z_crit}

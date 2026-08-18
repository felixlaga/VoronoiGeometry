"""MSD, diffusion coefficient, and the arrest density -- a consistency check.

The arrest density from ``D = A (eta_a - eta)^b`` is NOT the primary estimator
of this project (spec Sec. 4), and this module is written so it cannot be
quoted as if it were: :func:`fit_eta_a` always performs the fit-window scan
(landmine 7 -- the pilot's window systematic was 0.285 against a statistical
error of 0.012) and the reported uncertainty is the larger of the two.

The module is engine-agnostic -- it fits whatever ``(t, MSD)`` series it is
given -- but the label travels with the number: MSD from ``hsmc``/``hsmc2d``
is Monte Carlo (a Brownian proxy) and its ``eta_a`` must say so; MSD from
``edmd`` is Newtonian.  The two are not interchangeable.
"""

from __future__ import annotations

import numpy as np

__all__ = ["ETA_CLOSE_PACKED_3D", "PHI_CLOSE_PACKED_2D",
           "diffusion_coefficient", "fit_eta_a", "msd_exponent"]

ETA_CLOSE_PACKED_3D = 0.7404804896930611
PHI_CLOSE_PACKED_2D = 0.9068996821171089


def _loglog_slope(t: np.ndarray, m: np.ndarray) -> float:
    ok = (t > 0) & (m > 0)
    if np.count_nonzero(ok) < 3:
        return float("nan")
    return float(np.polyfit(np.log(t[ok]), np.log(m[ok]), 1)[0])


def msd_exponent(t, msd, *, tail: float = 0.5) -> float:
    """Log-log slope ``alpha`` over the last ``tail`` fraction in time.

    ``alpha = 1`` is diffusive; ``alpha -> 0`` is a plateau (caging).
    """
    t = np.asarray(t, float)
    m = np.asarray(msd, float)
    sel = t >= t.max() * (1.0 - tail)
    return _loglog_slope(t[sel], m[sel])


def diffusion_coefficient(t, msd, *, tail: float = 0.5, dim: int = 3) -> dict:
    """``D`` from the long-time MSD slope, ``MSD = 2 d D t``.

    The returned ``alpha`` says whether the fitted window is actually
    diffusive; a ``D`` read off a sub-diffusive window is an upper bound, not
    a measurement, and callers must honour the flag.
    """
    t = np.asarray(t, float)
    m = np.asarray(msd, float)
    sel = (t >= t.max() * (1.0 - tail)) & (t > 0)
    if np.count_nonzero(sel) < 3:
        return {"D": float("nan"), "alpha": float("nan"), "n_points": int(sel.sum())}
    slope, _ = np.polyfit(t[sel], m[sel], 1)
    alpha = _loglog_slope(t[sel], m[sel])
    return {
        "D": float(slope / (2.0 * dim)),
        "alpha": float(alpha),
        "diffusive": bool(alpha > 0.9),
        "n_points": int(sel.sum()),
        "t_min": float(t[sel].min()),
        "t_max": float(t[sel].max()),
    }


def _fit_power(eta: np.ndarray, D: np.ndarray, eta_a0: float):
    """Least squares for ``ln D = ln A + b ln(eta_a - eta)``."""
    from scipy.optimize import least_squares

    lnD = np.log(D)

    def resid(p):
        eta_a, lnA, b = p
        arg = eta_a - eta
        if np.any(arg <= 0):
            return np.full_like(lnD, 1e3)
        return lnA + b * np.log(arg) - lnD

    out = least_squares(resid, [eta_a0, 0.0, 3.5], method="lm", max_nfev=20000)
    dof = max(len(lnD) - 3, 1)
    s2 = float(out.cost * 2.0 / dof)
    try:
        cov = np.linalg.inv(out.jac.T @ out.jac) * s2
        err = float(np.sqrt(max(cov[0, 0], 0.0)))
    except np.linalg.LinAlgError:
        err = float("nan")
    return float(out.x[0]), float(out.x[2]), err, bool(out.success)


def fit_eta_a(
    eta,
    D,
    *,
    min_points: int = 4,
    max_eta_a: float | None = None,
    dim: int = 3,
    dynamics_label: str = "MC (Brownian proxy)",
) -> dict:
    """Arrest density from ``D = A (eta_a - eta)^b`` WITH the mandatory window scan.

    Every contiguous sub-window with at least ``min_points`` points is fitted.
    The returned ``eta_a`` is the widest converged window; ``err`` -- the only
    number that may be quoted -- is ``max(statistical, window spread)``.

    A window counts as converged only when its ``eta_a`` is physically
    possible: above close packing is a failed extrapolation, not a loose
    determination, and such windows are counted in ``n_rejected_unphysical``
    rather than allowed to inflate the systematic.

    Returns ``NaN`` throughout when nothing converges.  No fallback value.
    """
    if max_eta_a is None:
        max_eta_a = ETA_CLOSE_PACKED_3D if dim == 3 else PHI_CLOSE_PACKED_2D
    eta = np.asarray(eta, float)
    D = np.asarray(D, float)
    ok = np.isfinite(D) & (D > 0)
    eta, D = eta[ok], D[ok]
    order = np.argsort(eta)
    eta, D = eta[order], D[order]

    if len(eta) < min_points:
        return {
            "eta_a": float("nan"), "b": float("nan"), "err_stat": float("nan"),
            "err_window": float("nan"), "err": float("nan"),
            "windows": [], "n_windows": 0,
            "note": f"only {len(eta)} usable points, need {min_points}",
        }

    windows, n_rejected = [], 0
    for i in range(len(eta) - min_points + 1):
        for j in range(i + min_points, len(eta) + 1):
            e, d = eta[i:j], D[i:j]
            a, b, err, success = _fit_power(e, d, float(e.max()) + 0.02)
            if not (success and np.isfinite(a) and a > e.max()):
                continue
            if a > max_eta_a:
                n_rejected += 1
                continue
            windows.append(
                {"lo": float(e.min()), "hi": float(e.max()), "n": int(j - i),
                 "eta_a": a, "b": b, "err_stat": err}
            )
    if not windows:
        return {
            "eta_a": float("nan"), "b": float("nan"), "err_stat": float("nan"),
            "err_window": float("nan"), "err": float("nan"),
            "windows": [], "n_windows": 0, "n_rejected_unphysical": n_rejected,
            "note": f"no window converged to a physical eta_a <= {max_eta_a:.4f} "
                    f"({n_rejected} rejected as unphysical)",
        }

    full = [w for w in windows if w["n"] == len(eta)]
    ref = full[0] if full else max(windows, key=lambda w: w["n"])
    vals = np.array([w["eta_a"] for w in windows])
    err_window = float(0.5 * (vals.max() - vals.min()))
    err_stat = float(ref["err_stat"])
    return {
        "eta_a": float(ref["eta_a"]),
        "b": float(ref["b"]),
        "err_stat": err_stat,
        "err_window": err_window,
        "err": float(max(err_stat, err_window)) if np.isfinite(err_stat) else err_window,
        "window_min": float(vals.min()),
        "window_max": float(vals.max()),
        "n_windows": len(windows),
        "n_rejected_unphysical": n_rejected,
        "windows": windows,
        "note": f"{dynamics_label} arrest density; the engine label matters — "
                "MC and EDMD eta_a are not interchangeable",
    }


def threshold_sweep(eta, D, *, max_decades: int = 8) -> dict:
    """Apparent arrest density as a function of the observation-time criterion.

    Discriminator between a Kramers-type smooth barrier (Babu, arXiv:2607.19185:
    arrest is ``tau_alpha/tau_0`` crossing a conventional threshold, so the
    arrest density shifts systematically with the chosen threshold) and a
    geometric ground state (a structural feature pinned at a fixed density,
    indifferent to the criterion).  With ``D ~ 1/tau_alpha``, the criterion
    ladder is ``D_0 / 10^k``: for each decade ``k`` the sweep reports the
    density ``eta_x(k)`` where the measured ``D(eta)`` first falls below
    ``D_0 / 10^k`` (log-linear interpolation between bracketing grid points;
    no extrapolation -- decades not bracketed by the data are absent, never
    invented).  ``drift_per_decade`` is the slope of ``eta_x`` against ``k``:
    a smooth-barrier scenario predicts a nonzero, roughly constant drift; a
    pinned feature predicts the structural observables stay put regardless.

    ``D_0`` is ``D`` at the lowest density point, so the ladder is a relative
    criterion exactly like ``tau_alpha/tau_0``.  This never measures the
    barrier itself; it measures only how much the apparent arrest location
    depends on the convention -- which is the point in dispute.
    """
    eta = np.asarray(eta, float)
    D = np.asarray(D, float)
    ok = np.isfinite(eta) & np.isfinite(D) & (D > 0)
    eta, D = eta[ok], D[ok]
    order = np.argsort(eta)
    eta, D = eta[order], D[order]
    if len(eta) < 3:
        return {"D0": float("nan"), "crossings": [],
                "drift_per_decade": float("nan"), "drift_err": float("nan")}
    D0 = D[0]
    logD = np.log10(D)
    crossings = []
    for k in range(1, max_decades + 1):
        th = np.log10(D0) - k
        below = np.nonzero(logD <= th)[0]
        if len(below) == 0 or below[0] == 0:
            continue
        j = below[0]
        i = j - 1
        f = (logD[i] - th) / (logD[i] - logD[j])
        crossings.append({"decade": k, "threshold": float(10.0 ** th),
                          "eta_x": float(eta[i] + f * (eta[j] - eta[i]))})
    if len(crossings) >= 3:
        ks = np.array([c["decade"] for c in crossings], float)
        ex = np.array([c["eta_x"] for c in crossings], float)
        A = np.vstack([ks, np.ones_like(ks)]).T
        coef, res, *_ = np.linalg.lstsq(A, ex, rcond=None)
        dof = len(ks) - 2
        err = float(np.sqrt(res[0] / dof / np.sum((ks - ks.mean()) ** 2))) \
            if dof > 0 and len(res) else float("nan")
        drift, drift_err = float(coef[0]), err
    else:
        drift, drift_err = float("nan"), float("nan")
    return {"D0": float(D0), "crossings": crossings,
            "drift_per_decade": drift, "drift_err": drift_err}

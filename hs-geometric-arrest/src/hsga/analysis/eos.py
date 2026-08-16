"""Equation of state from the contact value of ``g(r)`` -- gates G1 and G1-2D.

3D virial route: ``Z = 1 + 4 eta g(sigma+)`` against Carnahan-Starling.
2D virial route: ``Z = 1 + 2 phi g(sigma+)`` against Henderson.

Everything hinges on extrapolating ``g(r)`` to contact.  Two hard-won details:

* A bin straddling ``sigma`` is only partly accessible and its ``g`` is
  depressed by the inaccessible volume fraction; a single such bin in the fit
  window biases ``g(sigma+)`` low by tens of percent.  A bin edge is therefore
  placed exactly at contact.
* ``g(r)`` has a cusp at contact and curves away from a straight line;
  the linear extrapolation degrades above ``eta ~ 0.4`` (4% at 0.45 in the
  pilot).  The fit here is a low-order polynomial in ``ln g`` over a narrow
  window, and the window is ALWAYS scanned: the spread over windows is
  reported as the systematic, per spec rule 6.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "carnahan_starling",
    "compressibility",
    "contact_value",
    "henderson_2d",
    "radial_distribution",
]


def carnahan_starling(eta) -> np.ndarray:
    """Carnahan-Starling ``Z(eta)`` for monodisperse 3D hard spheres."""
    eta = np.asarray(eta, float)
    return (1.0 + eta + eta**2 - eta**3) / (1.0 - eta) ** 3


def henderson_2d(phi) -> np.ndarray:
    """Henderson ``Z(phi)`` for monodisperse hard disks (accurate ~1% in the fluid)."""
    phi = np.asarray(phi, float)
    return (1.0 + phi**2 / 8.0) / (1.0 - phi) ** 2


def radial_distribution(
    frames, *, dim: int = 3, nbin: int = 600, rmax_over_sigma: float = 2.5,
    nbin_below: int = 10
) -> dict:
    """``g(r)`` of a monodisperse system, pooled over frames, bin edge at contact.

    Bins below contact are returned too: their counts must be zero, which is a
    free overlap audit (a violation raises).  Refuses polydisperse input --
    a mixture needs the partial ``g_ab(sigma_ab+)``, which is not implemented
    (DEBT.md).
    """
    frames = list(frames)
    _, rad0, _ = frames[0]
    if float(np.ptp(rad0)) > 1e-12:
        raise ValueError(
            "radial_distribution is monodisperse-only; a mixture needs partial "
            "contact values g_ab(sigma_ab+), which are not implemented"
        )
    sigma = 2.0 * float(rad0.mean())
    edges = np.concatenate(
        [
            np.linspace(0.9 * sigma, sigma, nbin_below + 1)[:-1],
            np.linspace(sigma, rmax_over_sigma * sigma, nbin + 1),
        ]
    )
    hist = np.zeros(len(edges) - 1)
    npair = 0.0
    n_below = 0.0
    for pos, rad, L in frames:
        n = len(pos)
        d = pos[:, None, :] - pos[None, :, :]
        d -= L * np.round(d / L)
        r = np.sqrt((d**2).sum(-1))[np.triu_indices(n, 1)]
        hist += np.histogram(r, bins=edges)[0]
        n_below += float(np.count_nonzero(r < sigma * (1.0 - 1e-12)))
        # N(N-1)/2 distinct pairs, not N^2/2: the difference is a 1/N bias in g
        npair += 0.5 * n * (n - 1) / L**dim
    mid = 0.5 * (edges[1:] + edges[:-1])
    if dim == 3:
        shell = (4.0 * np.pi / 3.0) * (edges[1:] ** 3 - edges[:-1] ** 3)
    else:
        shell = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    g = hist / (npair * shell)
    if n_below > 0:
        raise ValueError(
            f"{n_below:.0f} pair separations below contact: the configurations overlap"
        )
    return {
        "r": mid, "g": g, "edges": edges, "sigma": sigma,
        "n_frames": len(frames), "counts": hist, "n_below_contact": n_below,
    }


def contact_value(
    r, g, sigma: float,
    *,
    windows=((1.0, 1.05), (1.0, 1.08), (1.0, 1.12), (1.0, 1.15)),
    degree: int = 2,
    counts=None,
) -> dict:
    """Cusp-aware extrapolation of ``g(r)`` to contact, with a window scan.

    Fits ``ln g`` as a polynomial in ``(r/sigma - 1)`` over each window and
    evaluates at contact.  The first window is the reported one; the spread
    over all windows is the systematic and is never dropped.
    """
    r = np.asarray(r, float)
    g = np.asarray(g, float)
    x = r / sigma - 1.0

    results = []
    for lo, hi in windows:
        sel = (r >= lo * sigma) & (r <= hi * sigma) & (g > 0)
        if np.count_nonzero(sel) < degree + 2:
            continue
        coef = np.polyfit(x[sel], np.log(g[sel]), degree)
        results.append(
            {"window": (lo, hi), "g_contact": float(np.exp(np.polyval(coef, 0.0))),
             "n_bins": int(sel.sum())}
        )
    if not results:
        return {"g_contact": float("nan"), "err_window": float("nan"),
                "err_count": float("nan"), "err": float("nan"), "windows": []}

    vals = np.array([w["g_contact"] for w in results])
    err_window = float(0.5 * (vals.max() - vals.min()))

    err_count = float("nan")
    if counts is not None:
        counts = np.asarray(counts, float)
        lo, hi = results[0]["window"]
        sel = (r >= lo * sigma) & (r <= hi * sigma)
        ntot = counts[sel].sum()
        if ntot > 0:
            err_count = float(results[0]["g_contact"] / np.sqrt(ntot))

    err = err_window if not np.isfinite(err_count) else max(err_window, err_count)
    return {
        "g_contact": float(results[0]["g_contact"]),
        "err_window": err_window,
        "err_count": err_count,
        "err": err,
        "windows": results,
        "degree": degree,
    }


def compressibility(eta: float, g_contact: float, err: float = 0.0, *, dim: int = 3) -> dict:
    """Virial ``Z`` from the contact value, against the reference EOS."""
    coef = 4.0 if dim == 3 else 2.0
    Z = 1.0 + coef * eta * g_contact
    Z_ref = float(carnahan_starling(eta)) if dim == 3 else float(henderson_2d(eta))
    dZ = coef * eta * err
    return {
        "eta": float(eta),
        "dim": dim,
        "Z": float(Z),
        "Z_err": float(dZ),
        "Z_ref": Z_ref,
        "reference": "Carnahan-Starling" if dim == 3 else "Henderson",
        "rel_dev": float(abs(Z / Z_ref - 1.0)),
        "rel_dev_err": float(dZ / Z_ref),
    }

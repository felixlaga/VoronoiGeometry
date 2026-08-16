"""Radical Voronoi analysis of simulated configurations, 2D and 3D.

Per configuration: the isoperimetric quotient of every cell (``Q_iso`` in 3D,
``q`` in 2D), face counts, and the relative face gaps

    g_ij = (r_ij - sigma_ij) / sigma_ij,   sigma_ij = R_i + R_j,

from which ``f_c(delta)`` follows.  Also the large-particle ``q`` distribution
and its submode weight -- de Graaf's own primary structural signal in 2D.

A warning that is physics, not code: ``f_c`` at fixed tolerance is NOT a
usable headline observable.  Exact contact has measure zero in a thermal
hard-particle fluid, and the repaired ``f_c(delta)`` is degenerate with the
median gap because the gap distribution is single-scale
(:func:`single_scale_collapse` measures that degeneracy directly).  The
mechanistic observable is ``P_wrap^(k)`` on the pre-registered refscore
(``analysis/pwrap.py``); the primary campaign estimator is ``eps*``
(``analysis/percolation.py``).
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from ..geometry.tangential import Cell, PoolTooSmall, cell_from_neighbours

__all__ = [
    "config_cells",
    "config_observables",
    "f_c_of_delta",
    "q_submode_weight",
    "single_scale_collapse",
]


def config_cells(
    pos,
    rad,
    L: float,
    *,
    r_query: float | None = None,
    area_tol: float = 1e-9,
    contact_tol: float = 1e-7,
    indices=None,
) -> list[Cell]:
    """Radical Voronoi cells of a periodic configuration (dim from ``pos``).

    Each particle's neighbour pool is a ball, enlarged automatically until the
    cell is provably closed, so the result does not depend on the initial
    guess.  ``Cell.face_neighbours`` carries the particle indices, which is
    what ``pwrap`` builds cluster graphs from.
    """
    pos = np.asarray(pos, float)
    rad = np.asarray(rad, float)
    n, dim = pos.shape
    wrapped = pos - L * np.floor(pos / L)
    tree = cKDTree(wrapped, boxsize=L)

    sigma_max = 2.0 * float(rad.max())
    r0 = float(r_query) if r_query else 3.0 * sigma_max
    targets = range(n) if indices is None else list(indices)

    cells: list[Cell] = []
    for i in targets:
        r = r0
        for _ in range(6):
            idx = np.asarray(tree.query_ball_point(wrapped[i], r), dtype=int)
            idx = idx[idx != i]
            if len(idx) >= dim + 1:
                d = pos[idx] - pos[i]
                d -= L * np.round(d / L)
                try:
                    cells.append(
                        cell_from_neighbours(
                            d, rad[i], rad[idx], idx,
                            area_tol=area_tol, contact_tol=contact_tol,
                        )
                    )
                    break
                except PoolTooSmall:
                    pass
            r *= 1.5
            if r > 0.5 * L:
                raise PoolTooSmall(
                    f"cell of particle {i} not closed within half the box"
                )
        else:
            raise PoolTooSmall(f"cell of particle {i} did not close")
    return cells


def config_observables(pos, rad, L: float, **kw) -> dict:
    """Per-configuration Voronoi observables.

    Arrays of ``Q_iso`` and ``n_faces`` per particle, pooled face gaps (each
    pair counted from both sides), and the per-particle radii for downstream
    species splits.
    """
    cells = config_cells(pos, rad, L, **kw)
    Q = np.array([c.Q_iso for c in cells])
    nf = np.array([c.n_faces for c in cells])
    gaps = np.concatenate([c.face_gaps for c in cells])
    return {
        "Q_iso": Q,
        "n_faces": nf,
        "face_gaps": gaps,
        "median_gap": float(np.median(gaps)),
        "mean_Q_iso": float(Q.mean()),
        "mean_n_faces": float(nf.mean()),
        "n_cells": len(cells),
        "radii": np.asarray(rad, float).copy(),
        "cells": cells,
    }


def f_c_of_delta(gaps, deltas) -> np.ndarray:
    """Fraction of faces with gap at most ``delta``, for each ``delta``."""
    gaps = np.asarray(gaps, float)
    deltas = np.atleast_1d(np.asarray(deltas, float))
    return np.array([float(np.count_nonzero(gaps <= d) / len(gaps)) for d in deltas])


def single_scale_collapse(gaps, deltas_over_median=(0.01, 0.02, 0.05, 0.1)) -> dict:
    """Does the face-gap distribution carry information beyond its median?

    If ``P(g)`` collapses under ``g -> g/median`` then ``f_c(delta, eta)`` is
    a function of ``delta/median`` alone; the diagnostic ratio
    ``f_c(delta)/(delta/median)`` is then constant across densities.
    """
    gaps = np.asarray(gaps, float)
    med = float(np.median(gaps))
    x = np.asarray(deltas_over_median, float)
    fc = f_c_of_delta(gaps, x * med)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(x > 0, fc / x, np.nan)
    return {"median_gap": med, "delta_over_median": x, "f_c": fc, "ratio": ratio,
            "n_faces": len(gaps)}


def q_submode_weight(
    Q, radii, q_ref: float, *, half_width: float = 0.01, species: str = "large"
) -> dict:
    """Weight of the ``q`` distribution near a reference value, per species.

    De Graaf's primary structural signal in 2D is a submode (cusp) in the
    large-particle ``q`` distribution near the reference tile's ``q``.  This
    measures the fraction of the selected species' cells with
    ``|q - q_ref| < half_width``, plus the local density ratio against the
    two flanking windows of the same width -- a cusp shows up as a ratio
    above 1.  Histograms are the campaign deliverable; this scalar tracks the
    signal per state point.
    """
    Q = np.asarray(Q, float)
    radii = np.asarray(radii, float)
    med = float(np.median(radii))
    if species == "large":
        sel = radii > med + 1e-12
        if not sel.any():                       # monodisperse: use everything
            sel = np.ones(len(Q), bool)
    elif species == "small":
        sel = radii < med - 1e-12
        if not sel.any():
            sel = np.ones(len(Q), bool)
    else:
        sel = np.ones(len(Q), bool)
    q = Q[sel]
    inw = np.count_nonzero(np.abs(q - q_ref) < half_width)
    lo = np.count_nonzero(np.abs(q - (q_ref - 2 * half_width)) < half_width)
    hi = np.count_nonzero(np.abs(q - (q_ref + 2 * half_width)) < half_width)
    flank = 0.5 * (lo + hi)
    return {
        "q_ref": float(q_ref),
        "half_width": float(half_width),
        "n_species": int(sel.sum()),
        "weight": float(inw / max(1, sel.sum())),
        "flank_ratio": float(inw / flank) if flank > 0 else float("nan"),
    }

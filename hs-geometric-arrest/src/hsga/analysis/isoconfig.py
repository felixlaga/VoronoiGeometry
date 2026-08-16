"""Isoconfigurational propensity and the held-out predictive-power test (spec 2.5e).

From one equilibrated configuration, ``M`` independent MC move sequences are
launched (the engines' restart mode with different seeds); the per-particle
propensity is the mean displacement at each lag time over the ensemble.  The
question the module answers is de Graaf's own standard: do the structural
fields of this project -- refscores ``s_k``, marked-cluster size, topological
persistence -- predict propensity BEYOND a baseline of known-boring fields
{local packing fraction, tetrahedrality (3D) / bond-hexatic order (2D),
Voronoi anisotropy, cell volume} on held-out configurations?  If they add
nothing, the hypothesis fails that standard.

Definitions fixed here, before any propensity data existed (documented
choices; there is no unique convention):

* tetrahedrality (3D): Errington-Debenedetti orientational order
  ``q = 1 - 3/8 sum_{j<k in 4 nn} (cos psi_jik + 1/3)^2``;
* bond-orientational order (2D): ``|psi_6|`` over the Voronoi neighbours;
* Voronoi anisotropy: ``1 - r_in / R_circ`` of the cell;
* local packing fraction: ``v_i / V_cell_i``;
* displacements at lag are minimum-image (valid while displacement < L/2,
  which holds at the near-arrest densities this analysis targets).

The heavy part (running the M restarts) is orchestrated by
``scripts/run_isoconfig.py`` and has NOT been executed in this repository;
everything here is exercised by unit tests on synthetic data plus a smoke run.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "baseline_features",
    "held_out_r2",
    "propensity_from_frames",
    "psi6_2d",
    "structural_features",
    "tetrahedrality_3d",
]


# --------------------------------------------------------------------------- #
# propensity
# --------------------------------------------------------------------------- #
def propensity_from_frames(start, runs, L: float) -> np.ndarray:
    """Per-particle propensity at each lag: mean over runs of |dr(lag)|.

    ``start`` is the common initial position array; ``runs`` is a list of
    frame-lists (one per restart), each frame an ``(pos, rad, L)`` tuple at
    the same lag schedule.  Returns ``(n_lags, N)``.
    """
    start = np.asarray(start, float)
    n_lags = min(len(fr) for fr in runs)
    out = np.zeros((n_lags, len(start)))
    for fr in runs:
        for k in range(n_lags):
            pos = np.asarray(fr[k][0], float)
            d = pos - start
            d -= L * np.round(d / L)
            out[k] += np.linalg.norm(d, axis=1)
    return out / len(runs)


# --------------------------------------------------------------------------- #
# per-particle fields
# --------------------------------------------------------------------------- #
def tetrahedrality_3d(pos, L: float) -> np.ndarray:
    """Errington-Debenedetti ``q`` over the 4 nearest neighbours."""
    from scipy.spatial import cKDTree

    pos = np.asarray(pos, float)
    wrapped = pos - L * np.floor(pos / L)
    tree = cKDTree(wrapped, boxsize=L)
    _, nn = tree.query(wrapped, k=5)
    q = np.empty(len(pos))
    for i in range(len(pos)):
        vecs = pos[nn[i, 1:]] - pos[i]
        vecs -= L * np.round(vecs / L)
        vecs /= np.linalg.norm(vecs, axis=1)[:, None]
        s = 0.0
        for a in range(4):
            for b in range(a + 1, 4):
                s += (vecs[a] @ vecs[b] + 1.0 / 3.0) ** 2
        q[i] = 1.0 - 3.0 / 8.0 * s
    return q


def psi6_2d(pos, L: float, cells) -> np.ndarray:
    """``|psi_6|`` per particle over its Voronoi neighbours."""
    pos = np.asarray(pos, float)
    out = np.empty(len(pos))
    for i, c in enumerate(cells):
        acc = 0.0 + 0.0j
        for j in c.face_neighbours:
            d = pos[int(j)] - pos[i]
            d -= L * np.round(d / L)
            theta = np.arctan2(d[1], d[0])
            acc += np.exp(6j * theta)
        out[i] = abs(acc) / max(1, c.n_faces)
    return out


def baseline_features(cells, pos, rad, L: float) -> tuple[np.ndarray, list[str]]:
    """The baseline field set {local eta, tetrahedrality/psi6, anisotropy, V_cell}."""
    pos = np.asarray(pos, float)
    rad = np.asarray(rad, float)
    dim = pos.shape[1]
    V = np.array([c.V for c in cells])
    if dim == 3:
        v_i = (4.0 / 3.0) * np.pi * rad**3
        angular = tetrahedrality_3d(pos, L)
        angular_name = "tetrahedrality"
    else:
        v_i = np.pi * rad**2
        angular = psi6_2d(pos, L, cells)
        angular_name = "psi6"
    aniso = np.array([1.0 - c.r_in / c.R_circ for c in cells])
    X = np.column_stack([v_i / V, angular, aniso, V])
    return X, ["local_eta", angular_name, "anisotropy", "V_cell"]


def structural_features(cells, pos, rad, L: float) -> tuple[np.ndarray, list[str]]:
    """The candidate field set {s_k for every frozen reference, cluster size,
    persistence}, computed with the pre-registered weights and filters."""
    from .pwrap import marked_bonds
    from .percolation import UnionFind
    from .refscore import load_frozen, mark_cells, refscores
    from .topology import persistence as topo_persistence

    pos = np.asarray(pos, float)
    dim = pos.shape[1]
    scores, names = refscores(cells, rad, dim)

    # cluster size of the marked component the particle belongs to (0 if unmarked),
    # for the reference with the lowest median score (the most-populated class)
    best = int(np.argmin(np.median(scores, axis=0)))
    marked = mark_cells(scores, names, names[best])
    uf = UnionFind(len(cells), dim)
    for i, j, img in marked_bonds(cells, marked, pos, L):
        uf.union(i, j, img)
    comp: dict[int, int] = {}
    for i in np.flatnonzero(marked):
        r, _ = uf.find(int(i))
        comp[r] = comp.get(r, 0) + 1
    csize = np.zeros(len(cells))
    for i in np.flatnonzero(marked):
        r, _ = uf.find(int(i))
        csize[i] = comp[r]

    fz = load_frozen()["persistence_grid"]
    pers = np.array([topo_persistence(c, fz["g_grid"], fz["a_grid"]) for c in cells])

    X = np.column_stack([scores, csize, pers])
    return X, [f"s_{n}" for n in names] + [f"cluster_size[{names[best]}]", "persistence"]


# --------------------------------------------------------------------------- #
# held-out evaluation
# --------------------------------------------------------------------------- #
def _ridge_r2(X_tr, y_tr, X_te, y_te, lam: float = 1e-3) -> float:
    """R^2 on the held-out set of a ridge regression fit on the training set."""
    mu, sd = X_tr.mean(0), X_tr.std(0) + 1e-12
    Xt = (X_tr - mu) / sd
    Xe = (X_te - mu) / sd
    Xt = np.column_stack([np.ones(len(Xt)), Xt])
    Xe = np.column_stack([np.ones(len(Xe)), Xe])
    A = Xt.T @ Xt + lam * np.eye(Xt.shape[1])
    w = np.linalg.solve(A, Xt.T @ y_tr)
    resid = y_te - Xe @ w
    ss_res = float((resid**2).sum())
    ss_tot = float(((y_te - y_te.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def held_out_r2(feature_blocks, propensities, *, lam: float = 1e-3) -> dict:
    """Held-out predictive power: baseline vs baseline+structural, by configuration.

    ``feature_blocks`` is a list per configuration of ``(X_base, X_struct)``;
    ``propensities`` the matching per-particle propensity vectors.  Splitting
    is BY CONFIGURATION (never by particle): particles of one configuration
    share their environment and would leak across a per-particle split.
    Every leave-one-configuration-out fold is evaluated; the deliverable is
    ``delta_R2 = R2(base+struct) - R2(base)`` per fold, with its spread.
    """
    n_cfg = len(feature_blocks)
    if n_cfg < 2:
        raise ValueError("held-out evaluation needs >= 2 configurations")
    rows = []
    for hold in range(n_cfg):
        tr = [k for k in range(n_cfg) if k != hold]
        Xb_tr = np.vstack([feature_blocks[k][0] for k in tr])
        Xs_tr = np.vstack([np.hstack(feature_blocks[k]) for k in tr])
        y_tr = np.concatenate([propensities[k] for k in tr])
        Xb_te = feature_blocks[hold][0]
        Xs_te = np.hstack(feature_blocks[hold])
        y_te = propensities[hold]
        r2_base = _ridge_r2(Xb_tr, y_tr, Xb_te, y_te, lam)
        r2_full = _ridge_r2(Xs_tr, y_tr, Xs_te, y_te, lam)
        rows.append({"fold": hold, "r2_base": r2_base, "r2_full": r2_full,
                     "delta_r2": r2_full - r2_base})
    deltas = np.array([r["delta_r2"] for r in rows])
    return {
        "folds": rows,
        "delta_r2_mean": float(deltas.mean()),
        "delta_r2_spread": float(deltas.std(ddof=1)) if len(deltas) > 1 else float("nan"),
        "verdict_note": (
            "delta_r2 <= 0 within its spread means the structural fields add no "
            "held-out predictive power beyond the baseline: the hypothesis "
            "fails de Graaf's own standard"
        ),
    }

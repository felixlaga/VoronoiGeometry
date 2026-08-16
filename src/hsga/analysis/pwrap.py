"""``P_wrap^(k)`` -- the mechanistic observable (spec 2.5d).

Cells whose pre-registered refscore marks them as ``k``-like are connected
when they share a kept radical face; the marked-cell network is then tested
for a wrapping cluster with the same relative-displacement union-find as the
shell percolation.  This replaces the degenerate ``f_c`` observable: it probes
the *arrangement* of reference-like cells, not their number.

Also cluster statistics (size distribution, susceptibility, correlation
length via the radius of gyration) and the finite-size crossing analysis of
``P_wrap(eta; N)`` curves.
"""

from __future__ import annotations

import numpy as np

from ..geometry.tangential import Cell
from .percolation import UnionFind
from .refscore import load_frozen, mark_cells, refscores
from .topology import face_filter

__all__ = [
    "cluster_statistics",
    "finite_size_crossings",
    "marked_bonds",
    "p_wrap",
    "p_wrap_all_references",
]


def marked_bonds(cells, marked, pos, L):
    """Bonds between marked cells sharing a KEPT radical face.

    Face adjacency comes from ``Cell.face_neighbours`` after the frozen face
    filter -- the same filter the refscore itself used.  Bond image offsets are
    the minimum-image shifts of the particle positions.
    """
    fz = load_frozen()["face_filter"]
    pos = np.asarray(pos, float)
    marked_idx = set(np.flatnonzero(marked))
    bonds = []
    for i in marked_idx:
        c = face_filter(cells[i], fz["g_cut"], fz["a_cut"])
        for j in c.face_neighbours:
            j = int(j)
            if j in marked_idx and j > i:
                d = pos[i] - pos[j]
                img = np.round(d / L).astype(np.int64)
                bonds.append((int(i), j, img))
    return bonds


def p_wrap(cells, pos, rad, L, k_name: str, *, dim: int | None = None,
           scores=None, names=None) -> dict:
    """Does the ``k``-marked cell network wrap the box in all directions?

    Returns the wrap verdict, per-direction wraps, the marked fraction, and
    the cluster statistics.  ``scores``/``names`` can be passed to reuse a
    refscore evaluation across references.
    """
    pos = np.asarray(pos, float)
    if dim is None:
        dim = pos.shape[1]
    if scores is None:
        scores, names = refscores(cells, rad, dim)
    marked = mark_cells(scores, names, k_name)
    n = len(cells)
    uf = UnionFind(n, dim)
    bonds = marked_bonds(cells, marked, pos, L)
    for i, j, img in bonds:
        uf.union(i, j, img)
    stats = cluster_statistics(uf, marked, pos, L)
    return {
        "reference": k_name,
        "wraps": bool(uf.wrap.all()),
        "wrap_directions": uf.wrap.copy(),
        "marked_fraction": float(np.count_nonzero(marked) / n),
        "n_bonds": len(bonds),
        **stats,
    }


def p_wrap_all_references(cells, pos, rad, L, *, dim: int | None = None) -> dict:
    """``p_wrap`` for every frozen reference, sharing one refscore evaluation."""
    pos = np.asarray(pos, float)
    if dim is None:
        dim = pos.shape[1]
    scores, names = refscores(cells, rad, dim)
    return {k: p_wrap(cells, pos, rad, L, k, dim=dim, scores=scores, names=names)
            for k in names}


def cluster_statistics(uf: UnionFind, marked, pos, L) -> dict:
    """Cluster sizes, susceptibility and correlation length of the marked network.

    ``chi = sum_s s^2 n_s / sum_s s n_s`` and
    ``xi^2 = sum_s 2 R_g^2(s) s^2 n_s / sum_s s^2 n_s``, both excluding
    wrapping clusters (standard percolation practice; with none excluded the
    numbers diverge with the box instead of measuring correlation).
    """
    pos = np.asarray(pos, float)
    marked_idx = np.flatnonzero(marked)
    if not len(marked_idx):
        return {"cluster_sizes": [], "chi": 0.0, "xi": 0.0, "largest": 0}
    roots: dict[int, list[int]] = {}
    for i in marked_idx:
        r, _ = uf.find(int(i))
        roots.setdefault(r, []).append(int(i))

    sizes, weights, rg2 = [], [], []
    wrapping_roots = set()
    # a cluster wraps if any pair of members disagrees on unwrapped offsets --
    # the union-find already accumulated that in uf.wrap globally; per-cluster
    # wrap detection uses the offset spread:
    for r, members in roots.items():
        offs = np.array([uf.find(m)[1] for m in members])
        if len(members) > 1 and (np.ptp(offs, axis=0) * L > 0.6 * L).any():
            wrapping_roots.add(r)
    for r, members in roots.items():
        s = len(members)
        sizes.append(s)
        if r in wrapping_roots:
            continue
        # unfolded member positions: fold + image offset from the union-find
        offs = np.array([uf.find(m)[1] for m in members], dtype=float)
        unfolded = pos[members] + offs * L
        c = unfolded.mean(axis=0)
        rg2.append(float(((unfolded - c) ** 2).sum(axis=1).mean()))
        weights.append(s)

    sizes_arr = np.array(sizes)
    if weights:
        w = np.array(weights, dtype=float)
        chi = float((w**2).sum() / w.sum())
        r2 = np.array(rg2)
        xi = float(np.sqrt((2 * r2 * w**2).sum() / (w**2).sum()))
    else:
        chi, xi = 0.0, 0.0
    return {
        "cluster_sizes": sorted(sizes, reverse=True),
        "chi": chi,
        "xi": xi,
        "largest": int(sizes_arr.max()) if len(sizes_arr) else 0,
        "n_wrapping_clusters": len(wrapping_roots),
    }


def finite_size_crossings(curves: dict) -> list[dict]:
    """Pairwise crossings of ``P_wrap(eta)`` curves at different N.

    ``curves`` maps ``N -> (eta_array, P_array)``.  A size-independent crossing
    is the standard percolation-transition locator; the scatter of the pairwise
    crossings is its systematic and is returned, never averaged away.
    """
    out = []
    Ns = sorted(curves)
    for a in range(len(Ns)):
        for b in range(a + 1, len(Ns)):
            e1, p1 = (np.asarray(v, float) for v in curves[Ns[a]])
            e2, p2 = (np.asarray(v, float) for v in curves[Ns[b]])
            common = np.intersect1d(np.round(e1, 9), np.round(e2, 9))
            if len(common) < 2:
                continue
            d = (np.interp(common, e1, p1) - np.interp(common, e2, p2))
            sign = np.sign(d)
            for k in range(len(common)):
                if sign[k] == 0:                       # exact crossing on a grid point
                    out.append({"N_pair": (Ns[a], Ns[b]),
                                "eta_cross": float(common[k])})
                elif k + 1 < len(common) and sign[k + 1] != 0 and sign[k] != sign[k + 1]:
                    t = d[k] / (d[k] - d[k + 1])
                    out.append({
                        "N_pair": (Ns[a], Ns[b]),
                        "eta_cross": float(common[k] + t * (common[k + 1] - common[k])),
                    })
    return out

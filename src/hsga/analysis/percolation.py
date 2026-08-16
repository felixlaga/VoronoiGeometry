"""Shell percolation: the primary estimator of the campaign, in 2D and 3D.

Every particle is inflated by ``eps`` (in units of the mean diameter) and the
contact network -- bonds where ``r_ij < R_i + R_j + 2 eps <sigma>`` -- is
tested for a cluster that wraps the periodic box in **all** directions.
``eps_star`` is the smallest inflation at which that happens.

Wrapping is detected with union-find carrying **relative displacements**
(spec 2.5): each particle stores the integer image offset to its parent; when
a bond joins two particles already in the same cluster, a non-zero mismatch
around the loop means the cluster closes on a periodic image of itself, i.e.
wraps.  Merely touching both walls is not percolation -- a compact cluster can
span the box in extent without connecting around it.

Two estimators are provided and cross-checked in the tests:

* :func:`eps_star` -- exact: pairs are sorted by the inflation at which their
  bond appears and added in order; the threshold is the bond that completes
  wrapping in the last direction.  No bracketing resolution enters.
* :func:`eps_star_bisect` -- the literal bisection of the spec, kept as the
  independent cross-check (the union-find port was validated by recovering 2D
  RCP to 1%; see ``results/validation_pilot.md``).
"""

from __future__ import annotations

import itertools

import numpy as np

__all__ = [
    "UnionFind",
    "contact_bonds",
    "eps_star",
    "eps_star_bisect",
    "percolates",
    "wrapping_directions",
]


class UnionFind:
    """Union-find over periodic images with relative-displacement tracking.

    For a bond between ``i`` and ``j`` with minimum-image shift ``img``
    (``r_ij = pos_i - pos_j - L*img``), the unwrapped image indices satisfy
    ``n_i - n_j = -img``.  Accumulated around a loop that returns to the same
    root, a non-zero residual component is a wrap in that direction.
    """

    __slots__ = ("parent", "offset", "wrap", "dim")

    def __init__(self, n: int, dim: int = 3):
        self.dim = dim
        self.parent = np.arange(n)
        self.offset = np.zeros((n, dim), dtype=np.int64)   # n(self) - n(parent)
        self.wrap = np.zeros(dim, dtype=bool)

    def find(self, a: int):
        root = a
        disp = np.zeros(self.dim, dtype=np.int64)
        while self.parent[root] != root:
            disp += self.offset[root]
            root = self.parent[root]
        node, acc = a, disp.copy()
        while self.parent[node] != node:               # path compression
            nxt = self.parent[node]
            nxt_acc = acc - self.offset[node]
            self.parent[node] = root
            self.offset[node] = acc
            node, acc = nxt, nxt_acc
        return root, disp

    def union(self, a: int, b: int, img) -> None:
        img = np.asarray(img, dtype=np.int64)
        ra, da = self.find(a)
        rb, db = self.find(b)
        if ra == rb:
            residual = (da - db) + img
            self.wrap |= residual != 0
        else:
            self.parent[ra] = rb
            self.offset[ra] = db - da - img


def contact_bonds(pos: np.ndarray, rad: np.ndarray, L: float, eps_abs: float):
    """Yield ``(i, j, img)`` for every pair with ``r_ij < R_i + R_j + 2 eps_abs``."""
    pos = np.asarray(pos, float)
    rad = np.asarray(rad, float)
    n, dim = pos.shape
    cut = 2.0 * float(rad.max()) + 2.0 * eps_abs
    ncell = int(L / cut)

    if ncell < 3:
        for i in range(n):
            d = pos[i] - pos[i + 1:]
            img = np.round(d / L)
            d = d - L * img
            s = rad[i] + rad[i + 1:] + 2.0 * eps_abs
            for h in np.flatnonzero((d * d).sum(axis=1) < s * s):
                yield i, int(i + 1 + h), img[h].astype(np.int64)
        return

    cs = L / ncell
    idx = (pos / cs).astype(int) % ncell
    cells: dict[tuple, list[int]] = {}
    for i in range(n):
        cells.setdefault(tuple(idx[i]), []).append(i)
    hood = list(itertools.product((-1, 0, 1), repeat=dim))
    for i in range(n):
        cand: list[int] = []
        for off in hood:
            key = tuple((idx[i, k] + off[k]) % ncell for k in range(dim))
            cand.extend(cells.get(key, ()))
        cand = np.asarray([j for j in cand if j > i], dtype=int)
        if not len(cand):
            continue
        d = pos[i] - pos[cand]
        img = np.round(d / L)
        d = d - L * img
        s = rad[i] + rad[cand] + 2.0 * eps_abs
        for h in np.flatnonzero((d * d).sum(axis=1) < s * s):
            yield i, int(cand[h]), img[h].astype(np.int64)


def wrapping_directions(pos, rad, L: float, eps: float) -> np.ndarray:
    """Boolean per axis: does the ``eps``-inflated contact network wrap?"""
    pos = np.asarray(pos, float)
    rad = np.asarray(rad, float)
    eps_abs = float(eps) * 2.0 * float(rad.mean())
    uf = UnionFind(len(pos), pos.shape[1])
    for i, j, img in contact_bonds(pos, rad, L, eps_abs):
        uf.union(i, j, img)
        if uf.wrap.all():
            break
    return uf.wrap


def percolates(pos, rad, L: float, eps: float) -> bool:
    """True when the network wraps in ALL directions (spec 2.5)."""
    return bool(wrapping_directions(pos, rad, L, eps).all())


def eps_star_bisect(
    pos, rad, L: float, *, lo: float = 1e-5, hi: float = 0.2, iters: int = 18
) -> float:
    """Threshold by bisection in ``log eps`` -- the spec's literal algorithm.

    Monotone in ``eps`` (inflating only adds bonds), so bisection converges;
    returns ``NaN`` when the bracket fails rather than an endpoint.
    :func:`eps_star` computes the same number exactly and faster; this stays
    as its independent cross-check.
    """
    if percolates(pos, rad, L, lo):
        return float("nan")
    if not percolates(pos, rad, L, hi):
        return float("nan")
    a, b = np.log(lo), np.log(hi)
    for _ in range(iters):
        m = 0.5 * (a + b)
        if percolates(pos, rad, L, float(np.exp(m))):
            b = m
        else:
            a = m
    return float(np.exp(0.5 * (a + b)))


def eps_star(pos, rad, L: float, *, lo: float = 1e-5, hi: float = 0.2) -> float:
    """Smallest inflation at which the network wraps the box, exactly.

    Each pair carries the inflation at which its bond appears,
    ``eps_ij = (r_ij - R_i - R_j) / (2 <sigma>)``; adding bonds in that order
    makes the threshold the ``eps_ij`` of the bond completing the last wrap.
    Returns ``NaN`` outside the ``[lo, hi]`` bracket, matching
    :func:`eps_star_bisect`.
    """
    from scipy.spatial import cKDTree

    pos = np.asarray(pos, float)
    rad = np.asarray(rad, float)
    dim = pos.shape[1]
    sigma_mean = 2.0 * float(rad.mean())
    cut = 2.0 * float(rad.max()) + 2.0 * hi * sigma_mean

    wrapped = pos - L * np.floor(pos / L)
    tree = cKDTree(wrapped, boxsize=L)
    pairs = tree.query_pairs(r=cut, output_type="ndarray")
    if not len(pairs):
        return float("nan")

    i, j = pairs[:, 0], pairs[:, 1]
    d = pos[i] - pos[j]
    img = np.round(d / L)
    d = d - L * img
    r = np.linalg.norm(d, axis=1)
    eps_ij = (r - rad[i] - rad[j]) / (2.0 * sigma_mean)

    keep = eps_ij <= hi
    i, j, img, eps_ij = i[keep], j[keep], img[keep].astype(np.int64), eps_ij[keep]
    order = np.argsort(eps_ij, kind="stable")

    uf = UnionFind(len(pos), dim)
    for k in order:
        uf.union(int(i[k]), int(j[k]), img[k])
        if uf.wrap.all():
            e = float(eps_ij[k])
            return float("nan") if e < lo else e
    return float("nan")

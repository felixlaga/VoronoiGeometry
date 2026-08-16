"""Tangential Voronoi geometry: the exact analogue of ``phi = q`` in any dimension.

For a monodisperse packing of spheres (disks) of diameter ``sigma`` whose
Voronoi cells are *tangential* polytopes with the particle as insphere,
decomposing the cell into pyramids of height ``sigma/2`` erected on each face
gives ``V = sigma S / (2 d)``, i.e. the exact identity ``S = d V / R``, and
therefore

    3D:  eta = (pi/6) sigma^3 / V = 36 pi V^2 / S^3 = Q_iso = Psi^3
    2D:  phi = (pi/4) sigma^2 / A = 4 pi A / P^2    = q

(`paper.tex` Sec. II; verified as an exact symbolic identity by the sympy
certificates in ``tests/test_certificates.py``).  The hypothesis has a single
sharp reading: *every Voronoi neighbour must be a contact*.  If one face is
supported by a non-touching neighbour the pyramid decomposition fails and
``eta != Q_iso``.

This module builds radical (Laguerre) Voronoi cells by explicit half-space
intersection, in 2D and 3D through one code path.  The radical plane between
``i`` and ``j`` sits at

    d_ij = (|r_ij|^2 + R_i^2 - R_j^2) / (2 |r_ij|)

from ``i`` along ``r_ij``; with equal radii this is the ordinary Voronoi cell.
The strictly interior point that ``scipy.spatial.HalfspaceIntersection``
requires is the site itself, which works because every offset is positive for
a valid hard-particle configuration (landmine 6 of the task list).

Cells carry their full face combinatorics (vertex loops per face), which is
what ``analysis/topology.py`` builds the persistent-topology classifier on.

Invariants asserted in code rather than trusted:

* ``sum_k A_k h_k / d == V`` and ``sum_k A_k == S`` -- the face decomposition
  must rebuild the hull it came from; a violation is a bug in the face
  extraction, never a property of the structure.
* whenever ``is_tangential`` is true, ``abs(eta - Q_iso) < 1e-9`` -- the
  theorem; a violation means a bug, not a discovery.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

#: relative resolution at which two half-spaces count as the same plane
_PLANE_TOL = 1e-9

__all__ = [
    "Cell",
    "OverlapError",
    "PoolTooSmall",
    "analyse",
    "cell_from_neighbours",
    "is_tangential",
    "report",
    "voronoi_cell",
]


class PoolTooSmall(RuntimeError):
    """The supplied neighbour pool does not provably bound the cell."""


class OverlapError(ValueError):
    """The input configuration contains overlapping particles.

    Face gaps, contact counts and the tangentiality test are all defined
    relative to contact, so an overlap is bad input, not an edge case.
    """


# --------------------------------------------------------------------------- #
# the cell record
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Cell:
    """A single radical Voronoi cell in ``dim`` dimensions.

    ``V`` is volume (3D) or area (2D); ``S`` surface area or perimeter, summed
    over non-degenerate faces only.  ``Q_iso`` is the isoperimetric quotient of
    the dimension: ``36 pi V^2 / S^3`` in 3D (the cube of the Wadell
    sphericity) and ``4 pi V / S^2`` in 2D (de Graaf's ``q``).

    ``face_distances`` are the distances ``h_k`` from the site to each face
    plane -- for a tangential cell every entry equals ``sigma/2``.
    ``face_separations`` and ``face_sigmas`` are ``|r_ij|`` and ``R_i + R_j``
    for each face's generator, ``face_neighbours`` the generator's label.
    ``vertices`` are the merged cell vertices (site-centred coordinates) and
    ``face_vertices[k]`` the ordered vertex-index loop of face ``k`` (a vertex
    pair in 2D), from which ``analysis/topology.py`` builds p-vectors and the
    face-adjacency graph.
    """

    dim: int
    V: float
    S: float
    Q_iso: float
    n_faces: int
    face_distances: np.ndarray
    face_areas: np.ndarray
    face_separations: np.ndarray
    face_sigmas: np.ndarray
    face_neighbours: np.ndarray
    vertices: np.ndarray
    face_vertices: tuple
    r_in: float
    R_circ: float
    contact_tol: float = 1e-7

    @property
    def face_gaps(self) -> np.ndarray:
        """Relative face gaps ``g_ij = (r_ij - sigma_ij) / sigma_ij``."""
        return (self.face_separations - self.face_sigmas) / self.face_sigmas

    @property
    def n_contacts(self) -> int:
        """Number of faces whose generator touches at ``contact_tol``."""
        return int(np.count_nonzero(self.face_gaps <= self.contact_tol))

    @property
    def f_c(self) -> float:
        """Contact fraction of the faces."""
        return self.n_contacts / self.n_faces


def _iso_quotient(dim: int, V: float, S: float) -> float:
    if dim == 3:
        return 36.0 * np.pi * V**2 / S**3
    return 4.0 * np.pi * V / S**2


# --------------------------------------------------------------------------- #
# face extraction
# --------------------------------------------------------------------------- #
def _bounded(normals: np.ndarray) -> bool:
    """Do the half-spaces ``n_k . x <= h_k`` (all ``h_k > 0``) bound a polytope?

    They do exactly when the origin is interior to the convex hull of the
    normals.  Checking before calling Qhull matters: an unbounded region makes
    ``HalfspaceIntersection`` return points that are vertices of nothing, and
    the resulting volume is silently meaningless.
    """
    if len(normals) < normals.shape[1] + 1:
        return False
    try:
        hull = ConvexHull(normals)
    except Exception:
        return False
    return bool(np.all(hull.equations[:, -1] < -1e-12))


def _merge_vertices(verts: np.ndarray, tol: float) -> tuple[np.ndarray, np.ndarray]:
    """Merge near-duplicate vertices; returns (unique_vertices, index_map).

    Qhull duplicates intersection points wherever more than ``dim`` planes meet
    (every vertex of a cube or a rhombic dodecahedron); without merging, the
    polygon vertex counts -- and hence every p-vector -- come out wrong.
    """
    key = np.round(verts / tol).astype(np.int64)
    _, first, inverse = np.unique(key, axis=0, return_index=True, return_inverse=True)
    return verts[first], inverse


def _face_loop(pts_idx: np.ndarray, verts: np.ndarray, normal: np.ndarray, dim: int):
    """Order the vertices of one face into a loop; return (loop, measure)."""
    pts = verts[pts_idx]
    if dim == 2:
        if len(pts_idx) < 2:
            return None, 0.0
        # a 2D face is a segment; keep the two extreme points
        t = pts @ np.array([-normal[1], normal[0]])
        order = np.argsort(t)
        loop = pts_idx[[order[0], order[-1]]]
        return loop, float(t[order[-1]] - t[order[0]])
    if len(pts_idx) < 3:
        return None, 0.0
    centre = pts.mean(axis=0)
    rel = pts - centre
    seed = np.array([1.0, 0.0, 0.0])
    if abs(seed @ normal) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    u = np.cross(normal, seed)
    u /= np.linalg.norm(u)
    w = np.cross(normal, u)
    ang = np.arctan2(rel @ w, rel @ u)
    order = np.argsort(ang)
    loop = pts_idx[order]
    xs, ys = (rel @ u)[order], (rel @ w)[order]
    area = 0.5 * abs(np.dot(xs, np.roll(ys, -1)) - np.dot(ys, np.roll(xs, -1)))
    return loop, float(area)


def _build(normals: np.ndarray, offsets: np.ndarray, area_tol: float, dim: int):
    """Intersect ``n_k . x <= h_k`` (site at origin) and split into faces.

    Returns ``(V, S, live, areas, loops, verts, R_circ)``.
    """
    if not _bounded(normals):
        raise PoolTooSmall(
            f"the {len(normals)} candidate half-spaces do not bound a polytope; "
            "the neighbour pool is incomplete"
        )
    hs = np.hstack([normals, -offsets[:, None]])
    inter = HalfspaceIntersection(hs, np.zeros(dim))
    raw = inter.intersections
    R_circ = float(np.linalg.norm(raw, axis=1).max())
    scale = max(R_circ, 1e-30)
    verts, _ = _merge_vertices(raw, 1e-9 * scale)
    hull = ConvexHull(raw)

    ftol = 1e-9 * max(scale, 1.0) + 1e-12
    resid = verts @ normals.T - offsets[None, :]            # (n_vert, n_plane)
    live, areas, loops = [], [], []
    for k in range(len(offsets)):
        on = np.flatnonzero(np.abs(resid[:, k]) <= ftol)
        loop, a = _face_loop(on, verts, normals[k], dim)
        if loop is not None and a > area_tol:
            live.append(k)
            areas.append(a)
            loops.append(tuple(int(i) for i in loop))
    live = np.asarray(live, dtype=int)
    areas = np.asarray(areas, dtype=float)

    V = float(hull.volume)                 # 2D: area
    S_hull = float(hull.area)              # 2D: perimeter
    S = float(areas.sum())
    # Self-check: the pyramid decomposition about the site must rebuild the
    # hull.  Catches a mis-assigned or missed face immediately.
    V_pyr = float((areas * offsets[live]).sum() / dim)
    if abs(V_pyr - V) / max(abs(V), 1e-300) > 1e-6 or abs(S - S_hull) / max(S_hull, 1e-300) > 1e-6:
        raise RuntimeError(
            "face decomposition inconsistent with the hull: "
            f"V={V:.12g} V_pyramids={V_pyr:.12g} S={S:.12g} S_hull={S_hull:.12g}"
        )
    return V, S, live, areas, loops, verts, R_circ


def cell_from_neighbours(
    disp: np.ndarray,
    radius: float,
    neigh_radii: np.ndarray,
    neigh_ids: np.ndarray,
    *,
    area_tol: float = 1e-8,
    contact_tol: float = 1e-7,
    seed_planes: int = 48,
    check_overlap: bool = True,
    overlap_tol: float = 1e-9,
) -> Cell:
    """Radical Voronoi cell of a site from an explicit neighbour pool.

    ``disp`` are neighbour positions relative to the site (shape ``(n, dim)``),
    ``radius`` the site radius, ``neigh_radii`` the neighbour radii and
    ``neigh_ids`` an arbitrary per-neighbour label carried through to
    ``Cell.face_neighbours``.

    Half-spaces are added in order of increasing plane distance.  A plane at
    distance ``h`` can only cut a polytope of circumradius ``R_circ`` when
    ``h < R_circ``, so once the first excluded plane lies beyond the current
    circumradius the cell is final and exact -- the pool sufficiency is proved,
    not assumed; if the pool runs out first, :class:`PoolTooSmall` is raised
    and the caller enlarges it.
    """
    disp = np.asarray(disp, dtype=float)
    dim = disp.shape[1]
    dist = np.linalg.norm(disp, axis=1)
    keep = dist > 1e-12
    disp, dist = disp[keep], dist[keep]
    neigh_radii = np.asarray(neigh_radii, dtype=float)[keep]
    neigh_ids = np.asarray(neigh_ids)[keep]
    if len(dist) < dim + 1:
        raise PoolTooSmall(f"only {len(dist)} candidate neighbours supplied")

    if check_overlap:
        sig = radius + neigh_radii
        worst = float(np.min(dist - sig))
        if worst < -overlap_tol * float(np.max(sig)):
            raise OverlapError(
                f"closest pair overlaps by {-worst:.6g} (r_ij = {float(dist.min()):.6g}); "
                "not a valid hard-particle configuration"
            )

    normals = disp / dist[:, None]
    offsets = (dist**2 + radius**2 - neigh_radii**2) / (2.0 * dist)
    if offsets.min() <= 0.0:
        raise OverlapError(
            f"radical plane at h={offsets.min():.6g} <= 0; the site lies outside "
            "its own radical cell, which non-overlapping particles cannot do"
        )

    order = np.argsort(offsets)
    normals, offsets = normals[order], offsets[order]
    dist, neigh_radii, neigh_ids = dist[order], neigh_radii[order], neigh_ids[order]

    # coincident planes would each claim the same face and double-count it
    scale = float(offsets.min())
    key = np.round(
        np.hstack([normals, offsets[:, None] / scale]) / _PLANE_TOL
    ).astype(np.int64)
    _, first = np.unique(key, axis=0, return_index=True)
    first = np.sort(first)
    normals, offsets = normals[first], offsets[first]
    dist, neigh_radii, neigh_ids = dist[first], neigh_radii[first], neigh_ids[first]
    n = len(offsets)

    m = min(max(seed_planes, dim + 1), n)
    while True:
        V, S, live, areas, loops, verts, R_circ = _build(
            normals[:m], offsets[:m], area_tol, dim
        )
        if m == n:
            if offsets[-1] <= R_circ:
                raise PoolTooSmall(
                    f"outermost candidate plane h={offsets[-1]:.6g} is inside the "
                    f"circumradius {R_circ:.6g}: enlarge the neighbour pool"
                )
            break
        if offsets[m] > R_circ:
            break
        m = min(n, 2 * m)

    sig = radius + neigh_radii[live]
    return Cell(
        dim=dim,
        V=V,
        S=S,
        Q_iso=_iso_quotient(dim, V, S),
        n_faces=int(len(live)),
        face_distances=offsets[live].copy(),
        face_areas=areas.copy(),
        face_separations=dist[live].copy(),
        face_sigmas=sig,
        face_neighbours=neigh_ids[live].copy(),
        vertices=verts,
        face_vertices=tuple(loops),
        r_in=float(offsets[live].min()),
        R_circ=R_circ,
        contact_tol=contact_tol,
    )


# --------------------------------------------------------------------------- #
# lattice / basis entry point
# --------------------------------------------------------------------------- #
def _shifts(lattice: np.ndarray, nshell: int) -> np.ndarray:
    dim = lattice.shape[0]
    rng = range(-nshell, nshell + 1)
    n = np.array(list(itertools.product(rng, repeat=dim)), dtype=float)
    return n @ lattice


def _min_separation(lattice: np.ndarray, basis: np.ndarray) -> float:
    """Smallest non-zero distance between sites of the periodic structure."""
    shifts = _shifts(lattice, 2)
    best = np.inf
    for b in basis:
        d = np.linalg.norm(
            (basis[None, :, :] + shifts[:, None, :]).reshape(-1, basis.shape[1]) - b,
            axis=1,
        )
        d = d[d > 1e-9]
        if len(d):
            best = min(best, float(d.min()))
    if not np.isfinite(best):
        raise ValueError("structure has no pair of distinct sites")
    return best


def voronoi_cell(
    site,
    lattice,
    basis,
    *,
    radii=None,
    nshell: int = 3,
    area_tol: float = 1e-8,
    contact_tol: float = 1e-7,
) -> Cell:
    """Radical Voronoi cell of ``site`` in the periodic structure ``(lattice, basis)``.

    ``lattice`` holds the lattice vectors as rows (``d x d``); ``basis`` the
    Cartesian site positions of one period; ``site`` the position whose cell is
    wanted (it need not be an element of ``basis``, which the isolated-vacancy
    analysis needs).  ``radii=None`` means the monodisperse packing at contact:
    every radius is half the smallest inter-site distance.  ``nshell`` is only
    a starting value; it is enlarged until the pool provably bounds the cell.
    """
    lattice = np.asarray(lattice, dtype=float)
    dim = lattice.shape[0]
    lattice = lattice.reshape(dim, dim)
    basis = np.asarray(basis, dtype=float).reshape(-1, dim)
    site = np.asarray(site, dtype=float).reshape(dim)

    if radii is None:
        radii = np.full(len(basis), 0.5 * _min_separation(lattice, basis))
    else:
        radii = np.broadcast_to(np.asarray(radii, dtype=float), (len(basis),)).copy()

    # radius of the site itself: the basis image sitting on top of it, if any
    site_radius = None
    for shift in _shifts(lattice, 1):
        d = np.linalg.norm(basis + shift - site, axis=1)
        j = int(np.argmin(d))
        if d[j] < 1e-9:
            site_radius = float(radii[j])
            break
    if site_radius is None:                       # a vacancy site: use the mean
        site_radius = float(radii.mean())

    ids = np.arange(len(basis))
    for shells in range(nshell, nshell + 4):
        shifts = _shifts(lattice, shells)
        pts = (basis[None, :, :] + shifts[:, None, :]).reshape(-1, dim)
        pool_r = np.tile(radii, len(shifts))
        pool_id = np.tile(ids, len(shifts))
        try:
            return cell_from_neighbours(
                pts - site,
                site_radius,
                pool_r,
                pool_id,
                area_tol=area_tol,
                contact_tol=contact_tol,
            )
        except PoolTooSmall:
            continue
    raise PoolTooSmall(
        f"cell did not close within {nshell + 3} periodic shells; the lattice "
        "is probably degenerate"
    )


# --------------------------------------------------------------------------- #
# the admissibility criterion
# --------------------------------------------------------------------------- #
def is_tangential(cells, sigma: float, tol: float = 1e-7) -> bool:
    """True when every non-degenerate face of every cell is supported at ``sigma/2``.

    The admissibility criterion of ``paper.tex`` Sec. II: *every Voronoi
    neighbour is a contact*.
    """
    if isinstance(cells, Cell):
        cells = [cells]
    return all(
        bool(np.all(np.abs(c.face_distances - 0.5 * sigma) <= tol)) for c in cells
    )


# --------------------------------------------------------------------------- #
# structure-level analysis
# --------------------------------------------------------------------------- #
def _packing_fraction(dim: int, radii: np.ndarray, cell_volume: float) -> float:
    if dim == 3:
        return float((4.0 / 3.0) * np.pi * (radii**3).sum() / cell_volume)
    return float(np.pi * (radii**2).sum() / cell_volume)


def analyse(
    name: str,
    lattice,
    basis,
    radii=None,
    *,
    nshell: int = 3,
    tol: float = 1e-7,
    congruence_tol: float = 1e-8,
    sample=None,
) -> dict:
    """Full geometric analysis of one periodic structure (2D or 3D).

    Returns the packing fraction at contact, the per-site cells, the contact
    number, and the two booleans that decide admissibility: ``congruent`` (all
    cells identical in face count, volume and surface) and ``tangential``
    (every face supported at ``sigma/2``).  ``sample`` restricts the analysis
    to a subset of basis indices for large bases; congruence is then a sampled
    statement and is labelled as such in the returned record.

    Asserts the theorem: a congruent tangential structure must satisfy
    ``eta == Q_iso`` to 1e-9.
    """
    lattice = np.asarray(lattice, dtype=float)
    dim = lattice.shape[0]
    lattice = lattice.reshape(dim, dim)
    basis = np.asarray(basis, dtype=float).reshape(-1, dim)
    sigma = _min_separation(lattice, basis)
    if radii is None:
        radii = np.full(len(basis), 0.5 * sigma)
    radii = np.broadcast_to(np.asarray(radii, dtype=float), (len(basis),)).copy()

    Vcell = abs(float(np.linalg.det(lattice)))
    eta = _packing_fraction(dim, radii, Vcell)

    indices = list(range(len(basis))) if sample is None else list(sample)
    cells, congruent = [], True
    ref = None
    for i in indices:
        c = voronoi_cell(
            basis[i], lattice, basis, radii=radii, nshell=nshell, contact_tol=tol
        )
        cells.append(c)
        key = (c.n_faces, c.V, c.S)
        if ref is None:
            ref = key
        elif (
            key[0] != ref[0]
            or abs(key[1] - ref[1]) > congruence_tol
            or abs(key[2] - ref[2]) > congruence_tol
        ):
            congruent = False
            break

    tangential = congruent and is_tangential(cells, sigma, tol)
    Q = np.array([c.Q_iso for c in cells])
    z = np.array([c.n_contacts for c in cells])
    nf = np.array([c.n_faces for c in cells])

    if tangential and sample is None:
        # the theorem; a violation here is a bug in the cell builder
        assert abs(eta - float(Q.mean())) < 1e-9, (
            f"{name}: tangential but eta={eta!r} != Q_iso={Q.mean()!r}"
        )

    return {
        "name": name,
        "dim": dim,
        "sigma": sigma,
        "eta": eta,
        "Q": Q,
        "Q_iso": float(Q.mean()),
        "n_faces": nf,
        "z": z,
        "cells": cells,
        "congruent": congruent,
        "congruence_sampled": sample is not None,
        "tangential": tangential,
        "V_cell": Vcell,
    }


def report(res: dict) -> str:
    """One-line summary of :func:`analyse`, printed and returned."""
    nf, z = res["n_faces"], res["z"]
    line = (
        f"  {res['name']:36s} eta={res['eta']:.7f}  Q_iso={res['Q_iso']:.7f}  "
        f"faces={_uniq(nf):>5s}  z={_uniq(z):>5s}  "
        f"{'TANGENTIAL' if res['tangential'] else 'not tangential'}"
        f"{'' if res['congruent'] else '  (cells not congruent)'}"
        f"{'  [sampled]' if res.get('congruence_sampled') else ''}"
    )
    print(line)
    return line


def _uniq(a: np.ndarray) -> str:
    u = np.unique(a)
    return str(int(u[0])) if len(u) == 1 else "/".join(str(int(v)) for v in u)

"""Reference structures and the exact enumeration of tangential Bravais lattices.

By Voronoi's theorem the faces of a lattice Voronoi cell correspond to the seven
non-zero classes of ``L/2L``, and a class contributes a non-degenerate face
exactly when its minimum norm is attained by the antipodal pair alone.  Writing
the Selling parameters ``p_ij = -b_i . b_j`` of a reduced basis with
``b_4 = -(b_1 + b_2 + b_3)``, the seven class norms are *linear* functionals of
``p``:

    |b_i|^2       = sum_{j != i} p_ij
    |b_1 + b_2|^2 = p_13 + p_14 + p_23 + p_24        (and cyclic)

Tangentiality -- all live class norms equal -- is therefore a **linear system**,
not an optimisation.  Fixing ``sum p_ij = 1`` makes the shape space of all 3D
lattices the compact 5-simplex ``{p >= 0}``, so solving the system over every
live-face set and every zero-pattern of ``p`` is an exhaustive enumeration.

Landmine 3 of the task list, hit for real once: do NOT replace this with a
numerical minimisation.  A Nelder-Mead search over the shape space fails -- the
objective is non-smooth and flat-bottomed, and the simple cubic and simple
hexagonal solutions live on boundary strata where several ``p_ij`` vanish,
which the search walks straight past.

The module also carries the vacancy-shell computations behind the corrections
of ``paper.tex`` Sec. V, which are statements about named structures.
"""

from __future__ import annotations

import itertools

import numpy as np

from .tangential import analyse, voronoi_cell

__all__ = [
    "CLASS_NORMS",
    "FCC_GENERATORS",
    "NAMED",
    "NAMED_2D",
    "enumerate_tangential_lattices",
    "face_class_norms",
    "lattice_from_gram",
    "named_structure",
    "no_fourteen_faced_solution",
    "selling_gram",
    "tetrahedral_octahedral_honeycomb",
    "triangular_lattice",
    "vacancy_shells",
]

_S3 = np.sqrt(3.0)
_HCP_C = np.sqrt(8.0 / 3.0)

#: primitive FCC generators in integer coordinates, nearest-neighbour sqrt(2)
FCC_GENERATORS = np.array([[1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 0.0, 1.0]])

#: the Selling parameter pairs, in the order ``p12 p13 p14 p23 p24 p34``
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

#: the seven face-class norms as linear functionals of ``p``
#: rows: |b1|^2 |b2|^2 |b3|^2 |b4|^2 |b1+b2|^2 |b1+b3|^2 |b2+b3|^2
CLASS_NORMS = np.array(
    [
        [1, 1, 1, 0, 0, 0],
        [1, 0, 0, 1, 1, 0],
        [0, 1, 0, 1, 0, 1],
        [0, 0, 1, 0, 1, 1],
        [0, 1, 1, 1, 1, 0],
        [1, 0, 1, 1, 0, 1],
        [1, 1, 0, 0, 1, 1],
    ],
    dtype=float,
)

#: named 3D structures as ``(lattice_rows, basis)``
NAMED: dict[str, tuple[np.ndarray, np.ndarray]] = {
    "FCC": (
        np.eye(3),
        np.array([[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]], float),
    ),
    "HCP": (
        np.array([[1, 0, 0], [0.5, _S3 / 2, 0], [0, 0, _HCP_C]]),
        np.array([[0, 0, 0], [0.5, 1 / (2 * _S3), _HCP_C / 2]]),
    ),
    "BCC": (np.eye(3), np.array([[0, 0, 0], [0.5, 0.5, 0.5]])),
    "SC": (np.eye(3), np.array([[0.0, 0.0, 0.0]])),
    "simple hexagonal c=a": (
        np.array([[1, 0, 0], [0.5, _S3 / 2, 0], [0, 0, 1.0]]),
        np.array([[0.0, 0.0, 0.0]]),
    ),
    "diamond": (
        np.eye(3),
        np.array(
            [
                [0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0],
                [0.25, 0.25, 0.25], [0.25, 0.75, 0.75],
                [0.75, 0.25, 0.75], [0.75, 0.75, 0.25],
            ],
            float,
        ),
    ),
    "FCC-1in4": (
        np.eye(3),
        np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]], float),
    ),
}

#: the 2D triangular lattice (unit nearest-neighbour distance)
TRIANGULAR_GENERATORS = np.array([[1.0, 0.0], [0.5, _S3 / 2]])

NAMED_2D: dict[str, tuple[np.ndarray, np.ndarray]] = {
    "triangular": (TRIANGULAR_GENERATORS, np.array([[0.0, 0.0]])),
}

_FACE_LABEL = {
    6: "cube            -> simple cubic",
    8: "hexagonal prism -> simple hexagonal, c = a",
    12: "rhombic dodec.  -> FCC",
    14: "truncated octahedron",
}


def named_structure(name: str) -> tuple[np.ndarray, np.ndarray]:
    """Look up a named structure; returns ``(lattice, basis)`` copies."""
    table = NAMED if name in NAMED else NAMED_2D
    lat, bas = table[name]
    return np.array(lat, float), np.array(bas, float)


def triangular_lattice() -> tuple[np.ndarray, np.ndarray]:
    return named_structure("triangular")


# --------------------------------------------------------------------------- #
# Selling parameters
# --------------------------------------------------------------------------- #
def selling_gram(p) -> np.ndarray:
    """Gram matrix of the reduced basis with Selling parameters ``p``.

    ``p = (p12, p13, p14, p23, p24, p34)`` with ``p_ij = -b_i . b_j`` and
    ``b_4 = -(b_1 + b_2 + b_3)``: ``G_ii = sum_{j != i} p_ij``, ``G_ij = -p_ij``.
    """
    p = np.asarray(p, dtype=float).reshape(6)
    P = np.zeros((4, 4))
    for k, (i, j) in enumerate(PAIRS):
        P[i, j] = P[j, i] = p[k]
    G = np.zeros((3, 3))
    for i in range(3):
        G[i, i] = P[i].sum()
        for j in range(3):
            if i != j:
                G[i, j] = -P[i, j]
    return G


def face_class_norms(p) -> np.ndarray:
    """The seven Voronoi face-class squared norms; LINEAR in ``p``."""
    return CLASS_NORMS @ np.asarray(p, dtype=float).reshape(6)


def lattice_from_gram(G: np.ndarray) -> np.ndarray:
    """Lattice vectors (rows) realising a positive-definite Gram matrix."""
    return np.linalg.cholesky(np.asarray(G, dtype=float))


# --------------------------------------------------------------------------- #
# class-norm bookkeeping used by the enumeration
# --------------------------------------------------------------------------- #
_COEF = np.array(list(itertools.product(range(-2, 3), repeat=3)))
_COEF = _COEF[np.any(_COEF != 0, axis=1)]
_CLASS_OF = (_COEF[:, 0] % 2) + 2 * (_COEF[:, 1] % 2) + 4 * (_COEF[:, 2] % 2)
_CLASS_MASK = np.array([_CLASS_OF == c for c in range(1, 8)])


def _tangential_from_selling(p, *, tol: float = 1e-10):
    """``(n_faces, eta)`` if the lattice with Selling parameters ``p`` is tangential.

    A class is live when its minimum norm is attained by the antipodal pair
    alone; the cell is tangential when every live class has the same norm.
    """
    G = selling_gram(p)
    if np.linalg.eigvalsh(G).min() < 1e-11:
        return None
    norms = np.einsum("ij,jk,ik->i", _COEF, G, _COEF)
    live = []
    for mask in _CLASS_MASK:
        v = norms[mask]
        mn = v.min()
        if np.count_nonzero(v < mn * (1.0 + 1e-9)) == 2:
            live.append(mn)
    if len(live) < 3:
        return None
    if (max(live) - min(live)) / min(live) > tol:
        return None
    sigma = np.sqrt(min(live))
    V = np.sqrt(np.linalg.det(G))
    return 2 * len(live), float((np.pi / 6.0) * sigma**3 / V)


def enumerate_tangential_lattices(
    *, seed: int = 3, samples_per_null_dim: int = 200, verify: bool = True
) -> list[dict]:
    """Every tangential Bravais lattice in three dimensions.

    Enumerates all ``C(7,k)`` live-face sets for ``k = 3..7`` crossed with all
    zero-patterns of ``p`` up to three vanishing parameters, solves each linear
    system exactly, samples any non-trivial null space, and -- when ``verify``
    -- confirms every survivor by direct half-space cell construction.
    Returns one record per distinct solution, sorted by packing fraction.
    """
    rng = np.random.default_rng(seed)
    found: dict[tuple[int, float], np.ndarray] = {}
    n_systems = n_candidates = 0

    for nzero in range(4):
        for zeros in itertools.combinations(range(6), nzero):
            for k in range(3, 8):
                for live in itertools.combinations(range(7), k):
                    rows = [
                        CLASS_NORMS[live[i]] - CLASS_NORMS[live[i + 1]]
                        for i in range(k - 1)
                    ]
                    rows.append(np.ones(6))                   # sum p = 1
                    for s in zeros:
                        e = np.zeros(6)
                        e[s] = 1.0
                        rows.append(e)                        # p_s = 0
                    A = np.array(rows)
                    b = np.zeros(len(rows))
                    b[k - 1] = 1.0
                    n_systems += 1

                    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
                    if np.linalg.norm(A @ sol - b) > 1e-9:
                        continue                              # inconsistent
                    _, sv, Vt = np.linalg.svd(A)
                    null = Vt[np.count_nonzero(sv > 1e-9):]

                    cands = [sol]
                    for _ in range(samples_per_null_dim * len(null)):
                        cands.append(sol + null.T @ rng.normal(0.0, 0.5, len(null)))

                    for p in cands:
                        if p.min() < -1e-12 or abs(p.sum() - 1.0) > 1e-9:
                            continue
                        n_candidates += 1
                        res = _tangential_from_selling(np.clip(p, 0.0, None))
                        if res is not None:
                            found.setdefault(
                                (res[0], round(res[1], 9)), np.clip(p, 0.0, None)
                            )

    out = []
    for (n_faces, eta), p in sorted(found.items(), key=lambda kv: kv[0][1]):
        rec = {
            "n_faces": n_faces,
            "eta": eta,
            "selling": p,
            "label": _FACE_LABEL.get(n_faces, "?"),
            "n_systems": n_systems,
            "n_candidates": n_candidates,
        }
        if verify:
            lat = lattice_from_gram(selling_gram(p))
            a = analyse(f"selling{tuple(np.round(p, 6))}", lat, [[0.0, 0.0, 0.0]])
            rec["verified_tangential"] = bool(a["tangential"])
            rec["verified_eta"] = a["eta"]
            rec["verified_Q_iso"] = a["Q_iso"]
            rec["verified_n_faces"] = int(a["n_faces"][0])
            if (
                not a["tangential"]
                or abs(a["eta"] - eta) > 1e-9
                or int(a["n_faces"][0]) != n_faces
            ):
                raise RuntimeError(
                    "class-norm solution disagrees with the explicit cell: "
                    f"analytic (faces={n_faces}, eta={eta:.9f}) vs explicit "
                    f"(faces={a['n_faces'][0]}, eta={a['eta']:.9f}, "
                    f"tangential={a['tangential']})"
                )
        out.append(rec)
    return out


def no_fourteen_faced_solution() -> dict:
    """Exact exclusion of truncated-octahedral (BCC-type) tangential cells.

    Over the integers, each row block of :data:`CLASS_NORMS` sums to
    ``(2,2,2,2,2,2)``: the four ``|b_i|^2`` sum to ``2 sum(p)`` and so do the
    three ``|b_i + b_j|^2``.  All seven equal to ``t`` would give
    ``4t = 3t = 2 sum(p)``, hence ``t = 0`` and a degenerate lattice.  So no
    fourteen-faced tangential Bravais cell exists and the failure of BCC is
    structural.  At the symmetric point ``p_ij = 1/6`` the class norms are
    ``1/2`` (x4) and ``2/3`` (x3), the statement quoted in
    REFERENCE_VALUES.json.
    """
    ints = CLASS_NORMS.astype(np.int64)
    singles = ints[:4].sum(axis=0)
    pairs = ints[4:].sum(axis=0)
    two = np.full(6, 2, dtype=np.int64)
    ok = bool(np.array_equal(singles, two) and np.array_equal(pairs, two))
    p_sym = np.full(6, 1.0 / 6.0)
    return {
        "singles_coefficients": singles,
        "pairs_coefficients": pairs,
        "identities_hold": ok,
        "excluded": ok,
        "p_symmetric": p_sym,
        "class_norms_symmetric": face_class_norms(p_sym),
    }


# --------------------------------------------------------------------------- #
# paper.tex Sec. V: corrections to the original 3D discussion
# --------------------------------------------------------------------------- #
def vacancy_shells(structure: str, *, nsuper: int = 4, nshells: int = 2) -> dict:
    """Cells of the particles neighbouring an isolated vacancy (BCC or FCC).

    Builds an ``nsuper``-cubed supercell of the conventional cell, removes the
    central site, and computes the cell of one representative of each of the
    first ``nshells`` distance shells.  The redistributed vacancy volume
    ``sum_shells count * (V/V0 - 1)`` must equal exactly one ``V0``; the
    bookkeeping identity that fails with the neighbour counts quoted in the
    original work (REFERENCE_VALUES ``de_graaf_corrections``).
    """
    lat, basis = named_structure(structure)
    n = nsuper
    sites = np.array(
        [
            [i, j, k] + b
            for i in range(n)
            for j in range(n)
            for k in range(n)
            for b in basis
        ],
        dtype=float,
    )
    V0 = abs(float(np.linalg.det(lat))) / len(basis)
    sigma = {"BCC": _S3 / 2.0, "FCC": 1.0 / np.sqrt(2.0)}[structure]

    vac = np.full(3, n / 2.0)
    d_all = np.linalg.norm(sites - vac, axis=1)
    part = sites[d_all > 1e-9]
    super_lat = n * np.eye(3)
    radii = np.full(len(part), sigma / 2.0)

    d = np.linalg.norm(part - vac, axis=1)
    shells, total = [], 0.0
    for dist in np.unique(np.round(d, 6))[:nshells]:
        idx = np.flatnonzero(np.abs(d - dist) < 1e-6)
        c = voronoi_cell(part[idx[0]], super_lat, part, radii=radii, nshell=1)
        gain = c.V / V0 - 1.0
        total += len(idx) * gain
        shells.append(
            {
                "d": float(dist),
                "count": int(len(idx)),
                "V_over_V0": c.V / V0,
                "n_faces": c.n_faces,
                "eta_local": float((np.pi / 6.0) * sigma**3 / c.V),
                "Q_iso": c.Q_iso,
            }
        )
    return {
        "structure": structure,
        "V0": V0,
        "sigma": sigma,
        "shells": shells,
        "redistributed_volume": total,
        "balances": abs(total - 1.0) < 1e-9,
    }


def tetrahedral_octahedral_honeycomb() -> dict:
    """The corrected tet-oct honeycomb cell (REFERENCE ``de_graaf_corrections``).

    An octahedron of edge ``a`` has eight adjacent tetrahedra, each shared
    among four octahedra: ``V = sqrt2/3 a^3 + 2 a^3/(6 sqrt2) = a^3/sqrt2``,
    exactly the FCC rhombic dodecahedron.  The original takes two adjacent
    tetrahedra, fills 3/4 of space, and quotes ``eta = 0.537422``; corrected,
    ``eta = 2 sqrt3 pi/27 = 0.403067``.
    """
    a = 1.0
    V_oct = np.sqrt(2.0) / 3.0 * a**3
    V_tet = a**3 / (6.0 * np.sqrt(2.0))
    V_correct = V_oct + 8 * 0.25 * V_tet
    V_original = V_oct + 2 * 0.25 * V_tet
    sigma = 2.0 * a / np.sqrt(6.0)
    v_sphere = (np.pi / 6.0) * sigma**3

    lat, basis = named_structure("FCC")
    rd = voronoi_cell(basis[0], lat, basis)
    V_rd = rd.V * (a * np.sqrt(2.0)) ** 3   # rescale to nn distance a

    return {
        "V_cell_corrected": V_correct,
        "V_cell_original": V_original,
        "fill_fraction_original": V_original / V_correct,
        "V_rhombic_dodecahedron": V_rd,
        "reduces_to_rhombic_dodecahedron": abs(V_correct - V_rd) < 1e-12,
        "eta_corrected": v_sphere / V_correct,
        "eta_original": v_sphere / V_original,
        "eta_corrected_closed_form": 2.0 * np.sqrt(3.0) * np.pi / 27.0,
        "eta_original_closed_form": 8.0 * np.pi / (27.0 * np.sqrt(3.0)),
    }

"""The uniform-K vacancy theory: the ladder DERIVED, and its degeneracy probed.

Counting identity (paper.tex Sec. IV).  Let a contact lattice have coordination
``z0`` and remove a vacancy set that is *independent* (no two vacancies
adjacent) and *uniform*: every occupied site has exactly ``K`` vacant
neighbours.  Double-counting occupied-vacant edges gives ``z0 Nv = K No``,
hence

    eta_K = z0 / (z0 + K) * eta_cp .

Spectral bound.  A uniform-K independent vacancy set is an equitable
2-partition (perfect 2-colouring) of the contact graph with quotient matrix
``[[0, z0], [K, z0-K]]``, eigenvalues ``z0`` and ``-K``; on an infinite lattice
``-K`` must lie in the Bloch spectrum.  FCC band ``[-4, 12]`` gives ``K <= 4``;
triangular band ``[-3, 6]`` gives ``K <= 3`` -- which retracts the
once-claimed 2D ``phi = 0.453450`` rung (it would need ``K = 6``, ``z = 0``).

Independence is *forced*, not assumed: around an isolated vacancy every cell
stays tangential, while a divacancy produces a face-distance spread of exactly
``sqrt(2) - 1`` in the adjacent FCC cells (:func:`divacancy_spread`).

Degeneracy.  The identity needs no periodicity, so *all* uniform-K independent
sets on small tori are enumerated exhaustively (:func:`enumerate_uniform_K`,
a constraint-propagating DFS).  Golden counts (REFERENCE_VALUES
``degeneracy_counts``): triangular 4x4 K=2 -> 16 solutions / 4 orbits /
4 cosets; 6x6 K=3 -> 3/1/1; 6x6 K=2 -> 40/4/1, the non-coset orbits being
exactly tangential at ``Q = sqrt(3) pi / 8`` cell by cell.

Landmine 5, hit for real once: a torus carrying a modular rule must have side
a multiple of the modulus, or uniformity silently breaks.  The explicit
realisations below therefore never go through tori at all -- each modular rule
is converted to its exact vacancy *sublattice* (``depletion.sublattice_from_rule``)
and verified as an infinite periodic structure.
"""

from __future__ import annotations

import itertools

import numpy as np

from .depletion import depleted_structure, sublattice_from_rule
from .lattices import FCC_GENERATORS, TRIANGULAR_GENERATORS, named_structure
from .tangential import analyse, voronoi_cell

__all__ = [
    "MODULAR_RULES_2D",
    "MODULAR_RULES_3D",
    "bloch_band",
    "contact_vectors",
    "counting_eta",
    "divacancy_spread",
    "enumerate_uniform_K",
    "fcc_torus",
    "is_lattice_coset",
    "modular_vacancies",
    "modular_vacancies_2d",
    "realise_modular_K",
    "spectral_K_max",
    "triangular_torus",
    "verify_modular_K",
]

ETA_FCC = np.pi / np.sqrt(18.0)
PHI_TRI = np.pi / np.sqrt(12.0)


# --------------------------------------------------------------------------- #
# the counting identity
# --------------------------------------------------------------------------- #
def counting_eta(z0: int, K: int, eta_cp: float) -> float:
    """``eta_K = z0/(z0+K) * eta_cp`` -- the double-counting identity."""
    return z0 / (z0 + K) * eta_cp


# --------------------------------------------------------------------------- #
# Bloch bands and the spectral bound
# --------------------------------------------------------------------------- #
def contact_vectors(lattice, basis, tol: float = 1e-9) -> list[tuple[int, int, np.ndarray]]:
    """Contact bonds of a periodic structure: ``(a, b, R)`` with ``|R| = sigma``.

    Site ``basis[a]`` touches the image of ``basis[b]`` displaced by the
    Cartesian vector ``R``; every bond appears once per direction.
    """
    lattice = np.asarray(lattice, float)
    dim = lattice.shape[0]
    basis = np.asarray(basis, float).reshape(-1, dim)
    rng = range(-2, 3)
    shifts = np.array(list(itertools.product(rng, repeat=dim)), float) @ lattice
    # contact distance = global minimum separation
    sigma = np.inf
    for a in range(len(basis)):
        d = np.linalg.norm(
            (basis[None, :, :] + shifts[:, None, :]).reshape(-1, dim) - basis[a], axis=1
        )
        d = d[d > tol]
        sigma = min(sigma, float(d.min()))
    bonds = []
    for a in range(len(basis)):
        for b in range(len(basis)):
            for s in shifts:
                R = basis[b] + s - basis[a]
                r = np.linalg.norm(R)
                if tol < r < sigma * (1 + 1e-9):
                    bonds.append((a, b, R.copy()))
    return bonds


def bloch_band(lattice, basis, contacts=None, nq: int = 200_000, seed: int = 0):
    """Sampled Bloch band ``[lo, hi]`` of the contact-graph adjacency operator.

    ``H(q)_{ab} = sum_bonds(a->b) exp(i q . R)``; eigenvalues over ``nq``
    random ``q``.  Band edges at isolated quadratic minima are approached to
    ``O(nq^{-2/dim})``, which is why :func:`spectral_K_max` carries a 0.01
    tolerance.
    """
    lattice = np.asarray(lattice, float)
    dim = lattice.shape[0]
    basis = np.asarray(basis, float).reshape(-1, dim)
    if contacts is None:
        contacts = contact_vectors(lattice, basis)
    m = len(basis)
    rng = np.random.default_rng(seed)
    q = rng.uniform(-2 * np.pi, 2 * np.pi, (nq, dim))
    if m == 1:
        lam = np.zeros(nq)
        for _, _, R in contacts:
            lam += np.cos(q @ R)
        return float(lam.min()), float(lam.max())
    lo, hi = np.inf, -np.inf
    chunk = 20_000
    for start in range(0, nq, chunk):
        qq = q[start : start + chunk]
        H = np.zeros((len(qq), m, m), dtype=complex)
        for a, b, R in contacts:
            H[:, a, b] += np.exp(1j * (qq @ R))
        ev = np.linalg.eigvalsh(H)
        lo = min(lo, float(ev.min()))
        hi = max(hi, float(ev.max()))
    return lo, hi


def spectral_K_max(lattice, basis, contacts=None, nq: int = 200_000) -> int:
    """``floor(-band_min)``: the largest K the spectral condition allows.

    The 0.01 slack absorbs the sampling error at a quadratic band minimum; the
    exact band edges for FCC ``[-4,12]`` and triangular ``[-3,6]`` are integers,
    so the slack cannot change the answer for the gated cases.
    """
    lo, _ = bloch_band(lattice, basis, contacts, nq=nq)
    return int(np.floor(-lo + 0.01))


# --------------------------------------------------------------------------- #
# modular rules and their exact periodic realisations
# --------------------------------------------------------------------------- #
#: FCC vacancy rules on integer coordinates (x+y+z even), REFERENCE_VALUES
MODULAR_RULES_3D = {
    4: {"index": 4, "rule": lambda x, y, z: x % 2 == 0 and y % 2 == 0 and z % 2 == 0,
        "label": "x,y,z all even"},
    3: {"index": 5, "rule": lambda x, y, z: (y + 2 * z) % 5 == 0,
        "label": "y+2z=0 mod 5"},
    2: {"index": 7, "rule": lambda x, y, z: (x + 2 * y + 3 * z) % 7 == 0,
        "label": "x+2y+3z=0 mod 7"},
    1: {"index": 13, "rule": lambda x, y, z: (x + 3 * y + 4 * z) % 13 == 0,
        "label": "x+3y+4z=0 mod 13"},
}

#: triangular vacancy rules on lattice coordinates (i, j)
MODULAR_RULES_2D = {
    3: {"index": 3, "rule": lambda i, j: (i - j) % 3 == 0, "label": "i-j=0 mod 3",
        "name": "honeycomb"},
    2: {"index": 4, "rule": lambda i, j: i % 2 == 0 and j % 2 == 0,
        "label": "i,j both even", "name": "kagome"},
    1: {"index": 7, "rule": lambda i, j: (i + 3 * j) % 7 == 0, "label": "i+3j=0 mod 7",
        "name": "maple-leaf"},
}


def modular_vacancies(K: int):
    """The FCC vacancy rule for ``K`` deleted neighbours per particle."""
    return MODULAR_RULES_3D[K]["rule"]


def modular_vacancies_2d(K: int):
    """The triangular vacancy rule for ``K`` in 2D."""
    return MODULAR_RULES_2D[K]["rule"]


def realise_modular_K(K: int, dim: int = 3):
    """Exact periodic structure of the modular rule: ``(lattice, basis, H)``.

    The rule's kernel is a sublattice of the parent contact lattice; the
    depleted structure (parent minus that sublattice) is returned in the
    parent's integer coordinates.  No torus is involved, so landmine 5 cannot
    bite.
    """
    if dim == 3:
        spec, gen, rc = MODULAR_RULES_3D[K], FCC_GENERATORS, "cartesian"
    else:
        spec, gen, rc = MODULAR_RULES_2D[K], TRIANGULAR_GENERATORS, "lattice"
    H = sublattice_from_rule(spec["rule"], spec["index"], gen, rule_coords=rc)
    lat, basis = depleted_structure(H, gen)
    return lat, basis, H


def verify_modular_K(K: int, dim: int = 3) -> dict:
    """Verify one modular rule end to end: uniform-K, independent, tangential,
    congruent, correct p-vector, exact density.

    * independence: the shortest vacancy-sublattice vector exceeds the contact
      distance;
    * uniformity: every occupied coset representative has exactly ``K``
      nearest neighbours on the vacancy sublattice;
    * tangential + congruent: full cell analysis of the depleted structure;
    * density: against the counting identity ``z0/(z0+K) eta_cp``.
    """
    z0 = 12 if dim == 3 else 6
    eta_cp = ETA_FCC if dim == 3 else PHI_TRI
    gen = FCC_GENERATORS if dim == 3 else TRIANGULAR_GENERATORS
    sigma = np.sqrt(2.0) if dim == 3 else 1.0

    lat, basis, H = realise_modular_K(K, dim)

    # independence: shortest non-zero vacancy-lattice vector
    rng = range(-3, 4)
    combos = np.array(
        [c for c in itertools.product(rng, repeat=dim) if any(c)], dtype=float
    )
    vac_vecs = combos @ (H.astype(float) @ gen)
    shortest = float(np.linalg.norm(vac_vecs, axis=1).min())
    independent = shortest > sigma * (1 + 1e-9)

    # uniformity: K vacancy neighbours for every occupied representative
    nn = _nn_offsets(dim) @ gen if False else _nn_cartesian(dim)
    Linv = np.linalg.inv(lat)
    kcounts = []
    for b in basis:
        c = 0
        for v in nn:
            frac = (b + v) @ Linv
            # vacancy sites are the sublattice points: frac integer <=> on lattice
            if np.max(np.abs(frac - np.round(frac))) < 1e-9:
                c += 1
        kcounts.append(c)
    uniform = set(kcounts) == {K}

    res = analyse(f"modular K={K} ({dim}D)", lat, basis)
    pvec = sorted(len(loop) for loop in res["cells"][0].face_vertices)  # = topology.p_vector

    eta_exact = counting_eta(z0, K, eta_cp)
    return {
        "K": K,
        "dim": dim,
        "index": (MODULAR_RULES_3D if dim == 3 else MODULAR_RULES_2D)[K]["index"],
        "label": (MODULAR_RULES_3D if dim == 3 else MODULAR_RULES_2D)[K]["label"],
        "uniform": bool(uniform),
        "independent": bool(independent),
        "tangential": bool(res["tangential"]),
        "congruent": bool(res["congruent"]),
        "z": int(np.unique(res["z"])[0]),
        "n_faces": int(np.unique(res["n_faces"])[0]),
        "p_vector": pvec,
        "eta": res["eta"],
        "eta_counting": eta_exact,
        "eta_err": abs(res["eta"] - eta_exact),
        "kcounts": sorted(set(kcounts)),
    }


def _nn_cartesian(dim: int) -> np.ndarray:
    """Nearest-neighbour Cartesian vectors of the parent contact lattice."""
    if dim == 3:
        out = [
            v
            for v in itertools.product((-1, 0, 1), repeat=3)
            if sorted(map(abs, v)) == [0, 1, 1]
        ]
        return np.array(out, dtype=float)
    return np.array(
        [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)], dtype=float
    ) @ TRIANGULAR_GENERATORS


# --------------------------------------------------------------------------- #
# independence is forced: the divacancy counterexample
# --------------------------------------------------------------------------- #
def divacancy_spread(side: int = 8) -> dict:
    """Face-distance spread around a divacancy vs an isolated vacancy in FCC.

    Returns the maximum relative spread ``(h_max - h_min)/h_min`` over the
    cells adjacent to the defect.  Golden values: isolated vacancy -> 0
    exactly; divacancy -> ``sqrt(2) - 1 = 0.41421`` (a face supported at
    distance 1 against contact faces at ``1/sqrt(2)``).
    """
    assert side % 2 == 0, "FCC integer torus needs an even side"
    sites = np.array(
        [
            (x, y, z)
            for x in range(side)
            for y in range(side)
            for z in range(side)
            if (x + y + z) % 2 == 0
        ],
        dtype=float,
    )
    lat = side * np.eye(3)
    v1 = np.full(3, side // 2, dtype=float)
    if int(v1.sum()) % 2:
        v1[0] += 1.0
    v2 = v1 + np.array([1.0, 1.0, 0.0])

    def spread_with(vacancies):
        occ = np.array(
            [p for p in sites if not any(np.all(np.abs(p - v) < 1e-9) for v in vacancies)]
        )
        worst = 0.0
        for p in occ:
            if min(np.linalg.norm(p - v) for v in vacancies) < 1.5:
                c = voronoi_cell(p, lat, occ, nshell=1)
                worst = max(
                    worst,
                    float(
                        (c.face_distances.max() - c.face_distances.min())
                        / c.face_distances.min()
                    ),
                )
        return worst

    return {
        "isolated": spread_with([v1]),
        "divacancy": spread_with([v1, v2]),
        "exact": np.sqrt(2.0) - 1.0,
    }


# --------------------------------------------------------------------------- #
# tori and the exhaustive uniform-K enumeration
# --------------------------------------------------------------------------- #
class Torus:
    """A finite contact graph with its translation group.

    ``coords[i]`` are the integer lattice coordinates of site ``i``;
    ``nbr[i]`` the neighbour indices (with multiplicity, which tiny tori
    need); ``add(i, t)`` the site index of site ``i`` translated by the
    coordinates of site ``t``; ``cartesian[i]`` the embedding used for cell
    checks.
    """

    def __init__(self, coords, nbr, add_table, cartesian, lattice, z0, label):
        self.coords = coords
        self.nbr = nbr
        self.add_table = add_table          # (N, N) int array
        self.cartesian = cartesian
        self.lattice = lattice
        self.z0 = z0
        self.label = label
        self.N = len(coords)


def triangular_torus(N: int) -> Torus:
    """The N x N triangular torus (6 neighbours per site)."""
    coords = [(i, j) for i in range(N) for j in range(N)]
    pos = {p: k for k, p in enumerate(coords)}
    offs = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
    nbr = [
        [pos[((i + a) % N, (j + b) % N)] for a, b in offs] for i, j in coords
    ]
    add = np.empty((len(coords), len(coords)), dtype=np.int32)
    for s, (i, j) in enumerate(coords):
        for t, (a, b) in enumerate(coords):
            add[s, t] = pos[((i + a) % N, (j + b) % N)]
    cart = np.array([i * TRIANGULAR_GENERATORS[0] + j * TRIANGULAR_GENERATORS[1]
                     for i, j in coords])
    lat = np.array([N * TRIANGULAR_GENERATORS[0], N * TRIANGULAR_GENERATORS[1]])
    return Torus(coords, nbr, add, cart, lat, 6, f"triangular {N}x{N}")


def fcc_torus(side: int) -> Torus:
    """The FCC torus of even integer side (12 neighbours per site)."""
    assert side % 2 == 0, "FCC integer torus needs an even side (parity wraps)"
    coords = [
        (x, y, z)
        for x in range(side)
        for y in range(side)
        for z in range(side)
        if (x + y + z) % 2 == 0
    ]
    pos = {p: k for k, p in enumerate(coords)}
    offs = [
        v for v in itertools.product((-1, 0, 1), repeat=3)
        if sorted(map(abs, v)) == [0, 1, 1]
    ]
    nbr = [
        [pos[((x + a) % side, (y + b) % side, (z + c) % side)] for a, b, c in offs]
        for x, y, z in coords
    ]
    add = np.empty((len(coords), len(coords)), dtype=np.int32)
    for s, p in enumerate(coords):
        for t, q in enumerate(coords):
            add[s, t] = pos[tuple((p[k] + q[k]) % side for k in range(3))]
    cart = np.array(coords, dtype=float)
    lat = side * np.eye(3)
    return Torus(coords, nbr, add, cart, lat, 12, f"FCC side {side}")


def enumerate_uniform_K(
    torus: Torus,
    K: int,
    *,
    max_solutions: int | None = None,
    max_nodes: int | None = None,
) -> dict:
    """Exhaustively enumerate uniform-K independent vacancy sets on a torus.

    Constraint-propagating DFS over site states (occupied / vacant):
    vacancies pairwise non-adjacent; every occupied site ends with exactly
    ``K`` vacant neighbours; the total vacancy count is pinned to
    ``Nv = K N / (z0 + K)`` (a strong prune the reference prototype lacked).

    Returns solutions, translation-orbit count, and how many orbits are
    lattice cosets.  ``max_solutions`` / ``max_nodes`` cap runaway
    enumerations; if a cap is hit the result is flagged ``complete=False``
    and MUST NOT be quoted as a count.
    """
    N, z0 = torus.N, torus.z0
    if (K * N) % (z0 + K) != 0:
        return {
            "torus": torus.label, "K": K, "N": N, "Nv": None,
            "solutions": 0, "orbits": 0, "cosets": 0,
            "complete": True, "nodes": 0,
            "note": "Nv = K N/(z0+K) is not an integer; no solutions exist",
        }
    Nv = (K * N) // (z0 + K)
    nbr = torus.nbr
    state = np.full(N, -1, dtype=np.int8)      # -1 undecided, 0 occupied, 1 vacant
    vdec = np.zeros(N, dtype=np.int16)         # decided-vacant neighbours
    undec = np.array([len(nbr[i]) for i in range(N)], dtype=np.int16)
    sols: list[tuple] = []
    nodes = 0
    vac_used = 0
    aborted = False

    def feasible(i: int) -> bool:
        """Can site i (occupied, or undecided) still reach exactly K vacant nbrs?"""
        return vdec[i] <= K and vdec[i] + undec[i] >= K

    def dfs(s: int, vac_used: int):
        nonlocal nodes, aborted
        if aborted:
            return
        nodes += 1
        if max_nodes is not None and nodes > max_nodes:
            aborted = True
            return
        if s == N:
            if vac_used == Nv:
                sols.append(tuple(state))
                if max_solutions is not None and len(sols) >= max_solutions:
                    aborted = True
            return
        # count prune: enough undecided sites left to reach Nv?
        remaining = N - s
        if vac_used > Nv or vac_used + remaining < Nv:
            return
        for val in (0, 1):
            if val == 1:
                if vac_used == Nv:
                    continue
                if any(state[t] == 1 for t in nbr[s]):
                    continue                       # independence
            state[s] = val
            ok = True
            if val == 1:
                for t in nbr[s]:
                    vdec[t] += 1
            for t in nbr[s]:
                undec[t] -= 1
            # local consistency: s and its neighbours
            if state[s] == 0 and not feasible(s):
                ok = False
            if ok:
                for t in nbr[s]:
                    if state[t] == 0 and not feasible(t):
                        ok = False
                        break
                    if state[t] == -1 and val == 1:
                        pass
            if ok:
                dfs(s + 1, vac_used + (1 if val == 1 else 0))
            # undo
            for t in nbr[s]:
                undec[t] += 1
            if val == 1:
                for t in nbr[s]:
                    vdec[t] -= 1
            state[s] = -1
            if aborted:
                return

    dfs(0, 0)
    complete = not aborted

    # translation orbits
    add = torus.add_table
    orbit_reps: set[tuple] = set()
    for sol in sols:
        arr = np.array(sol, dtype=np.int8)
        best = min(tuple(arr[add[:, t]]) for t in range(N))
        orbit_reps.add(best)
    # coset test per orbit
    cosets = 0
    for rep in orbit_reps:
        vac = [i for i, v in enumerate(rep) if v == 1]
        if is_lattice_coset(vac, torus):
            cosets += 1
    return {
        "torus": torus.label, "K": K, "N": N, "Nv": Nv,
        "solutions": len(sols), "orbits": len(orbit_reps), "cosets": cosets,
        "orbit_representatives": sorted(orbit_reps),
        "complete": complete, "nodes": nodes,
    }


def is_lattice_coset(vacancy_indices, torus: Torus) -> bool:
    """Is the vacancy set a coset of a subgroup of the torus translations?

    Shift the set so one element sits at the origin; a coset is exactly a
    subgroup, i.e. closed under the torus addition.
    """
    if not len(vacancy_indices):
        return False
    add = torus.add_table
    coords = torus.coords
    N = torus.N
    pos = {p: k for k, p in enumerate(coords)}
    v0 = coords[vacancy_indices[0]]
    dim = len(v0)
    if dim == 2:
        n = int(round(torus.lattice[0, 0]))  # triangular: lattice = N * generators
        n = len({c[0] for c in coords})      # robust: side length
        shifted = {tuple((c[k] - v0[k]) % n for k in range(dim))
                   for c in (coords[i] for i in vacancy_indices)}
        return all(
            tuple((a[k] + b[k]) % n for k in range(dim)) in shifted
            for a in shifted for b in shifted
        )
    side = len({c[0] for c in coords})
    shifted = {tuple((c[k] - v0[k]) % side for k in range(dim))
               for c in (coords[i] for i in vacancy_indices)}
    return all(
        tuple((a[k] + b[k]) % side for k in range(dim)) in shifted
        for a in shifted for b in shifted
    )


def solution_cells(rep, torus: Torus, sample=None) -> list:
    """Voronoi cells of the occupied sites of one enumerated solution.

    Used to check the degeneracy claim cell by cell: every cell of every
    uniform-K solution -- crystalline or not -- is exactly tangential with the
    same ``Q``.
    """
    occ = np.array([torus.cartesian[i] for i, v in enumerate(rep) if v == 0])
    idx = range(len(occ)) if sample is None else sample
    return [voronoi_cell(occ[i], torus.lattice, occ, nshell=1) for i in idx]

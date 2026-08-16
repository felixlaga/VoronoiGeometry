"""The depleted-FCC ladder: FCC minus a periodic vacancy sublattice.

Removing a sublattice of FCC of index ``k`` and keeping the other ``k-1``
cosets leaves congruent tangential cells for exactly ``k in {4, 5, 7, 13}``:

    eta_k = (1 - 1/k) pi/sqrt(18),   z_k = 12 - 12/(k-1).

In the coloring language of ``geometry/coloring.py`` these are the uniform-K
states with ``K = 12/(k-1) in {4, 3, 2, 1}``; the scan here *observes* the
ladder from Hermite-normal-form enumeration alone, and gate G0b then checks
that the derivation (counting identity + spectral bound) reproduces it, so the
two routes stay mutually validating.

Landmines 1 and 2 of the task list, both hit for real once:

* HNF off-diagonal ranges are ``0 <= d, e < a`` and ``0 <= f < b``.  Wrong
  ranges silently enumerate an incomplete set of sublattices.
* Coset representatives **must be wrapped into the primitive cell** before any
  periodic cell computation.  Unwrapped representatives sit far outside the
  supercell, their periodic neighbourhoods are truncated by the finite image
  shell, and the scan silently returns "not tangential" for structures that
  are.
"""

from __future__ import annotations

import itertools

import numpy as np

from .lattices import FCC_GENERATORS, TRIANGULAR_GENERATORS
from .tangential import PoolTooSmall, analyse

__all__ = [
    "ETA_FCC",
    "PHI_TRIANGULAR",
    "coset_representatives",
    "depleted_fcc",
    "depleted_structure",
    "hermite_normal_forms",
    "hermite_normal_forms_2d",
    "ladder_closed_form",
    "scan_depletions",
    "sublattice_from_rule",
]

ETA_FCC = np.pi / np.sqrt(18.0)
PHI_TRIANGULAR = np.pi / np.sqrt(12.0)


def hermite_normal_forms(k: int) -> list[np.ndarray]:
    """All index-``k`` sublattices of ``Z^3`` in lower-triangular Hermite normal form.

    ``H = [[a,0,0],[d,b,0],[e,f,c]]`` with ``a b c = k``, ``0 <= d, e < a``,
    ``0 <= f < b``.  Each index-``k`` sublattice has exactly one such ``H``, so
    the enumeration is complete and duplicate-free.
    """
    out = []
    for a in range(1, k + 1):
        if k % a:
            continue
        for b in range(1, k // a + 1):
            if (k // a) % b:
                continue
            c = k // (a * b)
            for d in range(a):
                for e in range(a):
                    for f in range(b):
                        out.append(
                            np.array([[a, 0, 0], [d, b, 0], [e, f, c]], dtype=int)
                        )
    return out


def hermite_normal_forms_2d(k: int) -> list[np.ndarray]:
    """All index-``k`` sublattices of ``Z^2``: ``[[a,0],[d,b]]``, ``ab=k``, ``0<=d<a``."""
    out = []
    for a in range(1, k + 1):
        if k % a:
            continue
        b = k // a
        for d in range(a):
            out.append(np.array([[a, 0], [d, b]], dtype=int))
    return out


def coset_representatives(H: np.ndarray, generators: np.ndarray) -> np.ndarray:
    """Cartesian representatives of the ``k`` cosets of ``H`` in the parent lattice.

    The trivial coset (the vacancy orbit) comes first, at the origin.  Every
    representative is wrapped into the primitive cell of the sublattice
    (landmine 2).
    """
    H = np.asarray(H, dtype=np.int64)
    dim = H.shape[0]
    k = int(round(abs(np.linalg.det(H))))
    lat = H.astype(float) @ generators
    Linv = np.linalg.inv(lat)
    adj = np.round(np.linalg.det(H) * np.linalg.inv(H)).astype(np.int64)

    reps: list[np.ndarray] = []
    seen: set[tuple[int, ...]] = set()
    span = range(-k, k + 1)
    for n in itertools.chain(
        [tuple(0 for _ in range(dim))], itertools.product(span, repeat=dim)
    ):
        key = tuple(int(v) % k for v in (np.asarray(n, dtype=np.int64) @ adj))
        if key in seen:
            continue
        seen.add(key)
        x = np.asarray(n, dtype=float) @ generators
        frac = x @ Linv
        reps.append((frac - np.floor(frac + 1e-9)) @ lat)
        if len(reps) == k:
            break
    if len(reps) != k:
        raise RuntimeError(f"found {len(reps)} of {k} cosets; widen the search box")
    return np.array(reps)


def depleted_structure(H: np.ndarray, generators: np.ndarray):
    """``(lattice, basis)`` for "parent lattice minus the ``H`` sublattice orbit".

    Returns ``None`` when the coset layout is degenerate (duplicate
    representatives).
    """
    lat = np.asarray(H, dtype=float) @ generators
    reps = coset_representatives(H, generators)
    basis = reps[1:]
    if len(basis) > 1:
        d = np.linalg.norm(basis[:, None, :] - basis[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        if d.min() < 1e-8:
            return None
    return lat, basis


def scan_depletions(
    kmax: int = 20, *, generators: np.ndarray | None = None, dim: int = 3
) -> list[dict]:
    """Every vacancy sublattice of index ``k <= kmax`` with congruent tangential cells.

    3D scans FCC; ``dim=2`` scans the triangular lattice.  Returns one record
    per surviving ``k`` with packing fraction, contact number, face count and
    the number of distinct sublattices (HNFs) realising it.
    """
    if generators is None:
        generators = FCC_GENERATORS if dim == 3 else TRIANGULAR_GENERATORS
    hnfs = hermite_normal_forms if dim == 3 else hermite_normal_forms_2d
    hits: dict[int, dict] = {}
    eta_cp = ETA_FCC if dim == 3 else PHI_TRIANGULAR
    for k in range(2, kmax + 1):
        for H in hnfs(k):
            cand = depleted_structure(H, generators)
            if cand is None:
                continue
            lat, basis = cand
            try:
                res = analyse(f"dep-1in{k}", lat, basis)
            except (PoolTooSmall, RuntimeError, ValueError):
                continue
            if not (res["tangential"] and res["congruent"]):
                continue
            rec = hits.setdefault(
                k,
                {
                    "k": k,
                    "eta": res["eta"],
                    "Q_iso": res["Q_iso"],
                    "z": int(np.unique(res["z"])[0]),
                    "n_faces": int(np.unique(res["n_faces"])[0]),
                    "eta_closed_form": (1.0 - 1.0 / k) * eta_cp,
                    "n_sublattices": 0,
                    "hnf": H.copy(),
                },
            )
            rec["n_sublattices"] += 1
    return [hits[k] for k in sorted(hits)]


def depleted_fcc(k: int):
    """``(lattice, basis)`` of a tangential FCC depletion of index ``k``.

    Raises :class:`ValueError` for every ``k`` outside the ladder.
    """
    for H in hermite_normal_forms(k):
        cand = depleted_structure(H, FCC_GENERATORS)
        if cand is None:
            continue
        lat, basis = cand
        try:
            res = analyse(f"FCC-1in{k}", lat, basis)
        except (PoolTooSmall, RuntimeError, ValueError):
            continue
        if res["tangential"] and res["congruent"]:
            return lat, basis
    raise ValueError(f"no tangential FCC depletion of index k={k}")


def ladder_closed_form(k: int) -> dict:
    """Closed forms of the 3D ladder in the depletion index ``k`` (``K = 12/(k-1)``)."""
    return {
        "k": k,
        "K": 12 / (k - 1),
        "eta": (1.0 - 1.0 / k) * ETA_FCC,
        "z": 12 - 12 // (k - 1) if 12 % (k - 1) == 0 else None,
        "divides": 12 % (k - 1) == 0,
    }


# --------------------------------------------------------------------------- #
# modular rules -> explicit sublattices
# --------------------------------------------------------------------------- #
def sublattice_from_rule(
    rule, index: int, generators: np.ndarray, *, rule_coords: str = "cartesian"
) -> np.ndarray:
    """The HNF (w.r.t. ``generators``) of the vacancy sublattice defined by a rule.

    ``rule(coords) -> bool`` must be a linear congruence (true = vacancy); the
    vacancy set is then a sublattice of the parent, of the given ``index``.
    The unique HNF whose generators all satisfy the rule is returned:
    linearity makes the generator check sufficient, and HNF uniqueness makes
    the match unique.

    ``rule_coords`` names the coordinates the rule is written in --
    ``"cartesian"`` for the 3D FCC rules, which act on the integer Cartesian
    triples ``(x, y, z)``, and ``"lattice"`` for the 2D triangular rules,
    which act on the lattice indices ``(i, j)``.  Mixing these up silently
    matches nothing (or the wrong sublattice), so the caller must say which.
    """
    dim = generators.shape[0]
    hnfs = hermite_normal_forms(index) if dim == 3 else hermite_normal_forms_2d(index)
    matches = []
    for H in hnfs:
        if rule_coords == "cartesian":
            vecs = H.astype(float) @ generators
        elif rule_coords == "lattice":
            vecs = H.astype(float)
        else:
            raise ValueError(f"unknown rule_coords {rule_coords!r}")
        if all(rule(*np.round(v).astype(int)) for v in vecs):
            matches.append(H)
    if len(matches) != 1:
        raise ValueError(
            f"rule matched {len(matches)} sublattices of index {index}; "
            "expected exactly one"
        )
    return matches[0]

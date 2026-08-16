"""Per-particle distance to each reference cell (spec 2.5c) -- PRE-REGISTERED.

    s_k(i) = w_Q ((Q_i - eta_k)/eta_k)^2 + w_p d_p(pvec_i, pvec_k) + w_t eps_t(i)^2

with the tangentiality defect ``eps_t^2 = mean_j ((l_ij - R_i)/R_i)^2`` over
kept faces.  Weights, the face filter, the persistence grids and the marking
threshold live in ``refscore_frozen.json``, committed in the repository's
FIRST commit -- before any dynamical dataset existed.  Gate G5 verifies that
the file's git blob is unchanged since that commit.  Nothing in this module
may read a tunable from anywhere else, and nothing downstream may override
the frozen values (anti-pattern list: "tuning refscore weights after seeing
dynamical data").

Reference descriptors (``eta_k``, reference p-vectors) are computed at runtime
from the exact geometry modules; only the weights and thresholds are frozen.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np

from ..geometry.coloring import realise_modular_K
from ..geometry.lattices import named_structure
from ..geometry.tangential import Cell, voronoi_cell
from .topology import face_filter, p_vector

__all__ = [
    "frozen",
    "load_frozen",
    "mark_cells",
    "p_distance",
    "reference_descriptors",
    "refscores",
    "tangentiality_defect",
]

_REPO = Path(__file__).resolve().parents[3]


@lru_cache(maxsize=1)
def load_frozen() -> dict:
    return json.loads((_REPO / "refscore_frozen.json").read_text())


def frozen(key: str):
    return load_frozen()[key]


def p_distance(a, b) -> float:
    """``d_p(a,b) = sum |a'-b'| / (sum a' + sum b')`` on zero-padded sorted p-vectors.

    The exact formula pre-registered in ``refscore_frozen.json``: 0 means
    identical topology, 1 maximal mismatch.
    """
    a = sorted(a)
    b = sorted(b)
    n = max(len(a), len(b))
    a = [0] * (n - len(a)) + a
    b = [0] * (n - len(b)) + b
    tot = sum(a) + sum(b)
    if tot == 0:
        return 0.0
    return float(sum(abs(x - y) for x, y in zip(a, b)) / tot)


def tangentiality_defect(cell: Cell, R_i: float) -> float:
    """``eps_t(i)``: rms of ``(l_ij - R_i)/R_i`` over the (kept) faces."""
    return float(np.sqrt(np.mean(((cell.face_distances - R_i) / R_i) ** 2)))


@lru_cache(maxsize=4)
def reference_descriptors(dim: int = 3) -> dict:
    """``{name: {"eta": ..., "p_vector": ...}}`` from the exact geometry.

    The name list is the frozen one; the numbers are computed, not typed.
    """
    names = frozen("references_3d") if dim == 3 else frozen("references_2d")
    out = {}
    for name in names:
        if name.startswith("K") and dim == 3:
            K = int(name[1])
            lat, basis, _ = realise_modular_K(K, dim=3)
            cell = voronoi_cell(basis[0], lat, basis)
            eta = 12 / (12 + K) * np.pi / np.sqrt(18.0)
        elif dim == 2:
            K = int(name[1])
            lat, basis, _ = realise_modular_K(K, dim=2)
            cell = voronoi_cell(basis[0], lat, basis)
            eta = 6 / (6 + K) * np.pi / np.sqrt(12.0)
        else:
            key = {"simple_cubic": "SC", "simple_hexagonal_ca": "simple hexagonal c=a",
                   "fcc": "FCC"}[name]
            lat, basis = named_structure(key)
            cell = voronoi_cell(basis[0], lat, basis)
            eta = float((np.pi / 6.0) * (2 * cell.r_in) ** 3
                        / (abs(np.linalg.det(lat)) / len(basis)))
        out[name] = {"eta": float(eta), "p_vector": p_vector(cell)}
    return out


def refscores(cells, radii, dim: int = 3) -> tuple[np.ndarray, list[str]]:
    """``s_k`` for every cell against every reference: array ``(n_cells, n_refs)``.

    Cells are passed through the frozen face filter first; ``radii[i]`` is the
    particle's own radius (needed for the tangentiality defect).
    """
    fz = load_frozen()
    w = fz["weights"]
    ff = fz["face_filter"]
    refs = reference_descriptors(dim)
    names = list(refs)
    out = np.empty((len(cells), len(names)))
    for i, cell in enumerate(cells):
        c = face_filter(cell, ff["g_cut"], ff["a_cut"])
        pv = p_vector(c)
        et2 = tangentiality_defect(c, float(radii[i])) ** 2
        for k, name in enumerate(names):
            ref = refs[name]
            out[i, k] = (
                w["w_Q"] * ((c.Q_iso - ref["eta"]) / ref["eta"]) ** 2
                + w["w_p"] * p_distance(pv, ref["p_vector"])
                + w["w_t"] * et2
            )
    return out, names


def mark_cells(scores: np.ndarray, names: list[str], k_name: str) -> np.ndarray:
    """Boolean mask: cells marked as ``k``-like, ``s_k < mark_threshold`` (frozen)."""
    thr = frozen("mark_threshold")
    return scores[:, names.index(k_name)] < thr

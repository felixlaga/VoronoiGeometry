"""The depleted-FCC ladder, and the two landmines the task list names."""

import numpy as np
import pytest

from hsga.geometry.depletion import (
    ETA_FCC,
    coset_representatives,
    depleted_fcc,
    hermite_normal_forms,
    hermite_normal_forms_2d,
    ladder_closed_form,
    scan_depletions,
    sublattice_from_rule,
)
from hsga.geometry.lattices import FCC_GENERATORS, TRIANGULAR_GENERATORS
from hsga.geometry.tangential import analyse

LADDER_K_TO_INDEX = {4: 4, 3: 5, 2: 7, 1: 13}


def test_hermite_normal_form_ranges_and_completeness():
    """Landmine 1: 0 <= d,e < a and 0 <= f < b; wrong ranges under-enumerate."""
    for k in range(2, 13):
        hnfs = hermite_normal_forms(k)
        assert all(round(abs(np.linalg.det(H))) == k for H in hnfs)
        for H in hnfs:
            a, b = H[0, 0], H[1, 1]
            assert 0 <= H[1, 0] < a
            assert 0 <= H[2, 0] < a
            assert 0 <= H[2, 1] < b
        assert len({H.tobytes() for H in hnfs}) == len(hnfs)
        expected = sum(
            a * a * b
            for a in range(1, k + 1) if k % a == 0
            for b in range(1, k // a + 1) if (k // a) % b == 0
        )
        assert len(hnfs) == expected


def test_hermite_normal_forms_2d():
    for k in range(2, 10):
        hnfs = hermite_normal_forms_2d(k)
        assert all(round(abs(np.linalg.det(H))) == k for H in hnfs)
        # sigma_1(k) sublattices of index k in Z^2
        assert len(hnfs) == sum(d for d in range(1, k + 1) if k % d == 0)


def test_coset_representatives_are_wrapped():
    """Landmine 2: unwrapped representatives silently break the analysis."""
    H = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]])
    lat = H.astype(float) @ FCC_GENERATORS
    reps = coset_representatives(H, FCC_GENERATORS)
    assert len(reps) == 4
    assert np.allclose(reps[0], 0.0)
    frac = reps @ np.linalg.inv(lat)
    assert np.all(frac > -1e-9) and np.all(frac < 1 - 1e-9)


def test_unwrapped_representative_changes_the_verdict():
    """The documented failure mode must actually be a failure mode."""
    lat, basis = depleted_fcc(4)
    good = analyse("wrapped", lat, basis)
    assert good["tangential"] and good["congruent"]
    bad_basis = basis.copy()
    bad_basis[0] = bad_basis[0] + 9.0 * lat[0] + 7.0 * lat[1] - 8.0 * lat[2]
    try:
        bad = analyse("unwrapped", lat, bad_basis)
    except Exception:
        return                                   # failing loudly is also correct
    assert not (bad["tangential"] and bad["congruent"])


@pytest.mark.parametrize("k", [4, 5, 7, 13])
def test_ladder_members(k, golden):
    lat, basis = depleted_fcc(k)
    res = analyse(f"FCC-1in{k}", lat, basis)
    assert res["tangential"] and res["congruent"]
    K = 12 // (k - 1)
    ref = golden["fcc_depletion_ladder"][f"K{K}"]
    assert res["eta"] == pytest.approx(ref["eta"], abs=1e-12)
    assert res["Q_iso"] == pytest.approx(ref["eta"], abs=1e-9)
    assert int(np.unique(res["z"])[0]) == ref["z"]
    assert len(basis) == k - 1
    # p-vector against the golden block
    pvec = sorted(len(loop) for loop in res["cells"][0].face_vertices)
    assert pvec == ref["p_vector"]


@pytest.mark.parametrize("k", [2, 3, 6, 8, 9, 10, 11, 12])
def test_non_ladder_indices_have_no_tangential_depletion(k):
    with pytest.raises(ValueError):
        depleted_fcc(k)


@pytest.mark.slow
def test_scan_finds_exactly_the_ladder(golden):
    hits = scan_depletions(13)
    assert {h["k"] for h in hits} == {4, 5, 7, 13}
    for h in hits:
        assert h["eta"] == pytest.approx(h["eta_closed_form"], abs=1e-9)
        assert h["n_faces"] == h["z"]             # tangential: faces == contacts


def test_ladder_closed_forms(golden):
    for k, K in ((4, 4), (5, 3), (7, 2), (13, 1)):
        ref = golden["fcc_depletion_ladder"][f"K{K}"]
        got = ladder_closed_form(k)
        assert got["eta"] == pytest.approx(ref["eta"], abs=1e-12)
        assert got["K"] == K
        assert got["z"] == ref["z"]


def test_sublattice_from_rule_matches_known_hnf():
    H = sublattice_from_rule(
        lambda x, y, z: x % 2 == 0 and y % 2 == 0 and z % 2 == 0, 4, FCC_GENERATORS
    )
    # generators of H @ FCC_G must all be all-even
    vecs = np.round(H.astype(float) @ FCC_GENERATORS).astype(int)
    assert np.all(vecs % 2 == 0)
    with pytest.raises(ValueError):
        sublattice_from_rule(lambda x, y, z: False, 4, FCC_GENERATORS)


def test_sublattice_from_rule_2d_lattice_coords():
    H = sublattice_from_rule(
        lambda i, j: (i - j) % 3 == 0, 3, TRIANGULAR_GENERATORS, rule_coords="lattice"
    )
    assert np.all((H[:, 0] - H[:, 1]) % 3 == 0)

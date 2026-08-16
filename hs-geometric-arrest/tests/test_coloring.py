"""The derivation: counting identity, spectral bound, modular rules, degeneracy."""

import numpy as np
import pytest

from hsga.geometry.coloring import (
    ETA_FCC,
    PHI_TRI,
    counting_eta,
    divacancy_spread,
    enumerate_uniform_K,
    fcc_torus,
    is_lattice_coset,
    modular_vacancies,
    modular_vacancies_2d,
    solution_cells,
    spectral_K_max,
    triangular_torus,
    verify_modular_K,
)
from hsga.geometry.lattices import named_structure, triangular_lattice


def test_counting_identity_reproduces_the_ladder_to_1e12(golden):
    for K in (1, 2, 3, 4):
        ref = golden["fcc_depletion_ladder"][f"K{K}"]["eta"]
        assert abs(counting_eta(12, K, ETA_FCC) - ref) < 1e-12
        assert abs(counting_eta(12, K, ETA_FCC) - 2 * np.sqrt(2) * np.pi / (12 + K)) < 1e-15
    for K, key in ((3, "K3_honeycomb"), (2, "K2_kagome"), (1, "K1_maple_leaf")):
        ref = golden["family_2d"][key]["phi"]
        assert abs(counting_eta(6, K, PHI_TRI) - ref) < 1e-12


@pytest.mark.slow
def test_spectral_K_max():
    assert spectral_K_max(*named_structure("FCC"), nq=60_000) == 4
    assert spectral_K_max(*triangular_lattice(), nq=60_000) == 3


@pytest.mark.parametrize("K", [4, 3, 2, 1])
def test_modular_rules_3d(K, golden):
    ref = golden["fcc_depletion_ladder"][f"K{K}"]
    r = verify_modular_K(K, dim=3)
    assert r["uniform"] and r["independent"] and r["tangential"] and r["congruent"]
    assert r["z"] == ref["z"]
    assert r["p_vector"] == ref["p_vector"]
    assert abs(r["eta"] - ref["eta"]) < 1e-12
    assert r["eta_err"] < 1e-12


@pytest.mark.parametrize("K,key", [(3, "K3_honeycomb"), (2, "K2_kagome"), (1, "K1_maple_leaf")])
def test_modular_rules_2d(K, key, golden):
    ref = golden["family_2d"][key]
    r = verify_modular_K(K, dim=2)
    assert r["uniform"] and r["independent"] and r["tangential"] and r["congruent"]
    assert r["z"] == ref["z"]
    assert abs(r["eta"] - ref["phi"]) < 1e-12


def test_rules_callable_forms():
    assert modular_vacancies(4)(2, 4, 6) and not modular_vacancies(4)(1, 1, 0)
    assert modular_vacancies_2d(3)(3, 0) and not modular_vacancies_2d(3)(1, 0)


def test_divacancy_spread_is_sqrt2_minus_1(golden):
    d = divacancy_spread()
    assert d["isolated"] == pytest.approx(0.0, abs=1e-9)
    assert d["divacancy"] == pytest.approx(np.sqrt(2) - 1, abs=1e-9)


def test_degeneracy_counts_golden(golden):
    for key, torus, K in [
        ("tri_4x4_K2", triangular_torus(4), 2),
        ("tri_6x6_K3", triangular_torus(6), 3),
        ("tri_6x6_K2", triangular_torus(6), 2),
    ]:
        ref = golden["degeneracy_counts"][key]
        r = enumerate_uniform_K(torus, K)
        assert r["complete"]
        assert r["solutions"] == ref["solutions"]
        assert r["orbits"] == ref["orbits"]
        assert r["cosets"] == ref["cosets"]


def test_6x6_K2_noncoset_orbits_are_exactly_tangential(golden):
    """Every cell of every solution: tangential, Q = sqrt3 pi/8."""
    tor = triangular_torus(6)
    r = enumerate_uniform_K(tor, 2)
    target = np.sqrt(3) * np.pi / 8
    for rep in r["orbit_representatives"]:
        for c in solution_cells(rep, tor):
            assert abs(c.Q_iso - target) < 1e-9
            assert (c.face_distances.max() - c.face_distances.min()) < 1e-9


def test_non_integer_Nv_means_no_solutions():
    r = enumerate_uniform_K(triangular_torus(5), 2)   # Nv = 50/8, not integer
    assert r["solutions"] == 0 and r["complete"]
    assert "not an integer" in r["note"]


def test_caps_flag_incompleteness():
    r = enumerate_uniform_K(triangular_torus(6), 2, max_nodes=100)
    assert not r["complete"]


def test_fcc_torus_side4_K4_all_cosets():
    r = enumerate_uniform_K(fcc_torus(4), 4)
    assert r["complete"]
    assert r["solutions"] == 16 and r["orbits"] == 4 and r["cosets"] == 4


def test_is_lattice_coset_detects_non_cosets():
    tor = triangular_torus(6)
    r = enumerate_uniform_K(tor, 2)
    flags = [
        is_lattice_coset([i for i, v in enumerate(rep) if v == 1], tor)
        for rep in r["orbit_representatives"]
    ]
    assert sum(flags) == 1                       # only the kagome pattern

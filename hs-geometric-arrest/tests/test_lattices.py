"""Selling enumeration, the 14-face exclusion, and the Sec. V corrections."""

import numpy as np
import pytest

from hsga.geometry.lattices import (
    CLASS_NORMS,
    enumerate_tangential_lattices,
    face_class_norms,
    lattice_from_gram,
    no_fourteen_faced_solution,
    selling_gram,
    tetrahedral_octahedral_honeycomb,
    vacancy_shells,
)


def test_class_norms_are_the_real_squared_norms():
    rng = np.random.default_rng(7)
    for _ in range(25):
        p = rng.uniform(0.05, 1.0, 6)
        p /= p.sum()
        G = selling_gram(p)
        if np.linalg.eigvalsh(G).min() <= 0:
            continue
        B = lattice_from_gram(G)
        b4 = -(B[0] + B[1] + B[2])
        vectors = [B[0], B[1], B[2], b4, B[0] + B[1], B[0] + B[2], B[1] + B[2]]
        expected = np.array([v @ v for v in vectors])
        assert np.allclose(face_class_norms(p), expected, rtol=1e-12, atol=1e-12)


def test_class_norms_are_linear_in_p():
    rng = np.random.default_rng(11)
    p, q = rng.uniform(0, 1, 6), rng.uniform(0, 1, 6)
    assert np.allclose(
        face_class_norms(0.3 * p + 0.7 * q),
        0.3 * face_class_norms(p) + 0.7 * face_class_norms(q),
    )


@pytest.mark.slow
def test_exactly_three_tangential_bravais_lattices(golden):
    ref = golden["tangential_bravais_3d"]
    res = enumerate_tangential_lattices()
    assert len(res) == ref["count_is_complete"] == 3
    assert [r["n_faces"] for r in res] == [6, 8, 12]
    assert all(r["verified_tangential"] for r in res)
    expect = [ref["simple_cubic"]["eta"], ref["simple_hexagonal_ca"]["eta"], ref["fcc"]["eta"]]
    assert np.allclose([r["eta"] for r in res], expect, atol=1e-9)
    assert np.allclose([r["verified_eta"] for r in res], expect, atol=1e-9)


def test_no_fourteen_faced_tangential_cell():
    """4t = 3t = 2 sum(p) forces t = 0: exact, over the integers."""
    n = no_fourteen_faced_solution()
    assert n["identities_hold"] and n["excluded"]
    assert np.array_equal(CLASS_NORMS[:4].sum(axis=0).astype(int), np.full(6, 2))
    assert np.array_equal(CLASS_NORMS[4:].sum(axis=0).astype(int), np.full(6, 2))
    # the symmetric point p_ij = 1/6 quoted in REFERENCE_VALUES: norms 1/2, 2/3
    got = sorted(set(np.round(n["class_norms_symmetric"], 12)))
    assert got == [0.5, pytest.approx(2 / 3)]


@pytest.mark.parametrize(
    "structure,counts,ratios",
    [("BCC", [8, 6], [71 / 64, 49 / 48]), ("FCC", [12, 6], [13 / 12, 1.0])],
)
def test_vacancy_shell_counts_and_volume_balance(structure, counts, ratios):
    """A BCC vacancy has 8 nearest neighbours, not 4, and the volume balances."""
    v = vacancy_shells(structure)
    assert [s["count"] for s in v["shells"]] == counts
    assert np.allclose([s["V_over_V0"] for s in v["shells"]], ratios, rtol=1e-12)
    assert v["balances"]


def test_bcc_vacancy_eta_matches_reference(golden):
    ref = golden["de_graaf_corrections"]["bcc_vacancy_nn"]
    v = vacancy_shells("BCC")
    assert v["shells"][0]["count"] == ref["correct"] == 8
    assert v["shells"][0]["eta_local"] == pytest.approx(8 * np.sqrt(3) * np.pi / 71, rel=1e-9)
    # non-tangential either way: averaging them carries no geometric content
    for s in v["shells"]:
        assert abs(s["Q_iso"] - s["eta_local"]) > 1e-2


def test_tet_oct_honeycomb_correction(golden):
    ref = golden["de_graaf_corrections"]["tet_oct"]
    t = tetrahedral_octahedral_honeycomb()
    assert t["fill_fraction_original"] == pytest.approx(0.75, rel=1e-12)
    assert t["reduces_to_rhombic_dodecahedron"]
    assert t["eta_corrected"] == pytest.approx(2 * np.sqrt(3) * np.pi / 27, rel=1e-12)
    assert t["eta_original"] == pytest.approx(ref["paper_eta"], abs=1e-6)


def test_honeycomb_depletion_fraction(golden):
    """The 2D honeycomb is 1-in-3 removed, not 2-in-3."""
    assert golden["de_graaf_corrections"]["honeycomb_2d"]["correct"] == "remove 1 in 3"
    phi_quoted = np.pi / (3 * np.sqrt(3))
    assert (2 / 3) * np.pi / np.sqrt(12) == pytest.approx(phi_quoted, rel=1e-12)
    assert (1 / 3) * np.pi / np.sqrt(12) != pytest.approx(phi_quoted, rel=1e-3)

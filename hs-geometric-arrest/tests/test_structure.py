"""Persistent topology, refscore, P_wrap and the isoconfigurational machinery."""

import numpy as np
import pytest

from conftest import dense_configuration
from hsga.analysis.isoconfig import (
    baseline_features, held_out_r2, propensity_from_frames, tetrahedrality_3d,
)
from hsga.analysis.pwrap import finite_size_crossings, p_wrap
from hsga.analysis.refscore import (
    load_frozen, mark_cells, p_distance, reference_descriptors, refscores,
)
from hsga.analysis.topology import face_adjacency, face_filter, p_vector, persistence, topo_hash
from hsga.analysis.voronoi import config_cells
from hsga.geometry.coloring import realise_modular_K
from hsga.geometry.lattices import named_structure
from hsga.geometry.tangential import voronoi_cell


def _cell(name):
    if name.startswith("K"):
        lat, basis, _ = realise_modular_K(int(name[1]), dim=3)
    else:
        lat, basis = named_structure(name)
    return voronoi_cell(basis[0], lat, basis)


def test_p_vector_equals_loop_lengths_on_closed_cells():
    """On a closed cell each polygon edge is shared with exactly one face."""
    for name in ("FCC", "SC", "simple hexagonal c=a", "K4", "K1"):
        c = _cell(name)
        assert p_vector(c) == tuple(sorted(len(l) for l in c.face_vertices))


def test_face_filter_drops_only_marginal_faces():
    pos, rad, L = dense_configuration(nrep=3, eta=0.45, poly=0.08, seed=2)
    cells = config_cells(pos, rad, L)
    fz = load_frozen()["face_filter"]
    for c in cells[:10]:
        f = face_filter(c, fz["g_cut"], fz["a_cut"])
        assert f.n_faces <= c.n_faces
        # every dropped face was BOTH distant and small
        kept = set(map(tuple, np.array(f.face_vertices, dtype=object)))
        for k in range(c.n_faces):
            if tuple(c.face_vertices[k]) not in kept:
                assert c.face_gaps[k] > fz["g_cut"]
                assert c.face_areas[k] / c.S < fz["a_cut"]
    with pytest.raises(ValueError):
        face_filter(cells[0], -1.0, 2.0)          # would drop everything


def test_topo_hash_invariance_under_relabelling():
    """The hash is a graph invariant: rebuilding the same cell gives the same hash."""
    c1 = _cell("K3")
    c2 = _cell("K3")
    assert topo_hash(c1) == topo_hash(c2)
    assert topo_hash(c1) != topo_hash(_cell("K2"))


def test_persistence_is_low_for_fragile_topology():
    """A thermal cell's sliver faces come and go across the filter grid."""
    pos, rad, L = dense_configuration(nrep=3, eta=0.50, poly=0.10, seed=5)
    cells = config_cells(pos, rad, L)
    g = load_frozen()["persistence_grid"]
    vals = [persistence(c, g["g_grid"], g["a_grid"]) for c in cells[:20]]
    assert min(vals) < 1.0                        # at least one fragile cell
    assert all(0.0 < v <= 1.0 for v in vals)


def test_p_distance_properties():
    assert p_distance([4] * 12, [4] * 12) == 0.0
    assert p_distance([3] * 8, [4] * 12) > 0
    assert p_distance([3, 3, 4], [3, 4, 3]) == 0.0      # order-free
    assert 0.0 <= p_distance([3] * 8, [4] * 11) <= 1.0


def test_reference_descriptors_match_golden(golden):
    refs = reference_descriptors(3)
    for K in (1, 2, 3, 4):
        g = golden["fcc_depletion_ladder"][f"K{K}"]
        assert refs[f"K{K}"]["eta"] == pytest.approx(g["eta"], abs=1e-12)
        assert list(refs[f"K{K}"]["p_vector"]) == g["p_vector"]
    refs2 = reference_descriptors(2)
    assert refs2["K2_kagome"]["eta"] == pytest.approx(
        golden["family_2d"]["K2_kagome"]["phi"], abs=1e-12)


def test_refscore_discriminates_on_thermal_cells():
    """On a moderately dense fluid, cells are NOT marked as reference-like."""
    pos, rad, L = dense_configuration(nrep=3, eta=0.45, poly=0.10, seed=7)
    cells = config_cells(pos, rad, L)
    scores, names = refscores(cells, rad, 3)
    assert scores.shape == (len(cells), len(names))
    assert (scores >= 0).all()
    marked = mark_cells(scores, names, "fcc")
    assert marked.mean() < 0.5                    # a fluid is not mostly FCC


def test_p_wrap_null_on_a_fluid():
    pos, rad, L = dense_configuration(nrep=3, eta=0.45, poly=0.10, seed=9)
    cells = config_cells(pos, rad, L)
    pw = p_wrap(cells, pos, rad, L, "K1")
    assert not pw["wraps"]
    assert pw["marked_fraction"] < 0.5
    assert pw["chi"] >= 0.0 and pw["xi"] >= 0.0


def test_finite_size_crossings():
    eta = np.linspace(0.5, 0.7, 21)
    c1 = 1 / (1 + np.exp(-(eta - 0.60) / 0.02))       # N small: soft
    c2 = 1 / (1 + np.exp(-(eta - 0.60) / 0.008))      # N large: sharp
    out = finite_size_crossings({100: (eta, c1), 1000: (eta, c2)})
    assert len(out) >= 1
    assert abs(out[0]["eta_cross"] - 0.60) < 0.01


# ------------------------------------------------------------------ isoconfig
def test_propensity_from_frames_minimum_image():
    rng = np.random.default_rng(0)
    L = 10.0
    start = rng.uniform(0, L, (50, 3))
    runs = []
    for _ in range(4):
        fr = []
        for lag in (0.1, 0.2):
            step = rng.normal(0, lag, start.shape)
            fr.append(((start + step) % L, None, L))
            last = start + step
        runs.append(fr)
    prop = propensity_from_frames(start, runs, L)
    assert prop.shape == (2, 50)
    assert (prop >= 0).all()
    assert prop[1].mean() > prop[0].mean()        # more time, more displacement


def test_tetrahedrality_range():
    pos, rad, L = dense_configuration(nrep=3, eta=0.45, seed=3)
    q = tetrahedrality_3d(pos, L)
    assert q.shape == (len(pos),)
    assert (q <= 1.0 + 1e-9).all()


def test_baseline_features_shapes():
    pos, rad, L = dense_configuration(nrep=3, eta=0.45, poly=0.05, seed=4)
    cells = config_cells(pos, rad, L)
    X, names = baseline_features(cells, pos, rad, L)
    assert X.shape == (len(pos), 4)
    assert names[0] == "local_eta" and "V_cell" in names
    assert np.isfinite(X).all()


def test_held_out_r2_detects_a_real_signal_and_a_null():
    rng = np.random.default_rng(1)
    blocks, props = [], []
    for cfg in range(4):
        Xb = rng.normal(size=(200, 3))
        Xs = rng.normal(size=(200, 2))
        y = 1.0 * Xb[:, 0] + 2.0 * Xs[:, 0] + 0.1 * rng.normal(size=200)
        blocks.append((Xb, Xs))
        props.append(y)
    res = held_out_r2(blocks, props)
    assert res["delta_r2_mean"] > 0.2             # the structural block matters
    # null: propensity depends only on the baseline
    props0 = [1.0 * b[0][:, 0] + 0.1 * rng.normal(size=200) for b in blocks]
    res0 = held_out_r2(blocks, props0)
    assert abs(res0["delta_r2_mean"]) < 0.05      # nothing added
    with pytest.raises(ValueError):
        held_out_r2(blocks[:1], props[:1])

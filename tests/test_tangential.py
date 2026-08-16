"""The theorem of paper.tex Sec. II and the cell builder that implements it."""

import numpy as np
import pytest

from conftest import dense_configuration
from hsga.geometry.lattices import named_structure, triangular_lattice
from hsga.geometry.tangential import (
    OverlapError,
    PoolTooSmall,
    analyse,
    cell_from_neighbours,
    is_tangential,
    voronoi_cell,
)


def test_tangential_bravais_and_named_structures(golden):
    tol = golden["tolerance_exact_geometry"]
    ref3 = golden["tangential_bravais_3d"]
    names = {"SC": "simple_cubic", "simple hexagonal c=a": "simple_hexagonal_ca", "FCC": "fcc"}
    for name, key in names.items():
        res = analyse(name, *named_structure(name))
        assert res["tangential"] and res["congruent"]
        assert abs(res["eta"] - ref3[key]["eta"]) < tol
        assert abs(res["Q_iso"] - ref3[key]["eta"]) < tol
        assert int(np.unique(res["n_faces"])[0]) == ref3[key]["faces"]
        assert int(np.unique(res["z"])[0]) == ref3[key]["z"]
        # the theorem itself
        assert abs(res["eta"] - res["Q_iso"]) < 1e-9


def test_hcp_is_tangential_at_fcc_density(golden):
    res = analyse("HCP", *named_structure("HCP"))
    assert res["tangential"]
    assert abs(res["eta"] - golden["tangential_bravais_3d"]["fcc"]["eta"]) < 1e-8


def test_inadmissible_structures(golden):
    for name, key in [("BCC", "bcc"), ("diamond", "diamond")]:
        ref = golden["inadmissible"][key]
        res = analyse(name, *named_structure(name))
        assert not res["tangential"]
        assert abs(res["eta"] - ref["eta"]) < 1e-8
        assert abs(res["Q_iso"] - ref["Q_iso"]) < 1e-6
        assert int(np.unique(res["n_faces"])[0]) == ref["faces"]
        assert int(np.unique(res["z"])[0]) == ref["z"]
        # the whole point: eta != Q_iso once a face is not a contact
        assert abs(res["eta"] - res["Q_iso"]) > 1e-3


def test_bcc_contact_fraction_is_eight_fourteenths(golden):
    lat, bas = named_structure("BCC")
    sigma = np.sqrt(3.0) / 2.0
    cell = voronoi_cell(bas[0], lat, bas, radii=np.full(len(bas), sigma / 2))
    assert cell.n_faces == 14
    assert cell.n_contacts == 8
    assert cell.f_c == pytest.approx(golden["engine_gates"]["bcc_contact_fraction"])


def test_2d_triangular_lattice_theorem():
    """phi = q for the hexagonal cell of the triangular lattice."""
    res = analyse("triangular", *triangular_lattice())
    assert res["dim"] == 2
    assert res["tangential"]
    assert res["eta"] == pytest.approx(np.pi / np.sqrt(12.0), abs=1e-12)
    assert res["Q_iso"] == pytest.approx(res["eta"], abs=1e-9)
    assert int(np.unique(res["n_faces"])[0]) == 6


def test_face_decomposition_rebuilds_the_hull_3d():
    pos, rad, L = dense_configuration(nrep=3, eta=0.40, poly=0.08, seed=1)
    lat = np.eye(3) * L
    for i in range(5):
        c = voronoi_cell(pos[i], lat, pos, radii=rad)
        assert (c.face_areas * c.face_distances).sum() / 3.0 == pytest.approx(c.V, rel=1e-9)
        assert c.face_areas.sum() == pytest.approx(c.S, rel=1e-9)
        assert c.Q_iso <= 1.0 + 1e-12                 # isoperimetric inequality
        # face loops close: every loop has n_vertices >= 3 distinct vertices
        for loop in c.face_vertices:
            assert len(loop) >= 3
            assert len(set(loop)) == len(loop)


def test_face_decomposition_rebuilds_the_hull_2d():
    pos, rad, L = dense_configuration(nrep=5, eta=0.55, poly=0.05, seed=2, dim=2)
    lat = np.eye(2) * L
    for i in range(5):
        c = voronoi_cell(pos[i], lat, pos, radii=rad)
        assert c.dim == 2
        assert (c.face_areas * c.face_distances).sum() / 2.0 == pytest.approx(c.V, rel=1e-9)
        assert c.Q_iso <= 1.0 + 1e-12
        for loop in c.face_vertices:
            assert len(loop) == 2


def test_cells_tile_the_box():
    pos, rad, L = dense_configuration(nrep=3, eta=0.45, poly=0.1, seed=4)
    cells = [voronoi_cell(p, np.eye(3) * L, pos, radii=rad) for p in pos[:64]]
    # full tiling needs all cells; use a perfect lattice for the exact check
    lat, bas = named_structure("FCC")
    tot = sum(voronoi_cell(b, lat, bas).V for b in bas)
    assert tot == pytest.approx(1.0, rel=1e-9)        # unit conventional cube
    assert all(c.V > 0 for c in cells)


def test_radical_cell_reduces_to_voronoi_for_equal_radii():
    pos, rad, L = dense_configuration(nrep=3, eta=0.30, seed=0)
    lat = np.eye(3) * L
    r = float(rad.min())
    a = voronoi_cell(pos[0], lat, pos, radii=np.full(len(pos), 0.9 * r))
    b = voronoi_cell(pos[0], lat, pos, radii=np.full(len(pos), 0.2 * r))
    assert a.V == pytest.approx(b.V, rel=1e-12)
    assert a.S == pytest.approx(b.S, rel=1e-12)
    assert a.n_faces == b.n_faces


def test_overlapping_input_raises():
    disp = np.array([[0.5, 0, 0], [-0.5, 0, 0], [0, 0.5, 0],
                     [0, -0.5, 0], [0, 0, 0.5], [0, 0, -0.5]])
    with pytest.raises(OverlapError, match="overlaps by"):
        cell_from_neighbours(disp, 0.6, np.full(6, 0.6), np.arange(6))


def test_unbounded_half_space_set_is_reported_as_short_pool():
    disp = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0],
                     [1.0, 1.0, 0], [1.0, 0, 1.0], [0, 1.0, 1.0]])
    with pytest.raises(PoolTooSmall, match="do not bound"):
        cell_from_neighbours(disp, 0.1, np.full(6, 0.1), np.arange(6))


def test_is_tangential_is_sharp():
    lat, bas = named_structure("FCC")
    sigma = 1.0 / np.sqrt(2.0)
    c = voronoi_cell(bas[0], lat, bas, radii=np.full(len(bas), sigma / 2))
    assert is_tangential(c, sigma, tol=1e-9)
    assert not is_tangential(c, sigma * (1 + 1e-6), tol=1e-9)


def test_vertex_merging_gives_correct_polygon_sizes():
    """Qhull duplicates vertices where >3 planes meet; loops must not."""
    lat, bas = named_structure("FCC")
    c = voronoi_cell(bas[0], lat, bas)
    assert sorted(len(l) for l in c.face_vertices) == [4] * 12   # rhombic dodecahedron
    lat, bas = named_structure("SC")
    c = voronoi_cell(bas[0], lat, bas)
    assert sorted(len(l) for l in c.face_vertices) == [4] * 6    # cube
    lat, bas = named_structure("simple hexagonal c=a")
    c = voronoi_cell(bas[0], lat, bas)
    assert sorted(len(l) for l in c.face_vertices) == [4] * 6 + [6] * 2  # hex prism

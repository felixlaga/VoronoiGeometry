"""Percolation, EOS, Voronoi and dynamics analysis, 2D and 3D."""

import numpy as np
import pytest

from conftest import dense_configuration
from hsga.analysis.dynamics import diffusion_coefficient, fit_eta_a, msd_exponent
from hsga.analysis.eos import (
    carnahan_starling, compressibility, contact_value, henderson_2d, radial_distribution,
)
from hsga.analysis.percolation import (
    UnionFind, eps_star, eps_star_bisect, percolates, wrapping_directions,
)
from hsga.analysis.voronoi import (
    config_cells, config_observables, f_c_of_delta, q_submode_weight,
    single_scale_collapse,
)


# ------------------------------------------------------------------ percolation
def test_union_find_wrap_detection():
    uf = UnionFind(3, 3)
    uf.union(0, 1, [0, 0, 0])
    uf.union(1, 2, [0, 0, 0])
    uf.union(2, 0, [1, 0, 0])            # loop closes one image over in x
    assert uf.wrap.tolist() == [True, False, False]
    uf = UnionFind(3, 2)
    uf.union(0, 1, [0, 0])
    uf.union(1, 2, [0, 0])
    uf.union(2, 0, [0, 0])               # compact loop: no wrap
    assert not uf.wrap.any()


def test_chain_wraps_only_its_axis_3d():
    L, n = 10.0, 20
    pos = np.zeros((n, 3))
    pos[:, 0] = 5.0
    pos[:, 1] = 5.0
    pos[:, 2] = np.arange(n) * (L / n)
    rad = np.full(n, 0.5 * L / n)
    w = wrapping_directions(pos, rad, L, 1e-6)
    assert w.tolist() == [False, False, True]
    assert not percolates(pos, rad, L, 1e-6)   # all three axes required


@pytest.mark.engine
@pytest.mark.parametrize("dim", [2, 3])
def test_exact_and_bisect_agree_and_are_sharp(dim):
    """On real equilibrated fluids: a jittered lattice keeps near-contact pairs
    and percolates below the bracket, which is exactly why the estimator is
    validated on engine output (and, at scale, on the pilot data)."""
    from conftest import engine_run, has_compiler

    if not has_compiler():
        pytest.skip("no C compiler")
    r = engine_run(eta=0.55 if dim == 3 else 0.72, mode=0, dim=dim,
                   ncell=4 if dim == 3 else 8, eq=800, prod=1500, nsnap=2)
    pos, rad, L = r["frames"][-1]
    e = eps_star(pos, rad, L)
    b = eps_star_bisect(pos, rad, L, iters=22)
    assert np.isfinite(e) and e == pytest.approx(b, rel=1e-3)
    assert not percolates(pos, rad, L, e * 0.999)
    assert percolates(pos, rad, L, e * 1.001)


@pytest.mark.engine
def test_eps_star_monotone_with_density_2d():
    from conftest import engine_run, has_compiler

    if not has_compiler():
        pytest.skip("no C compiler")
    # widely spaced densities, averaged over frames: a single tiny-box
    # snapshot fluctuates by ~10% and adjacent densities can invert
    vals = []
    for phi in (0.55, 0.68, 0.79):
        r = engine_run(eta=phi, mode=0, dim=2, ncell=8, eq=800, prod=1500, nsnap=3)
        es = [eps_star(pos, rad, L) for pos, rad, L in r["frames"]]
        vals.append(float(np.mean(es)))
    assert all(np.isfinite(v) for v in vals)
    assert vals[0] > vals[1] > vals[2]


def test_failed_bracket_returns_nan():
    pos, rad, L = dense_configuration(nrep=4, eta=0.05, jitter=0.3, seed=23)
    assert np.isnan(eps_star(pos, rad, L, hi=0.01))
    assert np.isnan(eps_star_bisect(pos, rad, L, hi=0.01))


# ------------------------------------------------------------------ EOS
def test_reference_eos_values():
    assert float(carnahan_starling(0.30)) == pytest.approx(3.973764, rel=1e-6)
    for eta in (0.1, 0.25, 0.4):
        g = (1 - eta / 2) / (1 - eta) ** 3
        assert 1 + 4 * eta * g == pytest.approx(float(carnahan_starling(eta)), rel=1e-12)
    assert float(henderson_2d(0.0)) == pytest.approx(1.0)
    assert float(henderson_2d(0.5)) == pytest.approx((1 + 0.25 / 8) / 0.25, rel=1e-12)


def test_gr_bin_edge_exactly_at_contact():
    pos, rad, L = dense_configuration(nrep=4, eta=0.35, seed=0)
    rad = np.full_like(rad, rad.min())            # strictly monodisperse
    gr = radial_distribution([(pos, rad, L)])
    sigma = gr["sigma"]
    assert np.min(np.abs(gr["edges"] - sigma)) < 1e-12
    inside = (gr["edges"][:-1] < sigma - 1e-12) & (gr["edges"][1:] > sigma + 1e-12)
    assert not inside.any()
    assert gr["n_below_contact"] == 0


def test_gr_rejects_overlap_and_polydispersity():
    L = 10.0
    pos = np.array([[0.0, 0, 0], [0.5, 0, 0], [5.0, 5, 5]])
    with pytest.raises(ValueError, match="below contact"):
        radial_distribution([(pos, np.full(3, 0.5), L)], nbin=50)
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="monodisperse-only"):
        radial_distribution([(rng.uniform(0, L, (20, 3)), rng.uniform(0.4, 0.6, 20), L)])


def test_contact_value_recovers_known_curve_and_reports_systematic():
    sigma = 1.0
    r = np.linspace(sigma, 1.5 * sigma, 400)
    cv = contact_value(r, 3.25 * np.exp(-4.0 * (r - sigma)), sigma)
    assert cv["g_contact"] == pytest.approx(3.25, rel=1e-8)
    assert cv["err_window"] < 1e-6
    cv2 = contact_value(r, 3.0 + 12.0 * (r - sigma) ** 3, sigma)
    assert cv2["err_window"] > 0 and cv2["err"] >= cv2["err_window"]


def test_compressibility_2d_and_3d():
    c3 = compressibility(0.3, 2.4781, 0.01, dim=3)
    assert c3["Z"] == pytest.approx(1 + 4 * 0.3 * 2.4781)
    c2 = compressibility(0.6, 4.6, 0.01, dim=2)
    assert c2["Z"] == pytest.approx(1 + 2 * 0.6 * 4.6)
    assert c2["reference"] == "Henderson"


# ------------------------------------------------------------------ voronoi
def test_config_path_reproduces_lattice_exactly_3d():
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    reps = np.array([[i, j, k] for i in range(4) for j in range(4) for k in range(4)])
    pos = (reps[:, None, :] + basis[None, :, :]).reshape(-1, 3).astype(float)
    L = 4.0
    sigma = 1 / np.sqrt(2)
    o = config_observables(pos, np.full(len(pos), sigma / 2), L)
    assert o["mean_Q_iso"] == pytest.approx(np.pi / np.sqrt(18), abs=1e-9)
    assert o["mean_n_faces"] == pytest.approx(12.0)
    assert np.allclose(o["face_gaps"], 0.0, atol=1e-9)


def test_cells_tile_the_box_2d():
    pos, rad, L = dense_configuration(nrep=5, eta=0.60, poly=0.05, seed=4, dim=2)
    cells = config_cells(pos, rad, L)
    assert sum(c.V for c in cells) == pytest.approx(L**2, rel=1e-8)
    gaps = np.concatenate([c.face_gaps for c in cells])
    assert gaps.min() > 0                          # thermal fluid: no exact contact


def test_single_scale_distribution_collapses_and_two_scale_does_not():
    rng = np.random.default_rng(8)
    ratios = [single_scale_collapse(rng.exponential(s, 200000))["ratio"][0]
              for s in (0.25, 0.13, 0.07)]
    assert np.std(ratios) / np.mean(ratios) < 0.05
    ratios2 = []
    for w in (0.0, 0.5, 0.9):
        mix = np.where(rng.random(200000) < w,
                       rng.exponential(0.01, 200000), rng.exponential(0.3, 200000))
        ratios2.append(single_scale_collapse(mix)["ratio"][0])
    assert np.std(ratios2) / np.mean(ratios2) > 0.05


def test_f_c_of_delta_monotone():
    fc = f_c_of_delta(np.random.default_rng(6).exponential(0.1, 20000),
                      [1e-4, 1e-3, 1e-2, 1e-1])
    assert np.all(np.diff(fc) >= 0) and np.all((fc >= 0) & (fc <= 1))


def test_q_submode_weight_sees_a_planted_submode():
    rng = np.random.default_rng(11)
    radii = np.where(np.arange(4000) % 2 == 0, 0.5, 0.5 / 1.4)
    q = rng.uniform(0.6, 0.9, 4000)
    target = 0.7773425846718076
    # plant NET extra mass: move 10% of large-particle cells into the window
    large = radii > 0.4
    move = large & (rng.random(4000) < 0.10)
    q[move] = target + rng.normal(0, 0.003, int(move.sum()))
    r = q_submode_weight(q, radii, target)
    assert r["flank_ratio"] > 1.5
    flat = q_submode_weight(rng.uniform(0.6, 0.9, 4000), radii, target)
    assert 0.5 < flat["flank_ratio"] < 2.0


# ------------------------------------------------------------------ dynamics
def test_diffusion_and_exponent():
    t = np.linspace(1, 1000, 300)
    out3 = diffusion_coefficient(t, 6 * 3.5e-4 * t, dim=3)
    assert out3["D"] == pytest.approx(3.5e-4, rel=1e-9) and out3["diffusive"]
    out2 = diffusion_coefficient(t, 4 * 2e-4 * t, dim=2)
    assert out2["D"] == pytest.approx(2e-4, rel=1e-9)
    assert msd_exponent(t, np.full_like(t, 0.02)) == pytest.approx(0.0, abs=1e-6)


def test_fit_eta_a_window_scan_is_mandatory():
    eta = np.array([0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60])
    D = 2.5e-3 * (0.615 - eta) ** 3.8
    fit = fit_eta_a(eta, D)
    assert fit["eta_a"] == pytest.approx(0.615, abs=1e-4)
    assert fit["n_windows"] > 1
    assert fit["err"] == max(fit["err_stat"], fit["err_window"])
    assert "Brownian proxy" in fit["note"]
    fit2 = fit_eta_a(eta, D, dynamics_label="EDMD (Newtonian)")
    assert "EDMD" in fit2["note"]


def test_fit_eta_a_rejects_unphysical_windows():
    eta = np.array([0.40, 0.44, 0.48, 0.52, 0.56, 0.60])
    D = 1e-3 * np.exp(-30 * (eta - 0.40))          # extrapolates absurdly high
    fit = fit_eta_a(eta, D)
    for w in fit["windows"]:
        assert w["eta_a"] <= 0.7404804896930611
    assert "n_rejected_unphysical" in fit


def test_fit_eta_a_insufficient_data_returns_nan():
    fit = fit_eta_a([0.5, 0.55], [1e-4, 1e-5])
    assert np.isnan(fit["eta_a"]) and fit["n_windows"] == 0


def test_threshold_sweep_arrhenius_vs_pinned():
    """Constant drift for a smooth exponential slowdown; decelerating,
    bounded crossings when D vanishes at a fixed density."""
    from hsga.analysis.dynamics import threshold_sweep

    eta = np.linspace(0.58, 0.82, 60)
    r = threshold_sweep(eta, 10.0 ** (-(eta - 0.58) / 0.03))
    assert [c["decade"] for c in r["crossings"]] == list(range(1, 9))
    assert r["drift_per_decade"] == pytest.approx(0.03, abs=1e-6)

    with np.errstate(invalid="ignore"):
        D = np.where(eta < 0.78, np.clip(0.78 - eta, 0, None) ** 2.8, np.nan)
    r2 = threshold_sweep(eta, D)
    gaps = np.diff([c["eta_x"] for c in r2["crossings"]])
    assert (gaps > 0).all() and (np.diff(gaps) < 0).all()   # decelerating
    assert all(c["eta_x"] < 0.78 for c in r2["crossings"])  # bounded by eta_a

    # too little data: NaNs, never an invented drift
    r3 = threshold_sweep([0.6, 0.7], [1.0, 0.1])
    assert r3["crossings"] == [] and np.isnan(r3["drift_per_decade"])


def test_featuretest_calibration_and_detection():
    """Sideband/kink tests: unit-normal under a smooth null, detect a large
    planted kink, and the control geometry never lies about sensitivity."""
    from hsga.analysis.featuretest import (
        calibrate, empirical_p, kink_at_target, sideband_offset)

    rng = np.random.default_rng(11)
    t = 0.6801747615878316
    fine = np.arange(t - 0.006, t + 0.0061, 5e-4)
    coarse = np.array([0.64, 0.65, 0.66, 0.67, 0.69, 0.70, 0.71, 0.72])
    sb = np.array([t - 0.012, t - 0.0085, t - 0.007,
                   t + 0.007, t + 0.0085, t + 0.012])
    eta = np.sort(np.concatenate([coarse, fine, sb]))
    base = 0.04 - 0.185 * (eta - t) + 0.9 * (eta - t) ** 2
    sem = np.full(len(eta), 5e-5)

    null = calibrate(eta, base, sem, t, n_synth=300, degree=1, band=0.012)
    assert 0.8 < null["z_kink_null"].std() < 1.25   # unit-normal spread
    # curvature shifts the null centre; the centred p handles it

    y = base + rng.normal(0, sem)
    r = kink_at_target(eta, y, sem, t, degree=1, band=0.012)
    assert empirical_p(r["z_kink"], null["z_kink_null"]) > 0.05   # smooth: null

    y2 = base + 0.30 * 0.185 * np.clip(eta - t, 0, None) + rng.normal(0, sem)
    r2 = kink_at_target(eta, y2, sem, t, degree=1, band=0.012)
    assert empirical_p(r2["z_kink"], null["z_kink_null"]) < 0.01  # 30%: found

    off = sideband_offset(eta, y2 + 3e-4, sem, t)   # level shift
    assert off["tested"] and abs(off["z_offset"]) > 3

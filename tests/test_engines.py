"""The Monte Carlo engines and their driver, exercised at smoke scale."""

import numpy as np
import pytest

from conftest import engine_binary, engine_run, needs_engine
from hsga.engine.driver import RunSpec, parse_log, read_frames, run_one

pytestmark = [needs_engine, pytest.mark.engine]


def test_2d_engine_runs_audits_and_writes_frames():
    r = engine_run(eta=0.72, mode=0, dim=2, ncell=8, eq=500, prod=1000, nsnap=3)
    frames = r["frames"]
    assert len(frames) == 3
    pos, rad, L = frames[-1]
    assert pos.shape[1] == 2
    log = parse_log(r["prefix"] + ".log")
    assert log["final_overlap_audit"] == 0
    assert log["anneal_final_energy"] == 0.0
    assert log["radii_at_target"] == 1
    assert log["swap_moves_production"] == 0          # landmine 4, logged
    # no overlaps in the frame itself
    d = pos[:, None, :] - pos[None, :, :]
    d -= L * np.round(d / L)
    rr = np.sqrt((d**2).sum(-1))
    s = rad[:, None] + rad[None, :]
    np.fill_diagonal(rr, np.inf)
    assert (rr - s).min() > -1e-9


def test_2d_lattice_geometry_landmine():
    """nx = round(sqrt(3) ny): the box is near-square and the init fits phi=0.85."""
    import subprocess, tempfile

    d = tempfile.mkdtemp()
    exe = engine_binary("hsmc2d")
    p = subprocess.run(
        [exe, "--eta", "0.85", "--ncell", "8", "--mode", "1", "--seed", "1",
         "--eq", "50", "--prod", "100", "--nsnap", "1", "--melt", "50",
         "--prefix", f"{d}/hi"],
        capture_output=True, text=True,
    )
    assert p.returncode == 0                          # mono init survives 0.85
    log = parse_log(f"{d}/hi.log")
    assert log["nx"] == round(np.sqrt(3.0) * log["ny"])


def test_2d_restart_skips_init_and_preserves_radii():
    import tempfile

    base = engine_run(eta=0.75, mode=0, dim=2, ncell=8, eq=300, prod=600, nsnap=2)
    d = tempfile.mkdtemp()
    spec = RunSpec(eta=0.0, prefix=f"{d}/re", dim=2, ncell=8, mode=0, seed=77,
                   eq=100, prod=300, nsnap=2, infile=base["prefix"] + ".cfg")
    run_one(spec, engine_binary("hsmc2d"))
    log = parse_log(f"{d}/re.log")
    assert "restart_from" in log
    assert log["final_overlap_audit"] == 0
    frames = read_frames(f"{d}/re.cfg", dim=2)
    r0 = np.sort(base["frames"][-1][1])
    r1 = np.sort(frames[-1][1])
    assert np.allclose(r0, r1, atol=1e-8)             # same radius multiset


def test_2d_determinism_same_seed_same_output():
    import filecmp, subprocess, tempfile

    exe = engine_binary("hsmc2d")
    outs = []
    for tag in ("a", "b"):
        d = tempfile.mkdtemp()
        subprocess.run(
            [exe, "--eta", "0.70", "--ncell", "6", "--mode", "0", "--seed", "42",
             "--eq", "200", "--prod", "400", "--nsnap", "2", "--melt", "200",
             "--prefix", f"{d}/{tag}"],
            capture_output=True, check=True,
        )
        outs.append(f"{d}/{tag}")
    assert open(outs[0] + ".cfg").read() == open(outs[1] + ".cfg").read()
    assert open(outs[0] + ".msd").read() == open(outs[1] + ".msd").read()


def test_3d_engine_runs_audits_and_writes_frames():
    r = engine_run(eta=0.50, mode=0, dim=3, ncell=4, eq=500, prod=1000, nsnap=3)
    frames = r["frames"]
    assert len(frames) == 3
    pos, rad, L = frames[-1]
    assert pos.shape[1] == 3
    log = parse_log(r["prefix"] + ".log")
    assert log["final_overlap_audit"] == 0
    assert log["swap_moves_production"] == 0
    assert log["radii_at_target"] == 1


def test_3d_restart_mode():
    import tempfile

    base = engine_run(eta=0.50, mode=0, dim=3, ncell=4, eq=500, prod=1000, nsnap=3)
    d = tempfile.mkdtemp()
    spec = RunSpec(eta=0.0, prefix=f"{d}/re", dim=3, ncell=4, mode=0, seed=5,
                   eq=50, prod=200, nsnap=1, infile=base["prefix"] + ".cfg")
    run_one(spec, engine_binary("hsmc"))
    log = parse_log(f"{d}/re.log")
    assert "restart_from" in log and log["final_overlap_audit"] == 0


def test_unreachable_state_point_is_recorded_not_raised():
    """A binary 2D point above jamming exits 3; the driver records it."""
    import tempfile

    d = tempfile.mkdtemp()
    spec = RunSpec(eta=0.87, prefix=f"{d}/jam", dim=2, ncell=6, mode=0, seed=1,
                   eq=100, prod=100, nsnap=1, melt=200, anneal_max=2000)
    rec = run_one(spec, engine_binary("hsmc2d"), allow_anneal_failure=True)
    assert rec["unreachable"]
    import os

    assert not os.path.exists(f"{d}/jam.cfg")          # no data for a failed point
    with pytest.raises(RuntimeError):
        run_one(spec, engine_binary("hsmc2d"))         # without the flag it raises


def test_runspec_N():
    assert RunSpec(eta=0.5, prefix="x", dim=3, ncell=6).N == 864
    s = RunSpec(eta=0.7, prefix="x", dim=2, ncell=12)
    assert s.N == 2 * round(np.sqrt(3) * 12) * 12      # 504, the pilot's N
    assert s.N == 504

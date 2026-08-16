"""The event-driven MD engine: selftest, conservation, determinism.

EDMD exists for the stage-5 dynamical confirmation, which is downstream of
the T6 decision and has NOT been run; these tests validate the engine itself
at smoke scale (seconds).
"""

import subprocess

import numpy as np
import pytest

from conftest import engine_binary, engine_run, needs_engine
from hsga.engine.driver import EdmdSpec, parse_log, read_msd, run_edmd

pytestmark = [needs_engine, pytest.mark.engine]


def test_selftest_passes():
    exe = engine_binary("edmd")
    p = subprocess.run([exe, "--selftest"], capture_output=True, text=True)
    assert p.returncode == 0
    assert "FAIL" not in p.stdout
    assert p.stdout.count("PASS") >= 4


def test_conservation_audits_and_msd(tmp_path):
    base = engine_run(eta=0.45, mode=0, dim=3, ncell=4, eq=800, prod=1500, nsnap=2)
    spec = EdmdSpec(infile=base["prefix"] + ".cfg", prefix=str(tmp_path / "e"),
                    tmax=5.0, seed=7, nsample=50)
    rec = run_edmd(spec, engine_binary("edmd"))
    log = rec["log"]
    assert log["ke_drift_rel"] < 1e-9              # exact collision rule
    assert log["momentum_drift_rel"] < 1e-9
    assert log["final_overlap_audit_rel"] < 1e-7
    assert log["collisions"] > 0
    t, m = read_msd(f"{spec.prefix}.msd")
    assert t[0] == 0 and m[0] == 0
    assert m[-1] > 0 and (m >= -1e-15).all()


def test_determinism(tmp_path):
    base = engine_run(eta=0.45, mode=0, dim=3, ncell=4, eq=800, prod=1500, nsnap=2)
    exe = engine_binary("edmd")
    outs = []
    for tag in ("a", "b"):
        spec = EdmdSpec(infile=base["prefix"] + ".cfg", prefix=str(tmp_path / tag),
                        tmax=2.0, seed=99, nsample=20)
        run_edmd(spec, exe)
        outs.append(str(tmp_path / tag))
    assert open(outs[0] + ".msd").read() == open(outs[1] + ".msd").read()
    assert open(outs[0] + ".cfg").read() == open(outs[1] + ".cfg").read()


@pytest.mark.slow
def test_virial_eos_close_to_carnahan_starling(tmp_path):
    """Newtonian pressure route vs CS at eta=0.35, modest statistics (~3%)."""
    from hsga.analysis.eos import carnahan_starling

    base = engine_run(eta=0.35, mode=3, dim=3, ncell=4, eq=1500, prod=1500, nsnap=1)
    spec = EdmdSpec(infile=base["prefix"] + ".cfg", prefix=str(tmp_path / "z"),
                    tmax=80.0, seed=3, nsample=50)
    rec = run_edmd(spec, engine_binary("edmd"))
    Z = rec["log"]["Z_virial"]
    assert Z == pytest.approx(float(carnahan_starling(0.35)), rel=0.03)

import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))


@pytest.fixture(scope="session")
def golden():
    """REFERENCE_VALUES.json -- immutable golden numbers; tolerances are final."""
    return json.loads((REPO / "REFERENCE_VALUES.json").read_text())


@pytest.fixture(scope="session")
def frozen_refscore():
    return json.loads((REPO / "refscore_frozen.json").read_text())


# --------------------------------------------------------------------------- #
# synthetic dense configurations (no engine needed)
# --------------------------------------------------------------------------- #
def dense_configuration(nrep=5, eta=0.45, poly=0.0, jitter=0.35, seed=0, dim=3):
    """An overlap-free dense configuration in a cubic (square) box.

    Built by jittering a close-packed lattice, because that is the regime the
    analysis runs in; uniform random points at these densities overlap heavily
    and are not hard-particle configurations at all.  If jitter or
    polydispersity creates overlaps, all radii are shrunk by the single factor
    that just removes them, so the achieved packing fraction is at or below
    the requested one.
    """
    rng = np.random.default_rng(seed)
    if dim == 3:
        basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
        cells = np.array(
            [[i, j, k] for i in range(nrep) for j in range(nrep) for k in range(nrep)],
            dtype=float,
        )
        frac = (cells[:, None, :] + basis[None, :, :]).reshape(-1, 3) / nrep
    else:
        # centred-rectangular ~ triangular
        pts = []
        for i in range(nrep):
            for j in range(nrep):
                pts.append([i / nrep, j / nrep])
                pts.append([(i + 0.5) / nrep, (j + 0.5) / nrep])
        frac = np.array(pts)
    n = len(frac)
    scale = (
        np.clip(1.0 + poly * rng.standard_normal(n), 1 - 3 * poly, 1 + 3 * poly)
        if poly > 0
        else np.ones(n)
    )
    rad = 0.5 * scale
    if dim == 3:
        L = float(((4.0 / 3.0) * np.pi * (rad**3).sum() / eta) ** (1.0 / 3.0))
    else:
        L = float((np.pi * (rad**2).sum() / eta) ** 0.5)
    pos = frac * L
    nn = L / nrep / (np.sqrt(2.0) if dim == 3 else np.sqrt(2.0))
    pos = pos + rng.uniform(-0.5, 0.5, pos.shape) * jitter * nn
    pos -= L * np.floor(pos / L)

    from scipy.spatial import cKDTree

    tree = cKDTree(pos, boxsize=L)
    pairs = tree.query_pairs(r=2.0 * rad.max(), output_type="ndarray")
    if len(pairs):
        d = pos[pairs[:, 0]] - pos[pairs[:, 1]]
        d -= L * np.round(d / L)
        ratio = np.linalg.norm(d, axis=1) / (rad[pairs[:, 0]] + rad[pairs[:, 1]])
        rad = rad * min(1.0, float(ratio.min()) * (1.0 - 1e-9))
    return pos, rad, L


# --------------------------------------------------------------------------- #
# engine helpers (used by tests marked `engine`)
# --------------------------------------------------------------------------- #
_EXE_CACHE: dict = {}
_RUN_CACHE: dict = {}


def has_compiler():
    import shutil

    return bool(os.environ.get("CC") or shutil.which("cc") or shutil.which("gcc"))


needs_engine = pytest.mark.skipif(not has_compiler(), reason="no C compiler")


def engine_binary(name="hsmc"):
    from hsga.engine.driver import build_engine

    if name not in _EXE_CACHE:
        d = tempfile.mkdtemp(prefix="hsga-build-")
        _EXE_CACHE[name] = str(build_engine(d, name=name))
    return _EXE_CACHE[name]


def engine_run(eta=0.55, mode=0, seed=1, ncell=4, dim=3, eq=2000, prod=2000, nsnap=2):
    """A small cached MC run; returns {'frames': [...], 'prefix': str}."""
    from hsga.engine.driver import RunSpec, read_frames, run_one

    key = (dim, eta, mode, seed, ncell, eq, prod, nsnap)
    if key not in _RUN_CACHE:
        d = tempfile.mkdtemp(prefix="hsga-run-")
        name = "hsmc" if dim == 3 else "hsmc2d"
        spec = RunSpec(
            eta=eta, prefix=f"{d}/r", dim=dim, ncell=ncell, mode=mode, seed=seed,
            eq=eq, prod=prod, nsnap=nsnap, melt=min(eq, 3000),
        )
        run_one(spec, engine_binary(name))
        _RUN_CACHE[key] = {
            "frames": read_frames(f"{d}/r.cfg", dim=dim),
            "prefix": f"{d}/r",
        }
    return _RUN_CACHE[key]

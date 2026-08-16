"""Build, run and record the Monte Carlo engines (2D and 3D).

Owns the on-disk contract shared by ``hsmc.c`` and ``hsmc2d.c``:

* ``.cfg`` -- multi-frame; each frame has an ``N L`` header line, then ``N``
  lines of ``x y z r`` (3D) or ``x y r`` (2D);
* ``.msd`` -- ``sweep msd`` pairs from a single time origin;
* ``.log`` -- ``key=value`` lines, including the final overlap audit;
* exit codes -- 0 ok, 2 init failure, 3 anneal failed (physically unreachable
  state point), 4 overlap audit failed, 5 box too small for the cell list.

Every sweep records a provenance manifest: the exact command line of every
run, source hashes of both engines, compiler, library versions, and the parsed
log of each run.  Nothing downstream reads a ``.cfg`` whose run did not exit
zero.

Independence of configurations: a state point is sampled as ``n_replica``
independent runs times ``nsnap`` snapshots per run.  Only the replica axis is
independent by construction; both counts are recorded so a report can never
quote the product as the number of independent samples.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from itertools import islice
from pathlib import Path

import numpy as np

__all__ = [
    "ANNEAL_FAILED",
    "EdmdSpec",
    "MODES_2D",
    "MODES_3D",
    "RunSpec",
    "build_engine",
    "iter_frames",
    "parse_log",
    "read_frames",
    "read_msd",
    "run_edmd",
    "run_one",
    "run_sweep",
    "write_manifest",
]

MODES_3D = {
    0: "gaussian polydisperse (width --poly)",
    1: "binary 1:1 at R = 0.714",
    2: "trinary 1:1:1 at sigma, 0.8 sigma, 0.6 sigma",
    3: "monodisperse",
}
MODES_2D = {
    0: "binary 1:1 at R^-1 = 1.4 (the reference-study composition)",
    1: "monodisperse",
    2: "binary 1:1 at R^-1 = 1.7",
}

#: engine exit code for "could not remove the overlaps at this density"
ANNEAL_FAILED = 3

_CFLAGS = ["-O3", "-std=c11", "-Wall", "-Wextra"]


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #
def engine_source(name: str = "hsmc") -> Path:
    return Path(__file__).with_name(f"{name}.c")


def build_engine(build_dir: str | os.PathLike = "build", *, force: bool = False,
                 name: str = "hsmc") -> Path:
    """Compile one engine (``hsmc``, ``hsmc2d`` or ``edmd``); returns the path."""
    src = engine_source(name)
    out = Path(build_dir)
    out.mkdir(parents=True, exist_ok=True)
    exe = out / name
    if exe.exists() and not force and exe.stat().st_mtime > src.stat().st_mtime:
        return exe
    cc = os.environ.get("CC") or shutil.which("cc") or shutil.which("gcc")
    if cc is None:
        raise RuntimeError("no C compiler found; set CC")
    cmd = [cc, *_CFLAGS, "-o", str(exe), str(src), "-lm"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"engine build failed:\n{' '.join(cmd)}\n{proc.stderr}")
    return exe


# --------------------------------------------------------------------------- #
# run specification
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class RunSpec:
    """One Monte Carlo run, 2D or 3D."""

    eta: float                    # packing fraction (phi in 2D)
    prefix: str
    dim: int = 3
    ncell: int = 6                # nfcc (3D: N = 4 ncell^3) or ntri (2D)
    mode: int = 0
    seed: int = 1
    eq: int = 20000
    prod: int = 40000
    nsnap: int = 12
    melt: int = 20000
    poly: float = 0.10            # 3D mode 0 only
    swap: int = 1
    anneal_max: int = 200_000
    infile: str = ""              # restart: read this .cfg instead of initialising
    inframe: int = -1

    @property
    def N(self) -> int:
        if self.dim == 3:
            return 4 * self.ncell**3
        ny = self.ncell
        nx = int(round(np.sqrt(3.0) * ny))
        return 2 * nx * ny

    def argv(self, exe: str | os.PathLike) -> list[str]:
        args = [
            str(exe),
            "--eta", f"{self.eta:.6f}",
            "--ncell", str(self.ncell),
            "--mode", str(self.mode),
            "--seed", str(self.seed),
            "--eq", str(self.eq),
            "--prod", str(self.prod),
            "--nsnap", str(self.nsnap),
            "--melt", str(self.melt),
            "--swap", str(self.swap),
            "--anneal-max", str(self.anneal_max),
            "--prefix", self.prefix,
        ]
        if self.dim == 3:
            args += ["--poly", f"{self.poly:.6f}"]
        if self.infile:
            args += ["--in", self.infile, "--frame", str(self.inframe)]
        return args


@dataclass(frozen=True)
class EdmdSpec:
    """One event-driven MD run (3D; stage-5 confirmation only)."""

    infile: str
    prefix: str
    tmax: float = 100.0
    frame: int = -1
    seed: int = 1
    nsample: int = 200
    kT: float = 1.0
    equal_mass: int = 0

    def argv(self, exe: str | os.PathLike) -> list[str]:
        return [
            str(exe),
            "--in", self.infile,
            "--frame", str(self.frame),
            "--tmax", f"{self.tmax:.6f}",
            "--seed", str(self.seed),
            "--nsample", str(self.nsample),
            "--kT", f"{self.kT:.6f}",
            "--equal-mass", str(self.equal_mass),
            "--prefix", self.prefix,
        ]


def parse_log(path: str | os.PathLike) -> dict:
    """Parse a ``key=value`` engine log into floats/ints/strings."""
    out: dict = {}
    with open(path) as fh:
        for line in fh:
            if "=" not in line:
                continue
            k, _, v = line.strip().partition("=")
            v = v.split()[0] if v.split() else v
            try:
                out[k] = int(v)
            except ValueError:
                try:
                    out[k] = float(v)
                except ValueError:
                    out[k] = v
    return out


def run_one(spec: RunSpec, exe: str | os.PathLike, *,
            allow_anneal_failure: bool = False) -> dict:
    """Run one simulation.  Raises on non-zero exit or a non-zero overlap audit.

    ``allow_anneal_failure`` turns exit code 3 -- "overlaps could not be
    annealed away at this density" -- into a recorded outcome instead of an
    exception.  The failure is physical (a state point at or beyond its
    jamming density has no equilibrated fluid to sample); the run is marked
    ``unreachable``, produces no ``.cfg``, and nothing downstream can mistake
    it for data.  Every other failure still propagates.
    """
    Path(spec.prefix).parent.mkdir(parents=True, exist_ok=True)
    argv = spec.argv(exe)
    proc = subprocess.run(argv, capture_output=True, text=True)
    log_path = f"{spec.prefix}.log"
    log = parse_log(log_path) if os.path.exists(log_path) else {}
    rec = {"spec": asdict(spec), "argv": argv, "log": log, "unreachable": False}

    if proc.returncode == ANNEAL_FAILED and allow_anneal_failure:
        rec["unreachable"] = True
        rec["reason"] = (
            f"anneal did not reach zero overlap at eta={spec.eta} mode={spec.mode} "
            f"(residual energy {log.get('anneal_final_energy')})"
        )
        Path(f"{spec.prefix}.cfg").unlink(missing_ok=True)
        return rec
    if proc.returncode != 0:
        raise RuntimeError(
            f"engine exited {proc.returncode} for dim={spec.dim} eta={spec.eta} "
            f"mode={spec.mode} seed={spec.seed}\n  argv: {' '.join(argv)}\n"
            f"  log: {log}\n{proc.stderr}"
        )
    if log.get("final_overlap_audit", 1) != 0:
        raise RuntimeError(
            f"overlap audit non-zero for {spec.prefix}: {log.get('final_overlap_audit')}"
        )
    return rec


def run_edmd(spec: EdmdSpec, exe: str | os.PathLike) -> dict:
    """Run one EDMD simulation; raises on any non-zero exit."""
    Path(spec.prefix).parent.mkdir(parents=True, exist_ok=True)
    argv = spec.argv(exe)
    proc = subprocess.run(argv, capture_output=True, text=True)
    log_path = f"{spec.prefix}.log"
    log = parse_log(log_path) if os.path.exists(log_path) else {}
    if proc.returncode != 0:
        raise RuntimeError(
            f"edmd exited {proc.returncode} for {spec.infile}\n"
            f"  argv: {' '.join(argv)}\n  log: {log}\n{proc.stderr}"
        )
    return {"spec": asdict(spec), "argv": argv, "log": log}


def _run_one_star(args):
    spec, exe, allow = args
    return run_one(spec, exe, allow_anneal_failure=allow)


def run_sweep(specs, exe, *, nproc: int | None = None,
              allow_anneal_failure: bool = False) -> list[dict]:
    """Run many simulations in parallel.  Failures propagate (no skips), except
    anneal failures when ``allow_anneal_failure`` -- see :func:`run_one`."""
    specs = list(specs)
    if nproc is None:
        nproc = max(1, (os.cpu_count() or 2) - 1)
    if nproc == 1:
        return [run_one(s, exe, allow_anneal_failure=allow_anneal_failure) for s in specs]
    tasks = [(s, exe, allow_anneal_failure) for s in specs]
    with ProcessPoolExecutor(max_workers=nproc) as pool:
        return list(pool.map(_run_one_star, tasks))


# --------------------------------------------------------------------------- #
# provenance
# --------------------------------------------------------------------------- #
def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_rev() -> str | None:
    try:
        r = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parents[3],
            capture_output=True, text=True,
        )
        return r.stdout.strip() if r.returncode == 0 else None
    except OSError:
        return None


def write_manifest(path: str | os.PathLike, runs: list[dict], **extra) -> dict:
    """JSON provenance manifest next to the data it describes."""
    import scipy

    cc = os.environ.get("CC") or shutil.which("cc") or shutil.which("gcc")
    cc_version = ""
    if cc:
        v = subprocess.run([cc, "--version"], capture_output=True, text=True)
        cc_version = v.stdout.splitlines()[0] if v.stdout else ""
    man = {
        "engine_sources_sha256": {
            n: _sha256(engine_source(n)) for n in ("hsmc", "hsmc2d", "edmd")
        },
        "cflags": _CFLAGS,
        "cc": cc,
        "cc_version": cc_version,
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "git_rev": _git_rev(),
        "n_runs": len(runs),
        "n_unreachable": sum(1 for r in runs if r.get("unreachable")),
        "runs": runs,
        **extra,
    }
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(man, indent=2, default=str))
    return man


# --------------------------------------------------------------------------- #
# on-disk formats
# --------------------------------------------------------------------------- #
def iter_frames(path: str | os.PathLike, dim: int = 3):
    """Yield ``(positions, radii, L)`` for every frame of a ``.cfg`` file."""
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            header = header.split()
            if len(header) != 2:
                raise ValueError(f"{path}: malformed frame header {header!r}")
            n, L = int(header[0]), float(header[1])
            block = np.loadtxt(islice(fh, n)).reshape(n, dim + 1)
            yield block[:, :dim], block[:, dim], L


def read_frames(path: str | os.PathLike, dim: int = 3) -> list:
    return list(iter_frames(path, dim))


def read_msd(path: str | os.PathLike):
    d = np.loadtxt(path)
    if d.ndim == 1:
        d = d.reshape(1, -1)
    return d[:, 0], d[:, 1]

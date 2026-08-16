"""Gate G1 -- 3D equation of state against Carnahan-Starling.

Pass criterion (spec Sec. 3 / REFERENCE_VALUES ``engine_gates``): monodisperse
``Z = 1 + 4 eta g(sigma+)`` within 1% at ``eta = 0.30, 0.35, 0.40``.

The defaults ARE the statistics the gate needs -- with an order of magnitude
fewer pair counts the contact value scatters by several percent and the gate
fails on noise alone, which is a statement about sampling, not the engine.
At these defaults the gate is minutes of parallel compute, not hours.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from ..analysis.eos import compressibility, contact_value, radial_distribution
from ..engine.driver import RunSpec, build_engine, read_frames, run_sweep, write_manifest

REPO = Path(__file__).resolve().parents[3]
ETAS = (0.30, 0.35, 0.40)


def run(
    *,
    data_dir: str = "data/gate_s1",
    ncell: int = 6,
    eq: int = 20000,
    prod: int = 120000,
    nsnap: int = 120,
    seeds: tuple[int, ...] = (1, 2, 3, 4, 5, 6, 7, 8, 9),
    nproc: int | None = None,
    verbose: bool = True,
) -> dict:
    golden = json.loads((REPO / "REFERENCE_VALUES.json").read_text())
    tol = golden["engine_gates"]["eos_carnahan_starling_max_reldev_eta030"]
    exe = build_engine(name="hsmc")
    Path(data_dir).mkdir(parents=True, exist_ok=True)

    specs = [
        RunSpec(eta=eta, prefix=f"{data_dir}/mono_e{eta:.2f}_s{s}", dim=3,
                ncell=ncell, mode=3, seed=s, eq=eq, prod=prod, nsnap=nsnap, melt=eq)
        for eta in ETAS for s in seeds
    ]
    runs = run_sweep(specs, exe, nproc=nproc)
    write_manifest(f"{data_dir}/manifest.json", runs, gate="G1", etas=list(ETAS), tol=tol)

    rows, ok = [], True
    for eta in ETAS:
        frames = []
        for s in seeds:
            frames += read_frames(f"{data_dir}/mono_e{eta:.2f}_s{s}.cfg")
        gr = radial_distribution(frames)
        cv = contact_value(gr["r"], gr["g"], gr["sigma"], counts=gr["counts"])
        comp = compressibility(eta, cv["g_contact"], cv["err"], dim=3)
        good = comp["rel_dev"] < tol
        ok &= good
        rows.append({**comp, "g_contact": cv["g_contact"],
                     "g_err_window": cv["err_window"], "g_err_count": cv["err_count"],
                     "n_frames": len(frames), "passed": good})
        if verbose:
            print(f"  eta={eta:.2f}  g(sigma+)={cv['g_contact']:.5f} "
                  f"(window sys {cv['err_window']:.5f}, count {cv['err_count']:.5f})  "
                  f"Z={comp['Z']:.5f}  Z_CS={comp['Z_ref']:.5f}  "
                  f"dev={comp['rel_dev'] * 100:.3f}% +- {comp['rel_dev_err'] * 100:.3f}%  "
                  f"{'ok' if good else 'FAIL'}")
    if verbose:
        print(f"  GATE G1 {'PASSED' if ok else 'FAILED'}  (criterion < {tol:.0%})")
    return {"gate": "G1", "passed": bool(ok), "tol": tol, "rows": rows}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Gate G1: 3D EOS vs Carnahan-Starling")
    p.add_argument("--data-dir", default="data/gate_s1")
    p.add_argument("--ncell", type=int, default=6)
    p.add_argument("--eq", type=int, default=20000)
    p.add_argument("--prod", type=int, default=120000)
    p.add_argument("--nsnap", type=int, default=120)
    p.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5, 6, 7, 8, 9])
    p.add_argument("--nproc", type=int, default=None)
    a = p.parse_args(argv)
    print("GATE G1: 3D equation of state from the contact value of g(r)")
    res = run(data_dir=a.data_dir, ncell=a.ncell, eq=a.eq, prod=a.prod,
              nsnap=a.nsnap, seeds=tuple(a.seeds), nproc=a.nproc)
    return 0 if res["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

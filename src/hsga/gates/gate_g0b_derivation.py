"""Gate G0b -- the derivation of the ladder (geometry/coloring.py).

Pass criteria (task list T2 / spec Sec. 3):

* ``counting_eta`` reproduces the ladder to 1e-12;
* ``spectral_K_max(FCC) == 4`` and ``spectral_K_max(triangular) == 3``;
* the four 3D modular rules are uniform, independent, tangential and
  congruent, with the golden p-vectors;
* the kagome rung is exact at ``sqrt(3) pi / 8`` (with the 2D family);
* the divacancy spread is ``sqrt(2) - 1`` (isolated vacancy: 0);
* the ``degeneracy_counts`` block is reproduced exactly, and every cell of
  every 6x6 K=2 solution is exactly tangential at ``Q = sqrt(3) pi / 8``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from ..geometry.coloring import (
    ETA_FCC,
    PHI_TRI,
    counting_eta,
    divacancy_spread,
    enumerate_uniform_K,
    solution_cells,
    spectral_K_max,
    triangular_torus,
    verify_modular_K,
)
from ..geometry.lattices import named_structure, triangular_lattice

REPO = Path(__file__).resolve().parents[3]


def run(*, nq: int = 200_000, verbose: bool = True) -> dict:
    golden = json.loads((REPO / "REFERENCE_VALUES.json").read_text())
    checks: list[tuple[str, bool, str]] = []

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), detail))
        if verbose:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")

    # counting identity, 1e-12
    ok = all(
        abs(counting_eta(12, K, ETA_FCC) - golden["fcc_depletion_ladder"][f"K{K}"]["eta"]) < 1e-12
        for K in (1, 2, 3, 4)
    )
    ok2 = all(
        abs(counting_eta(6, K, PHI_TRI) - golden["family_2d"][key]["phi"]) < 1e-12
        for K, key in ((3, "K3_honeycomb"), (2, "K2_kagome"), (1, "K1_maple_leaf"))
    )
    check("counting_eta reproduces the 3D ladder to 1e-12", ok)
    check("counting_eta reproduces the 2D family to 1e-12", ok2)

    # spectral bound
    kf = spectral_K_max(*named_structure("FCC"), nq=nq)
    kt = spectral_K_max(*triangular_lattice(), nq=nq)
    check("spectral_K_max(FCC) == 4  (band [-4,12])", kf == 4, f"got {kf}")
    check("spectral_K_max(triangular) == 3  (band [-3,6]; 0.4534 rung retracted)",
          kt == 3, f"got {kt}")

    # the four 3D modular rules
    for K in (4, 3, 2, 1):
        ref = golden["fcc_depletion_ladder"][f"K{K}"]
        r = verify_modular_K(K, dim=3)
        ok = (
            r["uniform"] and r["independent"] and r["tangential"] and r["congruent"]
            and r["z"] == ref["z"] and r["p_vector"] == ref["p_vector"]
            and r["eta_err"] < 1e-12
        )
        check(f"3D K={K} [{r['label']}]: uniform/independent/tangential/congruent, "
              f"p={ref['p_vector']}", ok)

    # the 2D family, kagome exact
    for K, key in ((3, "K3_honeycomb"), (2, "K2_kagome"), (1, "K1_maple_leaf")):
        ref = golden["family_2d"][key]
        r = verify_modular_K(K, dim=2)
        ok = (
            r["uniform"] and r["independent"] and r["tangential"] and r["congruent"]
            and r["z"] == ref["z"] and abs(r["eta"] - ref["phi"]) < 1e-12
        )
        check(f"2D K={K} ({key}): phi={ref['phi']:.10f}, z={ref['z']}", ok,
              f"got phi={r['eta']:.10f}")

    # independence forced
    d = divacancy_spread()
    check("divacancy spread == sqrt2 - 1; isolated vacancy == 0",
          abs(d["divacancy"] - (np.sqrt(2) - 1)) < 1e-9 and d["isolated"] < 1e-9,
          f"got {d['divacancy']:.6f} / {d['isolated']:.2e}")

    # degeneracy counts, exact
    for key, N, K in [("tri_4x4_K2", 4, 2), ("tri_6x6_K3", 6, 3), ("tri_6x6_K2", 6, 2)]:
        ref = golden["degeneracy_counts"][key]
        r = enumerate_uniform_K(triangular_torus(N), K)
        ok = (
            r["complete"]
            and r["solutions"] == ref["solutions"]
            and r["orbits"] == ref["orbits"]
            and r["cosets"] == ref["cosets"]
        )
        check(f"{key}: solutions={ref['solutions']}, orbits={ref['orbits']}, "
              f"cosets={ref['cosets']}", ok,
              f"got {r['solutions']}/{r['orbits']}/{r['cosets']}")

    # every cell of every 6x6 K=2 solution exactly tangential at sqrt3 pi/8
    tor = triangular_torus(6)
    r = enumerate_uniform_K(tor, 2)
    target = np.sqrt(3) * np.pi / 8
    worstQ = worstS = 0.0
    for rep in r["orbit_representatives"]:
        for c in solution_cells(rep, tor):
            worstQ = max(worstQ, abs(c.Q_iso - target))
            worstS = max(worstS, float(c.face_distances.max() - c.face_distances.min()))
    check("every cell of every 6x6 K=2 solution tangential with Q = sqrt3 pi/8",
          worstQ < 1e-9 and worstS < 1e-9,
          f"max |Q-target|={worstQ:.1e}, max face spread={worstS:.1e}")

    passed = all(ok for _, ok, _ in checks)
    if verbose:
        print(f"  GATE G0b {'PASSED' if passed else 'FAILED'} "
              f"({sum(ok for _, ok, _ in checks)}/{len(checks)})")
    return {"gate": "G0b", "passed": passed, "checks": checks}


def main() -> int:
    print("GATE G0b: the ladder derivation (counting identity + spectral bound)")
    return 0 if run()["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

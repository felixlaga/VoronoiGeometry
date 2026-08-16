"""Gate G0 -- reproduce every geometric value in REFERENCE_VALUES.json.

Pass criteria (task list T1 / spec Sec. 3): all tabulated values reproduced to
``tolerance_exact_geometry``; tangential Bravais lattice count == 3; the
depletion ladder == {4, 5, 7, 13}.  Blocking; no tolerance may be widened.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from ..geometry.depletion import depleted_fcc, scan_depletions
from ..geometry.lattices import (
    enumerate_tangential_lattices,
    named_structure,
    no_fourteen_faced_solution,
    tetrahedral_octahedral_honeycomb,
    vacancy_shells,
)
from ..geometry.tangential import analyse

REPO = Path(__file__).resolve().parents[3]


def run(*, kmax: int = 13, verbose: bool = True) -> dict:
    golden = json.loads((REPO / "REFERENCE_VALUES.json").read_text())
    tol = golden["tolerance_exact_geometry"]
    checks: list[tuple[str, bool, str]] = []

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), detail))
        if verbose:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")

    # --- tangential Bravais lattices ------------------------------------- #
    ref3 = golden["tangential_bravais_3d"]
    for name, key in [("SC", "simple_cubic"), ("simple hexagonal c=a", "simple_hexagonal_ca"), ("FCC", "fcc")]:
        r = analyse(name, *named_structure(name))
        ok = (
            r["tangential"]
            and abs(r["eta"] - ref3[key]["eta"]) < tol
            and abs(r["Q_iso"] - ref3[key]["eta"]) < tol
            and int(np.unique(r["n_faces"])[0]) == ref3[key]["faces"]
            and int(np.unique(r["z"])[0]) == ref3[key]["z"]
        )
        check(f"{key}: eta=Q_iso={ref3[key]['eta']:.10f}, faces={ref3[key]['faces']}, z={ref3[key]['z']}", ok,
              f"got eta={r['eta']:.10f}")

    enum = enumerate_tangential_lattices()
    check(f"Bravais enumeration count == {ref3['count_is_complete']}",
          len(enum) == ref3["count_is_complete"],
          f"got {len(enum)} ({enum[0]['n_systems']} linear systems solved)")
    n14 = no_fourteen_faced_solution()
    got_norms = np.unique(np.round(n14["class_norms_symmetric"], 12))
    check("no 14-face solution (p_ij=1/6 gives norms 1/2, 2/3; 4t=3t=2*sum p forces t=0)",
          n14["excluded"] and len(got_norms) == 2
          and abs(got_norms[0] - 0.5) < 1e-9 and abs(got_norms[1] - 2 / 3) < 1e-9)

    # --- ladder ----------------------------------------------------------- #
    hits = scan_depletions(kmax)
    check("depletion ladder scan == {4,5,7,13}", {h["k"] for h in hits} == {4, 5, 7, 13},
          f"got {sorted(h['k'] for h in hits)}")
    for k, K in ((4, 4), (5, 3), (7, 2), (13, 1)):
        ref = golden["fcc_depletion_ladder"][f"K{K}"]
        lat, basis = depleted_fcc(k)
        r = analyse(f"FCC-1in{k}", lat, basis)
        pvec = sorted(len(loop) for loop in r["cells"][0].face_vertices)
        ok = (
            r["tangential"] and r["congruent"]
            and abs(r["eta"] - ref["eta"]) < tol
            and int(np.unique(r["z"])[0]) == ref["z"]
            and pvec == ref["p_vector"]
        )
        check(f"K{K} (k={k}): eta={ref['eta']:.10f}, z={ref['z']}, p-vector {ref['p_vector']}", ok,
              f"got eta={r['eta']:.10f}, p={pvec}")

    # exact V, S in integer coordinates for the two certified endpoints
    for K, k in ((4, 4), (1, 13)):
        ref = golden["fcc_depletion_ladder"][f"K{K}"]
        if "V_exact_intcoords" not in ref:
            continue
        lat, basis, _ = _integer_realisation(K)
        r = analyse(f"K{K} int", lat, basis)
        V_exp = eval(ref["V_exact_intcoords"].replace("sqrt2", "np.sqrt(2)"), {"np": np})
        S_exp = eval(ref["S_exact"].replace("sqrt2", "np.sqrt(2)"), {"np": np})
        c = r["cells"][0]
        check(f"K{K} exact cell V={ref['V_exact_intcoords']}, S={ref['S_exact']} (integer coords)",
              abs(c.V - V_exp) < 1e-10 and abs(c.S - S_exp) < 1e-10,
              f"got V={c.V:.12f} S={c.S:.12f}")

    # --- inadmissible ------------------------------------------------------ #
    for name, key in [("BCC", "bcc"), ("diamond", "diamond")]:
        ref = golden["inadmissible"][key]
        r = analyse(name, *named_structure(name))
        ok = (
            not r["tangential"]
            and abs(r["eta"] - ref["eta"]) < tol
            and abs(r["Q_iso"] - ref["Q_iso"]) < 1e-6
            and int(np.unique(r["n_faces"])[0]) == ref["faces"]
            and int(np.unique(r["z"])[0]) == ref["z"]
        )
        check(f"{key}: inadmissible, eta={ref['eta']:.7f} != Q_iso={ref['Q_iso']:.7f}", ok)

    # --- corrections ------------------------------------------------------- #
    v = vacancy_shells("BCC")
    check("BCC vacancy: 8 nearest neighbours (paper says 4); volume balances to 1 V0",
          v["shells"][0]["count"] == 8 and v["balances"],
          f"redistributed={v['redistributed_volume']:.9f}")
    check("BCC vacancy eta_correct = 8 sqrt3 pi/71 = 0.613115",
          abs(v["shells"][0]["eta_local"] - 8 * np.sqrt(3) * np.pi / 71) < 1e-9)
    t = tetrahedral_octahedral_honeycomb()
    check("tet-oct: 8 adjacent tetrahedra, cell = rhombic dodecahedron, eta = 2 sqrt3 pi/27",
          t["reduces_to_rhombic_dodecahedron"]
          and abs(t["eta_corrected"] - 2 * np.sqrt(3) * np.pi / 27) < 1e-12
          and abs(t["fill_fraction_original"] - 0.75) < 1e-12)

    passed = all(ok for _, ok, _ in checks)
    if verbose:
        print(f"  GATE G0 {'PASSED' if passed else 'FAILED'} ({sum(ok for _, ok, _ in checks)}/{len(checks)})")
    return {"gate": "G0", "passed": passed, "checks": checks}


def _integer_realisation(K: int):
    from ..geometry.coloring import realise_modular_K

    return realise_modular_K(K, dim=3)


def main() -> int:
    print("GATE G0: geometry vs REFERENCE_VALUES.json")
    return 0 if run()["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

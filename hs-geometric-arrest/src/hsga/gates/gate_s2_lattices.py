"""Gate G2 -- the radical-Voronoi pipeline on perfect lattices.

Pass criteria (spec Sec. 3, tolerances from REFERENCE_VALUES ``engine_gates``):
exact ``Q_iso`` (to ``pipeline_Q_max_abs_err_on_perfect_lattices``) on every
tangential structure, ``f_c = 1`` there, and ``f_c = 8/14`` for BCC.

This runs the SAME code path used on simulated configurations
(``analysis.voronoi.config_cells``), so it validates the pipeline, not merely
a lattice special case.  No simulation input; runs in seconds.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from ..analysis.voronoi import config_cells
from ..geometry.coloring import realise_modular_K
from ..geometry.lattices import named_structure
from ..geometry.tangential import voronoi_cell

REPO = Path(__file__).resolve().parents[3]


def _structures():
    for name in ("FCC", "HCP", "SC", "simple hexagonal c=a", "FCC-1in4", "BCC", "diamond"):
        lat, bas = named_structure(name)
        yield name, lat, bas
    for K in (3, 2, 1):
        lat, bas, _ = realise_modular_K(K, dim=3)
        yield f"K{K}", lat, bas


def _replicate(lattice, basis, reps):
    import itertools

    shifts = np.array(list(itertools.product(range(reps), repeat=3)), float) @ lattice
    return (basis[None, :, :] + shifts[:, None, :]).reshape(-1, 3), reps * lattice


def run(*, reps: int = 4, verbose: bool = True) -> dict:
    golden = json.loads((REPO / "REFERENCE_VALUES.json").read_text())
    tol_Q = golden["engine_gates"]["pipeline_Q_max_abs_err_on_perfect_lattices"]
    fc_bcc = golden["engine_gates"]["bcc_contact_fraction"]
    ref_eta = {
        "FCC": golden["tangential_bravais_3d"]["fcc"]["eta"],
        "HCP": golden["tangential_bravais_3d"]["fcc"]["eta"],
        "SC": golden["tangential_bravais_3d"]["simple_cubic"]["eta"],
        "simple hexagonal c=a": golden["tangential_bravais_3d"]["simple_hexagonal_ca"]["eta"],
        "FCC-1in4": golden["fcc_depletion_ladder"]["K4"]["eta"],
        "K3": golden["fcc_depletion_ladder"]["K3"]["eta"],
        "K2": golden["fcc_depletion_ladder"]["K2"]["eta"],
        "K1": golden["fcc_depletion_ladder"]["K1"]["eta"],
    }
    rows, ok = [], True
    for name, lat, bas in _structures():
        sites, super_lat = _replicate(lat, bas, reps)
        d = np.linalg.norm(sites - sites[0], axis=1)
        sigma = float(np.sort(d[d > 1e-9])[0])
        rad = np.full(len(sites), sigma / 2.0)

        if np.allclose(super_lat, np.diag(np.diag(super_lat))) and \
                np.ptp(np.diag(super_lat)) < 1e-12:
            L = float(super_lat[0, 0])
            pos = sites - L * np.floor(sites / L)
            cells = config_cells(pos, rad, L)
            path = "config_cells (cubic box)"
        else:
            cells = [voronoi_cell(b, lat, bas, radii=np.full(len(bas), sigma / 2))
                     for b in bas]
            path = "voronoi_cell (general lattice)"

        Q = float(np.mean([c.Q_iso for c in cells]))
        fc = float(np.mean([c.f_c for c in cells]))
        ref = ref_eta.get(name)
        if ref is None:                              # BCC, diamond
            expect_fc = fc_bcc if name == "BCC" else 4.0 / 16.0
            good = abs(fc - expect_fc) < 1e-12
            err = float("nan")
        else:
            err = abs(Q - ref)
            good = err < tol_Q and abs(fc - 1.0) < 1e-12
            expect_fc = 1.0
        ok &= good
        rows.append({"name": name, "Q_iso": Q, "reference": ref, "err": err,
                     "f_c": fc, "f_c_expected": expect_fc, "path": path,
                     "passed": good})
        if verbose:
            r = "        --" if ref is None else f"{ref:.9f}"
            e = "      --" if ref is None else f"{err:.1e}"
            print(f"  {name:22s} Q_iso={Q:.9f} ref={r} err={e} "
                  f"f_c={fc:.6f} (expect {expect_fc:.6f}) [{path}] "
                  f"{'ok' if good else 'FAIL'}")
    if verbose:
        print(f"  GATE G2 {'PASSED' if ok else 'FAILED'}")
    return {"gate": "G2", "passed": bool(ok), "rows": rows}


def main() -> int:
    print("GATE G2: radical-Voronoi pipeline on perfect lattices")
    return 0 if run()["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

"""Gate T8 -- the persistent-topology classifier and the refscore machinery.

Pass criteria (task list T8): the classifier reproduces every reference
p-vector on perfect lattices --

    3^8 (K=4), 3^2 4^7 (K=3), 4^10 (K=2), 4^11 (K=1),
    rhombic dodecahedron 4^12 (FCC), cube 4^6 (SC),
    hexagonal prism 4^6 6^2 (simple hexagonal)

-- with pairwise-distinct topo hashes, persistence 1.0 across the frozen
grid, refscore ~ 0 against a cell's own reference (and marked as such), and
``P_wrap^(k) = True`` on the perfect K-structures, where every cell is marked
and the marked network trivially wraps.
"""

from __future__ import annotations

import sys

import numpy as np

from ..analysis.pwrap import p_wrap
from ..analysis.refscore import load_frozen, mark_cells, refscores
from ..analysis.topology import p_vector, persistence, topo_hash
from ..geometry.coloring import realise_modular_K
from ..geometry.lattices import named_structure
from ..geometry.tangential import voronoi_cell

#: golden p-vectors on perfect lattices
TARGETS = {
    "K4": [3] * 8,
    "K3": [3, 3] + [4] * 7,
    "K2": [4] * 10,
    "K1": [4] * 11,
    "fcc": [4] * 12,
    "simple_cubic": [4] * 6,
    "simple_hexagonal_ca": [4] * 6 + [6] * 2,
}


def _reference_cells():
    out = {}
    for name in TARGETS:
        if name.startswith("K"):
            lat, basis, _ = realise_modular_K(int(name[1]), dim=3)
        else:
            key = {"simple_cubic": "SC", "simple_hexagonal_ca": "simple hexagonal c=a",
                   "fcc": "FCC"}[name]
            lat, basis = named_structure(key)
        out[name] = voronoi_cell(basis[0], lat, basis)
    return out


def run(*, verbose: bool = True) -> dict:
    checks = []

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), detail))
        if verbose:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")

    cells = _reference_cells()
    fz = load_frozen()

    # p-vectors
    for name, target in TARGETS.items():
        pv = list(p_vector(cells[name]))
        check(f"p-vector {name} == {target}", pv == target, f"got {pv}")

    # hashes pairwise distinct
    hashes = {name: topo_hash(cells[name]) for name in TARGETS}
    distinct = len(set(hashes.values())) == len(hashes)
    check("topo hashes pairwise distinct across the 7 reference classes", distinct,
          f"{len(set(hashes.values()))}/{len(hashes)} distinct")

    # persistence 1.0 on perfect cells over the frozen grid
    g = fz["persistence_grid"]
    worst = min(persistence(c, g["g_grid"], g["a_grid"]) for c in cells.values())
    check("persistence == 1.0 on every perfect reference cell", worst == 1.0,
          f"min {worst:.3f}")

    # refscore: each perfect K-cell scores ~0 against its own reference and is marked
    for K in (4, 3, 2, 1):
        lat, basis, _ = realise_modular_K(K, dim=3)
        cs = [voronoi_cell(b, lat, basis) for b in basis]
        radii = np.full(len(cs), cs[0].r_in)
        scores, names = refscores(cs, radii, 3)
        own = scores[:, names.index(f"K{K}")]
        marked = mark_cells(scores, names, f"K{K}")
        check(f"refscore(K{K} cells vs K{K}) ~ 0 and all marked",
              float(own.max()) < 1e-9 and marked.all(),
              f"max s = {own.max():.2e}")

    # P_wrap on a perfect K4 structure embedded on the FCC integer torus:
    # side 4 is even and a multiple of the mod-2 rule (landmine 5)
    side = 4
    sites = np.array([(x, y, z) for x in range(side) for y in range(side)
                      for z in range(side) if (x + y + z) % 2 == 0], float)
    rule = lambda x, y, z: x % 2 == 0 and y % 2 == 0 and z % 2 == 0
    occ = np.array([p for p in sites if not rule(int(p[0]), int(p[1]), int(p[2]))])
    from ..analysis.voronoi import config_cells

    radii = np.full(len(occ), np.sqrt(2.0) / 2.0)
    cs = config_cells(occ, radii, float(side))
    pw = p_wrap(cs, occ, radii, float(side), "K4")
    check("P_wrap^(K4) wraps on the perfect K4 structure (all cells marked)",
          pw["wraps"] and pw["marked_fraction"] == 1.0,
          f"marked={pw['marked_fraction']:.2f}, bonds={pw['n_bonds']}")

    passed = all(ok for _, ok, _ in checks)
    if verbose:
        print(f"  GATE T8 {'PASSED' if passed else 'FAILED'} "
              f"({sum(ok for _, ok, _ in checks)}/{len(checks)})")
    return {"gate": "T8", "passed": passed, "checks": checks}


def main() -> int:
    print("GATE T8: persistent-topology classifier + refscore machinery")
    return 0 if run()["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

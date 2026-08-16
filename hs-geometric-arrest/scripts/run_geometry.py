#!/usr/bin/env python3
"""Run the geometry gates G0 and G0b and write ``results/geometry.md``.

Both gates are blocking: any failure exits non-zero and the pipeline stops.
Everything here is theory -- no simulation input of any kind.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.gates import gate_g0_geometry, gate_g0b_derivation  # noqa: E402


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Gates G0 + G0b")
    p.add_argument("--kmax", type=int, default=13, help="depletion scan bound (20 for the DEBT claim)")
    p.add_argument("--nq", type=int, default=200_000, help="Bloch-band q samples")
    p.add_argument("--results", default=str(REPO / "results"))
    a = p.parse_args(argv)

    print("GATE G0: geometry vs REFERENCE_VALUES.json")
    g0 = gate_g0_geometry.run(kmax=a.kmax)
    print()
    print("GATE G0b: the ladder derivation")
    g0b = gate_g0b_derivation.run(nq=a.nq)

    out = Path(a.results)
    out.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Geometry gates G0 / G0b",
        "",
        "Every check compares against `REFERENCE_VALUES.json` (immutable; "
        "tolerances final). No simulation input.",
        "",
    ]
    for g in (g0, g0b):
        lines.append(f"## Gate {g['gate']} — {'PASSED' if g['passed'] else 'FAILED'}")
        lines.append("")
        lines.append("| check | result | detail |")
        lines.append("|---|---|---|")
        for name, ok, detail in g["checks"]:
            lines.append(f"| {name} | {'ok' if ok else 'FAIL'} | {detail} |")
        lines.append("")
    (out / "geometry.md").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {out / 'geometry.md'}")

    return 0 if (g0["passed"] and g0b["passed"]) else 1


if __name__ == "__main__":
    sys.exit(main())

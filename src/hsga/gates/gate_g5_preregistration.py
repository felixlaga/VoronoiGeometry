"""Gate G5 -- pre-registration of the refscore (spec Sec. 3).

Pass criteria: ``refscore_frozen.json`` (and ``REFERENCE_VALUES.json``) are
byte-identical to the versions in the repository's FIRST commit, and no later
commit has ever touched them.  Simulation data lives outside version control
(``data/`` is gitignored), so "committed before the first dynamical dataset it
is applied to" reduces to: frozen at commit zero, untouched since -- which is
exactly what this gate proves from the git history.

The gate also confirms that ``analysis/refscore.py`` actually reads the frozen
file (not a copy), so the proof covers the code path in use.
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
FILES = ("refscore_frozen.json", "REFERENCE_VALUES.json")


def _git(*args) -> str:
    r = subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed: {r.stderr.strip()}")
    return r.stdout


def run(*, verbose: bool = True) -> dict:
    checks = []

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), detail))
        if verbose:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")

    first = _git("rev-list", "--max-parents=0", "HEAD").strip().splitlines()[0]
    for name in FILES:
        blob_first = _git("show", f"{first}:{name}")
        current = (REPO / name).read_text()
        same = blob_first == current
        h = hashlib.sha256(current.encode()).hexdigest()[:16]
        check(f"{name} identical to the first commit ({first[:8]})", same,
              f"sha256 {h}")
        touching = _git("log", "--format=%H", "--", name).strip().splitlines()
        check(f"{name} touched by exactly one commit", len(touching) == 1,
              f"{len(touching)} commits")

    # the refscore module reads THIS file
    from ..analysis import refscore

    loaded = refscore.load_frozen()
    import json

    check("analysis.refscore reads the frozen file (values match on disk)",
          loaded == json.loads((REPO / "refscore_frozen.json").read_text()))

    passed = all(ok for _, ok, _ in checks)
    if verbose:
        print(f"  GATE G5 {'PASSED' if passed else 'FAILED'}")
    return {"gate": "G5", "passed": passed, "checks": checks, "first_commit": first}


def main() -> int:
    print("GATE G5: refscore pre-registration (frozen at commit zero, untouched since)")
    return 0 if run()["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

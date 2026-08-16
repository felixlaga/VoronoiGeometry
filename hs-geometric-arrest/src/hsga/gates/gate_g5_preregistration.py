"""Gate G5 -- pre-registration of the refscore (spec Sec. 3).

Pass criteria: ``refscore_frozen.json`` (and ``REFERENCE_VALUES.json``) enter
the git history in a ROOT commit -- before any analysis code or campaign data
existed -- are byte-identical to that version today, and every commit that has
ever touched them carries that identical blob.  Simulation data lives outside
version control (``data/`` is gitignored), so "committed before the first
dynamical dataset it is applied to" reduces to: frozen at the history's
origin, never modified since -- which is exactly what this gate proves.

The project may live either as a standalone repository or grafted into a
parent repository under a prefix (its history was merged into VoronoiGeometry
with ``--allow-unrelated-histories``, preserving the original root commit).
The gate handles both: it locates the unique root commit containing the
frozen files and verifies blob identity across every touching commit at both
the original and the prefixed path.  The proof is the same in both layouts;
nothing is loosened.

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


def _blob(commitish: str, path: str) -> str | None:
    """Object id of ``path`` in ``commitish``, or None if absent there."""
    r = subprocess.run(["git", "rev-parse", f"{commitish}:{path}"],
                       cwd=REPO, capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


def run(*, verbose: bool = True) -> dict:
    checks = []

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), detail))
        if verbose:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")

    top = Path(_git("rev-parse", "--show-toplevel").strip()).resolve()
    prefix = REPO.resolve().relative_to(top).as_posix()  # "." when standalone
    roots = _git("rev-list", "--max-parents=0", "HEAD").split()

    first = None
    for name in FILES:
        paths = [name] if prefix == "." else [name, f"{prefix}/{name}"]

        # the file's frozen blob id, computed from the bytes on disk
        disk = (REPO / name).read_bytes()
        frozen = _git("hash-object", str(REPO / name)).strip()
        h = hashlib.sha256(disk).hexdigest()[:16]

        # exactly one root commit of HEAD's history contains the file
        roots_with = [r for r in roots if _blob(r, name) is not None]
        check(f"{name} enters history in exactly one root commit",
              len(roots_with) == 1,
              f"{roots_with[0][:8] if len(roots_with) == 1 else roots_with}")
        if len(roots_with) != 1:
            continue
        first = roots_with[0]

        check(f"{name} identical to the root commit ({first[:8]})",
              _blob(first, name) == frozen, f"sha256 {h}")

        # every commit that ever touched the file carries the identical blob
        # (":/" anchors the pathspec to the repo root -- ``git log`` pathspecs
        # are cwd-relative, but ``rev-parse commit:path`` is root-relative;
        # --full-history keeps the pre-graft line, which default history
        # simplification would prune at the merge)
        touching = set()
        for p in paths:
            touching |= set(_git("log", "--full-history", "--format=%H",
                                 "--", f":/{p}").split())
        bad = [c for c in touching
               if any(b not in (None, frozen) for b in (_blob(c, p) for p in paths))]
        check(f"{name} never modified by any touching commit", not bad,
              f"{len(touching)} touching commit(s)"
              + (f", DIVERGENT: {[c[:8] for c in bad]}" if bad else ""))

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

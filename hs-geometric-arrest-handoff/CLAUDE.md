# CLAUDE.md — hs-geometric-arrest

You are building the repository specified in this bundle. Read, in order:
1. `CLAUDE_CODE_TASKS.md` — the ordered task list T0-T12 with gates and stop
   conditions. Execute it in order. Stop at every gate; never proceed on failure.
2. `IMPLEMENTATION_SPEC.md` — the normative module-by-module specification.
3. `REFERENCE_VALUES.json` — golden numbers. IMMUTABLE. Tolerances there are
   final; code adapts to this file, never the reverse.
4. `docs/paper.tex` Sections II-IV — the settled geometry you must reproduce,
   not re-derive.

## Non-negotiable rules
- No invented data, placeholders, or synthetic fallbacks. Failed fits are
  recorded as failures and propagate NaN.
- Gates are blocking. No --force flags, no widened tolerances, no skipped
  stages. A gate failure halts with a non-zero exit code.
- Never modify REFERENCE_VALUES.json or refscore_frozen.json after T0.
- Swap Monte Carlo moves are OFF in every production run, always, logged.
- Every reported eta_a carries the fit-window systematic, not only the
  statistical error.
- Repository documentation matches the repository exactly. Deferred work goes
  in DEBT.md, never in TODO comments.
- Check actively for bloat, shortcuts, task cheating, and half-completed work.

## Directory map of this bundle
- `prototypes/` — reference implementations from the research phase. They are
  CORRECT IN OUTPUT but carry known warts listed in CLAUDE_CODE_TASKS.md
  ("Status of the provided files"). Port them cleanly into src/hsga/; do not
  copy verbatim. `exact_certificates.py` goes into tests/ nearly as-is.
- `data/pilot3d/`, `data/pilot2d/` — the actual pilot configurations (N=864 3D,
  N=504 2D). Use them to validate ported analysis code against the numbers in
  docs/RESULTS_pilot.md before generating new data.
- `docs/` — paper, explainer, review verdict, pilot results, alternative-routes
  analysis. `docs/archive/` is superseded material kept for history.

## Landmines (each of these failed for real once; do not rediscover)
See "Landmines" in CLAUDE_CODE_TASKS.md: HNF ranges, coset wrapping, the
Nelder-Mead trap, torus-modulus matching, swap-off, interior-point requirement,
fit-window scans, revtex fallback.

## Supervision
Autonomous through T5. Stop and ask Felix at: end of T4, the T6 decision node,
and before any large compute at T9.

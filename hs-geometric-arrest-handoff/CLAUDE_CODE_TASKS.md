# Ordered task list for Claude Code — hs-geometric-arrest

Read this file first, then `IMPLEMENTATION_SPEC.md` (normative), then
`REFERENCE_VALUES.json` (golden numbers; tolerances there are final).
`SPEC_3D_hard_sphere_arrest.md` is SUPERSEDED by IMPLEMENTATION_SPEC.md —
where they conflict, IMPLEMENTATION_SPEC wins; treat the old file as history.

## Status of the provided files

(In this bundle, prototypes live under `prototypes/`, documents under `docs/`,
and the pilot data under `data/`.)

Prototypes, reference-only, known warts — port, don't copy:
- `voronoi_ground_states.py`, `enumerate_3d_ground_states.py`: correct, verified;
  clean API when porting to `src/hsga/geometry/`.
- `exact_certificates.py`: the sympy exact-arithmetic certificates. Port into
  the test suite; these are the strongest tests in the project.
- `verify_review.py`: correct results, messy code (dead branches in section B);
  reimplement, keep the printed reference outputs.
- `hsmc.c`: works, zero-overlap audited, but grew by patching; rewrite cleanly
  from the spec (Sec 2.4), keep the init/anneal/audit logic exactly.
- `hsmc2d.c`: works; contains dead variables from a fixed lattice-geometry bug
  (nx must be round(sqrt(3)*ny) for a near-square box). Clean up.
- `analysis.py`: contains a disabled loader line (`pass`) from a hotfix; the
  frame reader that works is the plain-Python one. Rewrite.
- `percolation.py`, `perc2d.py`: union-find with displacement tracking is
  correct and validated (recovers RCP to 1%); port as-is with tests.

## Landmines already hit once — do not rediscover them

1. HNF sublattice ranges are `0 <= d,e < a`, `0 <= f < b`; wrong ranges
   silently enumerate an incomplete set.
2. Coset representatives MUST be wrapped into the primitive cell before any
   periodic Voronoi computation (unwrapped reps silently return non-tangential).
3. The tangential-lattice search MUST be the linear Selling solve, never a
   numerical minimisation (Nelder-Mead misses the boundary strata where SC and
   simple-hex live — this failed for real).
4. Swap moves OFF in production, always; log the switch-off.
5. Modular-rule tori must have side = multiple of the modulus, or uniformity
   silently breaks (this failed for real at K=1 on a period-4 torus).
6. scipy HalfspaceIntersection needs a strictly interior point; the site itself
   works for tangential cells.
7. Never report eta_a without the fit-window scan; the window systematic was
   10x the statistical error.
8. revtex4-2 may be missing; paper_preview.tex shows the article-class fallback.

## Ordered tasks (stop at every gate; do not proceed on failure)

T0  Repo scaffold per IMPLEMENTATION_SPEC Sec 1; commit REFERENCE_VALUES.json
    and refscore_frozen.json (G5 pre-registration) in the FIRST commit.
T1  geometry/tangential.py + lattices.py + tests.
    GATE G0: all values in REFERENCE_VALUES.json reproduced; Bravais count == 3.
T2  geometry/coloring.py (counting_eta, bloch bands, spectral_K_max, modular
    rules, enumerate_uniform_K, coset test) + geometry/depletion.py.
    GATE G0b: ladder to 1e-12; K_max(FCC)=4, K_max(tri)=3; four modular rules
    uniform/independent/tangential/congruent; kagome exact; divacancy spread
    sqrt2-1; degeneracy_counts block reproduced exactly.
T3  Port exact_certificates.py into tests/ (sympy). GATE: all five certificates
    pass symbolically (identity == 0, not < tol).
T4  THEORY RUN (no simulation): enumerate_uniform_K on growing tori
    (2D: 6,9,12; larger if feasible; 3D K=4 small tori). Fit log(orbits) vs area.
    Deliverable: results/degeneracy.md. This can produce a standalone result.
T5  engine/hsmc2d.c clean rewrite + driver. GATE G1-2D: 2D EOS sanity + zero
    overlaps at phi=0.80.
T6  2D REPLICATION CAMPAIGN (the decision node): N=2048, >=350 configs/point,
    dense grids at 0.6046 / 0.6802 / 0.7773 per spec Sec 4. Analysis:
    eps*(phi) vs smooth null, P_wrap^(k), refscores, submode cusp.
    DECISION: feature at 0.7773 -> proceed to 3D. Nothing -> STOP 3D; write the
    negative-result section; the geometry paper + T4 stand alone.
T7  engine/hsmc.c rewrite + gates G1 (EOS < 1% at eta 0.30-0.40, cusp-aware
    g(sigma+) fit), G2 (pipeline on perfect lattices), G3 (equilibration),
    G4 (finite size).
T8  analysis/topology.py + refscore.py + pwrap.py + isoconfig.py per spec
    Sec 2.5b-e. GATE: classifier reproduces all reference p-vectors.
T9  3D campaign per spec Sec 4 (only if T6 passed).
T10 Baseline comparison (tetrahedrality etc., held-out). The hypothesis fails
    de Graaf's own standard if refscores add no predictive power.
T11 Update paper.tex Secs VI-VII with campaign results; rebuild.
T12 DEBT.md audit: every open item either closed with a commit hash or still
    listed. No silent disappearances.

## Supervision points for Felix

Claude Code should run T0-T5 autonomously. Stop for human review at: end of T4
(is the degeneracy result worth a standalone note?), the T6 decision, and
before committing cluster time at T9. Never let it widen a tolerance, touch
REFERENCE_VALUES.json, or edit refscore_frozen.json after T0.

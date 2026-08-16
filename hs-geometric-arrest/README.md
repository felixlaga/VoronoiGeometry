# hs-geometric-arrest

Tangential Voronoi cells as geometric ground states for hard spheres: the
exact geometry, the derivation of the depleted-lattice ladder, its degeneracy,
and the machinery to test whether the mechanism operates — built from the
handoff bundle (`../hs-geometric-arrest-handoff/`) along the task list
T0–T12, with every gate blocking and every golden number taken from the
immutable `REFERENCE_VALUES.json` committed in this repository's first commit.

The geometry is settled and verified three ways (float half-space cells, exact
sympy certificates, the derivation). The physics is deliberately **not**
claimed: the decision experiment (T6) is implemented but has not been run.

## The science in one page

**Theorem (Sec. II of the paper).** If every Voronoi cell of a monodisperse
packing is a tangential polytope with the particle as insphere — equivalently,
*every Voronoi neighbour is a contact* — then `S = d·V/R` exactly, so

- 3D: `eta = 36 pi V^2 / S^3 = Q_iso` (Wadell sphericity cubed),
- 2D: `phi = 4 pi A / P^2 = q`.

Proved as an exact symbolic identity in `tests/test_certificates.py`.

**Enumeration.** The seven Voronoi face-class norms are linear in the Selling
parameters, so tangentiality is a linear system (never an optimisation —
landmine 3): exactly **three** tangential Bravais lattices exist in 3D (SC
0.523599, simple hexagonal c=a 0.604600, FCC 0.740480), and no fourteen-faced
solution can exist (integer identity: `4t = 3t = 2 Σp` forces `t = 0`), so
BCC fails structurally.

**The ladder, derived.** Remove an independent vacancy set in which every
occupied site has exactly K vacant neighbours: double counting gives
`eta_K = z0/(z0+K) · eta_cp`, and the spectral bound (−K must lie in the
contact graph's Bloch band) caps K at 4 on FCC and 3 on the triangular
lattice. The 3D ladder is `K = 4,3,2,1` → 0.555360, 0.592384, 0.634698,
0.683520 (p-vectors 3⁸, 3²4⁷, 4¹⁰, 4¹¹); the 2D family is honeycomb (K=3,
0.604600), **kagome (K=2, 0.680175 — the novel falsifiable rung)** and
maple-leaf (K=1, 0.777343 = de Graaf's floret value). The once-claimed 2D
0.453450 rung is excluded by the bound (retracted). Independence is forced,
not assumed: a divacancy breaks tangentiality with a face-distance spread of
exactly `√2 − 1`.

**Degeneracy (T4, new here).** Exhaustive enumeration of all uniform-K states
on growing tori (`results/degeneracy.md`): the K=3 state is unique up to
translation at every commensurate size; K=1 is unique up to translation and
chirality; the K=2 class is degenerate with **mostly non-crystalline members —
every one exactly tangential at the ladder density** — but the count obeys
`S = 6·2^(side/2) − 8` exactly (side-16 predicted blind, then confirmed):
`ln S ~ √N`, a boundary/stacking law, **not** an extensive entropy. The same
formula matches both 3D FCC K=4 counts, where non-crystalline members also
exist. Consequence: the reference state is a degenerate family containing
disordered members, which weakens the configurational-entropy objection — but
sub-extensively, which does not answer it.

**Why no physics claim is made.** The 2D positive control (recorded in the
handoff, reproduced here by `scripts/validate_pilot.py`) showed the pipeline
at pilot statistics cannot resolve the known 2D transition, and calibrated the
requirement at ~350–500 configurations per density. The decision experiment is
therefore the **2D replication first** (T6): N=2016, fine 5·10⁻⁴ grids around
the three exact 2D targets, ~355 core-hours. A feature at 0.777343 validates
the pipeline; anything at the kagome rung is a new result either way; nothing
means the 3D campaign (T9, ~4800 core-hours) is not run at all.

## What has been executed in this repository (all of it on a laptop)

| stage | status | result |
|---|---|---|
| G0 geometry | **PASSED** 17/17 | every REFERENCE value reproduced; Bravais count 3; ladder {4,5,7,13}; exact K4/K1 cell V,S |
| G0b derivation | **PASSED** 16/16 | counting ladder to 1e-12; K_max(FCC)=4, K_max(tri)=3; 7 modular rules uniform/independent/tangential/congruent with golden p-vectors; divacancy √2−1; degeneracy counts exact; all 6×6-K2 cells at √3π/8 |
| T3 certificates | **PASSED** | five sympy certificates, identities `== 0`; exact V=8/3, S=8√2 (K4) and V=13/6, S=13√2/2 (K1) |
| T4 degeneracy | **executed** | `results/degeneracy.md`; boundary-law verdict above |
| G1-2D (T5) | **PASSED** | Z within 0.27–0.5% of Henderson at φ=0.55/0.60; binary R⁻¹=1.4 at φ=0.80 anneals to exactly zero overlap; ε* estimators agree |
| pilot validation | **PASSED** | every recorded ε* (2D control + 3D pilot) reproduced to all printed digits; RCP fit φ_m=0.8508, p=1.079 (sem-weighted, as recorded) |
| G2 pipeline | **PASSED** | Q exact to ≤5e-15 on ten structures (incl. K1–K3 via their modular rules); f_c=1 tangential, 8/14 BCC |
| G1 3D EOS | **PASSED** | 0.62% / 0.38% / 0.29% vs Carnahan–Starling at η=0.30/0.35/0.40 (criterion <1%), N=864, 9 seeds |
| G5 pre-registration | **PASSED** | `refscore_frozen.json` byte-identical to commit zero, touched once, and the module reads that file |
| T8 classifier | **PASSED** 14/14 | golden p-vectors 3⁸/3²4⁷/4¹⁰/4¹¹/4¹²/4⁶/4⁶6²; hashes distinct; persistence 1.0; refscore ~1e-28 on own class; P_wrap wraps on perfect K4 |
| EDMD engine | validated | selftest; energy/momentum drift <1e-9; determinism; virial Z within 3% of CS |

**Implemented but NOT executed** (this machine cannot carry the compute; every
entry has a make target and an honest cost estimate):
G3 (equilibration lengths to 1e5 sweeps), G4 (N to 10976),
**T6 the 2D replication campaign — the decision node** (~355 core-hours),
T9 the 3D campaign (~4800 core-hours; gated on T6 regardless),
T10 the isoconfigurational baseline test (needs campaign data; `--smoke`
exercises the full path), T11 the manuscript update (needs campaign results;
writing it without data would be invented data). Details in `DEBT.md`.

## Install and run

Python ≥ 3.10 with numpy, scipy, sympy, pytest; a C compiler.

```
pip install -e ".[dev]"
make engines          # hsmc (3D), hsmc2d (2D), edmd (event-driven, stage 5)
make test             # the test suite incl. the sympy certificates  (~1 min)
make gates-light      # G0, G0b, G2, G5, T8, G1-2D, G1 — everything a laptop can prove
make degeneracy       # regenerate results/degeneracy.md
make validate-pilot   # re-validate the analysis ports against the pilot data
make help             # the full map, light vs heavy
```

The heavy stages print their cost and refuse nothing — they are simply not run
here: `scripts/run_sweep.py --preset 2d-replication --dry-run`.

## Layout

```
REFERENCE_VALUES.json      golden numbers, immutable, first commit (G5-verified)
refscore_frozen.json       pre-registered refscore weights/filters (G5-verified)
src/hsga/
  geometry/tangential.py   the theorem; dim-generic radical cells with proved
                           pool sufficiency and face vertex loops
  geometry/lattices.py     named structures; Selling enumeration (linear solve);
                           integer 14-face exclusion; Sec.-V corrections
  geometry/depletion.py    HNF sublattices (landmine-1 ranges), wrapped cosets
                           (landmine 2), ladder scan, modular rule -> sublattice
  geometry/coloring.py     counting identity; Bloch bands + spectral K_max;
                           the 7 modular rules verified end to end; divacancy
                           lemma; exhaustive uniform-K enumeration on tori
  engine/hsmc.c, hsmc2d.c  NVT MC (3D/2D): melt, anneal-to-exactly-zero with
                           reheat + soft swaps, swap-off production (logged),
                           restart mode, audits, key=value logs
  engine/edmd.c            NVE event-driven MD (stage-5 confirmation only):
                           event calendar, exact elastic collisions, selftest
  engine/driver.py         RunSpec(dim), sweeps, provenance manifests,
                           unreachable-state-point recording, file formats
  analysis/percolation.py  eps* (exact sorted-bond + spec bisection), wrapping
                           union-find with relative displacements
  analysis/voronoi.py      per-config Q_iso, face gaps, single-scale test,
                           large-particle q submode weight
  analysis/eos.py          g(r) with a bin edge exactly at contact; cusp-aware
                           windowed contact value; CS (3D) and Henderson (2D)
  analysis/dynamics.py     MSD, D, eta_a with the mandatory window scan and
                           engine label
  analysis/topology.py     persistence filter, p-vectors, WL canonical hash
  analysis/refscore.py     the frozen s_k score; runtime-exact reference cells
  analysis/pwrap.py        P_wrap^(k), susceptibility, xi, finite-size crossings
  analysis/isoconfig.py    propensity, fixed baseline fields, held-out Delta-R2
  gates/                   G0, G0b, G1, G1-2D, G2, G3, G4, G5, T8
scripts/
  run_geometry.py          G0 + G0b -> results/geometry.md
  run_degeneracy.py        T4 -> results/degeneracy.md
  validate_pilot.py        ports vs recorded pilot numbers -> results/validation_pilot.md
  run_sweep.py             campaign generator (2d-replication / 3d-campaign / smokes)
  run_analysis.py          eps* + smooth null + look-elsewhere feature test +
                           P_wrap + submode + T6 decision text
  run_isoconfig.py         T10 driver (--smoke exercises the full path)
tests/                     the test suite; `-m "not slow"` skips the big scans
paper/                     manuscript + article-class preview (landmine 8) +
                           the checked-in compiled preview PDF
results/                   committed: geometry.md, degeneracy.md,
                           validation_pilot.md, smoke campaign artefacts
data/                      gitignored simulation output
```

## Design decisions worth knowing

- **Pool sufficiency is proved, not assumed**: half-spaces are added by
  distance and construction stops only when the next candidate provably cannot
  cut the cell; short pools raise, callers enlarge. The interior point for
  Qhull is the site itself (landmine 6).
- **Vertices are merged before topology**: Qhull duplicates intersection
  points wherever >d planes meet (every rhombic-dodecahedron vertex); without
  merging, every p-vector is wrong.
- **`eps*` is computed exactly** (sorted-bond insertion; the threshold is the
  bond completing the last wrap) and cross-checked against the spec's
  bisection; the port reproduced every recorded pilot value to all printed
  digits and recovers 2D RCP to 1%.
- **Feature claims must beat the smooth null**: the positive control showed a
  local-exponent "jump" that was entirely the `1/(phi_m - phi)` divergence.
  The campaign analysis fits the sem-weighted null, tests only standardised
  excursions inside the pre-declared fine windows, Bonferroni-corrects the
  look-elsewhere effect, bootstraps over replicas, and refuses to claim
  significance below 8 replicas. "Not tested" is reported as absent, never as
  negative.
- **Unreachable state points are data about jamming, not missing data**: runs
  that cannot anneal to zero overlap exit non-zero, write nothing, and appear
  in manifests and reports as unreachable; nothing interpolates over them.
- **Errors are over independent replicas only**; snapshots within a run are
  correlated and never counted as independent.
- **The refscore cannot be tuned**: weights, filter and threshold are frozen in
  the first commit and G5 proves it from the git history.

## Provenance

This project's git history was grafted into the `VoronoiGeometry` repository
(merge with `--allow-unrelated-histories`); the nine per-task commits and the
T0 root commit that freezes `REFERENCE_VALUES.json` / `refscore_frozen.json`
are preserved verbatim, and gate G5 verifies the pre-registration across the
graft (unique root commit, identical blob in every touching commit).

Every sweep writes a manifest with the exact command lines, engine source
hashes, compiler and library versions, and per-run logs. Prototype porting
notes, landmines and file status are in the handoff bundle's
`CLAUDE_CODE_TASKS.md`; everything listed there as a past failure mode has a
regression test here.

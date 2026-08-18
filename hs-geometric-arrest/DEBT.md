# DEBT — open items, deferred work, known limitations

Rules of this file (CLAUDE.md): deferred work lives here, never in TODO
comments. Every item is either open, or closed with the commit that closed it.
No silent disappearances.

## Theory

- **Tangential enumeration is complete for Bravais lattices only.** Non-lattice
  structures are covered only by the FCC sublattice scan to `k <= 20` and the
  uniform-K machinery; non-sublattice perfect coverings, HCP-derived depletions
  and general multi-orbit structures are not enumerated. A tangential structure
  outside these classes would not be found.
- ~~`(k-1) | 12` observed, not derived~~ **CLOSED at T2** (commit `0f6cc4d`):
  counting identity `eta_K = z0/(z0+K) eta_cp` + spectral bound
  (`geometry/coloring.py`), independence forced by tangentiality (divacancy
  spread exactly `sqrt2 - 1`). Gate G0b verifies all of it.
- **Whether uniform-K degeneracy is EXTENSIVE is now answered provisionally,
  not finally** (T4, `results/degeneracy.md`): at every accessible size the
  K=2-class count follows `S = 6·2^(side/2) - 8` exactly — a boundary/stacking
  law, `ln S ~ sqrt(N)`, carrying **no** finite entropy density. Open parts:
  the closed form is observed (side-16 predicted blind and confirmed), not
  derived; sizes stop at 16x16 (2D) and FCC side 6 (3D); the K=1 class is only
  checked at sides 7 and 14.
- **Existence of uniform-K states on HCP is open.** The spectral condition
  passes for all `K <= 4` (necessary, not sufficient). If they exist, density
  degenerates across stackings and the topology classifier
  (`analysis/topology.py`) becomes mandatory, not optional.
- **No closed-form 3D analogue of `q_r(n)` exists.** Any generalised face
  number is a documented choice; the default is to work with `Q_iso` directly.
- **Whether any of these structures is thermodynamically relevant at `T > 0` is
  open**, in 3D as in 2D. The T4 result (non-crystalline members exist, all
  exactly tangential, but the family is sub-extensive) weakens the
  configurational-entropy objection only marginally; it does not answer it.
- **REFERENCE_VALUES.json `family_2d.K1_maple_leaf.phi` decimal is off by
  ~4e-13 from its own exact string** `sqrt3*pi/7` (0.7773425846722069 quoted
  vs 0.77734258467180757...; found by the T3 sympy certificates, which prove
  the exact value symbolically). The file is immutable, so the decimal stays;
  every gate uses 1e-12 or looser and is unaffected. Fix it in the next
  revision of the golden file, outside this repository's rules.

## Observables

- **`eps*` is composition-independent but has NOT been shown to resolve the
  transition.** The 2D positive control was consistent with a featureless power
  law at 12 configurations per density; the calibrated requirement is ~350–500.
  This is why T6 exists and why it gates all 3D physics.
- **`f_c` at fixed tolerance is unmeasurable** (exact contact has measure zero)
  and its `f_c(delta)` repair is degenerate with the median gap
  (`analysis/voronoi.single_scale_collapse` measures the degeneracy). It must
  not return as a headline observable; the mechanistic observable is
  `P_wrap^(k)` on the pre-registered refscore.
- **MC dynamics is a Brownian proxy;** `eta_a` from it is not the EDMD `eta_a`.
  `analysis/dynamics.fit_eta_a` carries the engine label inside the returned
  record; the two are never interchangeable.
- **`g(sigma+)` linear extrapolation degrades above `eta ~ 0.4`** (4% at 0.45
  in the pilot). `analysis/eos.py` uses a cusp-aware `ln g` polynomial with a
  mandatory window scan; the spread is always reported as the systematic.
- **`eps* = C (eta_m - eta)^p` is the wrong functional form near contact
  percolation** (window spread 0.15 in the pilot). The campaign analysis fits
  it ONLY as the smooth null; `eta_m` is a nuisance parameter and is never
  quoted as a physical vanishing point. A free-volume-motivated form is needed
  before any vanishing point can be extracted.
- **Isoconfigurational field definitions are documented choices.**
  Tetrahedrality (Errington–Debenedetti over 4 nn), `psi6` over Voronoi
  neighbours, anisotropy `1 - r_in/R_circ`: fixed in `analysis/isoconfig.py`
  before any propensity data existed. None has been validated against real
  dynamical data (none exists yet — see Execution status).
- **The WL topo hash is an invariant, not a complete one.** Isomorphic face
  graphs always collide; non-isomorphic ones could in principle collide too,
  which would merge classes (conservative direction). The seven reference
  classes are pairwise separated (gate T8).

## Engine

- **The initialisation anneal extends the spec's literal wording.** Spec 2.4:
  set target radii, anneal `E = sum(sigma_ij - r_ij)^2` to exactly zero at
  decreasing `T`. A single monotonic cooling stalls above `eta ~ 0.6` (2D:
  `phi ~ 0.78`) with a few frustrated overlaps. Both engines therefore
  (a) reheat `T` to the current per-particle energy scale when the energy
  stalls and (b) include swap moves on the soft energy during the anneal. The
  end state is exactly what the spec requires — target radii (verified as a
  multiset, since swaps permute them), `E == 0` exactly, non-zero exit
  otherwise. Initialisation only; production is local-moves-only and logged.
- **The reachable density range is composition-dependent.** State points at or
  beyond their jamming density cannot be equilibrated by any initialisation
  protocol; such runs exit 3, write no configuration, and are recorded as
  `unreachable` in the sweep manifest. The analysis reports them as missing
  and never interpolates over them. Jamming ceilings have not been mapped for
  the campaign compositions.
- **Snapshot independence holds across replicas only.** Snapshots within a run
  are correlated; every error bar is taken over independent replicas, and
  `replicas x nsnap` is never quoted as an independent-sample count. No
  autocorrelation time is measured.
- **EDMD is Newtonian; the target experiments are Brownian.** The insensitivity
  of the arrest location to microscopic dynamics is assumed, not tested here.
  EDMD masses default to `m ∝ R^3`; `--equal-mass` changes `D` but should not
  move `eta_a` — an expectation, not a measurement.

## Execution status (the heavy stages are implemented, validated at smoke
scale, and deliberately NOT executed on this machine)

- **T0–T5 executed in full.** Gates G0 (17/17), G0b (16/16), the five T3 sympy
  certificates (exact identities), the T4 degeneracy run
  (`results/degeneracy.md`), gate G1-2D (Z within 0.27–0.5% of Henderson;
  zero-overlap audits at phi=0.80), and the pilot-data validation
  (`results/validation_pilot.md`: every recorded eps* reproduced to all
  printed digits; RCP fit phi_m=0.8508, p=1.079).
- **T7 partially executed:** G2 PASSED (Q exact to <=5e-15 on all ten
  structures incl. K1–K3; f_c = 8/14 for BCC), G1 PASSED at full statistics
  (0.62/0.38/0.29% vs Carnahan–Starling at eta=0.30/0.35/0.40). **G3 and G4
  implemented, NOT executed** (1e5-sweep equilibrations and N=10976 are hours
  of compute); `make gate-g3` / `make gate-g4` run them as specified.
- **T8 executed:** gate T8 PASSED 14/14. **G5 executed:** pre-registration
  proof from git history PASSED — re-verified after the project's history was
  grafted into the VoronoiGeometry repository (the gate is layout-aware: it
  finds the unique root commit containing the frozen files and demands the
  identical blob in every touching commit; same proof, standalone or nested).
- **T6 EXECUTED 2026-08-18** (4400/4400 runs, 7.1 h wall on the build
  machine, zero unreachable points, all audits zero;
  `results/campaign_2d/`). Outcome: at sem ~ 5e-5, eps* is locally LINEAR
  through all three fine windows, P_wrap on the frozen refscore is zero
  everywhere with flat chi and xi, and the apparent dynamical arrest
  drifts +0.041/decade with the observation-time criterion — no
  structural feature at 0.777343 or 0.680175. The pre-registered feature
  test MIS-FIRED (`report.md` says "FEATURE AT 0.777343"): its global
  power-law null is rejected at chi2/dof = 137 and the max-residual
  statistic tracks null misfit, proven by an equal-strength "feature" at
  the honeycomb window deep in the normal fluid. Full forensics in
  `results/campaign_2d/DIAGNOSTICS.md`. **Open: the smooth-null form must
  be replaced (free-volume-motivated) and the feature test re-specified
  before any re-test; the T6 decision (accept the negative / redesign /
  new observables) is Felix's supervision point. T9 stays blocked.**
- **T9 (3D campaign) NOT run** (~4800 core-hours) — gated on T6 in any case.
- **T10 NOT run** on real data (needs campaign output); the driver works end
  to end in `--smoke` mode.
- **T11 open:** `paper/paper.tex` Secs. VI–VII still describe the pilot;
  updating them requires campaign results that do not exist, and writing them
  without data would violate rule 1. `make paper` builds when a TeX
  installation is present (article-class fallback when revtex4-2 is missing —
  landmine 8); this machine has no TeX, so `paper/paper_preview.pdf` is the
  checked-in build from the handoff.
- **EDMD validated by selftest, conservation (drift < 1e-9), determinism and
  a virial-EOS check (3% at eta=0.35)**; the dynamical confirmation it exists
  for is downstream of T6/T9.

## Paper

- Affiliation is a placeholder; bibliography must be completed (add
  Donev–Torquato–Stillinger cond-mat/0408550 as prior art for the K=4
  structure — handoff README).
- The T4 degeneracy result (boundary-law count, observed closed form,
  non-crystalline 3D members) is not yet in the manuscript; it belongs in the
  degeneracy paragraph of Sec. IV and the discussion of the entropy objection.

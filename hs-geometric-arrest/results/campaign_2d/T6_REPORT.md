# T6 — the 2D replication campaign: full report

**Executed 2026-08-18 · repository `hs-geometric-arrest` · campaign commit
`e2f7817` (analysis pre-committed) · results commit `6fb10a6`**

---

## Executive summary

T6 was designed as the decision node of this program and upgraded, by the
appearance of Babu's cage theory (arXiv:2607.19185), into a decision
experiment between two live theories of two-dimensional arrest. It ran to
completion at full pre-registered statistics: 4400/4400 simulations, 100
densities, 44 independent replicas per density, N = 2016, zero unreachable
state points, every overlap audit zero, eps* measured to ~0.3% (absolute
sem ~ 5×10⁻⁵) at every density.

**Findings:**

1. **No structural feature exists at any predicted density** in the
   observables measured. eps*(φ) is locally a straight line through the
   maple-leaf window (0.777343, de Graaf's transition), through the kagome
   window (0.680175, this program's novel rung), and through the honeycomb
   window (0.604600): a weighted linear fit is fully consistent with the
   fine-grid data inside every window (chi²/dof 0.76–1.60; broken-line
   F-test p = 0.22–0.73).
2. **The mechanistic observable is silent.** P_wrap on the pre-registered
   frozen refscore is exactly zero in all 300 structural samples; marked
   fractions grow smoothly (0.13 → 0.31) with no acceleration at any
   target; marked-cluster susceptibility (χ ≈ 1.7) and size (ξ ≈ 0.7
   diameters) are flat across the entire density range. The submode weight
   — the fraction of large-particle cells whose local q sits at a target
   value — is identically zero at the kagome value everywhere and decays
   to zero at the maple-leaf value.
3. **The apparent dynamical arrest density is convention-dependent.** The
   pre-committed observation-time-threshold sweep finds the apparent
   arrest location at 0.7042 / 0.7561 / 0.7867 for criterion decades
   1/2/3, a drift of +0.0413 ± 0.0012 per decade (MC Brownian proxy).
4. **The pre-registered feature test mis-fired and its printed decision
   ("FEATURE AT 0.777343 → proceed to T9") must be discarded.** Its
   global power-law null is rejected by the data (chi²/dof = 137); the
   max-residual statistic tracks that misfit, demonstrated by an
   equal-strength "feature" (z = 8.4) at the honeycomb window deep in the
   normal fluid, where no theory predicts anything. Full forensics in
   `DIAGNOSTICS.md`.

**Reading:** on these observables, at this size and composition, the 2D
data favour a smooth, convention-dependent crossing (Babu) over a
structural geometric transition at a special density (the ground-state
reading of de Graaf's 0.7773, and this program's kagome prediction). The
negative is bounded and precise, but currently attackable on one front:
the failed test was the *pre-registered* one, and the local-linearity
tests above are post hoc. The follow-up program below exists to close
exactly that gap.

---

## 1. Background: what was at stake

**The geometric-ground-state picture** (this program, building on de
Graaf 2411.01199): tangential Voronoi packings define exact reference
densities; in 2D the family is honeycomb (K=3, φ = 0.604600), kagome
(K=2, φ = 0.680175 — the novel falsifiable rung derived here) and
maple-leaf (K=1, φ = 0.777343 = √3π/7, de Graaf's floret value). If
arrest is geometric, a structural feature sits pinned at a rung,
independent of composition and of any observation-time convention.

**The smooth-barrier picture** (Babu, arXiv:2607.19185): arrest is a
Kramers first-passage barrier crossing an observation-time threshold; no
structural feature exists at any special density; the apparent arrest
density shifts with the time criterion; 2D glass at φ ≈ 0.781.

The same published 2D data window (0.776–0.79) "corroborated" both — the
anti-numerology warning of this program demonstrated in the wild. T6 was
the mechanistic discriminator: fine grids of Δφ = 5×10⁻⁴ inside ±0.006
of each exact rung, statistics calibrated by the pilot positive control
(≥350 configurations per density), pre-registered observables
(composition-independent shell-percolation threshold eps*; P_wrap on a
refscore frozen at the repository's root commit), and a pre-committed
observation-time discriminator.

## 2. Campaign execution (all checks green)

| item | value |
|---|---|
| runs | 4400/4400 completed, none unreachable |
| grid | 100 densities in [0.580, 0.826]; Δφ = 5×10⁻⁴ inside the three windows |
| system | N = 2016, binary 1:1, R⁻¹ = 1.4 (reference-study composition) |
| per density | 44 independent replicas × 8 snapshots = 352 configurations |
| protocol | melt 20k → anneal to exactly zero overlap → eq 60k (swap on) → production 120k (swap OFF, logged) |
| audits | final overlap audit = 0 in all 4400 runs; radii multisets exact |
| wall time | 7.1 h on Apple M4, 7 workers (~2× faster than the conservative estimate) |
| provenance | manifest with command lines, source hashes, compiler, git rev: `data/2d-replication/manifest.json` |
| precision | eps* sem 2×10⁻⁵ – 2×10⁻⁴ per density (~0.3%) |

Notably, **no state point up to φ = 0.826 was unreachable**: this
composition equilibrates everywhere below its (pilot-fitted) RCP at
φ ≈ 0.85. The full 100-density eps*(φ) table is in `report.md`; the
per-replica data are on disk under `data/2d-replication/`.

## 3. Pre-registered outputs, verbatim

- Smooth-null fit: η_m = 0.8818, p = 1.424 (nuisance parameters).
- Feature test: honeycomb z_eff = 8.36 → "RESOLVABLE"; kagome −5.16 →
  "no feature"; maple-leaf 8.88 → "RESOLVABLE"; threshold 2.39.
- Printed decision: "FEATURE AT 0.777343: the pipeline resolves de
  Graaf's transition. T6 decision: proceed to the 3D campaign (T9). No
  kagome feature."
- MSD consistency: η_a = 0.864 ± 0.153 (MC proxy).
- Threshold sweep: drift +0.0413 ± 0.0012 per decade.

## 4. Why the printed decision is discarded

The feature test compares each fine window against a single global
three-parameter power law fitted to all 100 points. At pilot precision
that was adequate (the positive control validated it). At campaign
precision the null form itself is falsified: chi²/dof = 136.8, residuals
weaving between −63.7 and +14.6 sem across the grid — including at
coarse densities far from every window (+14.6 at φ = 0.58, −12.3 at
0.70, −63.7 at 0.82). The bootstrap normalisation in the test corrects
for statistical noise but cannot correct for deterministic misfit, so
z_eff grows without bound as sem shrinks. The decisive internal control:
the identical statistic reports an 8.4σ "feature" at the honeycomb
density, where the system is an ordinary liquid and **no theory on
either side predicts anything**. All three window verdicts are lobes of
one smooth misfit curve. (`DIAGNOSTICS.md` for the full forensics.)

This is a mis-specification the program had itself flagged in advance —
DEBT.md carries, from before the campaign: *"eps\* = C(η_m − η)^p is the
wrong functional form near contact percolation… a free-volume-motivated
form is needed."* The lesson is now quantitative: a global parametric
null cannot be used with data thirty times more precise than the form is
accurate.

## 5. What the data support

Exploratory analyses (labelled as such; see §7 for their pre-registered
successors):

- **Local linearity.** Inside every ±0.006 fine window (25–26 points), a
  sem-weighted straight line fits eps*(φ) within errors; adding a free
  breakpoint improves nothing (maple p = 0.22, kagome p = 0.69,
  honeycomb p = 0.73). Any kink at a rung is bounded to below the
  few-sem level at 5×10⁻⁴ density resolution.
- **No mechanistic signature.** P_wrap^(k) = 0 in every structural
  sample for all three reference classes; χ and ξ flat; marked fraction
  smooth. The frozen-refscore machinery works (gate T8: 14/14 on ideal
  states) — the campaign configurations simply never develop percolating
  regions of rung-like local geometry.
- **Convention dependence of arrest.** The apparent arrest density moves
  by ≈ 0.04 per decade of time criterion over the three accessible
  decades (gaps mildly decelerating: 0.052 → 0.031; three points cannot
  distinguish perpetual drift from eventual pinning above ~0.79).

Scope limits: one system size (N = 2016), one composition (R⁻¹ = 1.4),
structural observables limited to the eps* family and the frozen
refscore family, dynamics limited to an MC Brownian proxy over three
decades. These are exactly the axes the follow-up program attacks.

## 6. Implications

- **For de Graaf's 2D claim:** the replication does not reproduce a
  structural transition at 0.777343 at 30× the original statistics.
  Either the published signature was a statistical artefact of the same
  kind our pre-registered test just produced (a smooth-null failure), or
  it lives in an observable/size/composition outside this campaign's
  scope.
- **For this program's kagome rung:** unsupported. eps* is locally
  linear at 0.680175, the kagome submode weight is identically zero, and
  P_wrap^(kagome) never fires. The prediction was hedged ("a possible
  feature") and now carries an empirical bound.
- **For Babu:** the two observations his theory requires — no structural
  feature at any special density, arrest location shifting with the
  criterion — are both what we find, within scope limits. This is
  consistency, not confirmation: T6 measures 2D at one composition, and
  a drift of +0.041/decade over three decades is also compatible with
  other smooth scenarios.
- **For T9 (3D campaign):** remains blocked. Under the honest reading
  the applicable pre-registered branch is "NO FEATURE at any target…
  STOP — do not run the 3D campaign; write the negative result."
  Spending ~4800 core-hours on the mis-fired "proceed" branch is
  indefensible.
- **Theory standing:** every exact result of this program — the
  tangentiality theorem, the three-lattice enumeration, the depletion
  ladder and its spectral bound, the kagome rung's existence as an exact
  packing, the T4 degeneracy law — is untouched. What T6 bounds is the
  *dynamical relevance* of those states for this system class.

## 7. Follow-up program (design frozen before its data)

To make the outcome publishable in either direction, the analysis is
being redesigned so that a positive cannot be a null artefact and a
negative cannot be blamed on the method:

1. **Sideband hold-out null.** Per window, the background is measured
   from the data on both sides of the window alone (weighted local
   polynomial with propagated prediction covariance) — no global
   functional form anywhere. Statistic: standardised window offset.
2. **Kink-pinned-at-target test.** One extra degree of freedom: a slope
   discontinuity at exactly the predicted rung density. Sharp,
   theory-located, maximally powered.
3. **Empirical calibration.** The identical machinery run on (a) two new
   fine control windows at null densities (0.6300, 0.7250) with the full
   44-replica statistics — measuring the true false-positive rate at
   identical resolution — and (b) synthetic null curves with the
   campaign's noise structure.
4. **Quantified sensitivity.** Signal-injection: the minimum slope
   change and window offset detectable at 95% power, so a negative
   states "any transition-like feature larger than X was excluded".
5. **Equilibration defense.** Re-runs at the critical windows with 3×
   the equilibration, testing for drift in eps*.
6. **Independent kagome checks:** finite-size scan (N = 504 / 2016 /
   3520) for size-suppressed features and P_wrap crossings; a second
   composition (R⁻¹ = 1.7) for the composition-independence axis; and a
   direct stability test of the exact kagome packing itself
   (seeded-melting versus a hexagonal control at equal density) —
   probing the open DEBT question of thermodynamic relevance at T > 0.

Verdicts from these will be appended to this report as
`T6_REPORT_addendum` entries with their own commits.

---

*All raw data (2.5 GB), the manifest, per-run logs, and every analysis
input are on disk; every number in this report is reproducible from the
committed code at the stated commits. The pre-registered outputs are
preserved untouched in `report.md`; nothing in this report modifies
them.*

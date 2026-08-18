# T6 diagnostics — why the pre-registered decision text must not be taken at face value

Executed 2026-08-18. The campaign itself is impeccable: 4400/4400 runs, all
overlap audits zero, no unreachable state points, 352 configurations per
density, eps* sem ~ 5e-5 (0.3%) throughout. This file records the
post-run diagnostics of the ANALYSIS, separating the pre-registered output
(above, `report.md`, untouched) from what the data actually supports.
Everything below the first section is exploratory and labelled as such.

## 1. The pre-registered feature test is mis-specified at this precision

The test (committed at d8ec79a, before any campaign data) fits the smooth
null `eps* = C (eta_m - eta)^p` to ALL 100 densities and takes the maximum
standardised residual inside each fine window, normalised by its bootstrap
spread. Diagnostics on the campaign data:

- chi2/dof of the null over the full grid: **136.8**. The three-parameter
  form is falsified by its own fit at sem ~ 5e-5; residuals weave smoothly
  between -63.7 and +14.6 sem.
- The weave is global, not window-local. Coarse points far from every
  window carry |z| up to 14.6 (eta = 0.58) and -12.3 (eta = 0.70) — the
  same magnitude as the "features".
- The self-refuting tell: the test reports z_eff = 8.36 at the honeycomb
  window (eta = 0.6046), deep in the normal fluid where NO theory —
  ours or Babu's — predicts anything. All three window verdicts
  (maple +8.9, honeycomb +8.4, kagome -5.2) sit on lobes of one smooth
  misfit pattern.

This is the pre-declared DEBT weakness ("the wrong functional form...
a free-volume-motivated form is needed") biting at 30x the pilot
statistics. The bootstrap spread normalises statistical noise but not
deterministic misfit, so misfit inflates z_eff unboundedly as sem shrinks.
The positive control validated the test at pilot precision; at campaign
precision the null form itself is rejected, and the test loses meaning.

## 2. Exploratory local test: every window is featureless

Weighted straight line vs broken line (free breakpoint), fine windows only
(25-26 points each, span 0.012):

| window | chi2_line/dof | break improvement (F-test p) | verdict |
|---|---|---|---|
| honeycomb 0.6046 | 1.60 | p = 0.73 | locally linear |
| kagome 0.6802 | 1.14 | p = 0.69 | locally linear |
| maple-leaf 0.7773 | 0.76 | p = 0.22 | locally linear |

A plain line fits every window to within its errors. **No kink, no local
excursion, at 5e-5 absolute precision in eps*.**

## 3. The mechanistic observable is silent

P_wrap^(k) on the frozen refscore (mark threshold 0.10): **zero in all 300
structural rows** — no reference class ever percolates, at any density.
The K1 marked fraction grows smoothly 0.13 -> 0.31 across the whole range
with no acceleration near 0.7773; the marked-cluster susceptibility stays
chi ~ 1.7 and the cluster size xi ~ 0.7 particle diameters, flat
everywhere. Marked cells exist but never organise.

## 4. Dynamics (MC Brownian proxy, consistency only)

Observation-time-threshold sweep (pre-committed discriminator, e2f7817):
apparent arrest density 0.7042 / 0.7561 / 0.7867 at criterion decades
1/2/3, drift +0.0413 +- 0.0012 per decade. The arrest location moves with
the convention over the accessible three decades (gaps 0.052 -> 0.031,
mildly decelerating; three points cannot distinguish drift-forever from
eventual pinning far above 0.78). MSD extrapolation eta_a = 0.864 +- 0.153.

## 5. Honest reading

- The structural observables carry **no signature at 0.777343** (de
  Graaf's value), **none at the kagome rung 0.680175**, and none at the
  honeycomb density: eps* is locally linear through every window and
  P_wrap never fires.
- The apparent arrest density drifts with the observation-time criterion.
- Both facts favour the smooth-crossing picture (Babu, arXiv:2607.19185)
  over a structural geometric transition IN THIS OBSERVABLE SET, at
  N = 2016, binary 1:1 R^-1 = 1.4, MC dynamics.
- The mechanically printed decision ("FEATURE AT 0.777343... proceed to
  T9") is an artefact of the mis-specified null and must not gate T9.
  Under the honest reading the applicable branch is: "NO FEATURE at any
  target at these statistics... the negative result is publishable and
  must not be dressed up."
- Scope limits: one system size (no 2D finite-size scan), one
  composition, structural observables limited to eps* and the frozen
  refscore/P_wrap family, dynamics limited to three decades of a
  Brownian proxy. A different observable or the thermodynamic limit could
  still hide a feature; nothing here tests d = 3.

The T6 decision is a supervision point (CLAUDE.md): the call on how to
proceed — accept the negative, redesign the null and re-test on the same
data, or run new observables — belongs to Felix, not to this pipeline.

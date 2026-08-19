# T6 addendum 1 — redesigned test, calibration, and the follow-up battery

**Executed 2026-08-18/19 · design commit `1a7b921` (frozen before its data)
· 8054 additional runs (all groups 100% complete, zero unreachable)**

## 1. The redesigned feature test: final verdict

No global functional form anywhere: sideband hold-out offset (local
polynomial, full prediction covariance) and a kink test pinned at the
exact predicted density (linear + hinge over ±0.012 on the extended
near-sideband grid), with p-values from 2000 synthetic-noise calibrations
per window and Bonferroni over the three theory windows.

| window | z_offset | p_offset | z_kink | p_kink | kink excluded (95% power) | verdict |
|---|---|---|---|---|---|---|
| honeycomb 0.6046 (control) | −0.64 | 0.52 | 1.75 | 0.96 | ≥ 17.1% | null |
| **kagome 0.680175** | −0.48 | 0.50 | 0.55 | 0.95 | **≥ 14.4%** | **null** |
| **maple-leaf 0.777343** | 3.81 | 0.28 | −3.44 | 0.98 | **≥ 10.7%** | **null** |
| control 0.6300 (null density) | 0.05 | 0.87 | 0.80 | 0.78 | ≥ 16.0% | null |
| control 0.7250 (null density) | 0.74 | 0.56 | 0.20 | 0.87 | ≥ 12.5% | null |

**Controls clean.** The two control windows were measured at IDENTICAL
grid geometry and statistics (44 replicas, 33 densities each) at
densities where no theory predicts anything: the procedure does not
fire there, and it no longer fires at the honeycomb fluid window that
exposed the old test. A positive, had one existed, could not have been a
null artefact; the negative comes with quantified sensitivity: **any
transition-like slope change of ≥ 10.7% at de Graaf's density (≥ 14.4%
at the kagome rung) is excluded at 95% power.** A genuine arrest
transition implies order-unity slope changes.

## 2. Kagome: composition- and size-robustness of the null

Calibrated kink test on the kagome window re-measured independently:

| dataset | z_kink | p | excluded |
|---|---|---|---|
| binary R⁻¹ = 1.7, N = 2016 (44 reps) | 0.71 | 0.98 | ≥ 16.9% |
| binary R⁻¹ = 1.4, N = 504 (44 reps) | 1.49 | 0.87 | ≥ 27.1% |
| binary R⁻¹ = 1.4, N = 3520 (32 reps) | 0.43 | 0.79 | ≥ 15.5% |

No feature appears in a second composition and no feature grows with
system size. P_wrap remains zero at every size.

## 3. Equilibration defense — and a real equilibration limit

Sixteen replicas at 3× equilibration (180k sweeps) against the campaign:

| φ | eq×3 vs campaign | z |
|---|---|---|
| 0.680175 (kagome) | 0.03971 vs 0.03972 | −0.02 |
| 0.777343 (maple) | 0.01671 vs 0.01662 | +1.22 |
| 0.780000 | 0.01603 vs 0.01600 | +0.31 |
| **0.800000** | **0.01160 vs 0.01124** | **+4.94** |

The three theory windows are equilibrated. At φ = 0.800 the campaign
was NOT fully equilibrated (60k sweeps insufficient; eps* biased low by
~3%). Consequence: campaign points at φ ≳ 0.79 carry equilibration
bias; the kink tests (band ±0.012) never touch that region, the maple
offset test's outermost sidebands (band 0.045) graze it — one more
reason the offset statistic is secondary to the kink. All conclusions
at the three targets stand; the global eps*(φ) curve above 0.79 must
not be used quantitatively.

## 4. The mechanical story: kagome is real but marginal — and a new
prediction from the counting identity

Seeded-stability probe (`stability.json`): the exact kagome packing
(N = 2340, z = 4 verified) at its own density does not melt in 150k
sweeps (acceptance 0.09%, MSD → 7×10⁻⁵ σ²; K2-refscore 0.86) while a
hexagonal control at the SAME density melts completely.

Dilution ladder (`followup.json`), rigid ⇔ final MSD < 0.1 σ²:

| structure | 0.1% | 0.5% | 1% | 2% | 4% |
|---|---|---|---|---|---|
| kagome (z=4, tangential) | RIGID | melts | melts | melts | melts |
| square (z=4, NOT tangential) | RIGID | melts* | melts* | melts* | melts |
| hexagonal (z=6) | RIGID | RIGID | RIGID | RIGID | RIGID |

(*partial collapse, MSD ~0.7–0.8.)

Reading: kagome is an exact isostatic packing — mechanically rigid at
contact, with **zero thermal margin**, marginal exactly like the
non-tangential z = 4 control. It is mechanically real and
thermodynamically irrelevant, which *explains* the ensemble negative:
the equilibrium fluid has no reason to visit a measure-zero marginal
state. This closes the DEBT question ("thermodynamic relevance at
T > 0") for K = 2 in the negative.

**The rigidity ladder.** The counting identity fixes the contact number
of every ladder state: z = 6 − K. With the 2D rigidity threshold
z = 2d = 4 this predicts, with no free parameters: honeycomb (K = 3,
z = 3) **floppy**; kagome (K = 2, z = 4) **marginal** — measured above;
maple-leaf (K = 1, z = 5) **hyperstatic, rigid with a finite margin**.
If confirmed, de Graaf's density is mechanically meaningful after all —
not as an equilibrium transition (excluded in §1) but as the density of
a mechanically stable depleted solid. Measurement running
(`stab3`; verdict in addendum 2).

## 5. Standing decision state

- T6 ensemble verdict: **negative, calibrated, sensitivity-quantified,
  composition- and size-robust.** T9 (3D) remains blocked on the
  ensemble route.
- The mechanical route (rigidity ladder) is new, cheap, and follows
  from the exact mathematics; in 3D the same counting gives z = 12 − K
  ≥ 8 > 2d = 6 for every rung — all four 3D ladder states are
  hyperstatic, so the 3D mechanical question differs qualitatively
  from 2D. This is the natural impact direction for the paper.

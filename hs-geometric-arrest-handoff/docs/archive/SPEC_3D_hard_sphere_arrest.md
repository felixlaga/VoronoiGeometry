> **SUPERSEDED** by `IMPLEMENTATION_SPEC.md` + `CLAUDE_CODE_TASKS.md`. Kept for history only; where documents conflict, IMPLEMENTATION_SPEC.md wins.

# Geometric ground states and dynamical arrest in 3D hard spheres — study specification

Target question: **does a geometric/topological property of the Voronoi network
predict dynamical arrest in 3D hard-sphere systems**, in the sense established
for 2D disks by de Graaf, arXiv:2411.01199v2?

This document is the implementation spec. The geometry it rests on is closed and
verified in `enumerate_3d_ground_states.py`; the simulation campaign is not yet
run. Nothing here may be filled in with estimated numbers.

---

## 0. Settled geometry (do not re-derive; do re-run the script)

**Theorem (3D analogue of `phi = q`).** For a monodisperse packing in which every
Voronoi cell is a *tangential* polyhedron with the particle as its insphere,

```
V = sigma S / 6        =>        eta = (pi/6) sigma^3 / V = 36 pi V^2 / S^3 = Q_iso
```

where `Q_iso = 36 pi V^2 / S^3` is the 3D isoperimetric quotient (= Wadell
sphericity cubed). Proof: decompose the polyhedron into pyramids of height
`sigma/2` on each face.

**Admissibility criterion.** A candidate geometric ground state must satisfy
*every Voronoi neighbour is a contact*. If any Voronoi face is supported by a
non-touching neighbour then `eta != Q_iso` and the 2D argument does not transfer.
This is the single criterion that partitions the candidate list.

**Consequence: BCC and diamond are inadmissible.**

| structure | eta | Q_iso | faces | z | tangential |
|---|---|---|---|---|---|
| FCC | 0.740480 | 0.740480 | 12 | 12 | yes |
| HCP | 0.740480 | 0.740480 | 12 | 12 | yes |
| FCC − 1-in-13 | 0.683520 | 0.683520 | 11 | 11 | yes |
| FCC − 1-in-7 | 0.634698 | 0.634698 | 10 | 10 | yes |
| simple hexagonal, c=a | 0.604600 | 0.604600 | 8 | 8 | yes |
| FCC − 1-in-5 | 0.592384 | 0.592384 | 9 | 9 | yes |
| FCC − 1-in-4 | 0.555360 | 0.555360 | 8 | 8 | yes |
| simple cubic | 0.523599 | 0.523599 | 6 | 6 | yes |
| BCC | 0.680175 | **0.753367** | 14 | 8 | **no** |
| diamond | 0.340087 | **0.461735** | 16 | 4 | **no** |

**Complete enumeration of tangential Bravais lattices: exactly three** (SC,
simple hexagonal `c=a`, FCC). The seven Voronoi face-class norms are *linear* in
the Selling parameters `p_ij = -b_i.b_j`, so tangentiality is a linear system;
solving it over all live-face sets and all zero-patterns is exhaustive. No
14-faced solution exists (all seven norms equal forces `p_ij = 1/6`, giving
norms 1/2 and 2/3), so BCC-type cells can never be tangential.

**The depleted-FCC ladder.** Scanning all vacancy sublattices of FCC of index
`k <= 14` and keeping congruent + tangential cells gives exactly

```
k in {4, 5, 7, 13}        eta_k = (1 - 1/k) * pi/sqrt(18)
                          z_k   = 12 - 12/(k-1)
```

i.e. precisely the `k` with `(k-1) | 12`, each particle adjacent to `12/(k-1)`
vacancies. The 2D analogue on the triangular lattice is `(k-1) | 6`, giving
`k in {2,3,4,7}` and `phi_k = (1 - 1/k) pi/sqrt(12) =`
0.453450, 0.604600, 0.680175, 0.777343. de Graaf's floret pentagonal tiling is
`k = 7`; his honeycomb/caging state is `k = 3`. **Both are members of one
family**, and the family predicts two further 2D states he does not discuss.

### Errors found in arXiv:2411.01199v2 (all reproducible from the script)

1. **BCC vacancy neighbour count.** The paper reports "4 nearest neighbours for
   which the hexagonal face is extended, `eta = 4 sqrt3 pi/35 ≈ 0.621874`". A BCC
   vacancy has **8** nearest neighbours; their cells have `V/V0 = 71/64` and
   `eta = 8 sqrt3 pi/71 = 0.613115`. The 6 next-nearest values
   (`V/V0 = 49/48`, `eta = 6 sqrt3 pi/49 = 0.666294`) are correct. With 8 and 6
   the volume bookkeeping closes exactly (`8×0.109375 + 6×0.020833 = 1.000000 V0`);
   with the paper's counts it does not. Both cell types are non-tangential, so
   averaging them to reach ~0.644 has no geometric content.
2. **Tetrahedral-octahedral honeycomb.** The paper writes "each of the *two*
   tetrahedra to an octahedron contributes 1/4 of its volume". An octahedron
   there has 8 faces hence 8 adjacent tetrahedra, each shared by 4 octahedra —
   2 tetrahedra of volume in total, not 1/2. The paper's cell fills 3/4 of space.
   Corrected: `eta = 2 sqrt3 pi/27 ≈ 0.40307`, not `8 pi/(27 sqrt3) ≈ 0.537422`,
   and the cell reduces to the FCC rhombic dodecahedron. The claimed agreement
   with frictional RLP ≈ 0.536 does not survive.
3. **2D honeycomb depletion fraction.** Section VIII B says "two out of three
   particles are removed, with the remaining ones forming a regular honeycomb".
   Removing 2 in 3 leaves a sparser triangular lattice. The honeycomb is
   *1 in 3* removed, consistent with the paper's own
   `phi = pi/(3 sqrt3) = (2/3) pi/sqrt(12)`.
4. **Not an error, but unestablished in the paper:** the 1-in-13 FCC structure
   requires a perfect 1-error-correcting code on FCC (ball size 1+12=13). It
   exists — 96 generator triples found — and the resulting cells are congruent,
   11-faced and tangential with `eta = Q = 2 sqrt2 pi/13`.
5. The paper states a preliminary search found no candidate in the 3D glass
   range 0.58–0.64. Three fall in it: 0.592384, 0.604600, 0.634698.

### Honest reading of point 5

The candidate ladder is *dense* between 0.52 and 0.74. Matching a measured
`eta_a` to two significant figures against a literature range that itself spans
0.588–0.639 is therefore weak evidence, and much weaker than the corresponding
2D coincidence was. **The study must not be numerological.** It must test the
mechanism: the shape of `P_eta(Q_iso)`, the contact/Voronoi-neighbour ratio, and
independence of composition. Those are what would distinguish a real geometric
selection from a coincidence among many available numbers.

---

## 1. Falsifiable predictions

Ranked by discriminating power, strongest first.

- **P1 (mechanism).** As `eta -> eta_a`, the fraction of Voronoi neighbours that
  are contacts, `f_c = <z_contact / z_Voronoi>`, rises toward 1. Under the null
  hypothesis (no geometric selection) `f_c` has no feature at `eta_a`.
  This observable has no 2D precedent in the paper; it follows directly from the
  tangentiality criterion and is the cheapest decisive test.
- **P2 (composition independence).** `eta_a` is the same, within error, for
  binary, trinary and continuously polydisperse mixtures, over a wide range of
  size ratio — the 3D counterpart of de Graaf's near-constant `phi_a ≈ 0.777`
  over `R <~ 0.83` and his trinary check.
- **P3 (distributional).** `P_eta(Q_iso)` develops a resolvable feature (mode,
  submode cusp, or slope discontinuity) at `Q_iso = eta_a`, and `eta_a` coincides
  with one of the tabulated admissible values.
- **P4 (kinetic).** The neighbour-exchange rate between adjacent `n_Q` classes
  peaks near `eta_a`, as `r_{5<->6}` does in 2D (Appendix F of the paper).

Falsification: `eta_a` lands where no admissible structure exists (e.g. clearly
at 0.64 = frictionless RCP, or at BCC's 0.680175 which is inadmissible), **or**
`f_c` shows no feature, **or** `eta_a` drifts systematically with composition.
P2 failing alone kills the hypothesis regardless of P3.

---

## 2. System and protocol

**Engine.** Event-driven MD, Smallenburg (Eur. Phys. J. E 45, 22 (2022)), used
in 3D natively — no 2D port needed. Canonical NVT, thermalisation by redrawing
velocities from Maxwell-Boltzmann. Reduced time in units of the bare diffusion
time; verify `D -> 1` as `eta -> 0` before any production run.

**Compositions** (P2 is the point of this list; do not shorten it):
1. monodisperse (control, will crystallise — needed to locate the crystalline envelope)
2. binary 1:1, size ratio `R` on a grid in `[1/1.7, 1]`, mirroring the 2D study
3. trinary 1:1:1 at `sigma, 0.8 sigma, 0.6 sigma` (direct analogue of the paper's Fig. 9 validation)
4. continuous Gaussian polydispersity, 3% to 15% in steps of 2%

**Sizes.** `N = 10^4` production; `N = 10^5` at three `eta` per composition as a
finite-size check. Cubic box, periodic. Report both; do not average across `N`.

**State points.** 60+ non-equidistant `eta` per composition, denser above 0.55.
Upper bound `eta_m` from the 3D Desmond–Weeks construction (250 packings at the
production `N`); note in DEBT that this checks close packing, not disorder.

**Preparation and equilibration.**
- grow particles at dimensionless rate `1e-3` to the target `eta`
- equilibrate `1e4` reduced time units without growth before sampling
- **gate**: repeat at `1e3` and `1e5` for at least 4 `(eta, composition)` points
  near the expected transition and confirm structural observables change only
  nominally. This is the equilibration check the 2D paper performed; skipping it
  invalidates everything downstream.
- 350+ independent realisations per state point

**Do not** substitute soft potentials. The paper's own retrospective (Section I A)
attributes its earlier failure to place the peak correctly to an insufficiently
hard core.

---

## 3. Observables — exact definitions

Radical Voronoi tessellation via `voro++` (respects radius differences).

1. `Q_iso = 36 pi V^2 / S^3` per cell. Primary structural variable. **Not** the
   raw face count, **not** `W_{12,0}`.
2. `P_eta(Q_iso)`: PDF with non-linear binning; error bars = sqrt(counts).
   Compute for all particles, and separately for the largest size fraction —
   in 2D the signal was clearest there and vanished in the small-particle
   distribution.
3. `n_Q`: generalised face number. **Open definition problem.** In 2D,
   `q_r(n) = pi sin(2 pi/n) / (2n sin^2(pi/n))` is exact for regular n-gons.
   There is no such closed form in 3D. The reference family must be *chosen and
   documented*, not assumed. Options, in order of preference:
   (a) invert `Q_iso` against the admissible tangential cells in the table above
       (cube 6, hex prism 8, square bipyramid 8, ..., rhombic dodecahedron 12) —
       physically motivated but non-monotonic in face count;
   (b) invert against Catalan solids;
   (c) drop `n_Q` entirely and work with `Q_iso` directly.
   Option (c) is the safe default. Log the choice in DEBT.
4. `f_c = <z_contact / z_Voronoi>` with `z_contact` counted at separation
   `< sigma_ij (1 + 1e-6)`. For hard spheres exact contact is measure-zero, so
   report `f_c(delta)` for `delta` in `[1e-6, 1e-2]` and take the limit
   behaviour, exactly as the shell-percolation analysis handles `epsilon`.
5. `Q_6` (Steinhardt), for locating crystallisation only.
6. MSD `<|r(t)|^2>` over `1e4` time units, 150 independent curves.
7. `D(eta)` from fits with `alpha in (0.97, 1]`; `eta_a` from
   `D = A (eta_a - eta)^b`; inflection `eta_i` of `dD/deta` from a polynomial fit.
8. Shell percolation: inflate each sphere by `epsilon`, cluster, find the
   spanning threshold, locate the peak of
   `(d log[eta_m - eta]) / (d log epsilon)`. Independent cross-check on `eta_a`.
9. Neighbour-exchange rate between adjacent `Q_iso` classes, normalised by total
   particle displacements (P4).

---

## 4. Staging

Each stage must pass before the next begins. No stage may be started with a
placeholder for a previous stage's output.

- **S0** Reproduce `enumerate_3d_ground_states.py` output. Gate: table matches.
- **S1** EDMD harness, `D -> 1` at low `eta`, monodisperse equation of state
  against published hard-sphere data. Gate: agreement within error.
- **S2** Voronoi + `Q_iso` pipeline. Gate: applied to *perfect* FCC, HCP, SC,
  simple hexagonal and the four depleted-FCC lattices, it returns the tabulated
  `Q_iso` to 1e-9 and flags BCC/diamond as non-tangential.
- **S3** Binary sweep -> `eta_a(R)`, `eta_i(R)`. Gate: equilibration check passed.
- **S4** `f_c` (P1) and `P_eta(Q_iso)` (P3) at the S3 state points.
- **S5** Trinary and polydisperse (P2). **This is the decisive stage.**
- **S6** Neighbour-exchange rates (P4).
- **S7** Literature comparison table, structured as the paper's Appendix A.

## 5. Rules

- No invented parameters, no placeholder data, no synthetic fallbacks anywhere
  in the pipeline. If a fit does not converge, record the failure.
- Report failed attempts. The 2D paper's honesty about the SISF being unusable
  is the standard to match.
- Every claimed coincidence with a tabulated `eta_k` must state the fitted
  uncertainty on `eta_a` alongside it. Given the density of the candidate ladder,
  an agreement quoted without an error bar is worthless.
- Separate speculation from result explicitly in all write-ups.
- Repository documentation must match the repository exactly. Unfinished or
  deferred items go in `DEBT.md`.

## 6. Pilot outcome (see RESULTS_pilot.md)

S0, S1, S2 complete and gated. A pilot at N=864 with MC dynamics returned:
`eta_a ~ 0.60-0.63` with a fit-window systematic an order of magnitude larger
than the statistical error; **P1 null** (no feature in the Voronoi face-gap
distribution at arrest); and a composition-*dependent* maximum in `<Q_iso>`,
which is the wrong kind of feature. **Do not start S3-S7 until P1 is redefined
against a null model and `eta_a` is estimated by shell percolation rather than by
the power-law fit.**

## 7. Known gaps -> DEBT.md at repo creation

- Completeness of the tangential enumeration is established for **Bravais
  lattices only**. The depleted-FCC scan covers `k <= 14` and vacancy patterns
  that are *sublattices*; non-sublattice perfect coverings, HCP-derived
  depletions, and general multi-orbit structures are **not** enumerated.
- `k = 2, 3` are excluded by the scan rather than by proof; the argument that
  `(k-1) | 12` is necessary is observed, not derived.
- No 3D analogue of `q_r(n)` exists; the `n_Q` reference family is a choice.
- **`f_c` as originally defined is unmeasurable.** Exact contact has measure zero;
  a strict tolerance gives `f_c ~ 0.001` at all densities. The `f_c(delta)`
  reformulation is measurable but shows no feature. Needs a null model.
- MC dynamics is a Brownian proxy, not EDMD; `eta_a` from it is not the EDMD `eta_a`.
- The `g(sigma+)` linear extrapolation degrades above eta ~ 0.4 (4% at 0.45);
  replace with a cusp-aware fit before any EOS claim.
- Whether any of these structures is thermodynamically or kinetically relevant
  at `T > 0` is entirely open. The paper's own Section IX concedes the analogous
  2D objection (a single realisation, hence no configurational entropy) and never
  resolves it.

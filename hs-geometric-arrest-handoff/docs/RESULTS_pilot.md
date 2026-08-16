# Pilot results — 3D hard spheres, geometric ground states

**Status: preliminary and underpowered. These numbers do not answer the physics
question and must not be quoted as if they do.** They are recorded because the
protocol failures they expose are the useful output.

## Setup

N = 864, NVT hard-sphere Monte Carlo (`hsmc.c`), local displacement moves only in
production (swap moves used for equilibration only, then switched off). MC
dynamics is a proxy for Brownian dynamics, not EDMD. Two compositions: 10%
Gaussian polydisperse, and binary 1:1 at R = 0.714. eta from 0.40 to 0.62.
Equilibration 15k sweeps, production 40k sweeps, 12 snapshots. Every final
overlap audit returned 0.

Gates passed: EOS vs Carnahan-Starling 0.04% at eta=0.30 (4.0% at 0.45 — that
drift is the linear g(sigma+) extrapolation window, not the engine, and needs a
cusp fit). Radical-Voronoi pipeline reproduces the tabulated Q_iso of every
perfect tangential lattice to 1e-8 and flags BCC at f_c = 8/14.

## Measured

| eta | D_poly | alpha_poly | <Q>_poly | <n_f>_poly | D_bin | <Q>_bin |
|---|---|---|---|---|---|---|
|0.40|2.96e-04|1.03|0.67203|14.711|2.76e-04|0.66545|
|0.48|4.56e-05|1.05|0.69539|14.431|4.32e-05|0.68580|
|0.54|4.08e-06|0.93|0.71157|14.241|5.68e-06|0.70304|
|0.58|3.08e-07|0.52|0.72275|14.049|7.48e-07|0.71345|
|0.60|1.00e-07|0.00|0.72625|14.008|3.46e-07|0.71277|
|0.62|6.67e-08|0.08|0.72495|14.037|1.95e-07|0.70983|

Arrest fit D = A (eta_a - eta)^b on eta >= 0.48:
poly eta_a = 0.6122 +- 0.0117 (b=3.83); binary eta_a = 0.6257 +- 0.0121 (b=3.88).

## Three findings, two of them negative

**1. The quoted error on eta_a is fiction.** Varying the fit window gives
eta_a = 0.599 / 0.612 / 0.622 / 0.884 (poly) and 0.625 / 0.626 / 0.652 / 0.652
(binary). The window systematic dwarfs the statistical error by an order of
magnitude. The pilot constrains eta_a only to "roughly 0.60-0.63", which spans
the simple-hexagonal state (0.6046) and the k=7 depleted-FCC state (0.6347) and
therefore discriminates nothing. Any future report of eta_a must quote the
window systematic alongside the fit error.

**2. P1 is not measurable as written, and its repaired form shows no feature.**
With a strict contact tolerance the contact fraction of Voronoi faces is
f_c ~ 0.001 at every density: exact contact has measure zero in a thermal hard-
sphere fluid, so "every Voronoi neighbour is a contact" cannot be tested
directly. Reformulated as f_c(delta) with face gap (r_ij - sigma_ij)/sigma_ij:

|eta|1e-4|1e-3|3e-3|1e-2|3e-2|1e-1|median gap|
|---|---|---|---|---|---|---|---|
|0.44|0.0005|0.0032|0.0100|0.0304|0.0854|0.2551|0.24803|
|0.50|0.0004|0.0048|0.0146|0.0467|0.1316|0.3505|0.16868|
|0.54|0.0008|0.0074|0.0201|0.0612|0.1753|0.4388|0.12620|
|0.58|0.0009|0.0095|0.0283|0.0892|0.2341|0.5328|0.08808|
|0.60|0.0012|0.0120|0.0367|0.1095|0.2788|0.5941|0.07144|
|0.62|0.0017|0.0170|0.0517|0.1531|0.3598|0.6365|0.05445|

The gap distribution tightens smoothly and monotonically. d(ln median gap)/d eta
steepens monotonically (-6.4, -7.3, -9.0, -10.5, -13.6) with no inflection, cusp
or saturation anywhere near the arrest region. **At this system size and run
length there is no structural signature of approach to a tangential
configuration.** This is the cheapest decisive test in the spec and it came back
null. Either the effect needs far better statistics and larger N, or the
observable must be normalised against the trivial affine-compression trend, or
the mechanism does not operate in 3D.

**3. The one structural feature found is the wrong kind.** <Q_iso> rises with eta
and passes through a maximum — at eta = 0.60 (poly) and eta = 0.58 (binary).
This is the 3D analogue of de Graaf's max_phi(n_qbar). But in 2D that feature is
explicitly the *composition-dependent* one; the arrest density phi_a is the
composition-*independent* one. Here the maximum moves with composition, exactly
as expected, so it is not evidence for a selected geometric state.
<n_face> falls to ~14.0 and flattens where the dynamics arrests.

## Net

Geometry: solved (see SPEC section 0). Physics: unresolved, and the pilot's main
contribution is that it invalidates the planned P1 protocol and shows the eta_a
fit systematic is the binding constraint. The full campaign should not be run
until P1 is redefined against a null model and the eta_a estimator is replaced by
something less window-sensitive (shell percolation, which the spec already lists
as the independent cross-check, is the obvious candidate).

---

# Follow-up: is there a null that makes P1 sharp, and is there a better estimator?

## Shell percolation is the observable that works

Threshold `eps*(eta)` (in units of the mean diameter) at which the inflated
contact network first wraps the box in all three directions, by bisection,
union-find with relative-displacement tracking, 3 configurations per point:

| eta | poly | binary | rel. diff |
|---|---|---|---|
|0.40|0.02731|0.02811|2.9%|
|0.44|0.02154|0.01938|10.6%|
|0.48|0.01613|0.01565|3.0%|
|0.50|0.01248|0.01408|12.0%|
|0.52|0.01179|0.01177|0.2%|
|0.54|0.00999|0.01006|0.7%|
|0.56|0.00828|0.00843|1.8%|
|0.58|0.00679|0.00680|0.1%|
|0.60|0.00534|0.00496|7.4%|
|0.62|0.00372|0.00256|36.9%|

**`eps*(eta)` is composition-independent to 4.3% on average for eta <= 0.60**,
across a 10% polydisperse system and a binary 1:1 mixture at R = 0.714 — two very
different size distributions. In de Graaf's framework composition independence is
*the* signature of a geometric rather than a packing-detail origin, so this is the
most promising thread in the whole pilot. (The 0.62 binary point deviates by 37%
and is almost certainly under-equilibrated near jamming; it should be discarded,
not explained.)

The local slope `d ln eps*/d eta` = -5.9, -6.6, -11.0, -7.8, -5.6, -8.8, -9.7,
-11.0, -15.0, -18.1 steepens overall but wiggles at eta = 0.48-0.52 by more than
the quoted spread. With 3 configurations per point that wiggle is not resolvable
from noise. **No feature can be claimed; none can be excluded either.**

Fitting `eps* = C (eta_m - eta)^p` gives eta_m = 0.83 (poly) / 0.76 (binary) with
window spreads of 0.15 / 0.13 — the assumed power-law form is wrong near
contact percolation. Restricted windows converge to 0.675 / 0.633, bracketing RCP.
The vanishing point should be extracted with a form motivated by the free-volume
scaling, not a bare power law.

## Why P1 cannot be rescued in its current form

The face-gap distribution has a single scale. Rescaling `delta` by the median gap
collapses `f_c(delta, eta)` across all six densities: the small-delta ratio
`f_c / (delta/median)` is 0.79, 0.81, 0.93, 0.84, 0.86, 0.92 for
eta = 0.44...0.62 — constant within the ~15% scatter of two snapshots. If the
distribution is single-scale then `f_c(delta, eta)` contains **no information
beyond the median gap itself**, and the median gap has no feature (its log
derivative steepens monotonically). P1 is therefore not merely null here; as an
observable it is degenerate with a quantity that varies smoothly through arrest.

Any repaired version has to look at the *arrangement* of the near-contact faces
rather than their number — which is exactly what shell percolation does.

## Revised campaign (supersedes spec section 4 staging)

The expensive part of the original plan was the MSD-based `eta_a`, and it is the
part that failed: fit-window systematics of 0.285 (poly). The percolation
threshold needs no long-time dynamics at all — only equilibrated configurations.

1. Drop the MSD/`eta_a` power-law fit as the primary estimator. Keep MSD only as
   a consistency check.
2. Run **short** equilibrations at many eta (say 40 points, 0.40-0.64) and
   sample ~100 configurations each. This is roughly an order of magnitude cheaper
   than the original S3 and resolves the `eps*` local exponent.
3. Test for a peak in `d ln eps*/d eta`, and test whether its location is
   composition-independent across polydisperse / binary / trinary.
4. Only then commit to EDMD for the dynamical confirmation.

The decisive question is now narrow and cheap: **does the composition-independent
`eps*(eta)` curve have a resolvable feature, and if so, does it sit at one of
0.523599, 0.555360, 0.592384, 0.604600, 0.634698, 0.683520?**

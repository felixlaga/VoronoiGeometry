# Cheaper routes to the answer, and a positive control that changes the plan

The 3D campaign is expensive and, because the admissible ladder is dense, its
result would be ambiguous even if it succeeded. This document reports what I
found looking for cheaper routes. One of them overturns a recommendation I made
earlier, and I retract it explicitly below.

---

## 1. The positive control: run the pipeline where the answer is known

**Rationale.** A null result only means something if the instrument can detect
the effect. De Graaf reports arrest in 2D bidisperse disks at
$\phi_a \approx 0.777$, and the floret pentagonal value is $0.777343$. Running
*the same pipeline* on that system is the direct test of whether my observables
are sensitive at all. 2D is also far cheaper than 3D.

**Setup.** 2D hard-disk NVT Monte Carlo (`hsmc2d.c`), bidisperse 1:1 at
$R^{-1}=1.4$ matching the reference study, $N=504$, $\phi$ from 0.65 to 0.805,
4 independent seeds, 12 configurations per density. Zero overlaps in every audit.

**Result A — the machinery is correct.** Fitting the shell-percolation threshold
to $\varepsilon^* = C(\phi_m-\phi)^p$ returns $\phi_m = 0.8508$ against a
literature RCP of $\approx 0.84$ for this size ratio, with $p = 1.079\pm0.061$.
Recovering RCP to 1% from percolation alone is an independent validation of the
union-find wrapping code.

**Result B — no feature at 0.777, at this sensitivity.**

| $\phi$ | $\varepsilon^*$ | sem | residual from smooth law |
|---|---|---|---|
|0.650|0.04996|0.00082|+1.1|
|0.700|0.03418|0.00101|-1.8|
|0.740|0.02521|0.00064|-1.0|
|0.755|0.02118|0.00046|-2.0|
|0.765|0.01971|0.00032|+0.3|
|0.775|0.01798|0.00035|+2.3|
|0.785|0.01472|0.00045|-0.0|
|0.795|0.01255|0.00022|+1.0|
|0.805|0.00975|0.00019|-1.2|

Residuals are in units of the standard error. A single smooth power law
describes every point; the sign pattern `+---++-+-` shows no structure, and
nothing distinguishes the interval containing 0.777343.

The raw local exponent $\mathrm{d}\ln\varepsilon^*/\mathrm{d}\phi$ appears to
jump from $\approx-9$ to $-20$ exactly in that interval, which looked like a
detection. It is not: it is entirely accounted for by the $1/(\phi_m-\phi)$
divergence toward close packing. Running de Graaf's own diagnostic,
$\mathrm{d}\log(\phi_m-\phi)/\mathrm{d}\log\varepsilon$, which equals $1/p$ and
is flat for a pure power law, gives
$0.754, 1.012, 0.836, 1.536, 1.348, 0.706, 1.036, 0.781$ with errors of
$0.07$–$0.58$ against $1/p = 0.927$. **Every interval is consistent with flat.
No peak is resolvable.**

**Result C — this is a calibration, not a refutation.** The rms fractional
residual about the smooth law is 3.0% and the typical fractional standard error
is 2.2% per density with 12 configurations at $N=504$. To detect a 1% feature at
$3\sigma$ requires **~510 configurations per density**. De Graaf used 350
realizations at $N=2048$ and describes his feature as "small but measurable".

So the control does not show that the effect is absent. It shows that my pilot
was roughly forty times short of the statistics needed to see it, in 2D and by
extension in 3D.

## 2. Retraction

I previously recommended building the 3D campaign around $\varepsilon^*(\eta)$
because it was composition-independent to 4.3% while the diffusion-based
$\eta_a$ was dominated by a fit-window systematic. Composition independence
still holds and is still a point in its favour. But **I had not shown that
$\varepsilon^*$ resolves the transition at all**, and the control now shows it
does not at the statistics I had. That recommendation was premature and should
not be acted on in its earlier form. The corrected statement is: $\varepsilon^*$
is a candidate observable whose required sample size is now known
(~500 configurations per density), not an observable demonstrated to work.

## 3. The route this opens: replicate 2D first

This is the main practical finding.

The 3D campaign is expensive **and** its interpretation is ambiguous, because
the admissible values $0.5236, 0.5554, 0.5924, 0.6046, 0.6347, 0.6835$ are
spaced $\approx0.04$ apart across the range where the glass transition is
reported. The 2D system has neither problem: it is roughly thirty times cheaper
per configuration, and its target value $0.777343$ is isolated, with the nearest
family member at $0.680175$.

Most importantly: **de Graaf's central 2D claim has not been independently
replicated.** The entire 3D programme is downstream of it.

A 2D replication at $N=2048$ with 350 configurations per density across
$\phi \in [0.70, 0.82]$ is of order single-digit CPU-hours per density on one
core, i.e. a few days on one machine and hours on any cluster. It would test:

1. Does $\varepsilon^*(\phi)$ show a resolvable departure from the smooth
   close-packing law near 0.777?
2. Does the submode cusp in the large-particle $q$ distribution reproduce, which
   is de Graaf's own primary structural signal?
3. Does the family prediction hold — is there anything at the kagome rung
   $\phi = 0.680175$, which he does not discuss?

Item 3 (the kagome rung) is a genuinely new falsifiable test that costs nothing extra, since those
densities lie inside the range already being swept.

**Decision rule.** If the 2D replication produces a clean feature at 0.777, the
mechanism is real and the 3D search is worth its cost. If it does not, the 3D
campaign should not be run at all. Either way the 2D result is publishable.

## 4. A free analytical observation: isostaticity

The contact numbers of the ladder are $z_k = z_{\rm cp} - z_{\rm cp}/(k-1)$.
In 2D, $z_k = 6-6/(k-1)$ gives

| lattice | honeycomb | kagome | maple-leaf |
|---|---|---|---|
| $\phi_K$ | 0.604600 | 0.680175 | 0.777343 |
| $z_K$ | **3** | **4** | 5 |

(An earlier version of this table listed a $\phi=0.453450$ member with $z=0$;
the spectral bound $K\le3$ on the triangular lattice excludes it. Retracted.)

The 2D frictional isostatic number is $d+1 = 3$ and the frictionless one is
$2d = 4$. **The 2D ladder contains both exactly**, at $k=3$ (de Graaf's caging
state) and $k=4$ (the rung he does not discuss). His glass state, $k=7$, is
hyperstatic at $z=5$.

In 3D the corresponding numbers are $z_k = 12-12/(k-1) = 8, 9, 10, 11$ for
$k=4,5,7,13$ — the ladder contains **neither** isostatic number. The frictionless
value $2d=6$ is supplied instead by the simple cubic lattice
($\eta = \pi/6 = 0.523599$), and no tangential structure in the Bravais class has
$z = d+1 = 4$.

This asymmetry is a concrete, cost-free reason to expect 3D to behave differently
from 2D, and it sharpens what to look for: if rigidity rather than tiling is what
selects the density, the 3D candidate should be simple cubic, not a depleted-FCC
rung.

## 5. Routes considered and set aside

- **Reanalysis of published 2D data.** The family prediction at $\phi=0.680175$
  could in principle be tested against de Graaf's existing figures without any
  simulation. Set aside only because I do not have the underlying data; it should
  be requested, since it is the cheapest test in this document.
- **Free-volume / cell-model comparison of the candidates.** Compute the cell
  free energy of each admissible structure to ask which is thermodynamically
  competitive at $T>0$. Cheap, and it addresses the entropy objection directly
  rather than sidestepping it. Not attempted here; recommended.
- **Maximum-entropy null for $P(Q_{\rm iso})$.** Compare the measured
  distribution against a max-entropy prediction constrained by density and mean,
  looking for excess weight at the admissible values. Attractive but the null is
  not uniquely defined, which is the same problem that sank the contact-fraction
  observable.
- **Pure analytics on the arrest density.** No route found. The 2D argument is
  itself a geometric coincidence promoted to a mechanism; there is no derivation
  to extend.

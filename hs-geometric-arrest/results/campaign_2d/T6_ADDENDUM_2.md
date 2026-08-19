# T6 addendum 2 — the rigidity ladder confirmed

**Executed 2026-08-19 · 30 runs (stab3) · zero free parameters**

The counting identity fixes the contact number of every 2D ladder state:
z = 6 − K. Against the Maxwell threshold z = 2d = 4 this predicts the
mechanical character of each rung with no adjustable input. Measured by
seeded dilution (rigid ⇔ final MSD < 0.1 σ² over 150k sweeps, 3 seeds):

| rung | z | prediction | 0.1% | 0.5% | 1% | 2% | 4% |
|---|---|---|---|---|---|---|---|
| honeycomb K=3 (φ=0.6046) | 3 | floppy | melts (39) | melts | melts | melts | melts |
| kagome K=2 (φ=0.6802) | 4 | marginal | RIGID | melts (17) | melts | melts | melts |
| maple-leaf K=1 (φ=0.7773) | 5 | rigid, finite margin | RIGID | RIGID | RIGID | melts (1.0) | melts (4.0) |

All three rungs behave exactly as predicted: floppy / marginal /
hyperstatic with a rigidity margin an order of magnitude wider than the
isostatic rung's. The maple-leaf state is a mechanically stable depleted
solid at φ ≈ 0.77 — the density of de Graaf's claimed transition is
mechanically meaningful as the stability point of a hyperstatic depleted
lattice, even though the EQUILIBRIUM transition there is excluded at the
≥ 10.7%-kink level (addendum 1).

3D extension (untested, from the same identity): z = 12 − K ≥ 8 > 2d = 6
for every rung K ≤ 4 — ALL 3D ladder states are hyperstatic. The 2D
marginality of kagome does not carry over; the 3D mechanical question is
qualitatively different and is the natural next measurement (seeded
dilution in the existing 3D engine's restart mode).

Geometry caveats: seeds carry ~1% box strain (near-commensurate square
box), so contact densities sit 0.5–1% below the exact rungs; the
dilution grid is coarse (0.1/0.5/1/2/4%); "rigid" means no melting found
in 150k sweeps of local MC at N ≈ 450–2350 — a statement about
mechanical stability under thermal agitation, not a free-energy claim.

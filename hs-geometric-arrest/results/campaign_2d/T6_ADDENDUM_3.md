# T6 addendum 3 — the 3D rigidity ladder and literature positioning

**Executed 2026-08-19 · 60 runs (stab3d) · cubic tori exactly commensurate,
contact densities equal to the golden rung values to machine precision**

## 3D seeded dilution (rigid ⇔ final MSD < 0.1 σ², 150k sweeps, 2 seeds)

| structure | z | φ_contact | 0.1% | 0.5% | 1% | 2% | 4% |
|---|---|---|---|---|---|---|---|
| K4 (tangential) | 8 | 0.5554 | RIGID | RIGID | melts (3.2) | melts | melts |
| K3 | 9 | 0.5924 | RIGID | RIGID | RIGID | RIGID | melts (10.4) |
| K2 | 10 | 0.6347 | RIGID | RIGID | RIGID | RIGID | RIGID |
| K1 | 11 | 0.6835 | RIGID | RIGID | RIGID | RIGID | RIGID |
| FCC | 12 | 0.7405 | RIGID | RIGID | RIGID | RIGID | borderline (0.11) |
| BCC (NOT tangential) | 8 | 0.6802 | RIGID | RIGID | RIGID | RIGID | borderline (0.15) |

Prediction (z = 12 − K ≥ 8 > 2d = 6: all rungs hyperstatic with finite
margins) **confirmed**; the margin is monotone in z. Honest nuance: the BCC
control (z = 8, excluded from the tangential family) is as stable as K4 —
in 3D the mechanics follows the Maxwell count alone; tangentiality is the
density-selection principle (η = Q_iso), not the rigidity principle. The
borderline entries (MSD 0.11–0.15 at 4%) are onset-of-diffusion rattling,
not clean melts.

## Literature positioning (docs/LITERATURE_NOTES.md)

Established, now cited in Sec. VIII: kagome isostaticity and its marginal
mechanics (Kane & Lubensky, Nat. Phys. 10, 39 (2014)); mechanically stable
vacancy-depleted close-packed crystals in 3D (tunneled crystals, Torquato &
Stillinger, J. Appl. Phys. 102, 093511 (2007), strictly jammed down to
φ = 2π/9); the maple-leaf lattice as the 1/7-depleted triangular lattice of
frustrated magnetism. Appears new (no prior hits found): the z = z₀ − K
rigidity ladder as a corollary of the depletion counting identity tied to
the arrest-density ladder; the measured hyperstatic stability margin of the
maple-leaf *packing*; the calibrated equilibrium-route exclusion at the
rungs. Babu arXiv:2607.19185 remains v1/unrefereed; de Graaf 2411.01199
unchanged — re-check both before submission.

# Validation of ported analysis code against the recorded pilot numbers

Frame conventions follow the prototypes (2D: last 3 frames x 4 seeds;
3D: last 3 frames for eps*, last 6 for Voronoi). Criteria are stated per
block and reflect snapshot-selection and estimator differences; this is
consistency validation, not bit-identity.

## 2D shell percolation (positive control, ALTERNATIVE_ROUTES.md)

| phi | eps* recorded (sem) | eps* ported (sem) | n | pass |
|---|---|---|---|---|
| 0.650 | 0.04996 (0.00082) | 0.04996 (0.00085) | 12 | ok |
| 0.700 | 0.03418 (0.00101) | 0.03418 (0.00106) | 12 | ok |
| 0.740 | 0.02521 (0.00064) | 0.02521 (0.00067) | 12 | ok |
| 0.755 | 0.02118 (0.00046) | 0.02118 (0.00048) | 12 | ok |
| 0.765 | 0.01971 (0.00032) | 0.01971 (0.00033) | 12 | ok |
| 0.775 | 0.01798 (0.00035) | 0.01798 (0.00037) | 12 | ok |
| 0.785 | 0.01472 (0.00045) | 0.01472 (0.00047) | 12 | ok |
| 0.795 | 0.01255 (0.00022) | 0.01255 (0.00022) | 12 | ok |
| 0.805 | 0.00975 (0.00019) | 0.00975 (0.00019) | 12 | ok |

RCP recovery fit `eps* = C (phi_m - phi)^p`: phi_m = 0.8508 (recorded 0.8508; literature ~0.84 for this size ratio), p = 1.079 (recorded 1.079) — ok. This is the 1%
RCP validation of the union-find wrapping code.

## 3D shell percolation (RESULTS_pilot.md, 3 configurations per point)

| eta | poly recorded | poly ported | binary recorded | binary ported | pass |
|---|---|---|---|---|---|
| 0.40 | 0.02731 | 0.02732 | 0.02811 | 0.02811 | ok |
| 0.44 | 0.02154 | 0.02154 | 0.01938 | 0.01938 | ok |
| 0.48 | 0.01613 | 0.01613 | 0.01565 | 0.01565 | ok |
| 0.50 | 0.01248 | 0.01248 | 0.01408 | 0.01408 | ok |
| 0.52 | 0.01179 | 0.01179 | 0.01177 | 0.01177 | ok |
| 0.54 | 0.00999 | 0.00999 | 0.01006 | 0.01006 | ok |
| 0.56 | 0.00828 | 0.00828 | 0.00843 | 0.00843 | ok |
| 0.58 | 0.00679 | 0.00679 | 0.00680 | 0.00680 | ok |
| 0.60 | 0.00534 | 0.00534 | 0.00496 | 0.00496 | ok |
| 0.62 | 0.00372 | 0.00372 | 0.00256 | 0.00256 | ok |

## 3D Voronoi observables (spot values)

| eta | Q_poly rec | Q_poly port | n_f rec | n_f port | Q_bin rec | Q_bin port | pass |
|---|---|---|---|---|---|---|---|
| 0.40 | 0.67203 | 0.67156 | 14.711 | 14.722 | 0.66545 | 0.66506 | ok |
| 0.54 | 0.71157 | 0.71269 | 14.241 | 14.214 | 0.70304 | 0.70323 | ok |
| 0.62 | 0.72495 | 0.72504 | 14.037 | 14.026 | 0.70983 | 0.70940 | ok |

**Overall: PASSED**

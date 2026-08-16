# T4 — uniform-K degeneracy on growing tori (theory run, no simulation)

Exhaustive enumeration of all uniform-K independent vacancy sets
(`geometry/coloring.enumerate_uniform_K`, constraint-propagating DFS with
the vacancy count pinned to `Nv = K N/(z0+K)`). Every row below is a
complete enumeration unless flagged otherwise; capped runs are bounds,
not counts. Incommensurate sizes (non-integer `Nv`) are structural zeros
and are excluded from fits.

## Counts

| dim | lattice | side | N sites | K | Nv | solutions | orbits | cosets | complete | DFS nodes | s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | triangular | 3 | 9 | 3 | 3 | 3 | 1 | 1 | yes | 35 | 0.0 |
| 2 | triangular | 6 | 36 | 3 | 12 | 3 | 1 | 1 | yes | 347 | 0.0 |
| 2 | triangular | 9 | 81 | 3 | 27 | 3 | 1 | 1 | yes | 2049 | 0.0 |
| 2 | triangular | 12 | 144 | 3 | 48 | 3 | 1 | 1 | yes | 11846 | 0.0 |
| 2 | triangular | 15 | 225 | 3 | 75 | 3 | 1 | 1 | yes | 68370 | 0.2 |
| 2 | triangular | 4 | 16 | 2 | 4 | 16 | 4 | 4 | yes | 221 | 0.0 |
| 2 | triangular | 6 | 36 | 2 | 9 | 40 | 4 | 1 | yes | 1572 | 0.0 |
| 2 | triangular | 8 | 64 | 2 | 16 | 88 | 10 | 4 | yes | 8052 | 0.0 |
| 2 | triangular | 10 | 100 | 2 | 25 | 184 | 10 | 1 | yes | 36207 | 0.2 |
| 2 | triangular | 12 | 144 | 2 | 36 | 376 | 22 | 4 | yes | 154734 | 0.7 |
| 2 | triangular | 14 | 196 | 2 | 49 | 760 | 28 | 1 | yes | 656505 | 2.8 |
| 2 | triangular | 16 | 256 | 2 | 64 | 1528 | 58 | 4 | yes | 2815056 | 11.7 |
| 2 | triangular | 7 | 49 | 1 | 7 | 14 | 2 | 2 | yes | 939 | 0.0 |
| 2 | triangular | 14 | 196 | 1 | 28 | 14 | 2 | 2 | yes | 36597 | 0.2 |
| 3 | FCC | 4 | 32 | 4 | 8 | 16 | 4 | 4 | yes | 1334 | 0.0 |
| 3 | FCC | 6 | 108 | 4 | 27 | 40 | 4 | 1 | yes | 467377 | 2.6 |

## Scaling of the count

An extensive degeneracy (finite entropy density) means `ln(count) ~ N`;
a boundary-law degeneracy means `ln(count) ~ sqrt(N)`. Both fits below;
the better-fitting law is the honest headline, with the caveat that the
accessible sizes are small.

| dim | K | points | slope vs N | R²(N) | slope vs √N | R²(√N) |
|---|---|---|---|---|---|---|
| 2 | 3 | 5 | — | — | — | — |
| 2 | 2 | 7 | 0.0181 | 0.9545 | 0.3753 | 0.9979 |

## Reading

- **2D K=3**: count is constant (3) at every commensurate size — unique up to translation, zero degeneracy growth. No fit applies.
- **2D K=2**: ln(count) fits sqrt(N) (boundary law) better (R² 0.9979 vs 0.9545). 

## Observed closed form (K=2 class)

Every complete even-side triangular count obeys **S(side) = 6·2^(side/2) − 8** exactly: confirmed on all 7 sides. The side-16 count (1528) was predicted from the formula before being enumerated, and matched. This is an OBSERVED law, not a derived one (DEBT.md); its form — ln S ≈ (ln 2/2)·√N plus corrections — is exactly a boundary/stacking law: the degeneracy grows like independent layer choices along one axis, the same mechanism as ABC stacking freedom, not like a bulk entropy.

The two 3D FCC K=4 counts (16 at side 4, 40 at side 6) satisfy the same formula — two points prove nothing, but the coincidence is worth recording: the 3D degeneracy may be the same layered mechanism.

## Tangentiality of every enumerated state

- 2D 8x8 K=2, all orbits, all cells: worst deviation 4.44e-16 (exact within float)
- 3D FCC side-6 K=4, all orbits, sampled cells: worst deviation 5.55e-16 (exact within float)

## Caveats

- Sizes are small (up to 15x15 triangular, FCC side 6); both scaling laws
  are fitted through few points and the verdict is provisional, not a
  theorem.
- Orbit counts are per torus translation group only; point-group symmetry
  is not quotiented, so orbit counts overstate the number of genuinely
  distinct patterns by up to the point-group order.
- K=1 (maple-leaf class) has two orbits on the side-7 and side-14 tori —
  the two enantiomers of the chiral pattern — both lattice cosets: unique
  up to translation and chirality at accessible sizes.
- K=3 (honeycomb class) has exactly 3 solutions (1 orbit) at every
  commensurate size 3–15: unique up to translation. No degeneracy at all.

## Consequence for the configurational-entropy objection

The reference family does contain non-crystalline members (most orbits are
not lattice cosets, in 2D K=2 and in 3D K=4 side 6), and every one of them
is exactly tangential at the ladder density — so the geometric reference
state is a degenerate family, not a single crystal. But at the sizes
enumerated the degeneracy is a boundary law, not extensive: it carries no
finite entropy density. On this evidence the degeneracy weakens the
entropy objection only marginally; it does not answer it.

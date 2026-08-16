# Geometry gates G0 / G0b

Every check compares against `REFERENCE_VALUES.json` (immutable; tolerances final). No simulation input.

## Gate G0 — PASSED

| check | result | detail |
|---|---|---|
| simple_cubic: eta=Q_iso=0.5235987756, faces=6, z=6 | ok | got eta=0.5235987756 |
| simple_hexagonal_ca: eta=Q_iso=0.6045997881, faces=8, z=8 | ok | got eta=0.6045997881 |
| fcc: eta=Q_iso=0.7404804897, faces=12, z=12 | ok | got eta=0.7404804897 |
| Bravais enumeration count == 3 | ok | got 3 (4158 linear systems solved) |
| no 14-face solution (p_ij=1/6 gives norms 1/2, 2/3; 4t=3t=2*sum p forces t=0) | ok |  |
| depletion ladder scan == {4,5,7,13} | ok | got [4, 5, 7, 13] |
| K4 (k=4): eta=0.5553603673, z=8, p-vector [3, 3, 3, 3, 3, 3, 3, 3] | ok | got eta=0.5553603673, p=[3, 3, 3, 3, 3, 3, 3, 3] |
| K3 (k=5): eta=0.5923843918, z=9, p-vector [3, 3, 4, 4, 4, 4, 4, 4, 4] | ok | got eta=0.5923843918, p=[3, 3, 4, 4, 4, 4, 4, 4, 4] |
| K2 (k=7): eta=0.6346975626, z=10, p-vector [4, 4, 4, 4, 4, 4, 4, 4, 4, 4] | ok | got eta=0.6346975626, p=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4] |
| K1 (k=13): eta=0.6835204520, z=11, p-vector [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4] | ok | got eta=0.6835204520, p=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4] |
| K4 exact cell V=8/3, S=8*sqrt2 (integer coords) | ok | got V=2.666666666667 S=11.313708498985 |
| K1 exact cell V=13/6, S=13*sqrt2/2 (integer coords) | ok | got V=2.166666666667 S=9.192388155425 |
| bcc: inadmissible, eta=0.6801748 != Q_iso=0.7533666 | ok |  |
| diamond: inadmissible, eta=0.3400874 != Q_iso=0.4617350 | ok |  |
| BCC vacancy: 8 nearest neighbours (paper says 4); volume balances to 1 V0 | ok | redistributed=1.000000000 |
| BCC vacancy eta_correct = 8 sqrt3 pi/71 = 0.613115 | ok |  |
| tet-oct: 8 adjacent tetrahedra, cell = rhombic dodecahedron, eta = 2 sqrt3 pi/27 | ok |  |

## Gate G0b — PASSED

| check | result | detail |
|---|---|---|
| counting_eta reproduces the 3D ladder to 1e-12 | ok |  |
| counting_eta reproduces the 2D family to 1e-12 | ok |  |
| spectral_K_max(FCC) == 4  (band [-4,12]) | ok | got 4 |
| spectral_K_max(triangular) == 3  (band [-3,6]; 0.4534 rung retracted) | ok | got 3 |
| 3D K=4 [x,y,z all even]: uniform/independent/tangential/congruent, p=[3, 3, 3, 3, 3, 3, 3, 3] | ok |  |
| 3D K=3 [y+2z=0 mod 5]: uniform/independent/tangential/congruent, p=[3, 3, 4, 4, 4, 4, 4, 4, 4] | ok |  |
| 3D K=2 [x+2y+3z=0 mod 7]: uniform/independent/tangential/congruent, p=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4] | ok |  |
| 3D K=1 [x+3y+4z=0 mod 13]: uniform/independent/tangential/congruent, p=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4] | ok |  |
| 2D K=3 (K3_honeycomb): phi=0.6045997881, z=3 | ok | got phi=0.6045997881 |
| 2D K=2 (K2_kagome): phi=0.6801747616, z=4 | ok | got phi=0.6801747616 |
| 2D K=1 (K1_maple_leaf): phi=0.7773425847, z=5 | ok | got phi=0.7773425847 |
| divacancy spread == sqrt2 - 1; isolated vacancy == 0 | ok | got 0.414214 / 0.00e+00 |
| tri_4x4_K2: solutions=16, orbits=4, cosets=4 | ok | got 16/4/4 |
| tri_6x6_K3: solutions=3, orbits=1, cosets=1 | ok | got 3/1/1 |
| tri_6x6_K2: solutions=40, orbits=4, cosets=1 | ok | got 40/4/1 |
| every cell of every 6x6 K=2 solution tangential with Q = sqrt3 pi/8 | ok | max |Q-target|=4.4e-16, max face spread=3.3e-16 |


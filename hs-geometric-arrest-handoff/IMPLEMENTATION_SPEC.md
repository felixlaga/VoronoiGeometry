# Implementation specification — `hs-geometric-arrest`

Build specification for a coding agent. This document is normative: it defines
the repository, every module, every validation gate, and the acceptance criteria
for "done". It supersedes the earlier exploratory spec.

**Scientific context is in `paper.tex` and `EXPLAINER.md`. Read `paper.tex`
Sections II–IV before writing any code**; the geometry there is settled and must
be reproduced, not re-derived.

---

## 0. Non-negotiable rules

1. **No invented data.** No placeholder numbers, no synthetic fallbacks, no
   "example" values that could be mistaken for measurements. If a fit fails,
   record the failure and propagate `NaN`.
2. **No shortcuts.** No stage may be started with a stub standing in for a
   previous stage's output.
3. **Gates are blocking.** A stage that fails its gate halts the pipeline with a
   non-zero exit code. Do not add a `--force` flag.
4. **Repository documentation must match the repository exactly.** No generic
   README boilerplate, no files describing things that do not exist.
5. Unfinished or deferred work goes in `DEBT.md`, never in a comment saying
   "TODO later".
6. Every derived number in any report must carry the systematic that dominates
   it, not only the statistical error. The pilot showed the fit-window
   systematic exceeds the statistical error by an order of magnitude.
7. Production-ready means: README complete, install documented, every feature
   working, repository runs end-to-end without errors, all gates pass.

---

## 1. Repository layout

```
hs-geometric-arrest/
├── README.md
├── DEBT.md
├── LICENSE
├── pyproject.toml
├── Makefile
├── src/
│   └── hsga/
│       ├── __init__.py
│       ├── geometry/
│       │   ├── tangential.py        # eta = Q_iso theory, cell computation
│       │   ├── lattices.py          # reference structures + Selling enumeration
│       │   └── depletion.py         # FCC vacancy-sublattice ladder
│       ├── engine/
│       │   ├── hsmc.c               # NVT hard-sphere Monte Carlo
│       │   ├── edmd.c               # event-driven MD (stage 5 only)
│       │   └── driver.py            # parameter sweeps, provenance
│       ├── analysis/
│       │   ├── voronoi.py           # radical Voronoi, Q_iso, face gaps
│       │   ├── percolation.py       # shell percolation, union-find wrapping
│       │   ├── dynamics.py          # MSD, D(eta), alpha
│       │   └── eos.py               # g(r), contact value, Carnahan-Starling
│       ├── gates/
│       │   ├── gate_s1_eos.py
│       │   ├── gate_s2_lattices.py
│       │   ├── gate_s3_equilibration.py
│       │   └── gate_s4_finitesize.py
│       └── report/
│           └── tables.py
├── scripts/
│   ├── run_geometry.py
│   ├── run_sweep.py
│   └── run_analysis.py
├── tests/
├── data/            # gitignored; simulation output
└── results/         # committed; tables and figures
```

---

## 1b. External-data independence

The programme requires no data from de Graaf or anyone else. All targets are
generated in-repo: the 2D replication and the kagome test run from `hsmc2d.c`;
the degeneracy counts are pure combinatorics; the 3D campaign is gated on the 2D
outcome. Reanalysis of de Graaf's published data is an optional appendix item,
not a dependency.

## 2. Module specifications

### 2.1 `geometry/tangential.py`

```python
def voronoi_cell(site, lattice, basis, *, radii=None, nshell=3, area_tol=1e-8) -> Cell
```
Radical (Laguerre) Voronoi cell by half-space intersection. The radical plane
between `i` and `j` sits at `d_ij = (|r_ij|^2 + R_i^2 - R_j^2) / (2|r_ij|)` from
`i`. With `radii=None` this reduces to the ordinary Voronoi cell.

`Cell` carries: `V`, `S`, `Q_iso = 36*pi*V**2/S**3`, `n_faces` (non-degenerate
only), `face_distances`, `face_gaps`, `n_contacts`, `r_in`.

```python
def is_tangential(cells, sigma, tol=1e-7) -> bool
```
True when every non-degenerate face of every cell is supported at `sigma/2`.
This is the admissibility criterion: **every Voronoi neighbour is a contact**.

**Invariant to assert in code:** when `is_tangential` is True,
`abs(eta - Q_iso) < 1e-9`. This is the theorem of `paper.tex` Sec. II and a
violation means a bug, not a discovery.

### 2.2 `geometry/lattices.py`

Named structures as `(lattice, basis)`: FCC, HCP, BCC, SC, simple hexagonal at
`c=a`, diamond, FCC−1-in-4.

```python
def selling_gram(p) -> ndarray            # p = (p12,p13,p14,p23,p24,p34)
def face_class_norms(p) -> ndarray        # the 7 class norms, LINEAR in p
def enumerate_tangential_lattices() -> list
```
The seven Voronoi face-class norms (Voronoi's theorem: non-zero classes of
`L/2L`) are *linear* functionals of the Selling parameters. Tangentiality is
therefore a **linear system**, not an optimisation. Enumerate over all
`sum(k>=3) C(7,k)` live-face sets crossed with all zero-patterns of `p` up to
three vanishing parameters, solve each system exactly, sample any non-trivial
null space, and verify each candidate by direct cell computation.

Do **not** implement this as a numerical minimisation. A Nelder-Mead search over
the shape space fails: the objective is non-smooth and flat-bottomed, and it
misses the boundary strata where simple cubic and simple hexagonal live. This is
recorded because it was the actual failure mode.

### 2.2b `geometry/coloring.py`  (new)

The ladder is now DERIVED, and the derivation must be executable:

```python
def counting_eta(z0, K, eta_cp) -> float          # z0/(z0+K) * eta_cp
def bloch_band(lattice, basis, contacts, nq=200_000) -> (lo, hi)
def spectral_K_max(lattice, basis, contacts) -> int   # floor(-band_min)
def modular_vacancies(K) -> callable              # the mod-2/5/7/13 rules (FCC)
def enumerate_uniform_K(torus, K) -> list         # exhaustive DFS, small tori
def is_lattice_coset(vacancy_set, torus) -> bool
```

`enumerate_uniform_K` is the degeneracy probe: constraint-propagating DFS
(vacancies independent; every occupied site ends with exactly K vacant
neighbours), returning all solutions, their translation-orbits, and the coset
flag. Verified reference outputs: 6x6 triangular, K=2 -> 4 orbits, 1 coset;
K=3 -> 1 orbit. Extend to growing tori and fit log(count) vs N to test for
extensive degeneracy — this is the configurational-entropy question and is
PURE COMBINATORICS, no simulation needed.

Hard-won correctness notes: independence need not be imposed as a separate
axiom when testing candidate reference states geometrically (tangentiality
implies it — divacancy spread 0.414 in FCC); and the spectral condition is
necessary, not sufficient (HCP passes for all K<=4 but existence there is open —
`DEBT.md`).

### 2.3 `geometry/depletion.py`

```python
def hermite_normal_forms(k) -> list       # lower-triangular, 0<=d,e<a; 0<=f<b
def depleted_fcc(k) -> (lattice, basis)
def scan_depletions(kmax=20) -> list
```
Two mandatory correctness details, both of which were bugs on first attempt:
- HNF off-diagonal ranges are `0 <= d,e < a` and `0 <= f < b`. Wrong ranges
  silently enumerate an incomplete set of sublattices.
- Coset representatives **must be wrapped into the primitive cell** before the
  cell computation. Unwrapped representatives sit far outside the supercell and
  their periodic neighbourhoods are truncated, which silently returns
  "not tangential" for structures that are.

Expected output for `kmax >= 13`: exactly `k in {4,5,7,13}`.

### 2.4 `engine/hsmc.c`

NVT hard-sphere Monte Carlo, cell lists, deterministic seeded RNG.

- Initialisation: FCC at the target `eta` with all radii equal (overlap-free for
  `eta <= 0.7405`), melt, then set target radii and anneal the resulting overlaps
  to **exactly zero** with `E = sum (sigma_ij - r_ij)^2` at decreasing `T`.
  Abort with a non-zero exit code if `E > 0`.
- Equilibration: local displacement moves + swap moves.
- Production: **local moves only**. Swap moves destroy the physical dynamics.
  Log the switch-off explicitly.
- Final overlap audit; non-zero exit if any overlap.
- Output: `.cfg` multi-frame (`N L` header per frame, then `x y z r`), `.msd`,
  `.log`.

Local-move MC is a Brownian proxy, not EDMD. Every reported `eta_a` from it must
be labelled as such.

### 2.5 `analysis/percolation.py`

```python
def percolates(pos, rad, L, eps) -> bool
def eps_star(pos, rad, L, *, lo=1e-5, hi=0.2, iters=18) -> float
```
Union-find with **relative-displacement tracking**: when two particles already in
the same set are joined across a periodic image, the net displacement mismatch
identifies a wrapping cluster. Require wrapping in all three directions.
`eps_star` bisects in `log(eps)`; report in units of the mean diameter.

This is the primary estimator (Sec. 4). It requires no long-time dynamics.

### 2.5b `analysis/topology.py`  (new)

Raw face counts are brittle: tiny displacements create and destroy small faces.
Implement PERSISTENT radical-Voronoi topology:

```python
def face_filter(cell, g_cut, a_cut) -> Cell       # drop faces with gap>g_cut AND area/S<a_cut
def p_vector(cell) -> tuple                        # sorted polygon edge counts
def topo_hash(cell) -> str                         # canonical face-adjacency graph hash
def persistence(cell, g_grid, a_grid) -> float     # fraction of grid where topo_hash constant
```

Reference targets the classifier must reproduce on perfect lattices:
3^8 (K=4), 3^2 4^7 (K=3), 4^10 (K=2), 4^11 (K=1), rhombic dodecahedron (FCC),
cube (SC), hexagonal prism (simple hex). Use our own canonical graph hash
(colour-refinement / Weinberg-style); do not assume VoroTop is installable in
the sandbox network.

### 2.5c `analysis/refscore.py`  (new)

Per-particle distance to each reference cell C_k:
`s_k(i) = w_Q ((Q_i - eta_k)/eta_k)^2 + w_p d_p(pvec_i, pvec_k) + w_t eps_t(i)^2`
with the tangentiality defect `eps_t^2 = mean_j ((l_ij - R_i)/R_i)^2` over kept
faces. Weights and thresholds live in `refscore_frozen.json`, committed BEFORE
any dynamical analysis (gate G5 below verifies the commit predates the data).

### 2.5d `analysis/pwrap.py`  (new)

`P_wrap^(k)(eta, N)`: mark cells with `s_k < threshold`, connect marked cells
sharing a kept radical face, test wrapping with the existing union-find.
Also cluster susceptibility and correlation length. Finite-size crossing
analysis across N. This replaces the degenerate f_c observable as the
mechanistic test.

### 2.5e `analysis/isoconfig.py`  (new)

Isoconfigurational propensity: from one equilibrated configuration run M
independent MC move-sequences (different RNG streams), measure per-particle
mean displacement at several lag times. Evaluate whether {s_k, cluster size,
persistence} predict propensity BEYOND the baseline set {eta, tetrahedrality,
Voronoi anisotropy, V_cell} on held-out configurations. If they add nothing,
the hypothesis fails de Graaf's own standard.

### 2.6 `analysis/voronoi.py`, `dynamics.py`, `eos.py`

- `voronoi.py`: per-configuration `Q_iso`, `n_faces`, face gaps
  `g_ij = (r_ij - sigma_ij)/sigma_ij`, and `f_c(delta)`.
- `dynamics.py`: MSD, `D` from the long-time slope, exponent `alpha`, and
  `eta_a` from `D = A (eta_a - eta)^b` **with a mandatory fit-window scan**; the
  reported uncertainty is `max(statistical, window spread)`.
- `eos.py`: `g(r)`, contact value, `Z = 1 + 4 eta g(sigma+)`, Carnahan-Starling.
  The linear extrapolation of `g(sigma+)` degrades above `eta ~ 0.4`; implement a
  cusp-aware fit and document the window.

---

## 3. Validation gates

| Gate | Test | Pass criterion |
|---|---|---|
| **G0 geometry** | reproduce Sec. II–IV of `paper.tex` | all tabulated `Q_iso` to `1e-8`; tangential lattice count `== 3`; depletion ladder `== {4,5,7,13}` |
| **G1 EOS** | monodisperse `Z` vs Carnahan-Starling at `eta = 0.30, 0.35, 0.40` | `<1%` deviation |
| **G2 pipeline** | radical Voronoi on perfect lattices | exact `Q_iso`; `f_c = 1` for tangential; `f_c = 8/14` for BCC |
| **G3 equilibration** | vary equilibration length `1e3 / 1e4 / 1e5` at >= 4 `(eta, composition)` near arrest | structural observables change within error |
| **G4 finite size** | `N = 864 / 4000 / 10^4` at 3 densities | `eps*(eta)` consistent within error |
| **G0b derivation** | `geometry/coloring.py` | counting_eta reproduces the ladder to 1e-12; spectral_K_max(FCC)=4, (triangular)=3; the four modular rules give uniform/independent/tangential/congruent; kagome exact at sqrt(3)pi/8; divacancy spread 0.414 reproduced |
| **G5 pre-registration** | `refscore_frozen.json` | file hash committed before the first dynamical dataset it is applied to |

G0 and G2 already pass in the prototype and must keep passing.

---

## 4. The campaign

**Build the 2D replication first.** See `ALTERNATIVE_ROUTES.md`. A positive
control showed the pipeline cannot resolve a feature at the known 2D transition
with 12 configurations per density, and calibrated the requirement at
**~500 configurations per density**. The 2D system is ~30x cheaper per
configuration and its target value is isolated, whereas the 3D ladder is dense.
The 3D campaign below is gated on the 2D replication producing a feature. If it
does not, do not run stage 3D at all.

Add to the repository: `engine/hsmc2d.c`, and a `dim` parameter throughout
`analysis/` (the half-space cell construction and the union-find wrapping code
are dimension-generic already).

The primary observable is `eps*(eta)`, **not** `eta_a` from the MSD. Rationale in
`paper.tex` Sec. VI: the MSD fit carried a window systematic of 0.285, while
`eps*(eta)` is composition-independent to 4.3% and needs only equilibrated
configurations.

- **2D first (blocking).** Targets, all exact and in-repo:
  `0.604600` (honeycomb, K=3, z=3 frictional-isostatic),
  `0.680175` (kagome, K=2, z=4 frictionless-isostatic — the novel rung),
  `0.777343` (maple-leaf, K=1 — de Graaf's value). N=2048, >=350 configs/point,
  grid `Delta_phi = 5e-4` within `+-0.006` of each target, coarse elsewhere.
  Measure eps*(phi), P_wrap^(k), refscores, submode cusp of the large-particle
  q distribution. Success at 0.777 validates the pipeline; anything at 0.680
  is a new result either way.
- **Compositions (3D):** 10% Gaussian polydisperse; binary 1:1 at `R = 0.714`;
  trinary 1:1:1 at `sigma, 0.8 sigma, 0.6 sigma`; polydispersity series
  3,5,7,10,15%.
- **Densities (3D):** coarse over `[0.40, 0.70]`, plus `Delta_eta = 5e-4`
  within `+-0.006` of every admissible value
  `{0.523599, 0.555360, 0.592384, 0.604600, 0.634698, 0.683520}`. Note the
  literature landmarks phi_MCT ~ 0.592 and phi_0 ~ 0.635 (Berthier-Witten,
  PRE 80, 021502) sit on the K=3 and K=2 rungs to quoted precision; treat this
  as motivation for grid placement, never as evidence.
- **Statistics:** >= 500 independent configurations per point, calibrated from
  the 2D positive control (2.2% fractional sem at 12 configurations; a 1% feature
  needs ~510 at 3 sigma). The pilot had 3-12 and was ~40x short.
- **N:** 4000 production, `10^4` for G4.
- **Deliverable:** `eps*(eta)` per composition; local exponent
  `d ln eps* / d eta`; test for a resolvable peak; test whether the peak location
  is composition-independent.
- **Decision point:** only if a composition-independent peak exists, commit to
  `edmd.c` and the dynamical confirmation. If not, the project reports a negative
  result, which is a publishable outcome and must not be dressed up.

Comparison targets (exact, from `paper.tex` Table II):
`0.5235988, 0.5553604, 0.5923844, 0.6045998, 0.6346976, 0.6835205`.

---

## 5. Anti-patterns

The agent must actively check for and refuse:

- silently widening a tolerance to make a gate pass;
- reporting `eta_a` without the fit-window scan;
- quoting agreement with a tabulated `eta_k` without the uncertainty — the ladder
  is dense (spacing ~0.04), so unqualified agreement is meaningless;
- reintroducing `f_c` at fixed `delta` as a headline observable; it is degenerate
  with the median gap (Sec. VI C); the mechanistic observable is `P_wrap^(k)`
  built on the pre-registered refscore;
- quoting sub-1e-3 agreement with phi_MCT or phi_0 as evidence — both are
  two-to-three significant-figure fit parameters of avoided/extrapolated
  singularities and are system-dependent;
- classifying cells by raw face count instead of persistent topology;
- tuning refscore weights after seeing dynamical data (G5 exists for this);
- leaving swap moves on during production;
- using `pypdf`-style shortcuts, generic README text, or files that describe
  functionality that does not exist.

---

## 6. `DEBT.md` at repository creation

Seed with, at minimum:

- Tangential enumeration is complete for **Bravais lattices only**. Non-lattice
  structures are covered only by the FCC sublattice scan to `k <= 20`;
  non-sublattice perfect coverings, HCP-derived depletions and general
  multi-orbit structures are not enumerated.
- ~~`(k-1) | 12` observed, not derived~~ CLOSED: counting identity + spectral
  bound (`geometry/coloring.py`), with independence forced by tangentiality.
- Whether uniform-K degeneracy is EXTENSIVE (entropy density > 0) is open;
  attack by `enumerate_uniform_K` on growing tori. Highest-value theory item.
- Existence of uniform-K states on HCP is open (spectral condition passes for
  all K<=4); if they exist, density degenerates across stackings and topology
  classification becomes mandatory, not optional.
- No closed-form 3D analogue of `q_r(n)` exists; any generalised face number is a
  documented choice. Default: do not use one; work with `Q_iso` directly.
- `eps*` is composition-independent but has NOT been shown to resolve the
  transition; the 2D control was consistent with a featureless power law.
- `f_c` at fixed tolerance is unmeasurable (exact contact has measure zero) and
  its `f_c(delta)` repair is degenerate with the median gap.
- MC dynamics is a Brownian proxy; `eta_a` from it is not the EDMD `eta_a`.
- `g(sigma+)` linear extrapolation degrades above `eta ~ 0.4` (4% at 0.45).
- `eps* = C (eta_m - eta)^p` is the wrong functional form near contact
  percolation; window spread 0.15. A free-volume-motivated form is needed.
- Whether any of these structures is thermodynamically relevant at `T > 0` is
  entirely open, in 3D as in 2D.

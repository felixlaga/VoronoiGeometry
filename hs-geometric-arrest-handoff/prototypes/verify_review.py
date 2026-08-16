"""
Verify / refute the automated review's claims, and push the mathematics further.

 A. Counting theorem eta_K = z0/(z0+K) eta_cp  <->  my ladder (algebraic identity).
 B. The four modular constructions: uniform K, independence, tangential+congruent,
    face p-vectors, exact eta.  (Direct check, no trust.)
 C. NEW: spectral bound.  A uniform-K independent vacancy set is an equitable
    2-partition with quotient matrix [[0, z0],[K, z0-K]], whose second eigenvalue
    is -K.  For an infinite lattice, -K must lie in the Bloch spectrum of the
    contact graph.  FCC band = [-4, 12]  =>  K <= 4  (exact upper bound).
    Triangular band = [-3, 6]  =>  K <= 3  (2D).
 D. CONSEQUENCE: my earlier 2D claim of an index-2 rung (phi = 0.4535) is WRONG:
    it needs K = 6 > 3 and z = 0.  The 2D family is {honeycomb, kagome, maple-leaf}.
    Verify kagome (K=2) tangential with phi = sqrt(3) pi / 8 exactly, in 2D.
 E. NEW: independence is FORCED, not assumed: adjacent vacancies break
    tangentiality.  Numeric counterexample.
 F. HCP: compute the Bloch bands of the HCP contact graph and check which -K lie
    inside them (necessary condition for a uniform-K vacancy state on HCP).
 G. NEW: the counting theorem needs no periodicity.  Count ALL uniform-K
    independent sets on small tori (2D K=2; 3D K=4) to probe whether the states
    are unique-up-to-symmetry or extensive (entropy question).
"""
import itertools
import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection
from voronoi_ground_states import voronoi_cell, analyse

s3 = np.sqrt(3.0)

# ------------------------------------------------------------------ A
print("=" * 100)
print("A. COUNTING THEOREM vs MY LADDER (exact algebra)")
print("=" * 100)
eta_fcc = np.pi / np.sqrt(18)
for K in (4, 3, 2, 1):
    mine = (1 - 1 / (12 / K + 1)) * eta_fcc          # my index n = 12/K + 1
    theirs = 2 * np.sqrt(2) * np.pi / (12 + K)
    print(f"  K={K}: index n={12//K+1:2d}   mine={mine:.12f}   2sqrt2 pi/(12+K)={theirs:.12f}"
          f"   diff={abs(mine-theirs):.1e}")

# ------------------------------------------------------------------ B
print()
print("=" * 100)
print("B. MODULAR CONSTRUCTIONS (their table): uniform K, independent, tangential, congruent")
print("=" * 100)

def fcc_torus(n):
    pts = []
    for x in range(2 * n):
        for y in range(2 * n):
            for z in range(2 * n):
                if (x + y + z) % 2 == 0:
                    pts.append((x, y, z))
    return np.array(pts, float), 2.0 * n

RULES = {
    4: lambda x, y, z: x % 2 == 0 and y % 2 == 0 and z % 2 == 0,
    3: lambda x, y, z: (y + 2 * z) % 5 == 0,
    2: lambda x, y, z: (x + 2 * y + 3 * z) % 7 == 0,
    1: lambda x, y, z: (x + 3 * y + 4 * z) % 13 == 0,
}

def face_pvector(pos_i, pts, L):
    """p-vector (sorted polygon edge counts) of the periodic Voronoi cell."""
    sh = np.array(list(itertools.product((-1, 0, 1), repeat=3))) * L
    big = (pts[None] + sh[:, None]).reshape(-1, 3)
    d = big - pos_i
    dist = np.linalg.norm(d, axis=1)
    m = (dist > 1e-9) & (dist < 4.0)
    d, dist = d[m], dist[m]
    nr = d / dist[:, None]
    hs = np.hstack([nr, -(dist / 2)[:, None]])
    it = HalfspaceIntersection(hs, np.zeros(3))
    hull = ConvexHull(it.intersections)
    faces = {}
    for eq, simp in zip(hull.equations, hull.simplices):
        key = tuple(np.round(eq, 6))
        faces.setdefault(key, set()).update(
            tuple(np.round(it.intersections[v], 6)) for v in simp)
    pv = sorted(len(v) for v in faces.values() if len(v) >= 3)
    return pv, hull.volume

nn3 = [v for v in itertools.product((-1, 0, 1), repeat=3)
       if sorted(map(abs, v)) == [0, 1, 1]]

for K in (4, 3, 2, 1):
    n = {4: 2, 3: 5, 2: 7, 1: 13}[K]        # torus size = one period
    pts, L = fcc_torus(n if n <= 5 else 4) if False else fcc_torus(min(n, 4))
    pts, L = fcc_torus({4: 2, 3: 5, 2: 7, 1: 13}[K] // 1 and {4: 2, 3: 5, 2: 7, 1: 4}[K])
    rule = RULES[K]
    vac = np.array([rule(int(x) % ({4:2,3:5,2:7,1:13}[K]) if False else int(x),
                         int(y), int(z)) for x, y, z in pts])
    vac = np.array([rule(int(x), int(y), int(z)) for x, y, z in pts])
    occ = pts[~vac]; vpts = pts[vac]
    fv = len(vpts) / len(pts)
    # uniformity + independence on the torus (mod L)
    def wrap(v): return (v + L / 2) % L - L / 2
    kcounts = []
    for p in occ[:200]:
        c = 0
        for v in nn3:
            q = (p + v) % L
            if any(np.all(np.abs(wrap(q - w)) < 1e-9) for w in vpts):
                c += 1
        kcounts.append(c)
    indep = all(not any(np.all(np.abs(wrap((w + v) % L - w2)) < 1e-9)
                for v in nn3 for w2 in vpts) for w in vpts[:60])
    # tangential + congruent via periodic cells of a few occupied sites
    pv0, V0 = face_pvector(occ[0], occ, L)
    same = all(face_pvector(occ[i], occ, L)[0] == pv0 for i in
               np.random.default_rng(1).choice(len(occ), size=min(6, len(occ)), replace=False))
    eta = len(occ) * (np.pi / 6) / (L ** 3)   # sigma = 1 (nn distance sqrt2 scaled)
    sigma = np.sqrt(2.0)
    eta = len(occ) * (np.pi / 6) * sigma ** 3 / L ** 3
    ref = 2 * np.sqrt(2) * np.pi / (12 + K)
    print(f"  K={K}: f_v={fv:.4f} (=1/{round(1/fv)})  uniform-K={set(kcounts)=={K}}  "
          f"independent={indep}  p-vector={pv0}  congruent(sample)={same}  "
          f"eta={eta:.9f}  exact={ref:.9f}  diff={abs(eta-ref):.1e}")

# ------------------------------------------------------------------ C
print()
print("=" * 100)
print("C. SPECTRAL BOUND (new): -K must lie in the Bloch band of the contact graph")
print("=" * 100)
q = np.random.default_rng(0).uniform(-np.pi, np.pi, (200000, 3))
lam_fcc = 4 * (np.cos(q[:, 0]) * np.cos(q[:, 1]) + np.cos(q[:, 1]) * np.cos(q[:, 2])
               + np.cos(q[:, 2]) * np.cos(q[:, 0]))
print(f"  FCC band numeric: [{lam_fcc.min():.4f}, {lam_fcc.max():.4f}]  (exact [-4, 12])"
      f"   =>  K <= 4  ==> the 3D ladder {{4,3,2,1}} is COMPLETE within this class")
q2 = np.random.default_rng(1).uniform(-np.pi, np.pi, (200000, 2))
a1 = np.array([1, 0]); a2 = np.array([0.5, s3 / 2])
lam_tri = 2 * (np.cos(q2 @ a1) + np.cos(q2 @ a2) + np.cos(q2 @ (a1 - a2)))
print(f"  triangular band:  [{lam_tri.min():.4f}, {lam_tri.max():.4f}]  (exact [-3, 6])"
      f"   =>  K <= 3  ==> 2D family is {{K=3,2,1}} ONLY")
print("  => my earlier 2D 'index-2, phi=0.453450' rung REQUIRES K=6 and z=0: RETRACTED.")

# ------------------------------------------------------------------ D
print()
print("=" * 100)
print("D. 2D FAMILY VERIFIED DIRECTLY: honeycomb (K=3), kagome (K=2), maple-leaf (K=1)")
print("=" * 100)
A1 = np.array([1.0, 0.0]); A2 = np.array([0.5, s3 / 2])
def tri_torus(n):
    return np.array([i * A1 + j * A2 for i in range(n) for j in range(n)]), n
NN2 = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
RULES2 = {3: (3, lambda i, j: (i - j) % 3 == 0),
          2: (4, lambda i, j: i % 2 == 0 and j % 2 == 0),
          1: (7, lambda i, j: (i + 3 * j) % 7 == 0)}
NAME2 = {3: "honeycomb ", 2: "kagome    ", 1: "maple-leaf"}
def cell2d(p, occ, Lv1, Lv2):
    sh = np.array([a * Lv1 + b * Lv2 for a in (-1, 0, 1) for b in (-1, 0, 1)])
    big = (occ[None] + sh[:, None]).reshape(-1, 2)
    d = big - p; dist = np.linalg.norm(d, axis=1)
    m = (dist > 1e-9) & (dist < 4.0)
    d, dist = d[m], dist[m]
    nr = d / dist[:, None]
    it = HalfspaceIntersection(np.hstack([nr, -(dist / 2)[:, None]]), np.zeros(2))
    hull = ConvexHull(it.intersections)
    A, P = hull.volume, hull.area
    fdists = sorted(set(np.round(dist[np.unique(
        [np.argmin(np.abs(nr @ e[:2] - 1) + np.abs(dist / 2 + e[2]))
         for e in hull.equations])] / 2, 9)))
    return A, P, fdists
for K, (n, rule) in RULES2.items():
    N = 2 * n if K == 2 else n
    pts, _ = tri_torus(N)
    idx = np.array([[i, j] for i in range(N) for j in range(N)])
    vac = np.array([rule(i, j) for i, j in idx])
    occ = pts[~vac]
    Lv1, Lv2 = N * A1, N * A2
    A, P, fd = cell2d(occ[0], occ, Lv1, Lv2)
    qi = 4 * np.pi * A / P ** 2
    tang = max(fd) - min(fd) < 1e-9 and abs(min(fd) - 0.5) < 1e-9
    box = abs(np.cross(Lv1, Lv2))
    phi = len(occ) * np.pi * 0.25 / box
    ref = {3: np.pi / (3 * s3), 2: s3 * np.pi / 8, 1: s3 * np.pi / 7}[K]
    print(f"  K={K} ({NAME2[K]}): f_v=1/{round(len(pts)/vac.sum())}  tangential={tang}"
          f"  q={qi:.9f}  phi={phi:.9f}  exact={ref:.9f}  diff={abs(phi-ref):.1e}")

# ------------------------------------------------------------------ E
print()
print("=" * 100)
print("E. INDEPENDENCE IS FORCED: adjacent vacancies destroy tangentiality (FCC, numeric)")
print("=" * 100)
pts, L = fcc_torus(4)
vacpair = [np.array([4., 4., 4.]), np.array([5., 5., 4.])]     # nearest neighbours
occ = np.array([p for p in pts
                if not any(np.all(np.abs(p - v) < 1e-9) for v in vacpair)])
worst = 0.0
for p in occ:
    if min(np.linalg.norm(p - v) for v in vacpair) < 1.5:      # cells around the divacancy
        c = voronoi_cell(p, L * np.eye(3), occ, nshell=1)
        spread = (c["face_dists"].max() - c["face_dists"].min()) / c["face_dists"].min()
        worst = max(worst, spread)
print(f"  max face-distance spread around a divacancy: {worst:.4f}"
      f"   (tangential requires 0; single vacancy gives 0 exactly)")
print("  => the review's 'assumption (1)' is not an assumption: tangentiality implies it.")

# ------------------------------------------------------------------ F
print()
print("=" * 100)
print("F. HCP CONTACT GRAPH: Bloch bands and which K survive the necessary condition")
print("=" * 100)
h = np.sqrt(2.0 / 3.0)
hvec = np.array([[0, 1 / s3], [-0.5, -1 / (2 * s3)], [0.5, -1 / (2 * s3)]])
Q = np.random.default_rng(2).uniform(-2 * np.pi, 2 * np.pi, (400000, 3))
f = 2 * (np.cos(Q[:, :2] @ a1) + np.cos(Q[:, :2] @ a2) + np.cos(Q[:, :2] @ (a1 - a2)))
g = 2 * np.cos(Q[:, 2] * h) * np.abs(np.exp(1j * (Q[:, :2] @ hvec.T)).sum(axis=1))
lo, hi = f - g, f + g
allv = np.concatenate([lo, hi])
print(f"  HCP bands sampled: [{allv.min():.4f}, {allv.max():.4f}]")
for K in (1, 2, 3, 4):
    ok = np.any(np.abs(allv + K) < 0.02)
    print(f"    -K = {-K}: in band = {ok}")
print("  (a value in-band is necessary, not sufficient; existence on HCP stays open)")

# ------------------------------------------------------------------ G
print()
print("=" * 100)
print("G. ENTROPY PROBE: exhaustive count of uniform-K independent sets on small tori")
print("=" * 100)

def enumerate_2d(N, K):
    idx = [(i, j) for i in range(N) for j in range(N)]
    pos = {p: k for k, p in enumerate(idx)}
    nbr = [[pos[((i + a) % N, (j + b) % N)] for a, b in NN2] for i, j in idx]
    n = len(idx)
    state = [-1] * n            # -1 undecided, 0 occupied, 1 vacant
    sols = []
    def dfs(s):
        if s == n:
            sols.append(tuple(state)); return
        for val in (0, 1):
            state[s] = val
            ok = True
            if val == 1:
                if any(state[t] == 1 for t in nbr[s]): ok = False
            cells = set([s] + nbr[s])
            for c in cells:
                if not ok: break
                if state[c] == 0:
                    vdec = sum(1 for t in nbr[c] if state[t] == 1)
                    und = sum(1 for t in nbr[c] if state[t] == -1)
                    if vdec > K or vdec + und < K: ok = False
                if state[c] == 1:
                    if any(state[t] == 1 for t in nbr[c]): ok = False
            if ok: dfs(s + 1)
            state[s] = -1
    dfs(0)
    # orbits under torus translations
    def translate(sol, a, b):
        out = [0] * n
        for k, (i, j) in enumerate(idx):
            out[pos[((i + a) % N, (j + b) % N)]] = sol[k]
        return tuple(out)
    orbits = set()
    for sol in sols:
        orbits.add(min(translate(sol, a, b) for a in range(N) for b in range(N)))
    # coset test: vacancy differences closed under the torus group?
    coset = 0
    for o in orbits:
        vs = [idx[k] for k, v in enumerate(o) if v == 1]
        vset = set(vs)
        i0, j0 = vs[0]
        shifted = {((i - i0) % N, (j - j0) % N) for i, j in vs}
        closed = all((((x1 + x2) % N, (y1 + y2) % N) in shifted)
                     for x1, y1 in shifted for x2, y2 in shifted)
        coset += closed
    return len(sols), len(orbits), coset

for N, K in [(4, 2), (6, 3), (6, 2)]:
    ns, no, nc = enumerate_2d(N, K)
    print(f"  2D triangular {N}x{N} torus, K={K}:  solutions={ns:5d}  "
          f"translation-orbits={no:3d}  of which lattice-cosets={nc}")
print("  orbits > cosets  would mean non-crystalline uniform-K states exist")
print("  => bears directly on the configurational-entropy objection")

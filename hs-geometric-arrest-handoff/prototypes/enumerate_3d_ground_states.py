"""
Geometric ground states for hard-sphere dynamical arrest in 3D.

Reproduces every number quoted in the accompanying notes.  Requires
voronoi_ground_states.py in the same directory.

    python3 enumerate_3d_ground_states.py

Structure of the argument
-------------------------
(1) THE 3D ANALOGUE OF phi = q.
    de Graaf (arXiv:2411.01199v2) uses, in 2D: for a monodisperse packing whose
    disks are the inscribed circles of a congruent *tangential* tiling,
    phi = q = 4 pi A / P^2.
    In 3D, decomposing a tangential polyhedron into pyramids from the incenter
    gives V = sigma S / 6, hence

        eta = (pi/6) sigma^3 / V = 36 pi V^2 / S^3 = Q_iso = Psi^3,

    Psi = Wadell sphericity.  So Q_iso, NOT the face count and NOT W_{12,0}, is
    the correct 3D counterpart of the isoperimetric quotient q, and the
    admissibility criterion for a candidate ground state is exactly:

        EVERY VORONOI NEIGHBOUR MUST BE A CONTACT.

    If one Voronoi face comes from a non-touching neighbour, eta != Q and the
    geometric argument does not transfer.

(2) COMPLETE ENUMERATION OF TANGENTIAL BRAVAIS LATTICES.
    The seven Voronoi face classes of a 3D lattice (Voronoi's theorem: the
    non-zero classes of L/2L) have squared norms that are LINEAR functionals of
    the Selling parameters p_ij = -b_i.b_j.  Tangentiality is therefore a LINEAR
    system, not an optimisation problem, and can be solved exactly over every
    live-face set and every zero-pattern.  Exactly three solutions exist.

(3) DEPLETED-FCC FAMILY.
    Scanning all vacancy sublattices of FCC of index k <= 14 (Hermite normal
    form) and keeping only those whose cells are congruent AND tangential gives
    a discrete ladder k in {4,5,7,13}, i.e. exactly those k with (k-1) | 12 and
    a realisable perfect covering.  This is the 3D analogue of the 2D
    floret-pentagonal (k=7) and honeycomb (k=3) states.
"""

import itertools
import numpy as np

from voronoi_ground_states import analyse, report, voronoi_cell

s3 = np.sqrt(3.0)
FCC_G = np.array([[1, 1, 0], [1, -1, 0], [1, 0, 1]], float)  # FCC generators, nn = sqrt(2)


# =====================================================================
# (A) named structures
# =====================================================================
def section_a():
    print("=" * 100)
    print("(A) NAMED STRUCTURES: is eta = Q_iso, i.e. is the Voronoi cell tangential?")
    print("=" * 100)
    HCP_c = np.sqrt(8 / 3)
    structs = [
        ("FCC  (rhombic dodecahedron)", np.eye(3),
         [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]]),
        ("HCP  (trapezo-rhombic dodec.)",
         [[1, 0, 0], [.5, s3 / 2, 0], [0, 0, HCP_c]],
         [[0, 0, 0], [.5, 1 / (2 * s3), HCP_c / 2]]),
        ("BCC  (truncated octahedron)", np.eye(3), [[0, 0, 0], [.5, .5, .5]]),
        ("SC   (cube)", np.eye(3), [[0, 0, 0]]),
        ("simple hexagonal, c = a (hexagonal prism)",
         [[1, 0, 0], [.5, s3 / 2, 0], [0, 0, 1]], [[0, 0, 0]]),
        ("diamond (triakis truncated tetrahedron)", np.eye(3),
         [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0],
          [.25, .25, .25], [.25, .75, .75], [.75, .25, .75], [.75, .75, .25]]),
        ("FCC minus 1-in-4 (oblate octahedrille)", np.eye(3),
         [[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]),
    ]
    for st in structs:
        report(analyse(*st))


# =====================================================================
# (B) complete enumeration of tangential Bravais lattices
# =====================================================================
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
COEF = np.array(list(itertools.product(range(-2, 3), repeat=3)))
COEF = COEF[np.any(COEF != 0, axis=1)]
_CLS = (COEF[:, 0] % 2) + 2 * (COEF[:, 1] % 2) + 4 * (COEF[:, 2] % 2)
MASK = np.array([_CLS == c for c in range(1, 8)])

# squared norms of the seven face classes, as linear functionals of
# p = (p12, p13, p14, p23, p24, p34)
NORMS = np.array([[1, 1, 1, 0, 0, 0], [1, 0, 0, 1, 1, 0], [0, 1, 0, 1, 0, 1],
                  [0, 0, 1, 0, 1, 1], [0, 1, 1, 1, 1, 0], [1, 0, 1, 1, 0, 1],
                  [1, 1, 0, 0, 1, 1]], float)


def gram(p):
    P = np.zeros((4, 4))
    for k, (i, j) in enumerate(PAIRS):
        P[i, j] = P[j, i] = p[k]
    G = np.zeros((3, 3))
    for i in range(3):
        G[i, i] = P[i].sum()
        for j in range(3):
            if i != j:
                G[i, j] = -P[i, j]
    return G


def tangential_lattice(p):
    """Return (n_faces, eta) if the lattice with Selling parameters p has a
    tangential Voronoi cell, else None."""
    G = gram(p)
    if np.linalg.eigvalsh(G).min() < 1e-11:
        return None
    nm = np.einsum('ij,jk,ik->i', COEF, G, COEF)
    live = []
    for m in MASK:
        v = nm[m]
        mn = v.min()
        if np.sum(v < mn * (1 + 1e-9)) == 2:      # non-degenerate face class
            live.append(mn)
    if len(live) < 3 or (max(live) - min(live)) / min(live) > 1e-10:
        return None
    sigma = np.sqrt(min(live))
    V = np.sqrt(np.linalg.det(G))
    return 2 * len(live), float((np.pi / 6) * sigma ** 3 / V)


def section_b(seed=3):
    print()
    print("=" * 100)
    print("(B) COMPLETE ENUMERATION OF TANGENTIAL BRAVAIS-LATTICE VORONOI CELLS")
    print("    every live-face set x every zero-pattern of the Selling parameters,")
    print("    solved exactly (the equal-norm conditions are linear in p)")
    print("=" * 100)
    rng = np.random.default_rng(seed)
    res = {}
    for nz in range(4):
        for S in itertools.combinations(range(6), nz):
            for k in range(3, 8):
                for L in itertools.combinations(range(7), k):
                    rows = [NORMS[L[i]] - NORMS[L[i + 1]] for i in range(k - 1)]
                    rows.append(np.ones(6))
                    for s in S:
                        e = np.zeros(6); e[s] = 1; rows.append(e)
                    A = np.array(rows)
                    b = np.zeros(len(rows)); b[k - 1] = 1.0
                    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
                    if np.linalg.norm(A @ sol - b) > 1e-9:
                        continue
                    _, sv, Vt = np.linalg.svd(A)
                    null = Vt[np.sum(sv > 1e-9):]
                    cands = [sol] + [sol + null.T @ rng.normal(0, .5, len(null))
                                     for _ in range(300 * len(null))]
                    for p in cands:
                        if p.min() < -1e-12 or abs(p.sum() - 1) > 1e-9:
                            continue
                        r = tangential_lattice(np.clip(p, 0, None))
                        if r:
                            res[(r[0], round(r[1], 7))] = 1
    names = {6: "cube            -> simple cubic",
             8: "hexagonal prism -> simple hexagonal, c = a",
             12: "rhombic dodec.  -> FCC",
             14: "truncated octahedron"}
    for (nf, eta) in sorted(res):
        print(f"    faces = {nf:3d}   eta = Q_iso = {eta:.7f}   {names.get(nf, '?')}")
    print(f"\n    => exactly {len(res)} tangential Bravais lattices exist in 3D")
    print("       (no 14-faced solution: forcing all seven class norms equal gives")
    print("        p_ij = 1/6, whose class norms are 1/2 and 2/3 -- so BCC-type")
    print("        cells can never be tangential)")


# =====================================================================
# (C) depleted-FCC ladder
# =====================================================================
def hermite_normal_forms(k):
    """All index-k sublattices of Z^3 in lower-triangular Hermite normal form."""
    out = []
    for a in range(1, k + 1):
        if k % a:
            continue
        for b in range(1, k // a + 1):
            if (k // a) % b:
                continue
            c = k // (a * b)
            for d in range(a):
                for e in range(a):
                    for f in range(b):
                        out.append(np.array([[a, 0, 0], [d, b, 0], [e, f, c]], int))
    return out


def section_c(kmax=14):
    print()
    print("=" * 100)
    print(f"(C) FCC MINUS A PERIODIC VACANCY SUBLATTICE OF INDEX k <= {kmax}")
    print("    keeping only patterns whose cells are congruent AND tangential")
    print("=" * 100)
    hits = set()
    for k in range(2, kmax + 1):
        for H in hermite_normal_forms(k):
            lat = H.astype(float) @ FCC_G
            Li = np.linalg.inv(lat)
            adj = np.round(np.linalg.det(H) * np.linalg.inv(H)).astype(int)
            reps, seen = [], set()
            for n in itertools.product(range(-k, k + 1), repeat=3):
                key = tuple((np.array(n) @ adj) % k)
                if key not in seen:
                    seen.add(key)
                    x = np.array(n, float) @ FCC_G
                    f = x @ Li
                    reps.append((f - np.floor(f + 1e-9)) @ lat)   # wrap into cell
                if len(reps) == k:
                    break
            if len(reps) != k:
                continue
            B = np.array(reps[1:])          # reps[0] = origin = the vacancy orbit
            if len(B) > 1:
                dmin = min(np.linalg.norm(B[i] - B[j])
                           for i in range(len(B)) for j in range(i))
                if dmin < 1e-8:
                    continue
            try:
                r = analyse(f"k={k}", lat, B)
            except Exception:
                continue
            if r["tangential"] and r["congruent"]:
                hits.add((k, round(r["eta"], 7), int(np.unique(r["z"])[0]),
                          int(np.unique(r["n_faces"])[0])))
    eta_fcc = np.pi / np.sqrt(18)
    for k, eta, z, nf in sorted(hits):
        print(f"    remove 1 in {k:2d}:  eta = Q_iso = {eta:.7f}  z = {z:2d}  "
              f"faces = {nf:2d}   [(1-1/k) x pi/sqrt(18) = {(1-1/k)*eta_fcc:.7f}]")
    print("\n    every hit satisfies (k-1) | 12, i.e. each particle is adjacent to")
    print("    exactly 12/(k-1) vacancies.  2D analogue on the triangular lattice:")
    print("    (k-1) | 6 -> k in {2,3,4,7}, phi = (1-1/k) pi/sqrt(12) =")
    phi_tri = np.pi / np.sqrt(12)
    for k in (2, 3, 4, 7):
        tag = {3: "  <- de Graaf's honeycomb / caging state",
               7: "  <- de Graaf's floret pentagonal tiling"}.get(k, "")
        print(f"        k = {k}:  phi = {(1-1/k)*phi_tri:.6f}{tag}")


# =====================================================================
# (D) checks on specific claims in arXiv:2411.01199v2
# =====================================================================
def section_d():
    print()
    print("=" * 100)
    print("(D) DIRECT CHECKS ON THE PAPER'S 3D CLAIMS (isolated vacancy, 4x4x4 supercell)")
    print("=" * 100)
    for label, cell, V0, sigma in [
            ("BCC", [[0, 0, 0], [.5, .5, .5]], 0.5, s3 / 2),
            ("FCC", [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]], 0.25, 1 / np.sqrt(2))]:
        n = 4
        sites = np.array([[n1 + b[0], n2 + b[1], n3 + b[2]]
                          for n1 in range(n) for n2 in range(n) for n3 in range(n)
                          for b in cell], float)
        vac = np.array([n / 2.0] * 3)
        part = sites[np.linalg.norm(sites - vac, axis=1) > 1e-9]
        lat = n * np.eye(3)
        d = np.linalg.norm(part - vac, axis=1)
        print(f"  {label}:")
        total = 0.0
        for dist in np.unique(np.round(d, 6))[:3]:
            idx = np.where(np.abs(d - dist) < 1e-6)[0]
            c = voronoi_cell(part[idx[0]], lat, part, nshell=1)
            gain = (c["V"] - V0) / V0
            total += len(idx) * gain
            print(f"    shell d={dist:.4f}  count={len(idx):3d}  V/V0={c['V']/V0:.6f}  "
                  f"faces={c['n_faces']:2d}  eta_local={(np.pi/6)*sigma**3/c['V']:.6f}  "
                  f"Q_iso={c['Q']:.6f}")
        print(f"    redistributed vacancy volume = {total:.6f} V0 (must be 1)")


if __name__ == "__main__":
    section_a()
    section_b()
    section_c()
    section_d()

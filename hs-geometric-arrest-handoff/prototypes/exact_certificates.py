"""
EXACT-ARITHMETIC CERTIFICATES (independent of all previous float computations).

Idea: in integer lattice coordinates every bisector plane between sites p, q is
    (q - p) . x = (|q|^2 - |p|^2) / 2 ,
with rational coefficients.  Float geometry is used ONLY to locate which planes
and vertices are active; every claim is then verified in exact arithmetic
(sympy rationals and radicals):

  1. each vertex satisfies its defining planes exactly and violates no other
     plane (cell membership is exact);
  2. every active face plane is a CONTACT plane: distance from the site equals
     sigma/2 exactly  ->  tangentiality is an identity, not an approximation;
  3. V and S computed exactly by pyramid decomposition; the theorem
     eta = 36 pi V^2 / S^3 = (pi/6) sigma^3 / V  reduces to  S = 3 V / R,
     which is checked as an exact symbolic identity;
  4. the global density  eta = (pi/6) sigma^3 * (occupied sites / volume)
     equals the cell value exactly.

Cells certified: 3D FCC-depletion K = 4 and K = 1 (the two endpoints);
2D honeycomb (K=3), kagome (K=2), maple-leaf (K=1), with sqrt(3) kept symbolic.
"""
import itertools
import numpy as np
import sympy as sp
from scipy.spatial import HalfspaceIntersection, ConvexHull

R3 = sp.sqrt(3)


def certify(site, others, sigma2_exact, dim, label):
    """site, others: exact sympy coordinate tuples.  sigma2_exact: exact squared
    contact distance.  Returns after printing the certificate result."""
    sitef = np.array([float(c) for c in site])
    othf = np.array([[float(c) for c in o] for o in others])
    # float pass: which planes are active, where are the vertices
    d = othf - sitef
    dist = np.linalg.norm(d, axis=1)
    nr = d / dist[:, None]
    hs = np.hstack([nr, -(dist / 2)[:, None]])
    it = HalfspaceIntersection(hs, np.zeros(dim))
    hull = ConvexHull(it.intersections)
    vertsf = it.intersections
    # identify active planes per vertex (float, generous tol) then verify exact
    planes = []          # (normal_exact_vector, rhs_exact) in site-centred coords
    for o in others:
        n = sp.Matrix([oc - sc for oc, sc in zip(o, site)])
        rhs = sum((oc - sc) ** 2 for oc, sc in zip(o, site)) / sp.Integer(2)
        planes.append((n, sp.nsimplify(rhs)))
    # exact vertex coordinates: solve the defining plane triples symbolically,
    # seeded by the float vertex
    exact_verts = []
    for v in vertsf:
        resid = np.abs(np.array([[float(c) for c in p[0]] for p in planes]) @ v
                       - np.array([float(p[1]) for p in planes]))
        act = np.argsort(resid)[:dim]
        A = sp.Matrix([[planes[a][0][k] for k in range(dim)] for a in act])
        b = sp.Matrix([planes[a][1] for a in act])
        if A.det() == 0:
            # take next-best combination
            for combo in itertools.combinations(np.argsort(resid)[:dim + 3], dim):
                A = sp.Matrix([[planes[a][0][k] for k in range(dim)] for a in combo])
                b = sp.Matrix([planes[a][1] for a in combo])
                if A.det() != 0 and max(abs(float(x) - y) for x, y in
                                        zip(A.solve(b), v)) < 1e-6:
                    break
        ve = A.solve(b)
        assert max(abs(float(ve[k]) - v[k]) for k in range(dim)) < 1e-7, "vertex id failed"
        exact_verts.append(sp.Matrix(ve))
    # (1) membership: every exact vertex satisfies every plane n.x <= rhs EXACTLY
    for ve in exact_verts:
        for n, rhs in planes:
            val = sp.simplify((n.T * ve)[0] - rhs)
            assert sp.simplify(val) <= 0 or sp.simplify(val) == 0, "membership violated"
    # active faces = planes met with equality by >= dim vertices
    active = []
    for pi, (n, rhs) in enumerate(planes):
        onv = [ve for ve in exact_verts if sp.simplify((n.T * ve)[0] - rhs) == 0]
        if len(onv) >= dim:
            active.append((pi, n, rhs, onv))
    # (2) tangency: for each active plane, distance^2 from site = rhs^2/|n|^2
    # must equal sigma^2/4 exactly; equivalently |q-p|^2 = sigma^2 (contact)
    contact_ok = all(
        sp.simplify(sum((others[pi][k] - site[k]) ** 2 for k in range(dim))
                    - sigma2_exact) == 0
        for pi, *_ in active)
    # (3) exact V, S by pyramid/fan decomposition about the site (origin-shifted)
    Vtot = sp.Integer(0); Stot = sp.Integer(0)
    for pi, n, rhs, onv in active:
        pts = [sp.Matrix(ve) for ve in onv]
        if dim == 2:
            assert len(pts) == 2
            e = pts[1] - pts[0]
            face_meas = sp.sqrt(sp.simplify((e.T * e)[0]))
        else:
            # order polygon vertices by angle in the face plane (float order, exact pts)
            c = sum(pts, sp.zeros(3, 1)) / len(pts)
            nf = np.array([float(x) for x in n]); nf /= np.linalg.norm(nf)
            ref = np.array([float(x) for x in (pts[0] - c)])
            ref -= ref.dot(nf) * nf
            ref /= np.linalg.norm(ref)
            per = np.cross(nf, ref)
            ang = []
            for p in pts:
                w = np.array([float(x) for x in (p - c)])
                ang.append(np.arctan2(w.dot(per), w.dot(ref)))
            order = np.argsort(ang)
            pts = [pts[i] for i in order]
            face_meas = sp.Integer(0)
            for i in range(1, len(pts) - 1):
                u = pts[i] - pts[0]; w = pts[i + 1] - pts[0]
                cr = u.cross(w)
                face_meas += sp.sqrt(sp.simplify((cr.T * cr)[0])) / 2
        h = rhs / sp.sqrt(sp.simplify((n.T * n)[0]))       # exact distance site->plane
        Stot += face_meas
        Vtot += face_meas * h / dim
    Vtot = sp.simplify(sp.nsimplify(Vtot, [sp.sqrt(2), sp.sqrt(3)]))
    Stot = sp.simplify(sp.nsimplify(Stot, [sp.sqrt(2), sp.sqrt(3)]))
    R = sp.sqrt(sigma2_exact) / 2
    ident = sp.simplify(Stot - dim * Vtot / R)             # S = d V / R  <=> tangential
    if dim == 3:
        eta_cell = sp.simplify(sp.pi / 6 * sp.sqrt(sigma2_exact) ** 3 / Vtot)
    else:
        eta_cell = sp.simplify(sp.pi * sigma2_exact / 4 / Vtot)
    print(f"  {label}")
    print(f"    active faces: {len(active)}   all contact faces (exact): {contact_ok}")
    print(f"    V = {Vtot}    S = {Stot}")
    print(f"    identity S - {dim}V/R = {ident}   (0 <=> tangential, EXACT)")
    print(f"    eta = Q_iso (exact) = {eta_cell} = {sp.nsimplify(eta_cell)}"
          f" = {float(eta_cell):.12f}")
    return eta_cell


print("=" * 96)
print("EXACT CERTIFICATES, 3D (integer FCC coordinates, sigma^2 = 2)")
print("=" * 96)
# K=4: vacancies = all-even; site (1,1,0); neighbours from a 2-shell integer ball
def fcc_ball(center, r=3):
    out = []
    for dx in range(-r, r + 1):
        for dy in range(-r, r + 1):
            for dz in range(-r, r + 1):
                p = (center[0] + dx, center[1] + dy, center[2] + dz)
                if (p[0] + p[1] + p[2]) % 2 == 0 and (dx, dy, dz) != (0, 0, 0):
                    out.append(p)
    return out

site = (1, 1, 0)
others = [p for p in fcc_ball(site) if not (p[0] % 2 == 0 and p[1] % 2 == 0 and p[2] % 2 == 0)]
e4 = certify(tuple(sp.Integer(c) for c in site),
             [tuple(sp.Integer(c) for c in o) for o in others],
             sp.Integer(2), 3, "K=4 (square bipyramid), site (1,1,0), vacancies all-even")
print(f"    exact target pi/(4 sqrt 2): match = "
      f"{sp.simplify(e4 - sp.pi/(4*sp.sqrt(2))) == 0}")

site = (1, 1, 0)          # residue x+3y+4z = 4 mod 13 (occupied)
others = [p for p in fcc_ball(site)
          if (p[0] + 3 * p[1] + 4 * p[2]) % 13 != 0]
e1 = certify(tuple(sp.Integer(c) for c in site),
             [tuple(sp.Integer(c) for c in o) for o in others],
             sp.Integer(2), 3, "K=1 (11-faced cell), vacancies x+3y+4z=0 mod 13")
print(f"    exact target 2 sqrt2 pi/13: match = "
      f"{sp.simplify(e1 - 2*sp.sqrt(2)*sp.pi/13) == 0}")

print()
print("=" * 96)
print("EXACT CERTIFICATES, 2D (triangular coords x=i+j/2, y=j sqrt3/2, sigma^2 = 1)")
print("=" * 96)
def tri_pt(i, j):
    return (sp.Integer(i) + sp.Rational(j, 2), sp.Integer(j) * R3 / 2)

def tri_ball(ci, cj, r=4):
    return [(ci + di, cj + dj) for di in range(-r, r + 1) for dj in range(-r, r + 1)
            if (di, dj) != (0, 0)]

CASES2 = [
    ("K=3 honeycomb, vacancies i-j=0 mod 3", (1, 0),
     lambda i, j: (i - j) % 3 == 0, sp.pi / (3 * R3)),
    ("K=2 kagome, vacancies i,j both even", (1, 0),
     lambda i, j: i % 2 == 0 and j % 2 == 0, R3 * sp.pi / 8),
    ("K=1 maple-leaf, vacancies i+3j=0 mod 7", (1, 0),
     lambda i, j: (i + 3 * j) % 7 == 0, R3 * sp.pi / 7),
]
for label, (ci, cj), rule, target in CASES2:
    others = [tri_pt(i, j) for i, j in tri_ball(ci, cj) if not rule(i, j)]
    e = certify(tri_pt(ci, cj), others, sp.Integer(1), 2, label)
    print(f"    exact target {target}: match = {sp.simplify(e - target) == 0}")

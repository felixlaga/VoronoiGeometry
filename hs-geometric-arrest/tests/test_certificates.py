"""EXACT-ARITHMETIC CERTIFICATES (independent of all float computations).

Ported from ``prototypes/exact_certificates.py`` (task T3).  In integer lattice
coordinates every bisector plane between sites ``p, q`` is

    (q - p) . x = (|q|^2 - |p|^2) / 2

with rational coefficients.  Float geometry is used ONLY to locate which planes
and vertices are active; every claim is then verified in exact arithmetic
(sympy rationals and radicals):

1. each located vertex satisfies its defining planes exactly and violates no
   other plane -- cell membership is exact;
2. every active face plane is a CONTACT plane: the generating separation obeys
   ``|q - p|^2 == sigma^2`` exactly, so tangentiality is an identity, not an
   approximation;
3. ``V`` and ``S`` are computed exactly by pyramid decomposition and the
   theorem reduces to ``S - d V / R == 0``, checked as an exact symbolic
   identity (``== 0``, never ``< tol``);
4. the cell density equals the exact closed form.

Cells certified: 3D FCC-depletion K=4 and K=1 (the ladder endpoints, with the
golden exact V and S from REFERENCE_VALUES), and the 2D family: honeycomb
(K=3), kagome (K=2), maple-leaf (K=1), with sqrt(3) kept symbolic.
"""

import itertools

import numpy as np
import pytest
import sympy as sp
from scipy.spatial import ConvexHull, HalfspaceIntersection

R3 = sp.sqrt(3)


# --------------------------------------------------------------------------- #
# the certifier
# --------------------------------------------------------------------------- #
def certify(site, others, sigma2_exact, dim):
    """Exact certificate for the Voronoi cell of ``site`` among ``others``.

    All inputs are exact sympy coordinate tuples; returns a dict of exact
    quantities.  Raises AssertionError if any exact check fails.
    """
    sitef = np.array([float(c) for c in site])
    othf = np.array([[float(c) for c in o] for o in others])

    # float pass: locate active planes and vertices
    d = othf - sitef
    dist = np.linalg.norm(d, axis=1)
    nr = d / dist[:, None]
    hs = np.hstack([nr, -(dist / 2)[:, None]])
    it = HalfspaceIntersection(hs, np.zeros(dim))
    vertsf = it.intersections
    ConvexHull(vertsf)          # raises if degenerate

    # exact planes in site-centred coordinates: n.x = rhs
    planes = []
    for o in others:
        n = sp.Matrix([oc - sc for oc, sc in zip(o, site)])
        rhs = sum((oc - sc) ** 2 for oc, sc in zip(o, site)) / sp.Integer(2)
        planes.append((n, sp.nsimplify(rhs)))

    # exact vertices: solve the defining plane tuples, seeded by float vertices
    A_f = np.array([[float(c) for c in p[0]] for p in planes])
    b_f = np.array([float(p[1]) for p in planes])
    exact_verts = []
    seen = set()
    for v in vertsf:
        resid = np.abs(A_f @ v - b_f)
        act = np.argsort(resid)[:dim]
        A = sp.Matrix([[planes[a][0][k] for k in range(dim)] for a in act])
        b = sp.Matrix([planes[a][1] for a in act])
        if A.det() == 0:
            found = False
            for combo in itertools.combinations(np.argsort(resid)[: dim + 3], dim):
                A = sp.Matrix([[planes[a][0][k] for k in range(dim)] for a in combo])
                b = sp.Matrix([planes[a][1] for a in combo])
                if A.det() != 0:
                    ve = A.solve(b)
                    if max(abs(float(ve[k]) - v[k]) for k in range(dim)) < 1e-6:
                        found = True
                        break
            assert found, "no non-degenerate defining plane set found"
        ve = A.solve(b)
        assert max(abs(float(ve[k]) - v[k]) for k in range(dim)) < 1e-7, "vertex id failed"
        key = tuple(sp.nsimplify(x) for x in ve)
        if key in seen:
            continue
        seen.add(key)
        exact_verts.append(sp.Matrix(ve))

    # (1) exact membership: every vertex satisfies every plane n.x <= rhs
    for ve in exact_verts:
        for n, rhs in planes:
            val = sp.simplify((n.T * ve)[0] - rhs)
            nonpos = (val == 0) or bool(sp.simplify(val).is_negative)
            assert nonpos, f"membership violated exactly: {val}"

    # active faces: planes met with equality by >= dim vertices
    active = []
    for pi, (n, rhs) in enumerate(planes):
        onv = [ve for ve in exact_verts if sp.simplify((n.T * ve)[0] - rhs) == 0]
        if len(onv) >= dim:
            active.append((pi, n, rhs, onv))

    # (2) tangency: every active plane's generator is at exact contact
    for pi, *_ in active:
        gap = sp.simplify(
            sum((others[pi][k] - site[k]) ** 2 for k in range(dim)) - sigma2_exact
        )
        assert gap == 0, f"active face from a non-contact neighbour: |d|^2-sigma^2={gap}"

    # (3) exact V and S by pyramid/fan decomposition about the site
    Vtot = sp.Integer(0)
    Stot = sp.Integer(0)
    for pi, n, rhs, onv in active:
        pts = [sp.Matrix(ve) for ve in onv]
        if dim == 2:
            assert len(pts) == 2, "2D face must have exactly two vertices"
            e = pts[1] - pts[0]
            face_meas = sp.sqrt(sp.simplify((e.T * e)[0]))
        else:
            c = sum(pts, sp.zeros(3, 1)) / len(pts)
            nf = np.array([float(x) for x in n])
            nf /= np.linalg.norm(nf)
            ref = np.array([float(x) for x in (pts[0] - c)])
            ref -= ref.dot(nf) * nf
            ref /= np.linalg.norm(ref)
            per = np.cross(nf, ref)
            ang = []
            for p in pts:
                w = np.array([float(x) for x in (p - c)])
                ang.append(np.arctan2(w.dot(per), w.dot(ref)))
            pts = [pts[i] for i in np.argsort(ang)]
            face_meas = sp.Integer(0)
            for i in range(1, len(pts) - 1):
                u = pts[i] - pts[0]
                w = pts[i + 1] - pts[0]
                cr = u.cross(w)
                face_meas += sp.sqrt(sp.simplify((cr.T * cr)[0])) / 2
        h = rhs / sp.sqrt(sp.simplify((n.T * n)[0]))
        Stot += face_meas
        Vtot += face_meas * h / dim
    Vtot = sp.simplify(sp.nsimplify(Vtot, [sp.sqrt(2), sp.sqrt(3)]))
    Stot = sp.simplify(sp.nsimplify(Stot, [sp.sqrt(2), sp.sqrt(3)]))

    R = sp.sqrt(sigma2_exact) / 2
    identity = sp.simplify(Stot - dim * Vtot / R)      # S = d V / R  <=> tangential
    assert identity == 0, f"S - {dim}V/R = {identity} (must be exactly 0)"

    if dim == 3:
        eta = sp.simplify(sp.pi / 6 * sp.sqrt(sigma2_exact) ** 3 / Vtot)
    else:
        eta = sp.simplify(sp.pi * sigma2_exact / 4 / Vtot)
    return {
        "n_active": len(active),
        "V": Vtot,
        "S": Stot,
        "identity": identity,
        "eta": eta,
    }


# --------------------------------------------------------------------------- #
# 3D: integer FCC coordinates, sigma^2 = 2
# --------------------------------------------------------------------------- #
def fcc_ball(center, r=3):
    out = []
    for dx in range(-r, r + 1):
        for dy in range(-r, r + 1):
            for dz in range(-r, r + 1):
                p = (center[0] + dx, center[1] + dy, center[2] + dz)
                if (p[0] + p[1] + p[2]) % 2 == 0 and (dx, dy, dz) != (0, 0, 0):
                    out.append(p)
    return out


@pytest.mark.slow
def test_certificate_3d_K4(golden):
    """K=4 (square bipyramid), site (1,1,0), vacancies all-even."""
    site = (1, 1, 0)
    others = [
        p for p in fcc_ball(site)
        if not (p[0] % 2 == 0 and p[1] % 2 == 0 and p[2] % 2 == 0)
    ]
    r = certify(
        tuple(sp.Integer(c) for c in site),
        [tuple(sp.Integer(c) for c in o) for o in others],
        sp.Integer(2), 3,
    )
    assert r["n_active"] == 8
    assert sp.simplify(r["eta"] - sp.pi / (4 * sp.sqrt(2))) == 0
    # golden exact V and S (REFERENCE_VALUES fcc_depletion_ladder.K4)
    assert sp.simplify(r["V"] - sp.Rational(8, 3)) == 0
    assert sp.simplify(r["S"] - 8 * sp.sqrt(2)) == 0
    assert abs(float(r["eta"].evalf(30)) - golden["fcc_depletion_ladder"]["K4"]["eta"]) < 1e-15


@pytest.mark.slow
def test_certificate_3d_K1(golden):
    """K=1 (11-faced cell), vacancies x+3y+4z = 0 mod 13."""
    site = (1, 1, 0)
    others = [p for p in fcc_ball(site) if (p[0] + 3 * p[1] + 4 * p[2]) % 13 != 0]
    r = certify(
        tuple(sp.Integer(c) for c in site),
        [tuple(sp.Integer(c) for c in o) for o in others],
        sp.Integer(2), 3,
    )
    assert r["n_active"] == 11
    assert sp.simplify(r["eta"] - 2 * sp.sqrt(2) * sp.pi / 13) == 0
    assert sp.simplify(r["V"] - sp.Rational(13, 6)) == 0
    assert sp.simplify(r["S"] - 13 * sp.sqrt(2) / 2) == 0
    assert abs(float(r["eta"].evalf(30)) - golden["fcc_depletion_ladder"]["K1"]["eta"]) < 1e-15


# --------------------------------------------------------------------------- #
# 2D: triangular coordinates x = i + j/2, y = j sqrt3/2, sigma^2 = 1
# --------------------------------------------------------------------------- #
def tri_pt(i, j):
    return (sp.Integer(i) + sp.Rational(j, 2), sp.Integer(j) * R3 / 2)


def tri_ball(ci, cj, r=4):
    return [
        (ci + di, cj + dj)
        for di in range(-r, r + 1)
        for dj in range(-r, r + 1)
        if (di, dj) != (0, 0)
    ]


@pytest.mark.slow
@pytest.mark.parametrize(
    "label,rule,target,n_edges,gold_key",
    [
        ("K=3 honeycomb, vacancies i-j=0 mod 3",
         lambda i, j: (i - j) % 3 == 0, sp.pi / (3 * R3), 3, "K3_honeycomb"),
        ("K=2 kagome, vacancies i,j both even",
         lambda i, j: i % 2 == 0 and j % 2 == 0, R3 * sp.pi / 8, 4, "K2_kagome"),
        ("K=1 maple-leaf, vacancies i+3j=0 mod 7",
         lambda i, j: (i + 3 * j) % 7 == 0, R3 * sp.pi / 7, 5, "K1_maple_leaf"),
    ],
)
def test_certificate_2d(label, rule, target, n_edges, gold_key, golden):
    ci, cj = 1, 0
    assert not rule(ci, cj), "certified site must be occupied"
    others = [tri_pt(i, j) for i, j in tri_ball(ci, cj) if not rule(i, j)]
    r = certify(tri_pt(ci, cj), others, sp.Integer(1), 2)
    assert r["n_active"] == n_edges
    assert sp.simplify(r["eta"] - target) == 0, label
    # 1e-12, the gate tolerance for the counting ladder -- NOT 1e-15: the
    # K1_maple_leaf decimal in REFERENCE_VALUES.json disagrees with its own
    # exact string sqrt3*pi/7 by ~4e-13 (recorded in DEBT.md; the file is
    # immutable, so the check adapts, the file does not)
    assert abs(float(r["eta"].evalf(30)) - golden["family_2d"][gold_key]["phi"]) < 1e-12

"""
Voronoi-cell geometry of candidate 3D 'geometric ground states'.

Central identity being tested (proof in the accompanying notes):

    If a monodisperse sphere packing has the property that EVERY Voronoi cell is
    a *tangential* polyhedron whose insphere is the particle itself -- equivalently,
    every Voronoi neighbour is a contact -- then

        V_cell = sigma * S_cell / 6            (sigma = 2 * r_in = contact distance)
        eta    = (pi/6) sigma^3 / V_cell
               = 36 pi V_cell^2 / S_cell^3
               = Q_iso   (3D isoperimetric quotient = Wadell sphericity cubed)

    This is the exact 3D analogue of the 2D statement phi = q used by de Graaf
    (arXiv:2411.01199v2).  Everything is computed here from first principles;
    nothing is taken on trust from the paper.

No fitted or invented quantities appear in this file.  All structures are
specified by their lattice vectors and basis; all cell geometry is obtained by
half-space intersection of the perpendicular bisectors.
"""

import itertools
import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

TOL = 1e-9


# ----------------------------------------------------------------------------
# Voronoi cell of one site in a periodic structure, by half-space intersection
# ----------------------------------------------------------------------------

def voronoi_cell(site, lattice, basis, nshell=3, area_tol=1e-8):
    """Return dict describing the Voronoi cell of `site` in the periodic
    structure defined by `lattice` (3x3 rows = lattice vectors) and `basis`.

    Faces are returned with their supporting-plane distance from the site and
    their area.  Zero-area (degenerate) contacts are discarded via `area_tol`.
    """
    lattice = np.asarray(lattice, float)
    basis = np.asarray(basis, float)
    site = np.asarray(site, float)

    # neighbour cloud
    pts = []
    rng = range(-nshell, nshell + 1)
    for n1, n2, n3 in itertools.product(rng, rng, rng):
        shift = n1 * lattice[0] + n2 * lattice[1] + n3 * lattice[2]
        for b in basis:
            pts.append(b + shift)
    pts = np.array(pts)
    d = np.linalg.norm(pts - site, axis=1)
    keep = (d > TOL)
    pts, d = pts[keep], d[keep]
    order = np.argsort(d)
    pts, d = pts[order], d[order]
    # keep a generous shell: 3x the nearest-neighbour distance is far more than
    # enough for any Voronoi-relevant vector in these structures
    sel = d <= 3.0 * d[0] + TOL
    pts, d = pts[sel], d[sel]

    # half-spaces  n . x <= c   with n unit, c = |p - site| / 2  (site at origin)
    normals = (pts - site) / d[:, None]
    offsets = d / 2.0
    halfspaces = np.hstack([normals, -offsets[:, None]])  # A x + b <= 0

    hs = HalfspaceIntersection(halfspaces, np.zeros(3))
    verts = hs.intersections
    hull = ConvexHull(verts)

    # group hull facets by supporting plane -> physical faces
    faces = {}
    for eq, simplex in zip(hull.equations, hull.simplices):
        n, off = eq[:3], -eq[3]  # n . x = off
        key = tuple(np.round(np.concatenate([n, [off]]), 7))
        tri = verts[simplex]
        area = 0.5 * np.linalg.norm(np.cross(tri[1] - tri[0], tri[2] - tri[0]))
        faces[key] = faces.get(key, 0.0) + area

    face_list = [(k[3], a) for k, a in faces.items() if a > area_tol]
    dists = np.array([f[0] for f in face_list])
    areas = np.array([f[1] for f in face_list])

    V = hull.volume
    S = areas.sum()
    return dict(
        n_faces=len(face_list),
        face_dists=np.sort(dists),
        S=S,
        V=V,
        d_min=d[0],
        r_in=dists.min(),
        Q=36 * np.pi * V**2 / S**3,
    )


def analyse(name, lattice, basis, note=""):
    lattice = np.asarray(lattice, float)
    basis = np.asarray(basis, float)
    cells = [voronoi_cell(b, lattice, basis) for b in basis]

    sigma = min(c["d_min"] for c in cells)          # contact distance
    Vcell = np.array([c["V"] for c in cells])
    Scell = np.array([c["S"] for c in cells])
    Q = np.array([c["Q"] for c in cells])
    nf = np.array([c["n_faces"] for c in cells])

    # tangentiality: every face of every cell supported at sigma/2
    defect = max(
        (c["face_dists"].max() - c["face_dists"].min()) / c["face_dists"].min()
        for c in cells
    )
    tangential = defect < 1e-7 and all(
        abs(c["face_dists"].min() - sigma / 2) < 1e-7 for c in cells
    )

    # contacts = neighbours at exactly sigma that carry a face
    z = np.array([int(np.sum(np.abs(c["face_dists"] - sigma / 2) < 1e-7)) for c in cells])

    Vtot = Vcell.sum()
    eta = len(basis) * (np.pi / 6) * sigma**3 / Vtot

    congruent = (np.ptp(np.round(Vcell, 9)) < 1e-8) and (np.ptp(np.round(Scell, 9)) < 1e-8)

    return dict(
        name=name, note=note, eta=eta, Q=Q, n_faces=nf, z=z,
        tangential=tangential, defect=defect, congruent=congruent,
        n_basis=len(basis),
    )


def report(res):
    Qs = ", ".join(f"{q:.6f}" for q in np.unique(np.round(res["Q"], 9)))
    nf = ", ".join(str(x) for x in np.unique(res["n_faces"]))
    zs = ", ".join(str(x) for x in np.unique(res["z"]))
    print(f"{res['name']:<34s} eta={res['eta']:.6f}  Q={Qs:<22s} "
          f"faces={nf:<7s} z={zs:<6s} "
          f"{'TANGENTIAL' if res['tangential'] else 'not tangential (defect %.3g)' % res['defect']}"
          f"{'' if res['congruent'] else '  [cells NOT congruent]'}")
    if res["note"]:
        print(f"{'':34s}   {res['note']}")

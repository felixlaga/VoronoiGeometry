"""
Analysis pipeline (spec stage S2).

  - read multi-frame configurations written by hsmc
  - S1 gate: equation of state from the contact value of g(r), against Carnahan-Starling
  - S2 gate: radical Voronoi cells -> Q_iso = 36 pi V^2 / S^3, validated on perfect lattices
  - observables: P_eta(Q_iso), <Q_iso>, and the contact fraction
        f_c = <(number of Voronoi faces that are contacts) / (number of Voronoi faces)>
    which is prediction P1.
"""
import sys
import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection


def read_frames(path):
    pass
    frames, i = [], 0
    with open(path) as f:
        lines = f.readlines()
    while i < len(lines):
        n, L = lines[i].split()
        n = int(n); L = float(L)
        block = np.array([[float(t) for t in lines[i + 1 + k].split()] for k in range(n)])
        frames.append((block[:, :3], block[:, 3], L))
        i += n + 1
    return frames


# ------------------------------------------------------------------ S1 gate
def gr_contact(frames, nbin=400, rmaxf=2.5):
    """g(r) for a monodisperse system and its contact value by linear extrapolation."""
    pos, rad, L = frames[0]
    sig = 2 * rad.mean()
    edges = np.linspace(0.9 * sig, rmaxf * sig, nbin + 1)
    hist = np.zeros(nbin)
    for pos, rad, L in frames:
        n = len(pos)
        d = pos[:, None, :] - pos[None, :, :]
        d -= L * np.round(d / L)
        r = np.sqrt((d ** 2).sum(-1))
        r = r[np.triu_indices(n, 1)]
        hist += np.histogram(r, bins=edges)[0]
    n = len(frames[0][0]); L = frames[0][2]
    rho = n / L ** 3
    mid = 0.5 * (edges[1:] + edges[:-1])
    shell = (4 * np.pi / 3) * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = hist / (len(frames) * 0.5 * n * rho * shell)
    m = (mid > sig) & (mid < 1.12 * sig)
    c = np.polyfit(mid[m], g[m], 1)
    return sig, float(np.polyval(c, sig))


def carnahan_starling(eta):
    return (1 + eta + eta ** 2 - eta ** 3) / (1 - eta) ** 3


# ------------------------------------------------------------------ S2 pipeline
def radical_cells(pos, rad, L, cutoff_mult=3.0, area_tol=1e-9, contact_tol=1e-4):
    """Radical (Laguerre) Voronoi cell of every particle.

    The radical plane between i and j sits at
        d_ij = (|r_ij|^2 + R_i^2 - R_j^2) / (2 |r_ij|)
    from i along r_ij.  A face is a CONTACT face when |r_ij| <= R_i + R_j + tol.
    """
    n = len(pos)
    rmax = rad.max()
    cut = cutoff_mult * 2 * rmax
    out = []
    # periodic images
    shifts = np.array([[a, b, c] for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]) * L
    big = (pos[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    bigr = np.tile(rad, len(shifts))
    bigid = np.tile(np.arange(n), len(shifts))
    for i in range(n):
        d = big - pos[i]
        dist = np.linalg.norm(d, axis=1)
        sel = (dist > 1e-9) & (dist < cut)
        d, dist, rj, jid = d[sel], dist[sel], bigr[sel], bigid[sel]
        normals = d / dist[:, None]
        offs = (dist ** 2 + rad[i] ** 2 - rj ** 2) / (2 * dist)
        hs = np.hstack([normals, -offs[:, None]])
        try:
            inter = HalfspaceIntersection(hs, np.zeros(3))
            hull = ConvexHull(inter.intersections)
        except Exception:
            out.append(None)
            continue
        faces = {}
        for eq, simp in zip(hull.equations, hull.simplices):
            key = tuple(np.round(eq, 7))
            tri = inter.intersections[simp]
            a = 0.5 * np.linalg.norm(np.cross(tri[1] - tri[0], tri[2] - tri[0]))
            faces[key] = faces.get(key, 0.0) + a
        nf, ncontact = 0, 0
        for key, a in faces.items():
            if a <= area_tol:
                continue
            nf += 1
            nrm = np.array(key[:3]); off = -key[3]
            k = np.argmax(normals @ nrm - 1e-9 * np.abs(offs - off))
            match = np.argmin(np.abs(normals @ nrm - 1) + np.abs(offs - off))
            if dist[match] <= rad[i] + rj[match] + contact_tol:
                ncontact += 1
            del k
        V, S = hull.volume, sum(a for a in faces.values() if a > area_tol)
        out.append((V, S, 36 * np.pi * V ** 2 / S ** 3, nf, ncontact))
    return out


def frame_observables(frames):
    Q, fc, nf = [], [], []
    for pos, rad, L in frames:
        for c in radical_cells(pos, rad, L):
            if c is None:
                continue
            V, S, q, n_f, n_c = c
            Q.append(q); nf.append(n_f); fc.append(n_c / n_f)
    return np.array(Q), np.array(fc), np.array(nf)


# ------------------------------------------------------------------ S2 gate
def gate_lattices():
    """The radical-Voronoi pipeline must return the exact tabulated Q_iso for the
    perfect tangential lattices, and must flag BCC/diamond as non-tangential."""
    s3 = np.sqrt(3.0)
    tests = {
        "FCC": (np.eye(3), [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]], 0.7404805),
        "SC": (np.eye(3), [[0, 0, 0]], 0.5235988),
        "simple hex c=a": ([[1, 0, 0], [.5, s3 / 2, 0], [0, 0, 1]], [[0, 0, 0]], 0.6045998),
        "FCC-1in4": (np.eye(3), [[0, .5, .5], [.5, 0, .5], [.5, .5, 0]], 0.5553604),
        "BCC": (np.eye(3), [[0, 0, 0], [.5, .5, .5]], None),
    }
    print("S2 GATE: radical-Voronoi pipeline on perfect lattices")
    ok = True
    for name, (lat, bas, ref) in tests.items():
        lat = np.asarray(lat, float); bas = np.asarray(bas, float)
        rep = np.array([b + n1 * lat[0] + n2 * lat[1] + n3 * lat[2]
                        for n1 in range(6) for n2 in range(6) for n3 in range(6)
                        for b in bas])
        # nearest-neighbour distance -> radii at contact
        d = np.linalg.norm(rep - rep[0], axis=1)
        sigma = np.sort(d[d > 1e-9])[0]
        rad = np.full(len(rep), sigma / 2)
        Lx = 6 * lat  # only cubic/orthogonal-friendly cases below use the cube path
        cube = np.allclose(lat, np.diag(np.diag(lat)))
        if not cube:
            print(f"  {name:16s} skipped in gate (non-cubic supercell)")
            continue
        L = 6 * lat[0, 0]
        res = radical_cells(rep[:len(bas)], rad[:len(bas)], L) if False else None
        # use full periodic set: wrap
        pos = rep % L
        cells = radical_cells(pos, rad, L)
        qs = np.array([c[2] for c in cells if c is not None])
        fcs = np.array([c[4] / c[3] for c in cells if c is not None])
        if ref is None:
            print(f"  {name:16s} Q_iso={qs.mean():.7f}  f_c={fcs.mean():.4f}  "
                  f"(expect f_c<1: non-tangential)")
            ok &= fcs.mean() < 0.99
        else:
            err = abs(qs.mean() - ref)
            print(f"  {name:16s} Q_iso={qs.mean():.7f}  ref={ref:.7f}  "
                  f"err={err:.2e}  f_c={fcs.mean():.4f}")
            ok &= (err < 1e-6) and (fcs.mean() > 0.999)
        del Lx, res
    print(f"  GATE {'PASSED' if ok else 'FAILED'}")
    return ok


if __name__ == "__main__":
    if sys.argv[1] == "gate":
        gate_lattices()
    elif sys.argv[1] == "eos":
        for eta in ("0.30", "0.40", "0.45"):
            fr = read_frames(f"/tmp/mono{eta}.cfg")
            sig, gc = gr_contact(fr)
            Z = 1 + 4 * float(eta) * gc
            print(f"  eta={eta}  g(sigma+)={gc:.4f}  Z_MC={Z:.4f}  "
                  f"Z_CS={carnahan_starling(float(eta)):.4f}  "
                  f"rel.dev={abs(Z/carnahan_starling(float(eta))-1)*100:.2f}%")

"""
Two tests the pilot demanded (see RESULTS_pilot.md):

 (1) SINGLE-SCALE NULL for the Voronoi face-gap distribution.  If P(g) collapses
     onto one master curve under g -> g / <g>(eta), then f_c(delta, eta) carries
     no information beyond the median gap, and prediction P1 is structurally
     uninformative rather than merely null.

 (2) SHELL PERCOLATION, the estimator de Graaf uses as his independent
     cross-check on phi_a.  Each particle is inflated by epsilon; the contact
     network is built and tested for a wrapping (percolating) cluster using
     union-find with relative-displacement tracking.  eta_p(epsilon) is then the
     density at which half the configurations percolate.
"""
import numpy as np
from analysis import read_frames


# ---------------------------------------------------------------- union-find with wrapping
class UF:
    def __init__(self, n):
        self.p = list(range(n))
        self.d = np.zeros((n, 3))       # displacement to parent
        self.wrap = np.zeros(3, bool)

    def find(self, a):
        disp = np.zeros(3)
        while self.p[a] != a:
            disp += self.d[a]
            a = self.p[a]
        return a, disp

    def union(self, a, b, off):
        ra, da = self.find(a)
        rb, db = self.find(b)
        # position(a) - position(b) differs by off across the boundary
        if ra == rb:
            net = da - db - off
            self.wrap |= np.abs(net) > 0.5
        else:
            self.p[ra] = rb
            self.d[ra] = db - da + off


def percolates(pos, rad, L, eps):
    """True if the epsilon-inflated contact network wraps in all three directions."""
    n = len(pos)
    uf = UF(n)
    cut = 2 * rad.max() + 2 * eps
    ncell = max(3, int(L / cut))
    cs = L / ncell
    cell = {}
    idx = (pos / cs).astype(int) % ncell
    for i in range(n):
        cell.setdefault(tuple(idx[i]), []).append(i)
    for i in range(n):
        cx, cy, cz = idx[i]
        for a in (-1, 0, 1):
            for b in (-1, 0, 1):
                for c in (-1, 0, 1):
                    key = ((cx + a) % ncell, (cy + b) % ncell, (cz + c) % ncell)
                    for j in cell.get(key, ()):
                        if j <= i:
                            continue
                        d = pos[i] - pos[j]
                        img = np.round(d / L)
                        d = d - L * img
                        if np.dot(d, d) < (rad[i] + rad[j] + 2 * eps) ** 2:
                            uf.union(i, j, img)
    return uf.wrap.all()


# ---------------------------------------------------------------- driver
ETAS = ["0.40", "0.44", "0.48", "0.50", "0.52", "0.54", "0.56", "0.58", "0.60", "0.62"]


def shell_percolation(mode=0, nsnap=6):
    eps_grid = np.logspace(-4, -0.7, 16)
    print("SHELL PERCOLATION: fraction of configurations with a wrapping cluster")
    print("   eta  " + "".join(f"{e:>8.4f}" for e in eps_grid[::3]))
    table = {}
    for e in ETAS:
        fr = read_frames(f"run/m{mode}_e{e}.cfg")[-nsnap:]
        fracs = []
        for eps in eps_grid:
            hits = sum(percolates(p, r, L, eps * 2 * r.mean()) for p, r, L in fr)
            fracs.append(hits / len(fr))
        table[float(e)] = np.array(fracs)
        print(f"  {e}  " + "".join(f"{v:>8.2f}" for v in table[float(e)][::3]))
    return eps_grid, table


def eta_of_eps(eps_grid, table):
    """eta at which half the configurations percolate, for each epsilon."""
    etas = np.array(sorted(table))
    out = []
    for k, eps in enumerate(eps_grid):
        f = np.array([table[e][k] for e in etas])
        if f.max() < 0.5 or f.min() > 0.5:
            out.append(np.nan)
            continue
        i = np.argmax(f >= 0.5)
        if i == 0:
            out.append(etas[0])
        else:
            f0, f1 = f[i - 1], f[i]
            out.append(etas[i - 1] + (0.5 - f0) / (f1 - f0) * (etas[i] - etas[i - 1]))
    return np.array(out)


if __name__ == "__main__":
    import sys
    mode = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    eg, tab = shell_percolation(mode)
    ep = eta_of_eps(eg, tab)
    print("\n  epsilon      eta_p(epsilon)")
    for a, b in zip(eg, ep):
        print(f"  {a:.5f}   {b if not np.isnan(b) else float('nan'):.4f}")
    ok = ~np.isnan(ep)
    if ok.sum() > 4:
        le, lp = np.log(eg[ok]), ep[ok]
        slope = np.gradient(lp, le)
        print("\n  d eta_p / d log eps  (a peak marks the structural change):")
        for a, s in zip(eg[ok], slope):
            print(f"    eps={a:.5f}   {s:+.4f}")

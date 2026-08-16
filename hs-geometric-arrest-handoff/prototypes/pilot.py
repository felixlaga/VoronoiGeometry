import glob
import numpy as np
from analysis import read_frames, radical_cells

ETAS = ["0.40", "0.44", "0.48", "0.50", "0.52", "0.54", "0.56", "0.58", "0.60", "0.62"]
LABEL = {0: "10% polydisperse", 1: "binary 1:1, R=0.714"}


def diffusion(path):
    """Long-time slope of the MSD in MC sweeps.  Returns D (sigma^2 / sweep) and
    the effective exponent alpha over the last decade."""
    d = np.loadtxt(path)
    s, m = d[:, 0], d[:, 1]
    ok = (s > 0) & (m > 0)
    s, m = s[ok], m[ok]
    tail = s > 0.4 * s.max()
    a, b = np.polyfit(np.log(s[tail]), np.log(m[tail]), 1)
    D = m[tail][-1] / (6 * s[tail][-1])
    return D, a


def run(mode, nsnap=6):
    print(f"\n{'='*94}")
    print(f"PILOT  N=864  {LABEL[mode]}  (MC dynamics; swap moves OFF in production)")
    print(f"{'='*94}")
    print(f"{'eta':>6} {'D [sig^2/sweep]':>17} {'alpha':>7} {'<Q_iso>':>10} "
          f"{'sd(Q)':>8} {'f_c':>8} {'<n_face>':>9}")
    rows = []
    for e in ETAS:
        pre = f"run/m{mode}_e{e}"
        try:
            D, al = diffusion(pre + ".msd")
            fr = read_frames(pre + ".cfg")
        except Exception as ex:
            print(f"{e:>6}   missing ({ex})")
            continue
        fr = fr[-nsnap:]
        Q, FC, NF = [], [], []
        for pos, rad, L in fr:
            for c in radical_cells(pos, rad, L):
                if c is None:
                    continue
                V, S, q, nf, nc = c
                Q.append(q); FC.append(nc / nf); NF.append(nf)
        Q = np.array(Q); FC = np.array(FC); NF = np.array(NF)
        print(f"{e:>6} {D:>17.3e} {al:>7.3f} {Q.mean():>10.5f} {Q.std():>8.5f} "
              f"{FC.mean():>8.4f} {NF.mean():>9.3f}")
        rows.append((float(e), D, al, Q.mean(), Q.std(), FC.mean(), NF.mean(), Q))
    return rows


def arrest_fit(rows):
    """eta_a from D = A (eta_a - eta)^b, fitted on the decaying branch."""
    e = np.array([r[0] for r in rows]); D = np.array([r[1] for r in rows])
    m = (e >= 0.48) & (D > 0)
    from scipy.optimize import curve_fit
    f = lambda x, A, ea, b: A * np.clip(ea - x, 1e-9, None) ** b
    try:
        p, c = curve_fit(f, e[m], D[m], p0=[1e-3, 0.62, 1.5], maxfev=40000)
        err = np.sqrt(np.diag(c))
        return p[1], err[1]
    except Exception as ex:
        return None, str(ex)


if __name__ == "__main__":
    allrows = {}
    for mode in (0, 1):
        allrows[mode] = run(mode)
    print(f"\n{'='*94}")
    print("EXTRAPOLATED ARREST DENSITY   D = A (eta_a - eta)^b")
    print(f"{'='*94}")
    for mode in (0, 1):
        ea, err = arrest_fit(allrows[mode])
        if ea is None:
            print(f"  {LABEL[mode]:>22s}:  fit failed ({err})")
        else:
            print(f"  {LABEL[mode]:>22s}:  eta_a = {ea:.4f} +- {err:.4f}")
    print("\n  admissible tangential ground states:")
    for nm, v in [("FCC-1in4", 0.555360), ("FCC-1in5", 0.592384),
                  ("simple hex c=a", 0.604600), ("FCC-1in7", 0.634698),
                  ("simple cubic", 0.523599), ("FCC-1in13", 0.683520)]:
        print(f"     {nm:>16s}  eta = {v:.6f}")

#!/usr/bin/env python3
"""T6 follow-up simulation queue (design frozen before any of its data).

Groups, in execution order (cheap and decision-critical first):

  sb    near-sideband extension: +-{0.007, 0.0085, 0.010, 0.012} around each
        theory window, N=2016 mode 0, 44 replicas -- powers the primary
        kink test (deg 1 / band 0.012)
  ctrl  two control windows at NULL densities 0.6300 and 0.7250 with the
        IDENTICAL grid geometry (25 fine + 8 near-sideband points) and
        statistics -- empirical false-positive calibration; the redesigned
        procedure is invalid if it fires here
  eq    equilibration defense: 16 replicas at five critical densities with
        eq TRIPLED to 180k -- eps* must agree with the campaign values
  stab  seeded-stability probe: the exact (weakly strained) kagome packing
        and an equal-density hexagonal control, 3 seeds each, restart mode,
        plus a fluid baseline -- does kagome survive, convert, or melt?
  comp  kagome window + sidebands in a SECOND composition (mode 2, binary
        R^-1 = 1.7): a geometric feature must be composition-independent
  fs12  kagome window + sidebands at N=504  (ncell=12), 44 replicas
  fs32  kagome window + sidebands at N=3520 (ncell=32), 32 replicas
        -- finite-size axis; P_wrap crossings and size-dependence of any
        feature

All runs: melt 20k / eq 60k / prod 120k / nsnap 8 (except eq group and
stability probe), resume enabled -- the queue is relaunchable.
"""

from __future__ import annotations

import sys
import time
import zlib
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hsga.engine.driver import (  # noqa: E402
    RunSpec, build_engine, run_sweep, write_manifest,
)

DATA = REPO / "data" / "followup"
T_HC = 0.6045997880780726
T_KAG = 0.6801747615878316
T_MAP = 0.7773425846718076
SB_OFF = [0.007, 0.0085, 0.010, 0.012]
EQ_STD, PROD_STD, NSNAP_STD, MELT_STD = 60_000, 120_000, 8, 20_000


def seed_for(tag: str) -> int:
    return zlib.crc32(f"followup|{tag}".encode()) % (2**31)


def window_grid(center: float) -> list[float]:
    fine = [round(center + k * 5e-4, 6) for k in range(-12, 13)]
    sb = [round(center + s * o, 6) for o in SB_OFF for s in (-1, 1)]
    return sorted(fine + sb)


def spec(group, mode, eta, rep, *, ncell=24, eq=EQ_STD, prod=PROD_STD,
         nsnap=NSNAP_STD, infile=None):
    tag = f"{group}|{mode}|{eta:.6f}|{rep}|{ncell}|{eq}"
    return RunSpec(eta=float(eta), dim=2, ncell=ncell, mode=mode,
                   seed=seed_for(tag), eq=eq, prod=prod, nsnap=nsnap,
                   melt=MELT_STD if infile is None else 0, infile=infile,
                   inframe=0 if infile else None,
                   prefix=str(DATA / group / f"m{mode}" / f"e{eta:.6f}_r{rep}"))


def sb_specs():
    etas = [round(t + s * o, 6) for t in (T_HC, T_KAG, T_MAP)
            for o in SB_OFF for s in (-1, 1)]
    return [spec("sb", 0, e, r) for e in sorted(etas) for r in range(44)]


def ctrl_specs():
    return [spec("ctrl", 0, e, r) for c in (0.63, 0.725)
            for e in window_grid(c) for r in range(44)]


def eq_specs():
    dens = [0.680175, 0.725000, 0.777343, 0.780000, 0.800000]
    return [spec("eq", 0, e, r, eq=180_000) for e in dens for r in range(16)]


def comp_specs():
    return [spec("comp", 2, e, r) for e in window_grid(T_KAG) for r in range(44)]


def fs_specs(ncell, nrep):
    return [spec(f"fs{ncell}", 0, e, r, ncell=ncell)
            for e in window_grid(T_KAG) for r in range(nrep)]


# ---------------------------------------------------------------- stability
def build_stability_cfgs() -> list[tuple[str, str]]:
    """Exact kagome and an equal-density hexagonal control on a 52x60
    near-commensurate triangular torus (strain 7.5e-4) in the square box.

    Kagome = triangular lattice minus a TRIANGULAR superlattice of spacing
    2a.  In generator coordinates (g, h) the vacancy rule is g, h both
    even; in the row-offset embedding used here (row j, column i) that is
    ``j % 2 == 0 and (i - j // 2) % 2 == 0``.  Torus periodicity of the
    rule requires nx even and ny % 4 == 0 (landmine: a naive "i, j both
    even" removes a RECTANGULAR superlattice -- not kagome).  The builder
    verifies coordination exactly and refuses to write a wrong geometry.
    """
    out_dir = DATA / "stab"
    out_dir.mkdir(parents=True, exist_ok=True)
    nx, ny = 52, 60
    assert nx % 2 == 0 and ny % 4 == 0
    L = float(nx)                       # a = 1
    sites = []
    for j in range(ny):
        yy = (j + 0.5) * L / ny
        for i in range(nx):
            xx = ((i + 0.5 * (j % 2) + 0.25) % nx)
            vac = (j % 2 == 0) and ((i - j // 2) % 2 == 0)
            sites.append((xx, yy, vac))
    kag = [(x, y) for x, y, v in sites if not v]
    hexa = [(x, y) for x, y, _ in sites]
    assert len(kag) == 3 * len(hexa) // 4

    P = np.array(kag)
    def neighbour_counts(P, cut):
        n = len(P)
        counts = np.zeros(n, int)
        for dx in (-L, 0.0, L):
            for dy in (-L, 0.0, L):
                D = np.linalg.norm(P[None] + [dx, dy] - P[:, None], axis=-1)
                counts += ((D > 1e-9) & (D < cut)).sum(axis=1)
        return counts
    zc = neighbour_counts(P, 1.2)       # ~a with margin below sqrt(3) a
    if not np.all(zc == 4):
        raise RuntimeError(f"kagome coordination wrong: {np.bincount(zc)}")

    def min_dist(pts):
        Q = np.array(pts)
        best = np.inf
        for dx in (-L, 0.0, L):
            for dy in (-L, 0.0, L):
                D = np.linalg.norm(Q[None] + [dx, dy] - Q[:, None], axis=-1)
                if dx == 0 and dy == 0:
                    np.fill_diagonal(D, np.inf)
                best = min(best, float(D.min()))
        return best

    files = []
    r_k = 0.999 * min_dist(kag) / 2.0
    phi_k = len(kag) * np.pi * r_k**2 / L**2
    r_h = r_k * np.sqrt(len(kag) / len(hexa))   # hexagonal at the SAME phi
    for name, pts, r in (("kagome", kag, r_k), ("hex", hexa, r_h)):
        f = out_dir / f"seed_{name}.cfg"
        with open(f, "w") as fh:
            fh.write(f"{len(pts)} {L:.10f}\n")
            for x, y in pts:
                fh.write(f"{x:.8f} {y:.8f} {r:.8f}\n")
        files.append((name, str(f)))
    print(f"  stability seeds: kagome N={len(kag)} (z=4 verified) "
          f"phi={phi_k:.6f} (target {np.sqrt(3)*np.pi/8:.6f}); "
          f"hex N={len(hexa)} at the same phi")
    return files, phi_k


def stab_specs():
    files, phi_k = build_stability_cfgs()
    specs = []
    for name, f in files:
        for rep in range(3):
            tag = f"stab|{name}|{rep}"
            specs.append(RunSpec(
                eta=round(phi_k, 6), dim=2, ncell=24, mode=0,
                seed=seed_for(tag), eq=0, prod=150_000, nsnap=30, melt=0,
                infile=f, inframe=0,
                prefix=str(DATA / "stab" / f"{name}_r{rep}")))
    specs += [spec("stab", 1, round(phi_k, 6), r, ncell=12, eq=20_000,
                   prod=150_000, nsnap=30) for r in range(3)]
    return specs


def build_dilution_cfgs():
    """Dilution ladder for the stability margin: kagome and hexagonal each
    relative to their OWN contact density, plus a square-lattice control --
    z = 4 like kagome but NOT tangential (unbraced; mechanically unstable
    for disks).  If square collapses while kagome holds at equal gap, the
    rigidity is a property of the tangential geometry, not of z."""
    out_dir = DATA / "stab2"
    out_dir.mkdir(parents=True, exist_ok=True)
    files = []
    base = {}
    for name in ("kagome", "hex"):
        lines = (DATA / "stab" / f"seed_{name}.cfg").read_text().splitlines()
        n, L = lines[0].split()
        pts = [tuple(map(float, ln.split())) for ln in lines[1:]]
        base[name] = (float(L), [(x, y) for x, y, _ in pts])
    # square lattice 54x54, a = 1
    nsq = 54
    base["square"] = (float(nsq), [((i + 0.5), (j + 0.5))
                                   for j in range(nsq) for i in range(nsq)])
    for name, (L, pts) in base.items():
        P = np.array(pts)
        best = np.inf
        for dx in (-L, 0.0, L):
            for dy in (-L, 0.0, L):
                D = np.linalg.norm(P[None] + [dx, dy] - P[:, None], axis=-1)
                if dx == 0 and dy == 0:
                    np.fill_diagonal(D, np.inf)
                best = min(best, float(D.min()))
        for shrink in (0.999, 0.995, 0.99, 0.98, 0.96):
            r = shrink * best / 2.0
            f = out_dir / f"seed_{name}_s{shrink}.cfg"
            with open(f, "w") as fh:
                fh.write(f"{len(pts)} {L:.10f}\n")
                for x, y in pts:
                    fh.write(f"{x:.8f} {y:.8f} {r:.8f}\n")
            phi = len(pts) * np.pi * r * r / (L * L)
            files.append((name, shrink, str(f), phi))
    return files


def stab2_specs():
    specs = []
    for name, shrink, f, phi in build_dilution_cfgs():
        for rep in range(3):
            tag = f"stab2|{name}|{shrink}|{rep}"
            specs.append(RunSpec(
                eta=round(phi, 6), dim=2, ncell=24, mode=0,
                seed=seed_for(tag), eq=0, prod=150_000, nsnap=30, melt=0,
                infile=f, inframe=0,
                prefix=str(DATA / "stab2" / f"{name}_s{shrink}_r{rep}")))
    return specs


def build_ladder_cfgs():
    """Rigidity-ladder seeds: maple-leaf (K=1, z=5, hyperstatic) and
    honeycomb (K=3, z=3, floppy) on near-commensurate triangular tori.

    Prediction from the counting identity: occupied sites keep z = 6 - K
    contacts; 2D rigidity needs z >= 4, so the ladder is floppy (K=3) /
    marginal (K=2, measured) / rigid-with-margin (K=1).  Generator-coord
    rules mapped to the row-offset embedding via g = i - j//2, h = j;
    torus periodicity requires nx = 0 (mod m) and ny = 0 (mod 2m) for a
    modulus-m rule with even-row embedding.  Coordination is verified
    exactly; the builder refuses wrong geometry.  Box strain ~1%.
    """
    out_dir = DATA / "stab3"
    out_dir.mkdir(parents=True, exist_ok=True)
    lattices = {
        # name: (nx, ny, vacancy rule in generator coords, expected z)
        "maple": (49, 56, lambda g, h: (g + 3 * h) % 7 == 0, 5),
        "honey": (24, 28, lambda g, h: (g - h) % 3 == 0, 3),
    }
    files = []
    for name, (nx, ny, vacant, z_want) in lattices.items():
        L = float(nx)
        pts = []
        for j in range(ny):
            yy = (j + 0.5) * L / ny
            for i in range(nx):
                g, h = i - j // 2, j
                if not vacant(g, h):
                    pts.append((((i + 0.5 * (j % 2) + 0.25) % nx), yy))
        P = np.array(pts)
        best = np.inf
        counts = np.zeros(len(P), int)
        for dx in (-L, 0.0, L):
            for dy in (-L, 0.0, L):
                D = np.linalg.norm(P[None] + [dx, dy] - P[:, None], axis=-1)
                if dx == 0 and dy == 0:
                    np.fill_diagonal(D, np.inf)
                counts += (D < 1.25).sum(axis=1)
                best = min(best, float(D.min()))
        if not np.all(counts == z_want):
            raise RuntimeError(f"{name}: coordination {np.bincount(counts)} "
                               f"!= z={z_want}")
        for shrink in (0.999, 0.995, 0.99, 0.98, 0.96):
            r = shrink * best / 2.0
            f = out_dir / f"seed_{name}_s{shrink}.cfg"
            with open(f, "w") as fh:
                fh.write(f"{len(pts)} {L:.10f}\n")
                for x, y in pts:
                    fh.write(f"{x:.8f} {y:.8f} {r:.8f}\n")
            phi = len(pts) * np.pi * r * r / (L * L)
            files.append((name, shrink, str(f), phi))
        print(f"  {name}: N={len(pts)} z={z_want} verified, "
              f"phi(contact)={len(pts)*np.pi*(best/2)**2/L**2:.6f}")
    return files


def stab3_specs():
    specs = []
    for name, shrink, f, phi in build_ladder_cfgs():
        for rep in range(3):
            tag = f"stab3|{name}|{shrink}|{rep}"
            specs.append(RunSpec(
                eta=round(phi, 6), dim=2, ncell=24, mode=0,
                seed=seed_for(tag), eq=0, prod=150_000, nsnap=30, melt=0,
                infile=f, inframe=0,
                prefix=str(DATA / "stab3" / f"{name}_s{shrink}_r{rep}")))
    return specs


def build_3d_ladder_cfgs():
    """3D rigidity-ladder seeds: the four depleted-FCC rungs plus FCC and a
    BCC control, on cubic tori.

    FCC sites are the even-parity integer points; the modular vacancy rules
    (REFERENCE_VALUES / coloring.MODULAR_RULES_3D) act on those integer
    coordinates, so the cubic side must be a multiple of both 2 (FCC
    periodicity) and the rule modulus.  Occupied coordination z = 12 - K is
    verified exactly, as is the per-site vacancy count K.  BCC (z = 8, NOT
    tangential -- no 14-faced tangential cell exists, Sec. III of the
    paper) probes whether 3D rigidity needs tangentiality or only z: it is
    hyperstatic by counting yet shear-soft as a sphere packing.
    """
    out_dir = DATA / "stab3d"
    out_dir.mkdir(parents=True, exist_ok=True)
    rules = {
        "K4": (8, 2, lambda x, y, z: x % 2 == 0 and y % 2 == 0 and z % 2 == 0),
        "K3": (10, 5, lambda x, y, z: (y + 2 * z) % 5 == 0),
        "K2": (14, 7, lambda x, y, z: (x + 2 * y + 3 * z) % 7 == 0),
        "K1": (26, 13, lambda x, y, z: (x + 3 * y + 4 * z) % 13 == 0),
        "fcc": (8, 1, lambda x, y, z: False),
    }
    files = []
    for name, (L, mod, vacant) in rules.items():
        assert L % 2 == 0 and L % mod == 0
        occ, vac = [], set()
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    if (x + y + z) % 2:
                        continue
                    (vac.add if vacant(x, y, z) else
                     lambda q: occ.append(q))((x, y, z))
        P = np.array(occ, float)
        K_want = {"K4": 4, "K3": 3, "K2": 2, "K1": 1, "fcc": 0}[name]
        # coordination on the torus at nearest-neighbour distance sqrt(2)
        d = P[None] - P[:, None]
        d -= L * np.round(d / L)
        D = np.linalg.norm(d, axis=-1)
        np.fill_diagonal(D, np.inf)
        zc = (D < 1.5).sum(axis=1)
        if not np.all(zc == 12 - K_want):
            raise RuntimeError(f"{name}: z={np.bincount(zc)} != {12-K_want}")
        files.append(_write_seed(out_dir, name, P, L, np.sqrt(2.0)))
    # BCC control: two interpenetrating cubic grids, nn distance sqrt(3)/2
    Lb = 7
    pts = [(i, j, k) for i in range(Lb) for j in range(Lb) for k in range(Lb)]
    pts += [(i + .5, j + .5, k + .5) for i, j, k in pts]
    P = np.array(pts, float)
    files.append(_write_seed(out_dir, "bcc", P, Lb, np.sqrt(3.0) / 2))
    return files


def _write_seed(out_dir, name, P, L, contact):
    out = []
    for shrink in (0.999, 0.995, 0.99, 0.98, 0.96):
        r = shrink * contact / 2.0
        f = out_dir / f"seed_{name}_s{shrink}.cfg"
        with open(f, "w") as fh:
            fh.write(f"{len(P)} {float(L):.10f}\n")
            for x, y, z in P:
                fh.write(f"{x % L:.8f} {y % L:.8f} {z % L:.8f} {r:.8f}\n")
        phi = len(P) * 4 / 3 * np.pi * r**3 / L**3
        out.append((name, shrink, str(f), phi))
    print(f"  {name}: N={len(P)} phi(contact)="
          f"{len(P)*4/3*np.pi*(contact/2)**3/float(L)**3:.6f}")
    return out


def stab3d_specs():
    specs = []
    for group in build_3d_ladder_cfgs():
        for name, shrink, f, phi in group:
            for rep in range(2):
                tag = f"stab3d|{name}|{shrink}|{rep}"
                specs.append(RunSpec(
                    eta=round(phi, 6), dim=3, ncell=4, mode=0,
                    seed=seed_for(tag), eq=0, prod=150_000, nsnap=30, melt=0,
                    infile=f, inframe=0,
                    prefix=str(DATA / "stab3d" / f"{name}_s{shrink}_r{rep}")))
    return specs


GROUPS = [("sb", sb_specs), ("ctrl", ctrl_specs), ("eq", eq_specs),
          ("stab", stab_specs), ("comp", comp_specs),
          ("fs12", lambda: fs_specs(12, 44)), ("fs32", lambda: fs_specs(32, 32)),
          ("stab2", stab2_specs), ("stab3", stab3_specs),
          ("stab3d", stab3d_specs)]


def main(argv=None) -> int:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--nproc", type=int, default=9)
    p.add_argument("--groups", nargs="*", default=[g for g, _ in GROUPS])
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args(argv)

    total_ps = 0
    plan = []
    for gname, fn in GROUPS:
        if gname not in a.groups:
            continue
        specs = fn()
        ps = sum(s.N * (s.eq + s.prod + s.melt) for s in specs)
        total_ps += ps
        plan.append((gname, specs))
        print(f"group {gname:5s}: {len(specs):5d} runs  "
              f"{ps:.2e} particle-sweeps")
    print(f"TOTAL: {sum(len(s) for _, s in plan)} runs, {total_ps:.2e} "
          f"particle-sweeps (~{total_ps/8.9e7/3600:.1f} h at the measured "
          f"9-worker rate)")
    if a.dry_run:
        return 0
    for gname, specs in plan:
        t0 = time.time()
        exe = build_engine(name="hsmc" if specs[0].dim == 3 else "hsmc2d")
        print(f"\n=== running group {gname} ({len(specs)} runs, "
              f"dim={specs[0].dim}) ===", flush=True)
        runs = run_sweep(specs, exe, nproc=a.nproc, allow_anneal_failure=True,
                         resume=True)
        n_un = sum(1 for r in runs if r["unreachable"])
        n_res = sum(1 for r in runs if r.get("resumed"))
        write_manifest(DATA / gname / "manifest.json", runs, group=gname)
        print(f"group {gname}: {len(runs)-n_un}/{len(runs)} ok "
              f"({n_res} resumed, {n_un} unreachable) "
              f"in {(time.time()-t0)/60:.1f} min", flush=True)
    print("\nFOLLOWUP QUEUE COMPLETE")
    return 0


if __name__ == "__main__":
    sys.exit(main())

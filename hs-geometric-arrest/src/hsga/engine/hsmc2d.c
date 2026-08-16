/* hsmc2d.c -- NVT hard-disk Monte Carlo in 2D.
 *
 * The engine of the 2D-first campaign (spec Sec. 4): the positive control and
 * the replication of de Graaf's central 2D claim, with the kagome rung at
 * phi = sqrt(3) pi / 8 as the novel falsifiable prediction inside the same
 * density sweep.  Clean rewrite of the prototype, keeping its init / anneal /
 * audit logic; the lattice-geometry bug the prototype fixed is kept fixed:
 * for a near-square box the column count MUST be nx = round(sqrt(3) * ny)
 * (landmine: an equilateral triangular lattice in a square box needs cell
 * aspect by/bx = sqrt(3), i.e. nx/ny = sqrt(3)).
 *
 * Local single-particle displacement moves are a stochastic proxy for
 * Brownian dynamics.  Swap moves are used ONLY during initialisation and
 * equilibration; they destroy the physical dynamics and are switched off for
 * production, and the switch-off is logged (landmine 4).
 *
 * Initialisation: near-triangular lattice at the target phi with all radii
 * equal (overlap-free up to phi ~ 0.90), melt, impose the target radii, and
 * anneal the soft overlap energy E = sum (sigma_ij - r_ij)^2 to EXACTLY zero
 * at decreasing T -- with two additions that were needed above phi ~ 0.78 and
 * are documented in DEBT.md: T is reheated to the current per-particle energy
 * scale whenever the energy stalls, and swap moves act on the soft energy
 * during the anneal.  If E > 0 survives, the run aborts (exit 3) and writes
 * no configuration: a state point at or beyond its jamming density is
 * reported as unreachable, never approximated.
 *
 * Restart mode (--in cfg [--frame k]) skips initialisation entirely and reads
 * a frame of an earlier run: positions and radii.  Used by the
 * isoconfigurational analysis (spec 2.5e), which launches many production
 * runs from one equilibrated configuration with different seeds.
 *
 * Usage:
 *   hsmc2d --eta <phi> --ncell <ntri> --mode <0|1|2> --seed <u>
 *          --eq <k> --prod <k> --nsnap <k> --prefix <path>
 *          [--melt <k>] [--swap 0|1] [--anneal-max <k>] [--in cfg --frame k]
 *
 * modes: 0 binary 1:1 at R^-1 = 1.4 (the reference-study composition)
 *        1 monodisperse
 *        2 binary 1:1 at R^-1 = 1.7
 *
 * writes  <prefix>.cfg   multi-frame: "N L" header, then x y r per particle
 *         <prefix>.msd   sweep  msd
 *         <prefix>.log   key=value diagnostics and audits
 *
 * exit codes: 0 ok, 1 usage, 2 init/input failure, 3 anneal failed,
 *             4 final overlap audit failed, 5 box too small for cell list.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int N;
static double L, *x, *y, *xu, *yu, *rad, rmax;
static int *head, *lst, ncell;
static double cellsize;

/* ---- rng: same generator as hsmc.c ---- */
static unsigned long long S0, S1v;
static unsigned long long rotl(unsigned long long a, int k) { return (a << k) | (a >> (64 - k)); }
static unsigned long long nextr(void) {
    unsigned long long s0 = S0, s1 = S1v, r = s0 + s1;
    s1 ^= s0; S0 = rotl(s0, 55) ^ s1 ^ (s1 << 14); S1v = rotl(s1, 36); return r;
}
static double rnd(void) { return (nextr() >> 11) * (1.0 / 9007199254740992.0); }
static int randint(int n) { int i = (int)(rnd() * n); return i >= n ? n - 1 : i; }
static int cmp_double(const void *a, const void *b) {
    double u = *(const double *)a, v = *(const double *)b;
    return (u > v) - (u < v);
}

/* ---- cell list ---- */
static int build_cells(void) {
    ncell = (int)floor(L / (2.0 * rmax));
    if (ncell < 3) return 0;      /* stencil would miss pairs: refuse */
    cellsize = L / ncell;
    free(head); free(lst);
    head = malloc(sizeof(int) * (size_t)ncell * ncell);
    lst = malloc(sizeof(int) * (size_t)N);
    if (!head || !lst) return 0;
    for (int i = 0; i < ncell * ncell; i++) head[i] = -1;
    for (int i = 0; i < N; i++) {
        int cx = (int)(x[i] / cellsize), cy = (int)(y[i] / cellsize);
        if (cx >= ncell) cx = ncell - 1;
        if (cy >= ncell) cy = ncell - 1;
        int c = cx * ncell + cy;
        lst[i] = head[c]; head[c] = i;
    }
    return 1;
}
static inline int cellof(double a) {
    int c = (int)(a / cellsize);
    if (c >= ncell) c = ncell - 1;
    if (c < 0) c = 0;
    return c;
}
static void cell_remove(int i, double px, double py) {
    int c = cellof(px) * ncell + cellof(py);
    int p = head[c];
    if (p == i) { head[c] = lst[i]; return; }
    while (lst[p] != i) p = lst[p];
    lst[p] = lst[i];
}
static void cell_add(int i) {
    int c = cellof(x[i]) * ncell + cellof(y[i]);
    lst[i] = head[c]; head[c] = i;
}
static inline double pbc(double d) {
    if (d > 0.5 * L) d -= L; else if (d < -0.5 * L) d += L;
    return d;
}

/* overlap test (hard=1) or soft overlap energy (hard=0) */
static double local_energy(int i, double px, double py, double ri, int hard) {
    double e = 0.0;
    int cx = cellof(px), cy = cellof(py);
    for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) {
        int ix = (cx + a + ncell) % ncell, iy = (cy + b + ncell) % ncell;
        for (int j = head[ix * ncell + iy]; j >= 0; j = lst[j]) {
            if (j == i) continue;
            double dx = pbc(px - x[j]), dy = pbc(py - y[j]);
            double d2 = dx * dx + dy * dy, s = ri + rad[j];
            if (d2 < s * s) {
                if (hard) return 1.0;
                double o = s - sqrt(d2);
                e += o * o;
            }
        }
    }
    return e;
}
static double total_energy(void) {
    double e = 0;
    for (int i = 0; i < N; i++) e += local_energy(i, x[i], y[i], rad[i], 0);
    return 0.5 * e;
}

/* ---- arguments ---- */
static const char *argstr(int argc, char **argv, const char *key, const char *def) {
    for (int i = 1; i + 1 < argc; i++) if (!strcmp(argv[i], key)) return argv[i + 1];
    return def;
}
static double argf(int argc, char **argv, const char *key, double def) {
    const char *s = argstr(argc, argv, key, NULL);
    return s ? atof(s) : def;
}
static long argl(int argc, char **argv, const char *key, long def) {
    const char *s = argstr(argc, argv, key, NULL);
    return s ? atol(s) : def;
}

static int read_cfg_frame(const char *path, int want) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int n; double l; int got = -1; int allocated = 0;
    while (fscanf(f, "%d %lf", &n, &l) == 2) {
        if (!allocated) {
            N = n;
            x = malloc(N * sizeof(double)); y = malloc(N * sizeof(double));
            xu = malloc(N * sizeof(double)); yu = malloc(N * sizeof(double));
            rad = malloc(N * sizeof(double));
            allocated = 1;
        } else if (n != N) { fclose(f); return -2; }
        L = l;
        for (int i = 0; i < n; i++)
            if (fscanf(f, "%lf %lf %lf", &x[i], &y[i], &rad[i]) != 3) { fclose(f); return -3; }
        got++;
        if (want >= 0 && got == want) break;
    }
    fclose(f);
    if (got < 0) return -4;
    if (want >= 0 && got != want) return -5;
    for (int i = 0; i < N; i++) {
        x[i] -= L * floor(x[i] / L);
        y[i] -= L * floor(y[i] / L);
    }
    return got;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr,
            "usage: hsmc2d --eta <phi> --ncell <ntri> --mode <0|1|2> --seed <u> "
            "--eq <k> --prod <k> --nsnap <k> --prefix <path> "
            "[--melt <k>] [--swap 0|1] [--anneal-max <k>] [--in cfg --frame k]\n");
        return 1;
    }
    double phi = argf(argc, argv, "--eta", 0.70);
    int nt = (int)argl(argc, argv, "--ncell", 16);
    int mode = (int)argl(argc, argv, "--mode", 0);
    unsigned long long seed = strtoull(argstr(argc, argv, "--seed", "1"), 0, 10);
    long eq = argl(argc, argv, "--eq", 20000);
    long prod = argl(argc, argv, "--prod", 40000);
    long nsnap = argl(argc, argv, "--nsnap", 12);
    long melt = argl(argc, argv, "--melt", 20000);
    int use_swap = (int)argl(argc, argv, "--swap", 1);
    long anneal_max = argl(argc, argv, "--anneal-max", 200000);
    const char *pre = argstr(argc, argv, "--prefix", "run");
    const char *infile = argstr(argc, argv, "--in", NULL);
    int inframe = (int)argl(argc, argv, "--frame", -1);

    S0 = seed * 6364136223846793005ULL + 1442695040888963407ULL;
    S1v = seed ^ 0x9E3779B97F4A7C15ULL;
    for (int i = 0; i < 20; i++) nextr();

    char fn[1024];
    snprintf(fn, sizeof fn, "%s.log", pre);
    FILE *lg = fopen(fn, "w");
    if (!lg) { fprintf(stderr, "cannot write %s\n", fn); return 2; }

    long anneal = 0, reheats = 0;
    double dmax;

    if (infile) {
        /* ---- restart: skip initialisation entirely ---- */
        int got = read_cfg_frame(infile, inframe);
        if (got < 0) {
            fprintf(lg, "FATAL=cannot_read_input code=%d file=%s\n", got, infile);
            fclose(lg); return 2;
        }
        rmax = 0;
        for (int i = 0; i < N; i++) if (rad[i] > rmax) rmax = rad[i];
        if (!build_cells()) {
            fprintf(lg, "FATAL=box_too_small_for_cell_list\n"); fclose(lg); return 5;
        }
        double asum = 0;
        for (int i = 0; i < N; i++) asum += M_PI * rad[i] * rad[i];
        phi = asum / (L * L);
        fprintf(lg, "engine=hsmc2d\nrestart_from=%s\nframe_used=%d\nN=%d\nphi=%.6f\nL=%.8f\n",
                infile, got, N, phi, L);
        fprintf(lg, "seed=%llu\neq=%ld\nprod=%ld\nnsnap=%ld\nswap_equilibration=%d\n",
                seed, eq, prod, nsnap, use_swap);
        for (int i = 0; i < N; i++)
            if (local_energy(i, x[i], y[i], rad[i], 1) > 0.0) {
                fprintf(lg, "FATAL=input_overlaps\n"); fclose(lg); return 2;
            }
        dmax = 0.10 * rmax;
    } else {
        /* ---- fresh initialisation ---- */
        int ny = nt, nx = (int)lround(sqrt(3.0) * nt);   /* the landmine, fixed */
        N = 2 * nx * ny;
        x = malloc(N * sizeof(double)); y = malloc(N * sizeof(double));
        xu = malloc(N * sizeof(double)); yu = malloc(N * sizeof(double));
        rad = malloc(N * sizeof(double));
        head = 0; lst = 0;

        double *rt = malloc(N * sizeof(double)), asum = 0;
        for (int i = 0; i < N; i++) {
            if (mode == 0)      rt[i] = (i % 2) ? 0.5 : 0.5 / 1.4;
            else if (mode == 2) rt[i] = (i % 2) ? 0.5 : 0.5 / 1.7;
            else                rt[i] = 0.5;
            asum += M_PI * rt[i] * rt[i];
        }
        L = sqrt(asum / phi);

        /* centred-rectangular = near-equilateral triangular when nx/ny = sqrt3 */
        double bx = L / nx, by = L / ny;
        int k = 0;
        for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) {
            x[k] = i * bx;         y[k] = j * by;         k++;
            x[k] = (i + 0.5) * bx; y[k] = (j + 0.5) * by; k++;
        }
        double r0 = sqrt(phi * L * L / (N * M_PI));
        for (int i = 0; i < N; i++) rad[i] = r0;
        rmax = r0;
        if (!build_cells()) {
            fprintf(lg, "FATAL=box_too_small_for_cell_list\n"); fclose(lg); return 5;
        }
        double half_diag = 0.5 * sqrt(bx * bx + by * by);
        double nnmin = bx < by ? bx : by;
        if (half_diag < nnmin) nnmin = half_diag;
        fprintf(lg, "engine=hsmc2d\nN=%d\nphi=%.6f\nL=%.8f\nnx=%d\nny=%d\nr0=%.6f\nnn=%.6f\n",
                N, phi, L, nx, ny, r0, nnmin);
        fprintf(lg, "mode=%d\nseed=%llu\neq=%ld\nprod=%ld\nnsnap=%ld\nmelt=%ld\n"
                    "swap_equilibration=%d\nanneal_max=%ld\n",
                mode, seed, eq, prod, nsnap, melt, use_swap, anneal_max);
        if (2 * r0 > nnmin + 1e-12) {
            fprintf(lg, "FATAL=monodisperse_init_overlaps phi=%.4f\n", phi);
            fclose(lg); return 2;
        }

        /* melt the monodisperse crystal */
        dmax = 0.06 * 2 * r0;
        for (long s = 0; s < melt; s++) for (int t = 0; t < N; t++) {
            int i = randint(N);
            double px = x[i] + (2 * rnd() - 1) * dmax, py = y[i] + (2 * rnd() - 1) * dmax;
            px -= L * floor(px / L); py -= L * floor(py / L);
            if (local_energy(i, px, py, rad[i], 1) == 0.0) {
                cell_remove(i, x[i], y[i]); x[i] = px; y[i] = py; cell_add(i);
            }
        }

        /* impose target radii; anneal E to exactly zero (reheat + soft swaps) */
        for (int i = 0; i < N; i++) rad[i] = rt[i];
        rmax = 0;
        for (int i = 0; i < N; i++) if (rad[i] > rmax) rmax = rad[i];
        if (!build_cells()) {
            fprintf(lg, "FATAL=box_too_small_for_cell_list\n"); fclose(lg); return 5;
        }
        double E = total_energy();
        double E_ref = E, T = (E > 0 ? E / N : 0.0);
        long stall = 0;
        while (E > 0 && anneal < anneal_max) {
            for (int t = 0; t < N; t++) {
                int i = randint(N);
                double px = x[i] + (2 * rnd() - 1) * dmax, py = y[i] + (2 * rnd() - 1) * dmax;
                px -= L * floor(px / L); py -= L * floor(py / L);
                double eo = local_energy(i, x[i], y[i], rad[i], 0);
                double en = local_energy(i, px, py, rad[i], 0);
                if (en <= eo || (T > 0 && rnd() < exp(-(en - eo) / T))) {
                    cell_remove(i, x[i], y[i]); x[i] = px; y[i] = py; cell_add(i);
                }
            }
            for (int t = 0; t < N / 5; t++) {
                int i = randint(N), j = randint(N);
                if (i == j) continue;
                double ri = rad[i], rj = rad[j];
                if (fabs(ri - rj) < 1e-15) continue;
                double eo = local_energy(i, x[i], y[i], ri, 0) + local_energy(j, x[j], y[j], rj, 0);
                rad[i] = rj; rad[j] = ri;
                double en = local_energy(i, x[i], y[i], rad[i], 0) + local_energy(j, x[j], y[j], rad[j], 0);
                if (!(en <= eo || (T > 0 && rnd() < exp(-(en - eo) / T)))) { rad[i] = ri; rad[j] = rj; }
            }
            anneal++;
            if (anneal % 50 == 0) {
                T *= 0.85;
                E = total_energy();
                if (E < E_ref * (1.0 - 1e-6)) { E_ref = E; stall = 0; }
                else if (++stall >= 8) { T = (E > 0 ? E / N : 0.0); stall = 0; reheats++; }
            }
        }
        E = total_energy();
        /* the final radii must be the target multiset (swaps permute them) */
        int at_target = 1;
        {
            double *ra = malloc(N * sizeof(double)), *rb = malloc(N * sizeof(double));
            memcpy(ra, rad, N * sizeof(double));
            memcpy(rb, rt, N * sizeof(double));
            qsort(ra, N, sizeof(double), cmp_double);
            qsort(rb, N, sizeof(double), cmp_double);
            for (int i = 0; i < N; i++) if (ra[i] != rb[i]) { at_target = 0; break; }
            free(ra); free(rb);
        }
        fprintf(lg, "anneal_sweeps=%ld\nanneal_reheats=%ld\nanneal_final_energy=%.6e\n"
                    "radii_at_target=%d\n", anneal, reheats, E, at_target);
        if (E > 0 || !at_target) {
            fprintf(lg, "FATAL=anneal_failed phi=%.4f residual_energy=%.6e\n", phi, E);
            fclose(lg); return 3;
        }
        free(rt);
    }

    /* ---- equilibration: local moves (+ swap moves if enabled) ---- */
    long acc = 0, att = 0, swap_acc = 0, swap_att = 0;
    for (long s = 0; s < eq; s++) {
        for (int t = 0; t < N; t++) {
            int i = randint(N);
            double px = x[i] + (2 * rnd() - 1) * dmax, py = y[i] + (2 * rnd() - 1) * dmax;
            px -= L * floor(px / L); py -= L * floor(py / L);
            att++;
            if (local_energy(i, px, py, rad[i], 1) == 0.0) {
                cell_remove(i, x[i], y[i]); x[i] = px; y[i] = py; cell_add(i); acc++;
            }
        }
        if (use_swap) for (int t = 0; t < N / 5; t++) {
            int i = randint(N), j = randint(N);
            if (i == j) continue;
            double ri = rad[i], rj = rad[j];
            if (fabs(ri - rj) < 1e-15) continue;
            swap_att++;
            rad[i] = rj; rad[j] = ri;
            if (local_energy(i, x[i], y[i], rad[i], 1) > 0.0 ||
                local_energy(j, x[j], y[j], rad[j], 1) > 0.0) { rad[i] = ri; rad[j] = rj; }
            else swap_acc++;
        }
        if (s % 200 == 199) {
            double r = (double)acc / att;
            dmax *= (r > 0.3 ? 1.02 : 0.98);
            if (dmax > 0.3 * rmax) dmax = 0.3 * rmax;
            acc = att = 0;
        }
    }
    fprintf(lg, "equilibration_done=1\ndmax=%.8f\nswap_acceptance=%.6f\n",
            dmax, swap_att ? (double)swap_acc / swap_att : 0.0);
    fprintf(lg, "swap_moves_production=0\n");     /* landmine 4: explicit */

    /* ---- production: local moves only; MSD in unfolded coordinates ---- */
    for (int i = 0; i < N; i++) { xu[i] = x[i]; yu[i] = y[i]; }
    double *u0x = malloc(N * sizeof(double)), *u0y = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) { u0x[i] = xu[i]; u0y[i] = yu[i]; }
    char fm[1024];
    snprintf(fm, sizeof fm, "%s.msd", pre);
    FILE *fp = fopen(fm, "w");
    snprintf(fm, sizeof fm, "%s.cfg", pre);
    FILE *fc = fopen(fm, "w");
    if (!fp || !fc) { fprintf(lg, "FATAL=cannot_write_output\n"); fclose(lg); return 2; }

    long every = prod / 300; if (every < 1) every = 1;
    long nwritten = 0;
    acc = att = 0;
    for (long s = 0; s < prod; s++) {
        for (int t = 0; t < N; t++) {
            int i = randint(N);
            double ddx = (2 * rnd() - 1) * dmax, ddy = (2 * rnd() - 1) * dmax;
            double px = x[i] + ddx, py = y[i] + ddy;
            px -= L * floor(px / L); py -= L * floor(py / L);
            att++;
            if (local_energy(i, px, py, rad[i], 1) == 0.0) {
                cell_remove(i, x[i], y[i]); x[i] = px; y[i] = py; cell_add(i);
                xu[i] += ddx; yu[i] += ddy; acc++;
            }
        }
        if (s % every == 0) {
            double m = 0;
            for (int i = 0; i < N; i++) {
                double dx = xu[i] - u0x[i], dy = yu[i] - u0y[i];
                m += dx * dx + dy * dy;
            }
            fprintf(fp, "%ld %.8e\n", s, m / N);
        }
        /* snapshot k (0-based) at sweep ceil((k+1) prod/nsnap) - 1: exactly
         * nsnap frames, evenly spaced, last at the end of production */
        if (nsnap > 0 && nwritten < nsnap &&
            s == ((nwritten + 1) * prod + nsnap - 1) / nsnap - 1) {
            fprintf(fc, "%d %.10f\n", N, L);
            for (int i = 0; i < N; i++)
                fprintf(fc, "%.8f %.8f %.8f\n", x[i], y[i], rad[i]);
            nwritten++;
        }
    }
    fclose(fp); fclose(fc);
    fprintf(lg, "production_acceptance=%.6f\nframes_written=%ld\n",
            att ? (double)acc / att : 0.0, nwritten);

    int bad = 0;
    for (int i = 0; i < N; i++) if (local_energy(i, x[i], y[i], rad[i], 1) > 0.0) bad++;
    fprintf(lg, "final_overlap_audit=%d\nexit=%d\n", bad, bad ? 4 : 0);
    fclose(lg);
    return bad ? 4 : 0;
}

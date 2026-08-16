/* edmd.c -- NVE event-driven molecular dynamics for polydisperse hard spheres.
 *
 * Stage 5 of the spec: the engine for the dynamical confirmation.  The
 * dynamics is Newtonian and exact to machine precision -- free flight between
 * events, elastic hard-sphere collisions at events, no time step and no
 * potential softening.  The colloidal experiments the study targets are
 * Brownian; the location of structural arrest is conventionally taken to be
 * insensitive to the microscopic dynamics, and that assumption is recorded in
 * DEBT.md rather than silently relied on.  Every eta_a derived from this
 * engine must be labelled EDMD/Newtonian.
 *
 * Algorithm: an event calendar.  Events are pair collisions and cell-boundary
 * crossings, kept in a binary min-heap and invalidated lazily with
 * per-particle event counters: an event stores the counters of its particles
 * at scheduling time and is discarded on pop if either has since changed.
 * When a particle's trajectory changes (collision) or its neighbourhood
 * changes (cell crossing), its counter is bumped and its future is re-predicted
 * against the 27-cell neighbourhood.  Cell side >= max collision diameter, so
 * any pair is in adjacent cells before contact and the crossing that makes
 * them adjacent re-predicts the pair: no collision can be missed.
 *
 * Input is a configuration in the hsmc .cfg format ("N L" header, then N
 * lines of "x y z r"; --frame selects one frame, default the last).  Initial
 * velocities are Maxwell-Boltzmann at --kT with the centre-of-mass momentum
 * removed and the kinetic energy rescaled exactly; masses are m = (r/rbar)^3
 * (equal material density) unless --equal-mass 1.
 *
 * Usage:
 *   edmd --in run.cfg --tmax <T> --prefix out
 *        [--frame -1] [--seed 1] [--nsample 200] [--kT 1] [--equal-mass 0]
 *   edmd --selftest
 *
 * Writes  <prefix>.msd   "t msd" from a single origin at t = 0, real time
 *         <prefix>.cfg   final configuration (same format as the input)
 *         <prefix>.log   key=value diagnostics and audits
 *
 * The log carries the collisional-virial compressibility
 *   Z = 1 + sum_c J_c |r_ij|_c / (3 N kT_eff tmax),
 * with kT_eff = 2 KE / (3N - 3) (centre-of-mass momentum is fixed), and the
 * conservation audits: kinetic energy and momentum are conserved exactly by
 * the collision rule, so any drift is pure floating-point roundoff and is
 * reported, not assumed.
 *
 * exit codes: 0 ok, 1 usage, 2 bad/overlapping input, 4 final overlap audit
 *             failed, 5 box too small for the cell list, 9 out of memory.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define T_NONE 1e300
#define INIT_OVERLAP_MAX 1e-5    /* input files carry ~1e-8 truncation error */
#define FINAL_OVERLAP_MAX 1e-7   /* collisions resolve at ~1e-12 precision */

static int N;
static double L, t_end;
static double *rx, *ry, *rz, *ux, *uy, *uz, *u0x, *u0y, *u0z;
static double *vx, *vy, *vz, *rad, *mass, *tp;
static long *cnt;                /* per-particle event counter (invalidation) */
static int *cellx, *celly, *cellz, *cnext, *cprev, *chead;
static int ncell;
static double cs;
static double virial = 0.0;
static long n_coll = 0, n_cross = 0, n_proc = 0, n_stale = 0, n_rebuild = 0;

/* ---- rng: same generator as hsmc.c, deterministic per seed ---- */
static unsigned long long S0, S1v;
static unsigned long long rotl(unsigned long long a, int k) { return (a << k) | (a >> (64 - k)); }
static unsigned long long nextr(void) {
    unsigned long long s0 = S0, s1 = S1v, r = s0 + s1;
    s1 ^= s0; S0 = rotl(s0, 55) ^ s1 ^ (s1 << 14); S1v = rotl(s1, 36); return r;
}
static double rnd(void) { return (nextr() >> 11) * (1.0 / 9007199254740992.0); }
static double gauss(void) {
    double u1 = rnd(), u2 = rnd();
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
}
static void seed_rng(unsigned long long seed) {
    S0 = seed * 6364136223846793005ULL + 1442695040888963407ULL;
    S1v = seed ^ 0x9E3779B97F4A7C15ULL;
    for (int i = 0; i < 20; i++) nextr();
}

static inline double pbc(double d) {
    if (d > 0.5 * L) d -= L; else if (d < -0.5 * L) d += L;
    return d;
}

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) { fprintf(stderr, "edmd: out of memory\n"); exit(9); }
    return p;
}

static void alloc_particles(int n) {
    N = n;
    rx = xmalloc(n * sizeof(double)); ry = xmalloc(n * sizeof(double)); rz = xmalloc(n * sizeof(double));
    ux = xmalloc(n * sizeof(double)); uy = xmalloc(n * sizeof(double)); uz = xmalloc(n * sizeof(double));
    u0x = xmalloc(n * sizeof(double)); u0y = xmalloc(n * sizeof(double)); u0z = xmalloc(n * sizeof(double));
    vx = xmalloc(n * sizeof(double)); vy = xmalloc(n * sizeof(double)); vz = xmalloc(n * sizeof(double));
    rad = xmalloc(n * sizeof(double)); mass = xmalloc(n * sizeof(double)); tp = xmalloc(n * sizeof(double));
    cnt = xmalloc(n * sizeof(long));
    cellx = xmalloc(n * sizeof(int)); celly = xmalloc(n * sizeof(int)); cellz = xmalloc(n * sizeof(int));
    cnext = xmalloc(n * sizeof(int)); cprev = xmalloc(n * sizeof(int));
}

/* ---- event heap with lazy invalidation ---- */
typedef struct { double t; int i, j; long ci, cj; } Ev;   /* j >= 0 collision; j = -1-dir crossing */
static Ev *hp = 0;
static long hn = 0, hcap = 0;

static void hpush(double t, int i, int j, long ci, long cj) {
    if (t > t_end) return;                 /* never needed: beyond the run */
    if (hn == hcap) {
        hcap = hcap ? 2 * hcap : 4096;
        hp = realloc(hp, hcap * sizeof(Ev));
        if (!hp) { fprintf(stderr, "edmd: out of memory (heap)\n"); exit(9); }
    }
    long k = hn++;
    hp[k].t = t; hp[k].i = i; hp[k].j = j; hp[k].ci = ci; hp[k].cj = cj;
    while (k > 0) {
        long p = (k - 1) / 2;
        if (hp[p].t <= hp[k].t) break;
        Ev tmp = hp[p]; hp[p] = hp[k]; hp[k] = tmp; k = p;
    }
}
static Ev hpop(void) {
    Ev top = hp[0];
    hp[0] = hp[--hn];
    long k = 0;
    for (;;) {
        long a = 2 * k + 1, b = 2 * k + 2, m = k;
        if (a < hn && hp[a].t < hp[m].t) m = a;
        if (b < hn && hp[b].t < hp[m].t) m = b;
        if (m == k) break;
        Ev tmp = hp[m]; hp[m] = hp[k]; hp[k] = tmp; k = m;
    }
    return top;
}

/* ---- cell list (doubly linked, O(1) removal) ---- */
static inline int cidx(int cx, int cy, int cz) { return (cx * ncell + cy) * ncell + cz; }
static void cell_insert(int i) {
    int c = cidx(cellx[i], celly[i], cellz[i]);
    cprev[i] = -1; cnext[i] = chead[c];
    if (chead[c] >= 0) cprev[chead[c]] = i;
    chead[c] = i;
}
static void cell_delete(int i) {
    int c = cidx(cellx[i], celly[i], cellz[i]);
    if (cprev[i] >= 0) cnext[cprev[i]] = cnext[i]; else chead[c] = cnext[i];
    if (cnext[i] >= 0) cprev[cnext[i]] = cprev[i];
}

/* ---- exact free flight ---- */
static inline void advance(int i, double t) {
    double dt = t - tp[i];
    if (dt > 0) {
        rx[i] += vx[i] * dt; ry[i] += vy[i] * dt; rz[i] += vz[i] * dt;
        ux[i] += vx[i] * dt; uy[i] += vy[i] * dt; uz[i] += vz[i] * dt;
    }
    tp[i] = t;
}

/* Collision time of i and j, both already advanced to time t.
 * Solve |d + w s| = sigma:  (w.w) s^2 + 2 (d.w) s + d.d - sigma^2 = 0.
 * A collision needs d.w < 0 (approaching) and a real root; the smaller root is
 * the contact time.  A roundoff overlap (d.d < sigma^2) with d.w < 0 collides
 * immediately, which is what resolves file-precision input overlaps at t=0. */
static double pair_time(int i, int j, double t) {
    double dx = pbc(rx[i] - rx[j]), dy = pbc(ry[i] - ry[j]), dz = pbc(rz[i] - rz[j]);
    double wx = vx[i] - vx[j], wy = vy[i] - vy[j], wz = vz[i] - vz[j];
    double b = dx * wx + dy * wy + dz * wz;
    if (b >= 0.0) return T_NONE;
    double a = wx * wx + wy * wy + wz * wz;
    if (a <= 0.0) return T_NONE;
    double sig = rad[i] + rad[j];
    double c2 = dx * dx + dy * dy + dz * dz - sig * sig;
    if (c2 < 0.0) return t;
    double disc = b * b - a * c2;
    if (disc <= 0.0) return T_NONE;
    return t + (-b - sqrt(disc)) / a;
}

static void predict_pairs(int i, double t, int only_higher) {
    for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) for (int c = -1; c <= 1; c++) {
        int cx = (cellx[i] + a + ncell) % ncell;
        int cy = (celly[i] + b + ncell) % ncell;
        int cz = (cellz[i] + c + ncell) % ncell;
        for (int j = chead[cidx(cx, cy, cz)]; j >= 0; j = cnext[j]) {
            if (j == i) continue;
            if (only_higher && j < i) continue;
            advance(j, t);                    /* exact: trajectory unchanged */
            double tc = pair_time(i, j, t);
            if (tc < T_NONE) hpush(tc, i, j, cnt[i], cnt[j]);
        }
    }
}

static void predict_cross(int i, double t) {
    double best = T_NONE;
    int dir = -1;
    double pos[3] = { rx[i], ry[i], rz[i] }, vel[3] = { vx[i], vy[i], vz[i] };
    int cc[3] = { cellx[i], celly[i], cellz[i] };
    for (int ax = 0; ax < 3; ax++) {
        if (vel[ax] > 0.0) {
            double dtc = ((cc[ax] + 1) * cs - pos[ax]) / vel[ax];
            if (dtc < 0) dtc = 0;
            if (dtc < best) { best = dtc; dir = 2 * ax + 1; }
        } else if (vel[ax] < 0.0) {
            double dtc = (cc[ax] * cs - pos[ax]) / vel[ax];
            if (dtc < 0) dtc = 0;
            if (dtc < best) { best = dtc; dir = 2 * ax; }
        }
    }
    if (dir >= 0) hpush(t + best, i, -1 - dir, cnt[i], 0);
}

/* Elastic hard-sphere collision along the line of centres.  Conserves momentum
 * and kinetic energy exactly and reverses the normal relative velocity.
 * Returns the pre-collision normal relative velocity (negative = approaching);
 * accumulates the collisional virial J |r_ij|. */
static double collide_pair(int i, int j) {
    double dx = pbc(rx[i] - rx[j]), dy = pbc(ry[i] - ry[j]), dz = pbc(rz[i] - rz[j]);
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    double nx = dx / dist, ny = dy / dist, nz = dz / dist;
    double bn = (vx[i] - vx[j]) * nx + (vy[i] - vy[j]) * ny + (vz[i] - vz[j]) * nz;
    if (bn < 0.0) {
        double mu = mass[i] * mass[j] / (mass[i] + mass[j]);
        double J = -2.0 * mu * bn;            /* impulse on i, along +n */
        vx[i] += J * nx / mass[i]; vy[i] += J * ny / mass[i]; vz[i] += J * nz / mass[i];
        vx[j] -= J * nx / mass[j]; vy[j] -= J * ny / mass[j]; vz[j] -= J * nz / mass[j];
        virial += J * dist;
    }
    return bn;
}

static void do_collide(int i, int j, double te) {
    advance(i, te);
    advance(j, te);
    if (collide_pair(i, j) < 0.0) n_coll++;
    /* a grazing pop (bn >= 0 from roundoff) still re-predicts, harmlessly */
    cnt[i]++; cnt[j]++;
    predict_pairs(i, te, 0); predict_cross(i, te);
    predict_pairs(j, te, 0); predict_cross(j, te);
}

static void do_cross(int i, double te, int dir) {
    advance(i, te);
    cell_delete(i);
    int ax = dir / 2, up = dir & 1;
    int *cc = ax == 0 ? cellx : ax == 1 ? celly : cellz;
    double *pp = ax == 0 ? rx : ax == 1 ? ry : rz;
    if (up) { cc[i]++; if (cc[i] == ncell) { cc[i] = 0;         pp[i] -= L; } }
    else    { cc[i]--; if (cc[i] < 0)      { cc[i] = ncell - 1; pp[i] += L; } }
    cell_insert(i);
    cnt[i]++;         /* old predictions may miss the new neighbourhood */
    n_cross++;
    predict_pairs(i, te, 0);
    predict_cross(i, te);
}

/* rebuild the calendar from scratch (heap hygiene; never affects trajectories) */
static void rebuild(double t) {
    hn = 0;
    for (int i = 0; i < N; i++) { advance(i, t); cnt[i]++; }
    for (int i = 0; i < N; i++) { predict_cross(i, t); predict_pairs(i, t, 1); }
    n_rebuild++;
}

/* max signed relative overlap (sig - dist)/sig over all near pairs */
static double max_overlap_rel(void) {
    double worst = -T_NONE;
    for (int i = 0; i < N; i++) {
        for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) for (int c = -1; c <= 1; c++) {
            int cx = (cellx[i] + a + ncell) % ncell;
            int cy = (celly[i] + b + ncell) % ncell;
            int cz = (cellz[i] + c + ncell) % ncell;
            for (int j = chead[cidx(cx, cy, cz)]; j >= 0; j = cnext[j]) {
                if (j <= i) continue;
                double dx = pbc(rx[i] - rx[j]), dy = pbc(ry[i] - ry[j]), dz = pbc(rz[i] - rz[j]);
                double sig = rad[i] + rad[j];
                double ov = (sig - sqrt(dx * dx + dy * dy + dz * dz)) / sig;
                if (ov > worst) worst = ov;
            }
        }
    }
    return worst;
}

/* ---- input ---- */
static int read_cfg(const char *path, int want) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int n; double l; int got = -1; int allocated = 0;
    while (fscanf(f, "%d %lf", &n, &l) == 2) {
        if (!allocated) { alloc_particles(n); allocated = 1; }
        else if (n != N) { fclose(f); return -2; }
        L = l;
        for (int i = 0; i < n; i++) {
            if (fscanf(f, "%lf %lf %lf %lf", &rx[i], &ry[i], &rz[i], &rad[i]) != 4) {
                fclose(f); return -3;
            }
        }
        got++;
        if (want >= 0 && got == want) break;
    }
    fclose(f);
    if (got < 0) return -4;
    if (want >= 0 && got != want) return -5;
    for (int i = 0; i < N; i++) {          /* enforce folded positions */
        rx[i] -= L * floor(rx[i] / L);
        ry[i] -= L * floor(ry[i] / L);
        rz[i] -= L * floor(rz[i] / L);
    }
    return got;
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
static int has_flag(int argc, char **argv, const char *key) {
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], key)) return 1;
    return 0;
}

/* ---- selftest: exercises the actual production routines ---- */
static int selftest(void) {
    seed_rng(12345ULL);
    L = 1e6;                                /* wraps never fire */
    t_end = T_NONE;                          /* selftest pushes are unbounded */
    alloc_particles(2);
    int fails = 0;

    /* 1. elastic collision: conservation and normal-velocity reversal */
    int bad = 0;
    for (int k = 0; k < 2000; k++) {
        rad[0] = 0.3 + rnd(); rad[1] = 0.3 + rnd();
        mass[0] = 0.2 + 3.0 * rnd(); mass[1] = 0.2 + 3.0 * rnd();
        double sig = rad[0] + rad[1];
        double nx = gauss(), ny = gauss(), nz = gauss();
        double nn = sqrt(nx * nx + ny * ny + nz * nz);
        nx /= nn; ny /= nn; nz /= nn;
        rx[1] = 500.0; ry[1] = 500.0; rz[1] = 500.0;
        rx[0] = rx[1] + sig * nx; ry[0] = ry[1] + sig * ny; rz[0] = rz[1] + sig * nz;
        for (int i = 0; i < 2; i++) { vx[i] = 2 * rnd() - 1; vy[i] = 2 * rnd() - 1; vz[i] = 2 * rnd() - 1; tp[i] = 0; }
        double bn0 = (vx[0] - vx[1]) * nx + (vy[0] - vy[1]) * ny + (vz[0] - vz[1]) * nz;
        if (fabs(bn0) < 1e-3) continue;
        if (bn0 > 0) {                       /* make it approaching: swap velocities */
            double t;
            t = vx[0]; vx[0] = vx[1]; vx[1] = t;
            t = vy[0]; vy[0] = vy[1]; vy[1] = t;
            t = vz[0]; vz[0] = vz[1]; vz[1] = t;
            bn0 = -bn0;
        }
        double P0[3] = { mass[0] * vx[0] + mass[1] * vx[1],
                         mass[0] * vy[0] + mass[1] * vy[1],
                         mass[0] * vz[0] + mass[1] * vz[1] };
        double KE0 = 0.5 * mass[0] * (vx[0] * vx[0] + vy[0] * vy[0] + vz[0] * vz[0])
                   + 0.5 * mass[1] * (vx[1] * vx[1] + vy[1] * vy[1] + vz[1] * vz[1]);
        collide_pair(0, 1);
        double P1[3] = { mass[0] * vx[0] + mass[1] * vx[1],
                         mass[0] * vy[0] + mass[1] * vy[1],
                         mass[0] * vz[0] + mass[1] * vz[1] };
        double KE1 = 0.5 * mass[0] * (vx[0] * vx[0] + vy[0] * vy[0] + vz[0] * vz[0])
                   + 0.5 * mass[1] * (vx[1] * vx[1] + vy[1] * vy[1] + vz[1] * vz[1]);
        double bn1 = (vx[0] - vx[1]) * nx + (vy[0] - vy[1]) * ny + (vz[0] - vz[1]) * nz;
        double dp = fabs(P1[0] - P0[0]) + fabs(P1[1] - P0[1]) + fabs(P1[2] - P0[2]);
        if (dp > 1e-12 || fabs(KE1 - KE0) / KE0 > 1e-12 || fabs(bn1 + bn0) > 1e-12) bad++;
    }
    printf("selftest collision_conservation: %s (%d bad)\n", bad ? "FAIL" : "PASS", bad);
    fails += bad != 0;

    /* 2. predicted collision time lands exactly on contact */
    bad = 0;
    for (int k = 0; k < 2000; k++) {
        rad[0] = 0.3 + rnd(); rad[1] = 0.3 + rnd();
        double sig = rad[0] + rad[1];
        double nx = gauss(), ny = gauss(), nz = gauss();
        double nn = sqrt(nx * nx + ny * ny + nz * nz);
        nx /= nn; ny /= nn; nz /= nn;
        rx[1] = 500.0; ry[1] = 500.0; rz[1] = 500.0;
        double d0 = sig * (1.2 + rnd());
        rx[0] = rx[1] + d0 * nx; ry[0] = ry[1] + d0 * ny; rz[0] = rz[1] + d0 * nz;
        for (int i = 0; i < 2; i++) { vx[i] = 2 * rnd() - 1; vy[i] = 2 * rnd() - 1; vz[i] = 2 * rnd() - 1; tp[i] = 0; }
        double tc = pair_time(0, 1, 0.0);
        if (tc >= T_NONE) continue;
        double wx = vx[0] - vx[1], wy = vy[0] - vy[1], wz = vz[0] - vz[1];
        double ex = (rx[0] - rx[1]) + wx * tc, ey = (ry[0] - ry[1]) + wy * tc, ez = (rz[0] - rz[1]) + wz * tc;
        double dc = sqrt(ex * ex + ey * ey + ez * ez);
        ex = (rx[0] - rx[1]) + wx * tc * 0.999; ey = (ry[0] - ry[1]) + wy * tc * 0.999; ez = (rz[0] - rz[1]) + wz * tc * 0.999;
        double db = sqrt(ex * ex + ey * ey + ez * ez);
        if (fabs(dc - sig) / sig > 1e-9 || db <= sig) bad++;
    }
    printf("selftest collision_time: %s (%d bad)\n", bad ? "FAIL" : "PASS", bad);
    fails += bad != 0;

    /* 3. heap orders events */
    bad = 0;
    hn = 0;
    for (int k = 0; k < 5000; k++) hpush(rnd() * 1000.0, 0, -1, 0, 0);
    double last = -1.0;
    while (hn > 0) {
        Ev e = hpop();
        if (e.t < last) bad++;
        last = e.t;
    }
    printf("selftest heap_order: %s (%d bad)\n", bad ? "FAIL" : "PASS", bad);
    fails += bad != 0;

    printf("selftest %s\n", fails ? "FAIL" : "PASS");
    return fails ? 1 : 0;
}

int main(int argc, char **argv) {
    if (has_flag(argc, argv, "--selftest")) return selftest();

    const char *infile = argstr(argc, argv, "--in", NULL);
    const char *pre = argstr(argc, argv, "--prefix", NULL);
    double tmax = argf(argc, argv, "--tmax", 0.0);
    if (!infile || !pre || tmax <= 0.0) {
        fprintf(stderr,
            "usage: edmd --in <cfg> --tmax <T> --prefix <path>\n"
            "            [--frame -1] [--seed 1] [--nsample 200] [--kT 1] [--equal-mass 0]\n"
            "       edmd --selftest\n");
        return 1;
    }
    int frame = (int)argl(argc, argv, "--frame", -1);
    unsigned long long seed = strtoull(argstr(argc, argv, "--seed", "1"), 0, 10);
    long nsample = argl(argc, argv, "--nsample", 200);
    double kT = argf(argc, argv, "--kT", 1.0);
    int equal_mass = (int)argl(argc, argv, "--equal-mass", 0);
    if (nsample < 1) nsample = 1;
    seed_rng(seed);
    t_end = tmax * (1.0 + 1e-9);

    int frame_used = read_cfg(infile, frame);
    if (frame_used < 0) {
        fprintf(stderr, "edmd: cannot read frame %d from %s (code %d)\n", frame, infile, frame_used);
        return 2;
    }

    /* masses and velocities */
    double rbar = 0.0;
    for (int i = 0; i < N; i++) rbar += rad[i];
    rbar /= N;
    for (int i = 0; i < N; i++)
        mass[i] = equal_mass ? 1.0 : (rad[i] / rbar) * (rad[i] / rbar) * (rad[i] / rbar);
    for (int i = 0; i < N; i++) {
        double sd = sqrt(kT / mass[i]);
        vx[i] = sd * gauss(); vy[i] = sd * gauss(); vz[i] = sd * gauss();
    }
    double Msum = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0;
    for (int i = 0; i < N; i++) {
        Msum += mass[i];
        Px += mass[i] * vx[i]; Py += mass[i] * vy[i]; Pz += mass[i] * vz[i];
    }
    for (int i = 0; i < N; i++) { vx[i] -= Px / Msum; vy[i] -= Py / Msum; vz[i] -= Pz / Msum; }
    double KE = 0.0;
    for (int i = 0; i < N; i++)
        KE += 0.5 * mass[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    double scale = sqrt(1.5 * N * kT / KE);
    for (int i = 0; i < N; i++) { vx[i] *= scale; vy[i] *= scale; vz[i] *= scale; }
    double KE0 = 1.5 * N * kT;
    double kT_eff = 2.0 * KE0 / (3.0 * N - 3.0);   /* COM momentum is fixed */

    /* cells */
    double rmax = 0.0, vol = 0.0;
    for (int i = 0; i < N; i++) {
        if (rad[i] > rmax) rmax = rad[i];
        vol += (4.0 / 3.0) * M_PI * rad[i] * rad[i] * rad[i];
    }
    double eta = vol / (L * L * L);
    ncell = (int)floor(L / (2.0 * rmax));
    if (ncell < 3) {
        fprintf(stderr, "edmd: box too small for the cell list (L=%.6f rmax=%.6f)\n", L, rmax);
        return 5;
    }
    cs = L / ncell;
    chead = xmalloc((size_t)ncell * ncell * ncell * sizeof(int));
    for (int c = 0; c < ncell * ncell * ncell; c++) chead[c] = -1;
    for (int i = 0; i < N; i++) {
        cellx[i] = (int)(rx[i] / cs); if (cellx[i] >= ncell) cellx[i] = ncell - 1;
        celly[i] = (int)(ry[i] / cs); if (celly[i] >= ncell) celly[i] = ncell - 1;
        cellz[i] = (int)(rz[i] / cs); if (cellz[i] >= ncell) cellz[i] = ncell - 1;
        cell_insert(i);
    }

    double ov0 = max_overlap_rel();
    if (ov0 > INIT_OVERLAP_MAX) {
        fprintf(stderr, "edmd: input overlaps by %.3e relative (tolerance %.0e); "
                        "not a valid hard-sphere configuration\n", ov0, INIT_OVERLAP_MAX);
        return 2;
    }

    for (int i = 0; i < N; i++) {
        ux[i] = rx[i]; uy[i] = ry[i]; uz[i] = rz[i];
        u0x[i] = rx[i]; u0y[i] = ry[i]; u0z[i] = rz[i];
        tp[i] = 0.0; cnt[i] = 0;
    }
    for (int i = 0; i < N; i++) { predict_cross(i, 0.0); predict_pairs(i, 0.0, 1); }

    /* output files */
    char fn[1024];
    snprintf(fn, sizeof fn, "%s.msd", pre);
    FILE *fmsd = fopen(fn, "w");
    snprintf(fn, sizeof fn, "%s.log", pre);
    FILE *lg = fopen(fn, "w");
    if (!fmsd || !lg) { fprintf(stderr, "edmd: cannot write output for prefix %s\n", pre); return 1; }
    fprintf(fmsd, "%.8f %.10e\n", 0.0, 0.0);

    clock_t wall0 = clock();
    double dt_s = tmax / nsample;
    for (long s = 1; s <= nsample; s++) {
        double ts = s * dt_s;
        if (s == nsample) ts = tmax;
        while (hn > 0 && hp[0].t <= ts) {
            Ev e = hpop();
            n_proc++;
            if (e.j >= 0) {
                if (e.ci != cnt[e.i] || e.cj != cnt[e.j]) { n_stale++; continue; }
                do_collide(e.i, e.j, e.t);
            } else {
                if (e.ci != cnt[e.i]) { n_stale++; continue; }
                do_cross(e.i, e.t, -1 - e.j);
            }
            if (hn > 200L * N + 10000) rebuild(e.t);
        }
        double m = 0.0;
        for (int i = 0; i < N; i++) {
            advance(i, ts);
            double dx = ux[i] - u0x[i], dy = uy[i] - u0y[i], dz = uz[i] - u0z[i];
            m += dx * dx + dy * dy + dz * dz;
        }
        fprintf(fmsd, "%.8f %.10e\n", ts, m / N);
    }
    fclose(fmsd);
    double wall = (double)(clock() - wall0) / CLOCKS_PER_SEC;

    /* audits: conservation is exact in the collision rule, so drift is roundoff */
    double KE1 = 0.0, P1x = 0.0, P1y = 0.0, P1z = 0.0;
    for (int i = 0; i < N; i++) {
        KE1 += 0.5 * mass[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
        P1x += mass[i] * vx[i]; P1y += mass[i] * vy[i]; P1z += mass[i] * vz[i];
    }
    double ke_drift = fabs(KE1 - KE0) / KE0;
    double p_drift = sqrt(P1x * P1x + P1y * P1y + P1z * P1z) / (N * sqrt((Msum / N) * kT));
    double ov1 = max_overlap_rel();
    double Z = 1.0 + virial / (3.0 * N * kT_eff * tmax);

    snprintf(fn, sizeof fn, "%s.cfg", pre);
    FILE *fc = fopen(fn, "w");
    if (fc) {
        fprintf(fc, "%d %.10f\n", N, L);
        for (int i = 0; i < N; i++)
            fprintf(fc, "%.12f %.12f %.12f %.12f\n", rx[i], ry[i], rz[i], rad[i]);
        fclose(fc);
    }

    fprintf(lg, "engine=edmd\ninput=%s\nframe_used=%d\n", infile, frame_used);
    fprintf(lg, "N=%d\nL=%.8f\neta=%.6f\nsigma_mean=%.8f\n", N, L, eta, 2.0 * rbar);
    fprintf(lg, "kT_input=%.6f\nkT_effective=%.8f\nmass_mode=%s\n",
            kT, kT_eff, equal_mass ? "equal" : "density");
    fprintf(lg, "seed=%llu\ntmax=%.6f\nnsample=%ld\n", seed, tmax, nsample);
    fprintf(lg, "collisions=%ld\ncrossings=%ld\nevents_processed=%ld\nevents_stale=%ld\n",
            n_coll, n_cross, n_proc, n_stale);
    fprintf(lg, "heap_rebuilds=%ld\n", n_rebuild);
    fprintf(lg, "collision_rate_per_particle=%.6f\n", 2.0 * n_coll / ((double)N * tmax));
    fprintf(lg, "Z_virial=%.8f\n", Z);
    fprintf(lg, "ke_drift_rel=%.3e\nmomentum_drift_rel=%.3e\n", ke_drift, p_drift);
    fprintf(lg, "max_overlap_rel_initial=%.3e\nfinal_overlap_audit_rel=%.3e\n", ov0, ov1);
    fprintf(lg, "wall_seconds=%.3f\n", wall);
    int bad = ov1 > FINAL_OVERLAP_MAX;
    fprintf(lg, "exit=%d\n", bad ? 4 : 0);
    fclose(lg);
    return bad ? 4 : 0;
}

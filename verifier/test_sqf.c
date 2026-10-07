// test_sqf.c -- unit tests for sqf.h (fast certified squarefree test) against ref.h (slow reference)
// and for the deterministic Miller-Rabin primality test against sieves.
//
// usage: test_sqf [nrandom_per_class=2000000] [threads=10]
#include "common.h"
#include "sqf.h"
#include "ref.h"
#include <pthread.h>
#include <stdatomic.h>

static int NTHR = 10;
static long NRAND = 2000000;
static atomic_long g_fail;

static void report_mismatch(const char *what, u64 m, int got, int want) {
    long f = atomic_fetch_add(&g_fail, 1);
    if (f < 50) printf("MISMATCH [%s] m=%llu got=%d want=%d\n", what, (unsigned long long)m, got, want);
}

static u64 rand_prime_in(u64 *s, u64 lo, u64 hi) {  // random (MR-)prime in [lo,hi)
    for (;;) { u64 x = lo + rng_below(s, hi - lo); x |= 1; if (x < hi && sqf_is_prime(x)) return x; }
}

// ------------------------------------------------------------------ 1. primality vs sieves
static void test_prime_small(u32 L) {
    unsigned char *c = xcalloc(L + 1, 1);
    c[0] = c[1] = 1;
    for (u64 i = 2; i * i <= L; i++) if (!c[i]) for (u64 j = i * i; j <= L; j += i) c[j] = 1;
    long bad = 0;
    for (u64 n = 0; n <= L; n++) if (sqf_is_prime(n) != !c[n]) { if (bad++ < 10) printf("PRIME MISMATCH n=%llu\n", (unsigned long long)n); }
    printf("[prime] all n <= %u vs sieve: mismatches=%ld\n", L, bad);
    if (bad) atomic_fetch_add(&g_fail, bad);
    free(c);
}

// segmented sieve check of a window [X, X+W) using base primes up to sqrt(X+W)
static u32 *bp; static u32 nbp; static u64 bpmax;
static void make_base_primes(u64 L) {  // odd-only byte sieve; idx i <-> 2i+1
    u64 H = L / 2 + 1;
    unsigned char *c = xcalloc(H, 1);
    for (u64 p = 3; p * p <= L; p += 2) if (!c[p / 2]) for (u64 j = p * p; j <= L; j += 2 * p) c[j / 2] = 1;
    u64 n = 1; for (u64 i = 1; i < H && 2 * i + 1 <= L; i++) if (!c[i]) n++;
    bp = xmalloc(4 * n); nbp = 0; bp[nbp++] = 2;
    for (u64 i = 1; i < H && 2 * i + 1 <= L; i++) if (!c[i]) bp[nbp++] = (u32)(2 * i + 1);
    bpmax = L; free(c);
}
static long window_check(u64 X, u64 W) {
    unsigned char *c = xcalloc(W, 1);
    u64 hi = X + W - 1; u64 r = isqrt64(hi);
    if (r > bpmax) die("base primes too small");
    for (u32 i = 0; i < nbp && bp[i] <= r; i++) {
        u64 p = bp[i]; u64 st = X + (p - X % p) % p;   /* first multiple of p >= X (overflow-safe) */
        if (st < X) continue;
        if (st < p * p) st = p * p;
        for (u64 j = st; j <= hi && j >= st; j += p) c[j - X] = 1;
    }
    long bad = 0, np = 0;
    for (u64 k = 0; k < W; k++) {
        u64 n = X + k; int sv = (n >= 2) && !c[k]; np += sv;
        if (sqf_is_prime(n) != sv) { if (bad++ < 10) printf("PRIME MISMATCH (window) n=%llu sieve=%d\n", (unsigned long long)n, sv); }
    }
    free(c);
    return bad;
}

// ------------------------------------------------------------------ 2. generators
enum { G_UNIFORM, G_UV1, G_UV1_SMALL, G_P2Q, G_P2, G_PQCLOSE, G_PQR, G_P3, G_PQ_MED, G_P2QR, G_SQ, G_KP2, G_P2Q_ROUGH, G_PQ_ROUGH, NGEN };
static const char *gname[NGEN] = {"uniform64", "u*v+1(<2^32)", "u*v+1(<1e9)", "p^2*q", "p^2", "p*q close", "p*q*r ~2^21",
                                  "p^3,p^3*q", "p*q (p in 2^12..2^22)", "p^2*q*r", "s^2", "k*p^2", "p^2*q rough", "p*q rough balanced"};

static u64 gen(int g, u64 *s) {
    for (;;) {
        u64 m = 0;
        switch (g) {
        case G_UNIFORM: m = splitmix64(s); break;
        case G_UV1: { u64 u = 1 + rng_below(s, 0xFFFFFFFFULL), v = 1 + rng_below(s, 0xFFFFFFFFULL); m = u * v + 1; break; }
        case G_UV1_SMALL: { u64 u = 1 + rng_below(s, 1000000000ULL), v = 1 + rng_below(s, 1000000000ULL); m = u * v + 1; break; }
        case G_P2Q: {  // p in [2^k, 2^(k+1)), k random 1..31; q random integer cofactor
            int k = 1 + (int)rng_below(s, 31); u64 p = rand_prime_in(s, 1ULL << k, 1ULL << (k + 1));
            u64 p2 = p * p; u64 qmax = UINT64_MAX / p2; if (qmax < 1) continue;
            u64 q = 1 + rng_below(s, qmax); m = p2 * q; break; }
        case G_P2: { int k = 8 + (int)rng_below(s, 24); u64 p = rand_prime_in(s, 1ULL << k, (1ULL << (k + 1)) - 1); m = p * p; break; }
        case G_PQCLOSE: { int k = 20 + (int)rng_below(s, 12); u64 p = rand_prime_in(s, 1ULL << k, (1ULL << (k + 1)) - 200000);
            u64 q = p + 2; while (!sqf_is_prime(q)) q += 2; if ((u128)p * q >= ((u128)1 << 64)) continue; m = p * q; break; }
        case G_PQR: { u64 a = rand_prime_in(s, 4096, 1ULL << 21), b = rand_prime_in(s, 4096, 1ULL << 21), c = rand_prime_in(s, 4096, 1ULL << 22);
            u128 t = (u128)a * b * c; if (t >> 64) continue; m = (u64)t; break; }
        case G_P3: { u64 p = rand_prime_in(s, 3, 2642245); u128 t = (u128)p * p * p; if (t >> 64) continue;
            u64 q = 1 + rng_below(s, UINT64_MAX / (u64)t); if (rng_below(s, 2)) q = 1; m = (u64)t * q; break; }
        case G_PQ_MED: { u64 p = rand_prime_in(s, 4096, 1ULL << 22); u64 qmax = UINT64_MAX / p;
            u64 q = rand_prime_in(s, qmax / 2 + 1, qmax); m = p * q; break; }
        case G_P2QR: { u64 p = rand_prime_in(s, 4096, 1ULL << 20); u64 rest = UINT64_MAX / (p * p); if (rest < 16) continue;
            u64 r1 = isqrt64(rest); u64 q = rand_prime_in(s, 2, r1), r = rand_prime_in(s, 2, rest / q);
            u128 t = (u128)p * p * q * r; if (t >> 64) continue; m = (u64)t; break; }
        case G_SQ: { u64 x = 2 + rng_below(s, 0xFFFFFFFEULL); m = x * x; break; }
        case G_KP2: { u64 p = rand_prime_in(s, 50, 1ULL << 31); u64 k = 1 + rng_below(s, 1000); u128 t = (u128)k * p * p; if (t >> 64) continue; m = (u64)t; break; }
        case G_P2Q_ROUGH: {  // p and q both rough (> 4096): exercises rho / fallback trial division
            u64 p = rand_prime_in(s, 4099, 1ULL << 26); u64 qmax = UINT64_MAX / (p * p); if (qmax < 5000) continue;
            u64 q = rand_prime_in(s, 4099, qmax); m = p * p * q; break; }
        case G_PQ_ROUGH: {  // balanced rough semiprime, both factors > 2^21.3 (needs full search to certify)
            u64 p = rand_prime_in(s, 2700000, 1ULL << 32); u64 qmax = UINT64_MAX / p; if (qmax < p) continue;
            u64 q = rand_prime_in(s, p, qmax < (1ULL << 40) ? qmax : (1ULL << 40)); m = p * q; break; }
        }
        if (m == 0) continue;
        return m;
    }
}

typedef struct { int g; long n; u64 seed; long yes, no, unk; double sec; sqf_stats_t st; } job_t;
static void *worker(void *arg) {
    job_t *J = arg; u64 s = J->seed; double t0 = now_sec();
    memset(&sqf_st, 0, sizeof sqf_st);
    for (long i = 0; i < J->n; i++) {
        u64 m = gen(J->g, &s);
        int a = sqf_cert(m);
        int b = ref_sf(m);
        if (a == SF_UNK) { J->unk++; continue; }
        if (a != b) report_mismatch(gname[J->g], m, a, b);
        if (a == SF_YES) J->yes++; else J->no++;
        if (sqf_quick_no(m) && b) report_mismatch("quick_no", m, 0, b);
    }
    J->sec = now_sec() - t0; J->st = sqf_st;
    return NULL;
}

static void run_gen(int g, long n, u64 seed) {
    pthread_t th[64]; job_t J[64];
    for (int t = 0; t < NTHR; t++) { memset(&J[t], 0, sizeof J[t]); J[t].g = g; J[t].n = n / NTHR; J[t].seed = seed * 1000003 + t * 7919 + g; pthread_create(&th[t], NULL, worker, &J[t]); }
    long y = 0, no = 0, u = 0; double sec = 0; sqf_stats_t S; memset(&S, 0, sizeof S);
    for (int t = 0; t < NTHR; t++) {
        pthread_join(th[t], NULL); y += J[t].yes; no += J[t].no; u += J[t].unk; if (J[t].sec > sec) sec = J[t].sec;
        u64 *a = (u64 *)&S, *b = (u64 *)&J[t].st; for (size_t k = 0; k < sizeof S / 8; k++) a[k] += b[k];
    }
    printf("[%-22s] n=%ld yes=%ld no=%ld unk=%ld  (%.1fs)  paths: td_no=%llu early=%llu cofsmall=%llu sq=%llu ltP3=%llu mr=%llu prime=%llu cube=%llu hard=%llu rho=%llu/%llu td2=%llu(found %llu, none %llu)\n",
           gname[g], y + no + u, y, no, u, sec, (unsigned long long)S.td_no, (unsigned long long)S.td_early_yes, (unsigned long long)S.cof_small_yes,
           (unsigned long long)S.sq_no, (unsigned long long)S.lt_p3_yes, (unsigned long long)S.mr_calls, (unsigned long long)S.mr_prime,
           (unsigned long long)S.cube_no, (unsigned long long)S.hard_calls, (unsigned long long)S.rho_found, (unsigned long long)S.rho_calls,
           (unsigned long long)S.td2_calls, (unsigned long long)S.td2_found, (unsigned long long)S.td2_none_yes);
    fflush(stdout);
    if (u) atomic_fetch_add(&g_fail, u);
}

// exhaustive over small m
typedef struct { u64 lo, hi; long bad; } rjob_t;
static void *rworker(void *arg) {
    rjob_t *J = arg;
    for (u64 m = J->lo; m < J->hi; m++) {
        int a = sqf_cert(m), b = ref_sf(m);
        if (a != b) { if (J->bad++ < 5) printf("MISMATCH [exhaustive] m=%llu got=%d want=%d\n", (unsigned long long)m, a, b); }
        if (sqf_quick_no(m) && b) { if (J->bad++ < 5) printf("MISMATCH [quick exhaustive] m=%llu\n", (unsigned long long)m); }
    }
    return NULL;
}
static void exhaustive(u64 lo, u64 hi) {
    pthread_t th[64]; rjob_t J[64]; u64 step = (hi - lo) / NTHR + 1;
    for (int t = 0; t < NTHR; t++) {
        u64 off = (u64)t * step; if (off > hi - lo) off = hi - lo;
        J[t].lo = lo + off; J[t].hi = (hi - J[t].lo > step) ? J[t].lo + step : hi; J[t].bad = 0;
        pthread_create(&th[t], NULL, rworker, &J[t]);
    }
    long bad = 0; for (int t = 0; t < NTHR; t++) { pthread_join(th[t], NULL); bad += J[t].bad; }
    printf("[exhaustive] all m in [%llu,%llu): mismatches=%ld\n", (unsigned long long)lo, (unsigned long long)hi, bad);
    if (bad) atomic_fetch_add(&g_fail, bad);
}

int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    if (argc > 1) NRAND = atol(argv[1]);
    if (argc > 2) NTHR = atoi(argv[2]);
    double T0 = now_sec();
    sqf_init(4096); ref_init();
    printf("sqf: B1=4096 (%u odd primes), table primes up to %u (%u primes)\n", sqf_nB1, sqf_primes[sqf_nprimes - 1], sqf_nprimes);

    // --- primality
    test_prime_small(1u << 26);
    make_base_primes((1ULL << 32) - 1);
    {
        u64 s = 12345; long bad = 0;
        int mags[] = {33, 40, 48, 56, 60, 62};
        for (int i = 0; i < 6; i++) for (int rep = 0; rep < 2; rep++) {
            u64 X = (1ULL << mags[i]) + rng_below(&s, 1ULL << mags[i]);
            if (mags[i] == 62) X = (1ULL << 62) + rng_below(&s, 3ULL << 62);  // anywhere up to 2^64
            if (X > UINT64_MAX - (1 << 20)) X = UINT64_MAX - (1 << 20) - 5;
            bad += window_check(X, 1 << 18);
        }
        bad += window_check(UINT64_MAX - (1 << 18) + 1, (1 << 18) - 1);
        printf("[prime] 13 sieve windows of 2^18 at magnitudes 2^33..2^64: mismatches=%ld\n", bad);
        if (bad) atomic_fetch_add(&g_fail, bad);
        // known strong pseudoprimes (composite!)
        const u64 spsp[] = {2047, 3277, 4033, 4681, 8321, 15841, 29341, 42799, 49141, 52633, 65281, 74665, 80581, 85489, 88357, 90751,
                            1373653, 25326001, 3215031751ULL, 2152302898747ULL, 3474749660383ULL, 341550071728321ULL, 3825123056546413051ULL,
                            4759123141ULL, 1122004669633ULL, 9080191, 561, 1105, 1729, 2465, 2821, 6601, 8911, 41041, 825265, 321197185,
                            5394826801ULL, 232250619601ULL, 9746347772161ULL, 1436697831295441ULL, 60977817398996785ULL,
                            7156857700403137441ULL, 1791562810662585601ULL, 18446744073709551557ULL /* prime: 2^64-59 */};
        int nb = 0;
        for (size_t i = 0; i < sizeof spsp / 8; i++) {
            u64 n = spsp[i]; int isp = sqf_is_prime(n);
            // independent: trial division up to min(sqrt, 2^31) using base primes
            int td = 1; u64 r = isqrt64(n);
            for (u32 j = 0; j < nbp && bp[j] <= r; j++) if (n % bp[j] == 0) { td = 0; break; }
            if (r > bpmax && td) td = -1;  // undecided
            if (td >= 0 && td != isp) { nb++; printf("PRIME MISMATCH special n=%llu mr=%d td=%d\n", (unsigned long long)n, isp, td); }
            if (td < 0 && n == 18446744073709551557ULL && !isp) { nb++; printf("2^64-59 not prime?\n"); }
        }
        printf("[prime] %zu special (pseudo)primes: mismatches=%d\n", sizeof spsp / 8, nb);
        if (nb) atomic_fetch_add(&g_fail, nb);
    }

    // --- squarefree: exhaustive small range
    exhaustive(0, 30000000);
    exhaustive(UINT64_MAX - 2000000, UINT64_MAX);
    exhaustive((1ULL << 32) - 1000000, (1ULL << 32) + 1000000);

    // --- squarefree: random classes, default parameters
    printf("--- random classes, rho budget %llu x %d tries\n", (unsigned long long)sqf_rho_budget, sqf_rho_tries);
    for (int g = 0; g < NGEN; g++) run_gen(g, (g == G_PQ_ROUGH) ? NRAND / 4 : NRAND, 1);
    // --- same with rho disabled (forces deterministic fallback trial division)
    sqf_rho_tries = 0;
    printf("--- random classes, rho DISABLED (fallback trial division only)\n");
    for (int g = 0; g < NGEN; g++) run_gen(g, (g == G_PQ_ROUGH) ? NRAND / 20 : NRAND / 5, 2);
    // --- tiny rho budget (rho often fails -> mixture)
    sqf_rho_tries = 3; sqf_rho_budget = 64;
    printf("--- random classes, rho budget 64 x 3 tries\n");
    for (int g = 0; g < NGEN; g++) run_gen(g, (g == G_PQ_ROUGH) ? NRAND / 20 : NRAND / 5, 3);
    sqf_rho_tries = 2; sqf_rho_budget = 1 << 12;

    // --- different stage-1 bounds
    long tot = atomic_load(&g_fail);
    printf("\nTOTAL FAILURES: %ld   (%.1fs)\n", tot, now_sec() - T0);
    return tot ? 1 : 0;
}

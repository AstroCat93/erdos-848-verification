// brute_cmp.c -- brute-force cross-checks for small N (default 1e5), independent of the squarefree
// certificate code: the graph G_N is built by an EDGE SIEVE (for every prime p <= N and every a in V,
// the partners b = -a^{-1} mod p^2 are adjacent), exactly as in ../clique_check.c.
//
//  (1) for ALL pairs a < b in V_N:  adjacent(a,b)  <=>  sqf_cert(ab+1) == NO  <=>  ref_sf(ab+1) == 0
//      (also every a in V has a loop: a^2+1 not squarefree)
//  (2) replay of colour's move log: after EVERY step n, every vertex placed/moved in that step is
//      non-adjacent (edge sieve) to every other member of its final class, class k contains 25k-18,
//      #classes = t(n), all of V_n placed  -> the colouring is proper for every n <= N
//  (3) (optional, "clique") exact induction check of omega(G_n) <= t(n) for every n <= N via
//      branch-and-bound max-clique on earlier neighbourhoods (algorithm of ../clique_check.c).
// usage: brute_cmp N logfile [clique]
#include "common.h"
#include "sqf.h"
#include "ref.h"
#include <pthread.h>

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) { u64 r = 1 % m; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static int64_t inv_mod(int64_t a, int64_t m) {
    int64_t old_r = ((a % m) + m) % m, r = m, old_s = 1, s = 0;
    while (r) { int64_t q = old_r / r, t; t = old_r - q * r; old_r = r; r = t; t = old_s - q * s; old_s = s; s = t; }
    if (old_r != 1) return -1;
    old_s %= m; if (old_s < 0) old_s += m; return old_s;
}
static u64 sqrtm1(u64 p) { for (u64 c = 2;; c++) if (powmod(c, (p - 1) / 2, p) == p - 1) return powmod(c, (p - 1) / 4, p); }

static int NMAX, nv, W; static int *vid, *verts; static u64 *adj;
#define BIT(a, i) ((a)[(i) >> 6] >> ((i) & 63) & 1)
#define SETB(a, i) ((a)[(i) >> 6] |= 1ULL << ((i) & 63))

// ---------------- (3) clique decision, from ../clique_check.c
static int popc(const u64 *a, int w) { int c = 0; for (int i = 0; i < w; i++) c += __builtin_popcountll(a[i]); return c; }
static long long nodes;
static int expand(u64 *cand, int k, u64 *scratch) {
    nodes++;
    int cnt = popc(cand, W); if (cnt < k) return 0; if (k <= 0) return 1;
    int *ord = malloc(sizeof(int) * cnt), *bnd = malloc(sizeof(int) * cnt);
    u64 *Q = scratch, *U = scratch + W;
    memcpy(U, cand, sizeof(u64) * W);
    int color = 0, m = 0;
    while (m < cnt) {
        color++; memcpy(Q, U, sizeof(u64) * W);
        for (int wi = 0; wi < W; wi++) while (Q[wi]) {
            int b = __builtin_ctzll(Q[wi]); int v = wi * 64 + b;
            Q[wi] &= Q[wi] - 1; U[wi] &= ~(1ULL << b);
            const u64 *av = adj + (size_t)v * W;
            for (int j = wi; j < W; j++) Q[j] &= ~av[j];
            ord[m] = v; bnd[m] = color; m++;
        }
    }
    int found = 0;
    u64 *C = malloc(sizeof(u64) * W); memcpy(C, cand, sizeof(u64) * W);
    u64 *next = malloc(sizeof(u64) * W);
    for (int i = m - 1; i >= 0 && !found; i--) {
        if (bnd[i] < k) break;
        int v = ord[i]; const u64 *av = adj + (size_t)v * W;
        for (int j = 0; j < W; j++) next[j] = C[j] & av[j];
        if (expand(next, k - 1, scratch + 2 * W)) found = 1;
        C[v >> 6] &= ~(1ULL << (v & 63));
    }
    free(C); free(next); free(ord); free(bnd);
    return found;
}

// ---------------- (1) all-pairs comparison
typedef struct { int i0, i1; u64 pairs, edges, bad_fast, bad_ref; } pj_t;
static void *pworker(void *a) {
    pj_t *J = a;
    for (int i = J->i0; i < J->i1; i++) for (int j = i + 1; j < nv; j++) {
        u64 x = verts[i], y = verts[j], m = x * y + 1;
        int e = BIT(adj + (size_t)i * W, j);
        int f = sqf_cert(m), r = ref_sf(m);
        J->pairs++; J->edges += e;
        if (f == SF_UNK || (f == SF_NO) != e) { if (!J->bad_fast++) printf("MISMATCH sqf_cert: %llu*%llu+1 edge=%d cert=%d\n", (unsigned long long)x, (unsigned long long)y, e, f); }
        if ((r == 0) != e) { if (!J->bad_ref++) printf("MISMATCH ref_sf: %llu*%llu+1 edge=%d ref=%d\n", (unsigned long long)x, (unsigned long long)y, e, r); }
    }
    return NULL;
}

int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    NMAX = argc > 1 ? (int)strtod(argv[1], 0) : 100000;
    const char *logp = argc > 2 ? argv[2] : NULL;
    int do_clique = argc > 3 && !strcmp(argv[3], "clique");
    double t0 = now_sec();
    sqf_init(4096); ref_init();
    char *comp = calloc(NMAX + 2, 1); int *pr = malloc(sizeof(int) * (NMAX + 2)); int np = 0;
    for (int i = 2; i <= NMAX; i++) if (!comp[i]) { pr[np++] = i; for (long long j = (long long)i * i; j <= NMAX; j += i) comp[j] = 1; }
    int *isV = calloc(NMAX + 1, sizeof(int));
    for (int k = 0; k < np; k++) {
        u64 p = pr[k]; if (p % 4 != 1) continue; u64 p2 = p * p;
        u64 r = sqrtm1(p), f = (r * r + 1) % p2, q = f / p;
        u64 t = mulmod((p - q % p) % p, (u64)inv_mod((int64_t)(2 * r % p), (int64_t)p), p);
        u64 R = (r + t * p) % p2; if ((mulmod(R, R, p2) + 1) % p2) die("lift fail");
        u64 roots[2] = {R, p2 - R};
        for (int s = 0; s < 2; s++) for (u64 n = roots[s]; n <= (u64)NMAX; n += p2) if (n >= 1) isV[n] = 1;
    }
    vid = malloc(sizeof(int) * (NMAX + 1)); verts = malloc(sizeof(int) * (NMAX + 1)); nv = 0;
    for (int n = 1; n <= NMAX; n++) { vid[n] = -1; if (isV[n]) { vid[n] = nv; verts[nv++] = n; } }
    W = (nv + 63) / 64;
    adj = calloc((size_t)nv * W, sizeof(u64));
    for (int k = 0; k < np; k++) {
        int64_t p = pr[k], p2 = p * p;
        for (int i = 0; i < nv; i++) {
            int64_t a = verts[i]; if (a % p == 0) continue;
            int64_t c = (p2 - inv_mod(a % p2, p2)) % p2; if (c == 0) c = p2;
            for (int64_t b = c; b <= NMAX; b += p2) { int j = vid[b]; if (j >= 0) SETB(adj + (size_t)i * W, j); }
        }
    }
    long loops = 0, asym = 0;
    for (int i = 0; i < nv; i++) { loops += BIT(adj + (size_t)i * W, i); for (int j = 0; j < nv; j++) if (BIT(adj + (size_t)i * W, j) != BIT(adj + (size_t)j * W, i)) asym++; }
    printf("edge sieve: N=%d |V|=%d, loops present for %ld/%d vertices, asymmetric entries=%ld (%.1fs)\n", NMAX, nv, loops, nv, asym, now_sec() - t0);
    for (int i = 0; i < nv; i++) (adj + (size_t)i * W)[i >> 6] &= ~(1ULL << (i & 63));

    // (1) all pairs
    double t1 = now_sec(); int NT = 10; pthread_t th[16]; pj_t J[16];
    // balance: thread t takes rows i with i % NT == t would be better; use interleaved blocks
    for (int t = 0; t < NT; t++) {
        // rows near the start are longer; split by sqrt to balance
        double a = 1.0 - sqrt(1.0 - (double)t / NT), b = 1.0 - sqrt(1.0 - (double)(t + 1) / NT);
        J[t].i0 = (int)(a * nv); J[t].i1 = (t == NT - 1) ? nv : (int)(b * nv); J[t].pairs = J[t].edges = J[t].bad_fast = J[t].bad_ref = 0;
        pthread_create(&th[t], NULL, pworker, &J[t]);
    }
    u64 pairs = 0, edges = 0, bf = 0, brf = 0;
    for (int t = 0; t < NT; t++) { pthread_join(th[t], NULL); pairs += J[t].pairs; edges += J[t].edges; bf += J[t].bad_fast; brf += J[t].bad_ref; }
    printf("[all pairs] %llu pairs a<b in V_N, %llu edges; mismatches edge-sieve vs sqf_cert=%llu, vs ref_sf=%llu (%.1fs)\n",
           (unsigned long long)pairs, (unsigned long long)edges, (unsigned long long)bf, (unsigned long long)brf, now_sec() - t1);
    int ok = (loops == nv && !asym && !bf && !brf);

    // (2) replay
    if (logp) {
        FILE *f = fopen(logp, "rb"); if (!f) die("cannot open %s", logp);
        int *where = calloc(NMAX + 1, sizeof(int));
        int TM = (NMAX + 18) / 25; int **cls = malloc(sizeof(int *) * (TM + 2)); int *cnt = calloc(TM + 2, sizeof(int));
        for (int k = 0; k <= TM + 1; k++) cls[k] = malloc(sizeof(int) * 64);
        int vi = 0, T = 0, step[256], ns = 0; u64 checks = 0, errs = 0, steps = 0; u32 r[2];
        while (vi < nv && fread(r, 4, 2, f) == 2) {
            int v = (int)r[0], k = (int)r[1];
            if (v == verts[vi]) {
                if (v % 25 == 7) { if (k != T + 1 || k != (v + 18) / 25) errs++; T = k; }
                else if (k < 1 || k > T || T != (v + 18) / 25) errs++;
                if (k >= 1 && k <= TM) { cls[k][cnt[k]++] = v; where[v] = k; }
                step[ns++] = v;
                for (int s = 0; s < ns; s++) {
                    int u = step[s], kk = where[u];
                    if (cnt[kk] == 0 || cls[kk][0] != 25 * kk - 18) errs++;
                    for (int j = 0; j < cnt[kk]; j++) { int w = cls[kk][j]; if (w == u) continue; checks++;
                        if (BIT(adj + (size_t)vid[u] * W, vid[w])) { if (!errs++) printf("REPLAY FAIL: %d and %d adjacent but same class after step %d\n", u, w, v); } }
                }
                ns = 0; vi++; steps++;
            } else if (v < verts[vi] && where[v] && v % 25 != 7 && k >= 1 && k <= T) {
                int a = where[v]; for (int j = 0; j < cnt[a]; j++) if (cls[a][j] == v) { cls[a][j] = cls[a][--cnt[a]]; break; }
                cls[k][cnt[k]++] = v; where[v] = k; step[ns++] = v;
            } else { printf("REPLAY FAIL: unexpected record (%d,%d)\n", v, k); errs++; break; }
        }
        fclose(f);
        int rok = (vi == nv && errs == 0 && T == (NMAX + 18) / 25);
        printf("[replay vs edge sieve] steps=%llu (|V|=%d), same-class pairs checked=%llu, errors=%llu => %s\n",
               (unsigned long long)steps, nv, (unsigned long long)checks, (unsigned long long)errs, rok ? "PASS: proper t(n)-colouring for every n <= N" : "FAIL");
        ok = ok && rok;
    }

    // (3) clique induction
    if (do_clique) {
        double t3 = now_sec();
        u64 *cand = calloc(W, sizeof(u64)); u64 *scratch = calloc((size_t)W * 2 * (NMAX / 10 + 10), sizeof(u64));
        long long tot = 0; int fails = 0; int maxt = 0;
        for (int i = 0; i < nv; i++) {
            int n = verts[i]; int t = (n + 18) / 25; maxt = t;
            const u64 *an = adj + (size_t)i * W;
            for (int j = 0; j < W; j++) cand[j] = an[j];
            for (int j = i; j < nv; j++) cand[j >> 6] &= ~(1ULL << (j & 63));
            nodes = 0;
            if (expand(cand, t, scratch)) { fails++; printf("CLIQUE VIOLATION at n=%d\n", n); }
            tot += nodes;
            if (i % 2000 == 0) printf("  clique check progress: n=%d (%.0fs)\n", n, now_sec() - t3);
        }
        printf("[clique induction] for every vertex n <= %d: earlier neighbourhood has no clique of size t(n) -> omega(G_n) <= t(n) for all n <= %d; violations=%d, B&B nodes=%lld (%.1fs)\n",
               NMAX, NMAX, fails, tot, now_sec() - t3);
        ok = ok && !fails;
    }
    printf("BRUTE_CMP %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

// Exact verification of Erdos #848 for all N <= NMAX.
// Claim: if A subset [N] and ab+1 is non-squarefree for all a,b in A (a=b allowed),
// then |A| <= t(N) := #{n <= N : n = 7 mod 25}.
//
// Induction: w(G_n) <= max(w(G_{n-1}), 1 + w(G[N^-(n)])) where N^-(n) = earlier neighbours.
// So it suffices to show, for every vertex n, that G[N^-(n)] has no clique of size t(n).
//
// Vertices: V = {n <= N : n^2+1 not squarefree}. Edge a~b iff ab+1 not squarefree.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) { u64 r = 1 % m; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static int64_t inv_mod(int64_t a, int64_t m) { // a coprime to m
    int64_t g = m, x = 0, x1 = 1, a1 = a % m; if (a1 < 0) a1 += m;
    int64_t b = a1; // extended Euclid on (b, g)
    int64_t old_r = b, r = g, old_s = 1, s = 0;
    while (r) { int64_t q = old_r / r, t;
        t = old_r - q * r; old_r = r; r = t;
        t = old_s - q * s; old_s = s; s = t; }
    (void)x; (void)x1;
    if (old_r != 1) return -1;
    old_s %= m; if (old_s < 0) old_s += m; return old_s;
}

// sqrt of -1 mod p for p = 1 mod 4
static u64 sqrtm1(u64 p) {
    for (u64 c = 2;; c++) { if (powmod(c, (p - 1) / 2, p) == p - 1) return powmod(c, (p - 1) / 4, p); }
}

int NMAX; int *isV; int *vid; int *verts; int nv; int W;
u64 *adj; // nv x W bitset

#define BIT(a, i) ((a)[(i) >> 6] >> ((i) & 63) & 1)
#define SETB(a, i) ((a)[(i) >> 6] |= 1ULL << ((i) & 63))

static int popc(const u64 *a, int w) { int c = 0; for (int i = 0; i < w; i++) c += __builtin_popcountll(a[i]); return c; }

// ---- max clique decision: does G[cand] contain a clique of size >= k ? ----
static long long nodes;
static int *colorOrder, *colorBound; // scratch per depth allocated dynamically

static int expand(u64 *cand, int k, int depth, u64 *scratch) {
    // greedy colouring of cand, produce vertex order with colour bounds
    nodes++;
    int cnt = popc(cand, W); if (cnt < k) return 0; if (k <= 0) return 1;
    int *ord = (int *)malloc(sizeof(int) * cnt), *bnd = (int *)malloc(sizeof(int) * cnt);
    u64 *Q = scratch, *U = scratch + W; // scratch has room for 2W + next level
    memcpy(U, cand, sizeof(u64) * W);
    int color = 0, m = 0;
    while (m < cnt) {
        color++; memcpy(Q, U, sizeof(u64) * W);
        for (int wi = 0; wi < W; wi++) {
            while (Q[wi]) {
                int b = __builtin_ctzll(Q[wi]); int v = wi * 64 + b;
                Q[wi] &= Q[wi] - 1; U[wi] &= ~(1ULL << b);
                const u64 *av = adj + (size_t)v * W;
                for (int j = wi; j < W; j++) Q[j] &= ~av[j];
                ord[m] = v; bnd[m] = color; m++;
            }
        }
    }
    int found = 0;
    u64 *C = (u64 *)malloc(sizeof(u64) * W); memcpy(C, cand, sizeof(u64) * W);
    u64 *next = (u64 *)malloc(sizeof(u64) * W);
    for (int i = m - 1; i >= 0 && !found; i--) {
        if (bnd[i] < k) break; // colour bound: cannot reach k
        int v = ord[i];
        const u64 *av = adj + (size_t)v * W;
        for (int j = 0; j < W; j++) next[j] = C[j] & av[j];
        if (expand(next, k - 1, depth + 1, scratch + 2 * W)) found = 1;
        C[v >> 6] &= ~(1ULL << (v & 63));
    }
    free(C); free(next); free(ord); free(bnd);
    return found;
}

int main(int argc, char **argv) {
    NMAX = argc > 1 ? atoi(argv[1]) : 5000;
    int verbose = argc > 2 ? atoi(argv[2]) : 0;
    // primes up to NMAX (p^2 | ab+1 <= NMAX^2+1 implies p <= NMAX)
    char *comp = calloc(NMAX + 2, 1); int *pr = malloc(sizeof(int) * (NMAX + 2)); int np = 0;
    for (int i = 2; i <= NMAX; i++) { if (!comp[i]) { pr[np++] = i; for (long long j = (long long)i * i; j <= NMAX; j += i) comp[j] = 1; } }
    // vertex set
    isV = calloc(NMAX + 1, sizeof(int));
    for (int k = 0; k < np; k++) {
        u64 p = pr[k]; if (p % 4 != 1) continue; u64 p2 = p * p;
        u64 r = sqrtm1(p); // lift to mod p^2: (r + t p)^2 = -1 mod p^2
        u64 f = (r * r + 1) % p2; // divisible by p
        u64 q = f / p; // need 2 r t + q = 0 mod p
        u64 t = mulmod((p - q % p) % p, (u64)inv_mod((int64_t)(2 * r % p), (int64_t)p), p);
        u64 R = (r + t * p) % p2; if ((mulmod(R, R, p2) + 1) % p2) { fprintf(stderr, "lift fail p=%llu\n", (unsigned long long)p); return 1; }
        u64 roots[2] = {R, p2 - R};
        for (int s = 0; s < 2; s++) for (u64 n = roots[s]; n <= (u64)NMAX; n += p2) if (n >= 1) isV[n] = 1;
    }
    vid = malloc(sizeof(int) * (NMAX + 1)); verts = malloc(sizeof(int) * (NMAX + 1)); nv = 0;
    for (int n = 1; n <= NMAX; n++) { vid[n] = -1; if (isV[n]) { vid[n] = nv; verts[nv++] = n; } }
    W = (nv + 63) / 64;
    fprintf(stderr, "N=%d |V|=%d W=%d\n", NMAX, nv, W);
    adj = calloc((size_t)nv * W, sizeof(u64));
    // edges: for each prime p, each a in V, b = -a^{-1} mod p^2
    for (int k = 0; k < np; k++) {
        int64_t p = pr[k], p2 = p * p;
        for (int i = 0; i < nv; i++) {
            int64_t a = verts[i]; if (a % p == 0) continue;
            int64_t c = (p2 - inv_mod(a % p2, p2)) % p2; if (c == 0) c = p2;
            for (int64_t b = c; b <= NMAX; b += p2) { int j = vid[b]; if (j >= 0) { SETB(adj + (size_t)i * W, j); } }
        }
    }
    // sanity: symmetric, loops present
    for (int i = 0; i < nv; i++) if (!BIT(adj + (size_t)i * W, i)) { fprintf(stderr, "missing loop %d\n", verts[i]); return 1; }
    // remove loops for clique search
    for (int i = 0; i < nv; i++) (adj + (size_t)i * W)[i >> 6] &= ~(1ULL << (i & 63));
    u64 *cand = calloc(W, sizeof(u64));
    u64 *scratch = calloc((size_t)W * 2 * (NMAX / 10 + 10), sizeof(u64));
    long long totnodes = 0; int fails = 0;
    for (int i = 0; i < nv; i++) {
        int n = verts[i]; int t = (n >= 7) ? (n - 7) / 25 + 1 : 0;
        // candidates: earlier neighbours of n
        memset(cand, 0, sizeof(u64) * W);
        const u64 *an = adj + (size_t)i * W;
        for (int j = 0; j < W; j++) cand[j] = an[j];
        // keep only indices < i
        for (int j = i; j < nv; j++) cand[j >> 6] &= ~(1ULL << (j & 63));
        nodes = 0;
        int has = expand(cand, t, 0, scratch);
        totnodes += nodes;
        if (has) { fails++; printf("POSSIBLE VIOLATION at n=%d: earlier-neighbourhood has clique of size %d\n", n, t); }
        if (verbose && (i % 500 == 0)) fprintf(stderr, "n=%d t=%d nodes=%lld tot=%lld\n", n, t, nodes, totnodes);
    }
    printf("NMAX=%d |V|=%d checked all vertices, violations=%d, total B&B nodes=%lld\n", NMAX, nv, fails, totnodes);
    return 0;
}

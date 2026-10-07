// vsieve2.h -- INDEPENDENT construction of V_N = {n <= N : n^2+1 not squarefree} (N < 2^32),
// sharing no code with vsieve.h.  Different algorithm:
//   * primes p = 1 mod 4 up to N from a segmented sieve whose base primes come from trial division;
//   * a square root r of -1 modulo p (NOT p^2) by Tonelli-Shanks;
//   * every candidate n = +-r (mod p), n <= N, is tested DIRECTLY: (n*n+1) % (p*p) == 0  (64-bit exact).
// Correctness: q^2 | n^2+1 forces q | n^2+1, i.e. n = +-r (mod q); q = 1 mod 4, q <= N (see vsieve.h).
#ifndef VSIEVE2_H
#define VSIEVE2_H
#include "common.h"
#include <pthread.h>
#include <stdatomic.h>

static u64 v2_mulmod(u64 a, u64 b, u64 p) { return (u64)((u128)a * b % p); }
static u64 v2_pow(u64 a, u64 e, u64 p) { u64 r = 1; a %= p; while (e) { if (e & 1) r = v2_mulmod(r, a, p); a = v2_mulmod(a, a, p); e >>= 1; } return r; }
// Tonelli-Shanks: x with x^2 = a mod p (a a quadratic residue, p odd prime)
static u64 v2_tonelli(u64 a, u64 p) {
    u64 Q = p - 1; int S = 0; while (!(Q & 1)) { Q >>= 1; S++; }
    u64 z = 2; while (v2_pow(z, (p - 1) / 2, p) != p - 1) z++;
    int M = S; u64 c = v2_pow(z, Q, p), t = v2_pow(a, Q, p), R = v2_pow(a, (Q + 1) / 2, p);
    while (t != 1) {
        int i = 0; u64 tt = t; while (tt != 1) { tt = v2_mulmod(tt, tt, p); i++; if (i == M) die("tonelli: non-residue"); }
        u64 b = c; for (int j = 0; j < M - i - 1; j++) b = v2_mulmod(b, b, p);
        M = i; c = v2_mulmod(b, b, p); t = v2_mulmod(t, c, p); R = v2_mulmod(R, b, p);
    }
    return R;
}

typedef struct { u64 *V; u64 N; u32 *bp; u32 nbp; atomic_ullong next; u64 blk; atomic_ullong nprimes; } v2_ctx_t;
static void *v2_worker(void *arg) {
    v2_ctx_t *C = arg;
    unsigned char *seg = xmalloc(C->blk);
    u64 np = 0;
    for (;;) {
        u64 b = atomic_fetch_add(&C->next, 1);
        u64 lo = b * C->blk; if (lo > C->N) break;
        u64 hi = lo + C->blk - 1; if (hi > C->N) hi = C->N;  // [lo, hi]
        memset(seg, 1, hi - lo + 1);                        // 1 = prime candidate
        for (u32 i = 0; i < C->nbp; i++) {
            u64 p = C->bp[i]; if (p * p > hi) break;
            u64 st = lo / p * p; if (st < lo) st += p; if (st < p * p) st = p * p;
            for (u64 j = st; j <= hi; j += p) seg[j - lo] = 0;
        }
        for (u64 x = lo; x <= hi; x++) {
            if (!seg[x - lo] || x < 2 || x % 4 != 1) continue;
            u64 p = x; np++;
            u64 r = v2_tonelli(p - 1, p);
            if (v2_mulmod(r, r, p) != p - 1) die("v2: sqrt(-1) mod %llu wrong", (unsigned long long)p);
            u64 p2 = p * p, rr[2] = {r, p - r};
            for (int s = 0; s < 2; s++)
                for (u64 n = rr[s]; n <= C->N; n += p)
                    if ((n * n + 1) % p2 == 0) __atomic_fetch_or(&C->V[n >> 6], 1ULL << (n & 63), __ATOMIC_RELAXED);
        }
    }
    atomic_fetch_add(&C->nprimes, np);
    free(seg);
    return NULL;
}
static u64 *vsieve2_build(u64 N, int nthreads, u64 *nprimes_out) {
    if (N >= (1ULL << 32)) die("vsieve2: N too large");
    u64 *V = xcalloc((N >> 6) + 2, 8);
    // base primes up to 65536 by trial division (independent of any sieve code)
    u32 *bp = xmalloc(sizeof(u32) * 7000); u32 nbp = 0;
    for (u32 x = 2; x <= 65536; x++) { int pr = 1; for (u32 i = 0; i < nbp && bp[i] * bp[i] <= x; i++) if (x % bp[i] == 0) { pr = 0; break; } if (pr) bp[nbp++] = x; }
    v2_ctx_t C; C.V = V; C.N = N; C.bp = bp; C.nbp = nbp; atomic_init(&C.next, 0); atomic_init(&C.nprimes, 0); C.blk = 1 << 20;
    pthread_t th[64]; if (nthreads < 1) nthreads = 1; if (nthreads > 64) nthreads = 64;
    for (int t = 0; t < nthreads; t++) pthread_create(&th[t], NULL, v2_worker, &C);
    for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);
    if (nprimes_out) *nprimes_out = atomic_load(&C.nprimes);
    free(bp);
    return V;
}
#endif

// vsieve.h -- main vertex sieve: bitset of V_N = {1 <= n <= N : n^2+1 not squarefree}, N < 2^32.
//
// n in V  <=>  q^2 | n^2+1 for some prime q.  q = 2 is impossible (n^2+1 = 1,2 mod 4), q = 3 mod 4 is
// impossible (-1 is a non-residue), and q^2 <= n^2+1 forces q <= N.  So for every prime q = 1 mod 4,
// q <= N, we compute the square roots +-R of -1 modulo q^2 and mark n = +-R mod q^2.
//   sqrt(-1) mod q:  x = a^((q-1)/4) mod q for a quadratic non-residue a (Jacobi symbol), check x^2 = -1.
//   Hensel lift:     R = x + t q,  t = -((x^2+1)/q) * (2x)^-1 mod q;  check R^2+1 = 0 mod q^2 in 128-bit.
// Parallel over blocks of q (segmented Eratosthenes per block); bits are set with atomic OR.
#ifndef VSIEVE_H
#define VSIEVE_H
#include "common.h"
#include <pthread.h>
#include <stdatomic.h>

static inline int vbit(const u64 *V, u64 n) { return (V[n >> 6] >> (n & 63)) & 1; }

static u64 vs_powmod32(u64 a, u64 e, u64 q) {  // q < 2^32
    u64 r = 1 % q; a %= q;
    while (e) { if (e & 1) r = r * a % q; a = a * a % q; e >>= 1; }
    return r;
}
static int vs_jacobi(u64 a, u64 n) {  // n odd positive
    int s = 1; a %= n;
    while (a) {
        while (!(a & 1)) { a >>= 1; u64 r = n & 7; if (r == 3 || r == 5) s = -s; }
        u64 t = a; a = n; n = t;
        if ((a & 3) == 3 && (n & 3) == 3) s = -s;
        a %= n;
    }
    return n == 1 ? s : 0;
}
static u64 vs_invmod(u64 a, u64 m) {  // a coprime to m, m < 2^62
    int64_t t = 0, nt = 1, r = (int64_t)m, nr = (int64_t)(a % m);
    while (nr) { int64_t q = r / nr, tmp;
        tmp = t - q * nt; t = nt; nt = tmp;
        tmp = r - q * nr; r = nr; nr = tmp; }
    if (r != 1) die("vs_invmod: not invertible");
    if (t < 0) t += (int64_t)m;
    return (u64)t;
}
// square root of -1 modulo q^2 (q prime = 1 mod 4, q < 2^32); returns R in (0, q^2)
static u64 vs_sqrtm1_q2(u64 q) {
    u64 a;
    if ((q & 7) == 5) a = 2;
    else { a = 3; while (vs_jacobi(a, q) != -1) a++; }
    u64 x = vs_powmod32(a, (q - 1) / 4, q);
    if (x * x % q != q - 1) die("sqrt(-1) mod %llu failed", (unsigned long long)q);
    u64 f = x * x + 1;  // x < 2^32 so no overflow
    if (f % q) die("Hensel: f not divisible");
    u64 k = (f / q) % q;
    u64 t = (q - k) % q * vs_invmod((2 * x) % q, q) % q;
    u64 R = x + t * q, q2 = q * q;
    if ((u64)(((u128)R * R + 1) % q2) != 0) die("Hensel lift failed for q=%llu", (unsigned long long)q);
    return R;
}

typedef struct {
    u64 *V; u64 N; u32 *bp; u32 nbp; atomic_ullong next; u64 blk; atomic_ullong nq;
} vs_ctx_t;

static void *vs_worker(void *arg) {
    vs_ctx_t *C = arg;
    u64 blk = C->blk;
    unsigned char *seg = xmalloc(blk / 2 + 1);
    u64 nq = 0;
    for (;;) {
        u64 b = atomic_fetch_add(&C->next, 1);
        u64 lo = b * blk; if (lo > C->N) break;
        u64 hi = lo + blk; if (hi > C->N + 1) hi = C->N + 1;   // q in [lo, hi)
        // odd numbers lo|1 .. hi-1 : index (n - lo0)/2 where lo0 = lo | 1
        u64 lo0 = lo | 1; if (lo0 >= hi) continue;
        u64 cnt = (hi - lo0 + 1) / 2;
        memset(seg, 0, cnt);
        for (u32 i = 1; i < C->nbp; i++) {  // odd base primes
            u64 p = C->bp[i]; if (p * p >= hi) break;
            u64 st = (lo0 + p - 1) / p * p; if (st < p * p) st = p * p; if (!(st & 1)) st += p;
            for (u64 j = st; j < hi; j += 2 * p) seg[(j - lo0) / 2] = 1;
        }
        for (u64 i = 0; i < cnt; i++) {
            if (seg[i]) continue;
            u64 q = lo0 + 2 * i; if (q < 3) continue;
            if ((q & 3) != 1) continue;
            nq++;
            u64 R = vs_sqrtm1_q2(q), q2 = q * q;
            u64 rr[2] = {R, q2 - R};
            for (int s = 0; s < 2; s++)
                for (u64 n = rr[s]; n <= C->N; n += q2)
                    __atomic_fetch_or(&C->V[n >> 6], 1ULL << (n & 63), __ATOMIC_RELAXED);
        }
    }
    atomic_fetch_add(&C->nq, nq);
    free(seg);
    return NULL;
}

// returns bitset with bit n set iff n in V (0 <= n <= N); *nq_out = number of primes q = 1 mod 4 used
static u64 *vsieve_build(u64 N, int nthreads, u64 *nq_out) {
    if (N >= (1ULL << 32)) die("vsieve: N must be < 2^32");
    u64 words = (N >> 6) + 2;
    u64 *V = (u64 *)xcalloc(words, 8);
    // base primes up to sqrt(N)+1 (<= 65536)
    u64 L = 65537;
    unsigned char *c = xcalloc(L + 1, 1); u32 nbp = 0;
    u32 *bp = xmalloc(sizeof(u32) * 7000);
    for (u64 i = 2; i <= L; i++) if (!c[i]) { bp[nbp++] = (u32)i; for (u64 j = i * i; j <= L; j += i) c[j] = 1; }
    free(c);
    vs_ctx_t C; C.V = V; C.N = N; C.bp = bp; C.nbp = nbp; atomic_init(&C.next, 0); atomic_init(&C.nq, 0);
    C.blk = 1ULL << 21;
    if (nthreads < 1) nthreads = 1;
    pthread_t th[64]; if (nthreads > 64) nthreads = 64;
    for (int t = 0; t < nthreads; t++) pthread_create(&th[t], NULL, vs_worker, &C);
    for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);
    if (nq_out) *nq_out = atomic_load(&C.nq);
    free(bp);
    V[0] &= ~1ULL;  // n = 0 is not a vertex (never marked anyway)
    return V;
}

static u64 vcount_upto(const u64 *V, u64 n) {  // #{1 <= x <= n : x in V}
    u64 c = 0, w = n >> 6;
    for (u64 i = 0; i < w; i++) c += __builtin_popcountll(V[i]);
    u64 rem = (n & 63) + 1; u64 mask = rem == 64 ? ~0ULL : ((1ULL << rem) - 1);
    c += __builtin_popcountll(V[w] & mask);
    return c;
}
#endif

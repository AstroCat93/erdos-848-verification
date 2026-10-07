// ref.h -- slow but obviously-correct squarefree test, sharing NO code with sqf.h.
//
// ref_sf(m), m < 2^64:  divide out every prime p (in increasing order) while p^3 <= (current m),
// using plain hardware '%' and '/'.  If some p divides twice -> not squarefree.
// When the loop stops, every prime factor of the remaining m is >= p with p^3 > m, so m has at most
// two prime factors: m in {1, q, q*r, q^2}.  Hence squarefree iff m is not a perfect square > 1.
// The integer square root is computed by plain binary search (no floating point).
// Returns 1 if squarefree, 0 if not.
#ifndef REF_H
#define REF_H
#include "common.h"

#define REF_PMAX 2700000u  /* 2700000^3 > 2^64, so the p^3 > m stop always triggers inside the table */
static u32 *ref_primes; static u32 ref_np;

static void ref_init(void) {
    if (ref_primes) return;
    unsigned char *c = (unsigned char *)xcalloc(REF_PMAX + 1, 1);
    for (u32 i = 2; (u64)i * i <= REF_PMAX; i++)
        if (!c[i]) for (u32 j = i * i; j <= REF_PMAX; j += i) c[j] = 1;
    u32 n = 0;
    for (u32 i = 2; i <= REF_PMAX; i++) if (!c[i]) n++;
    ref_primes = (u32 *)xmalloc(sizeof(u32) * n);
    n = 0;
    for (u32 i = 2; i <= REF_PMAX; i++) if (!c[i]) ref_primes[n++] = i;
    ref_np = n;
    free(c);
}

static u64 ref_isqrt(u64 m) {  // largest s with s*s <= m, by binary search
    u64 lo = 0, hi = 1ULL << 32;  // invariant: lo^2 <= m < hi^2
    while (hi - lo > 1) {
        u64 mid = lo + (hi - lo) / 2;
        if ((u128)mid * mid <= m) lo = mid; else hi = mid;
    }
    return lo;
}

static int ref_sf(u64 m) {
    if (m == 0) return 0;
    u32 i;
    for (i = 0; i < ref_np; i++) {
        u64 p = ref_primes[i];
        if ((u128)p * p * p > m) break;
        if (m % p == 0) {
            m /= p;
            if (m % p == 0) return 0;
        }
    }
    if (i == ref_np) { fprintf(stderr, "ref_sf: prime table exhausted\n"); exit(4); }
    if (m == 1) return 1;
    u64 s = ref_isqrt(m);
    return s * s != m;
}

#endif

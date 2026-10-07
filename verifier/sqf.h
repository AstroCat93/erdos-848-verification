// sqf.h -- rigorous squarefree certification for integers 1 <= m < 2^64 (fast path).
//
// sqf_cert(m) returns
//   SF_YES  only if m has been PROVEN squarefree,
//   SF_NO   only if m has been proven NOT squarefree (p^2 | m exhibited, or a perfect square/cube > 1 divides m),
//   SF_UNK  on any internal inconsistency (never expected).  Callers must treat SF_UNK as "incompatible".
//
// Method (all steps exact integer arithmetic):
//  stage 1: trial division by every prime p <= B1 using Granlund-Montgomery multiply-by-inverse tests
//           (odd p:  p | m  <=>  m*inv(p) mod 2^64 <= floor((2^64-1)/p); then m/p = m*inv(p) mod 2^64).
//           p | m twice  =>  NO.
//  stage 2: the cofactor c has every prime factor >= P (P = least prime > B1).  sqf_rough(c):
//     c == 1 or c < P^2               -> c is 1 or prime          -> YES
//     c a perfect square               ->                           NO
//     c < P^3 (and not a square)      -> Omega(c)<=2, c != p^2     -> YES
//     c prime (deterministic MR)      ->                           YES
//     c a perfect cube                 ->                           NO
//     otherwise (composite, >= P^3): Pollard-Brent rho; split c = f*g:
//         gcd(f,g) > 1 -> NO ; else YES iff sqf_rough(f) and sqf_rough(g)
//     if rho fails within its budget: deterministic trial division of c by all primes in [P, icbrt(c)]
//         (table goes to 2.7e6 > 2^(64/3)); first factor p found -> p^2|c ? NO : sqf_rough(c/p);
//         no factor <= icbrt(c) -> every prime factor > icbrt(c) -> Omega(c) <= 2, c not a square -> YES.
//  Deterministic Miller-Rabin base sets (all bases < n whenever used, so no "a = 0 mod n" subtleties):
//     n < 1,373,653: {2,3};  n < 4,759,123,141: {2,7,61};  n < 1,122,004,669,633: {2,13,23,1662803};
//     otherwise (n < 2^64): {2,325,9375,28178,450775,9780504,1795265022}  (Sinclair 2011).
#ifndef SQF_H
#define SQF_H
#include "common.h"
#include <math.h>

#define SF_NO 0
#define SF_YES 1
#define SF_UNK 2

#define SQF_PMAX 2700000u /* > cbrt(2^64) = 2642245.95..., with room for the next prime */

typedef struct { u64 inv, lim; } sqf_divtab_t;

static u32 *sqf_primes;          // odd primes ascending; sqf_primes[0] = 3
static sqf_divtab_t *sqf_div;    // divisibility data per odd prime
static u64 *sqf_psq_next;        // sqf_psq_next[i] = sqf_primes[i+1]^2
static u32 sqf_nprimes;
static u32 sqf_nB1;              // number of odd primes <= B1 (stage 1)
static u64 sqf_rho_budget = 1u << 12;   // rho iterations per attempt
static int sqf_rho_tries = 2;           // rho attempts before deterministic fallback

#define SQF_NQ 14 /* quick filter: odd primes 3..47 */
static u64 sqf_qinv[SQF_NQ], sqf_qlim[SQF_NQ];

typedef struct {
    u64 calls, td_no, td_early_yes, cof_small_yes, sq_no, lt_p3_yes, mr_calls, mr_prime;
    u64 cube_no, hard_calls, rho_calls, rho_found, rho_iters, td2_calls, td2_found, td2_none_yes, unk;
    u64 quick_calls, quick_no;
} sqf_stats_t;
static _Thread_local sqf_stats_t sqf_st;

static inline u64 inv64(u64 a) { u64 x = a; for (int i = 0; i < 5; i++) x *= 2 - a * x; return x; }

static void sqf_init(u64 B1) {
    if (sqf_primes) return;
    u32 L = SQF_PMAX;
    unsigned char *comp = (unsigned char *)xcalloc(L + 1, 1);
    for (u64 i = 2; i * i <= L; i++) if (!comp[i]) for (u64 j = i * i; j <= L; j += i) comp[j] = 1;
    u32 cnt = 0;
    for (u32 i = 3; i <= L; i += 2) if (!comp[i]) cnt++;
    sqf_primes = (u32 *)xmalloc(sizeof(u32) * cnt);
    sqf_div = (sqf_divtab_t *)xmalloc(sizeof(sqf_divtab_t) * cnt);
    sqf_psq_next = (u64 *)xmalloc(sizeof(u64) * cnt);
    cnt = 0;
    for (u32 i = 3; i <= L; i += 2) if (!comp[i]) sqf_primes[cnt++] = i;
    sqf_nprimes = cnt;
    for (u32 i = 0; i < cnt; i++) {
        u64 p = sqf_primes[i];
        sqf_div[i].inv = inv64(p);
        sqf_div[i].lim = UINT64_MAX / p;
        if (sqf_div[i].inv * p != 1) die("inv64 failed for p=%llu", (unsigned long long)p);
        sqf_psq_next[i] = (i + 1 < cnt) ? (u64)sqf_primes[i + 1] * sqf_primes[i + 1] : UINT64_MAX;
    }
    free(comp);
    sqf_nB1 = 0;
    while (sqf_nB1 < cnt && sqf_primes[sqf_nB1] <= B1) sqf_nB1++;
    if (sqf_nB1 >= cnt - 1) die("B1 too large");
    for (int i = 0; i < SQF_NQ; i++) {
        u64 p = sqf_primes[i], p2 = p * p;
        sqf_qinv[i] = inv64(p2); sqf_qlim[i] = UINT64_MAX / p2;
    }
    if (sqf_primes[SQF_NQ - 1] != 47) die("quick table mismatch");
}

// ---------------- exact integer roots ----------------
static inline u64 isqrt64(u64 n) {
    u64 r = (u64)sqrt((double)n);
    while ((u128)r * r > n) r--;
    while ((u128)(r + 1) * (r + 1) <= n) r++;
    return r;
}
static inline int is_square64(u64 n) { u64 r = isqrt64(n); return (u128)r * r == n; }
static inline u64 icbrt64(u64 n) {
    u64 r = (u64)cbrt((double)n);
    while (r > 0 && (u128)r * r * r > n) r--;
    while ((u128)(r + 1) * (r + 1) * (r + 1) <= n) r++;
    return r;
}
static inline int is_cube64(u64 n) { u64 r = icbrt64(n); return (u128)r * r * r == n; }

// ---------------- Montgomery arithmetic (any odd n < 2^64) ----------------
typedef struct { u64 n, ni, one, r2; } mont_t;  // ni = n^-1 mod 2^64, one = 2^64 mod n, r2 = 2^128 mod n
static inline void mont_init(mont_t *M, u64 n) {
    M->n = n; M->ni = inv64(n);
    M->one = (u64)(0 - n) % n;
    M->r2 = (u64)(((u128)M->one * M->one) % n);
}
// returns a*b*2^-64 mod n, for a,b < n.  t - m*n has zero low word; result = hi(t) - hi(m*n) in (-n, n).
static inline u64 mont_mul(u64 a, u64 b, const mont_t *M) {
    u128 t = (u128)a * b;
    u64 lo = (u64)t, hi = (u64)(t >> 64);
    u64 m = lo * M->ni;
    u64 mh = (u64)(((u128)m * M->n) >> 64);
    return hi >= mh ? hi - mh : hi - mh + M->n;
}
static inline u64 mont_to(u64 a, const mont_t *M) { return mont_mul(a % M->n, M->r2, M); }
static inline u64 mont_pow(u64 b, u64 e, const mont_t *M) {
    u64 r = M->one;
    while (e) { if (e & 1) r = mont_mul(r, b, M); b = mont_mul(b, b, M); e >>= 1; }
    return r;
}
static inline u64 addmod64(u64 a, u64 b, u64 n) { u64 t = a + b; if (t < a || t >= n) t -= n; return t; }

static inline int mr_sprp(u64 a, const mont_t *M, u64 d, int s) {
    u64 n = M->n;
    a %= n; if (a == 0) return 1;
    u64 x = mont_pow(mont_to(a, M), d, M);
    u64 one = M->one, mone = n - M->one;
    if (x == one || x == mone) return 1;
    for (int i = 1; i < s; i++) {
        x = mont_mul(x, x, M);
        if (x == mone) return 1;
        if (x == one) return 0;
    }
    return 0;
}

static int sqf_is_prime(u64 n) {
    static const u32 sp[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return 0;
    for (int i = 0; i < 12; i++) { if (n == sp[i]) return 1; if (n % sp[i] == 0) return 0; }
    if (n < 1681) return 1;  // no prime factor <= 37 and n < 41^2
    mont_t M; mont_init(&M, n);
    u64 d = n - 1; int s = __builtin_ctzll(d); d >>= s;
    static const u64 b1[] = {2, 3}, b2[] = {2, 7, 61}, b3[] = {2, 13, 23, 1662803},
                     b4[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    const u64 *b; int nb;
    if (n < 1373653ULL) { b = b1; nb = 2; }
    else if (n < 4759123141ULL) { b = b2; nb = 3; }
    else if (n < 1122004669633ULL) { b = b3; nb = 4; }
    else { b = b4; nb = 7; }
    for (int i = 0; i < nb; i++) if (!mr_sprp(b[i], &M, d, s)) return 0;
    return 1;
}

static inline u64 gcd64(u64 a, u64 b) {
    if (!a) return b;
    if (!b) return a;
    int sh = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do {
        b >>= __builtin_ctzll(b);
        if (a > b) { u64 t = a; a = b; b = t; }
        b -= a;
    } while (b);
    return a << sh;
}

// Pollard-Brent rho on odd composite n.  Returns a nontrivial factor or 0.
static u64 sqf_rho(u64 n, u64 cadd, u64 budget) {
    mont_t M; mont_init(&M, n);
    u64 c = mont_to(cadd, &M);
    u64 y = mont_to(2, &M), x = y, ys = y, q = M.one, g = 1;
    u64 r = 1, it = 0;
    const u64 mblk = 128;
#define SQF_F(v) addmod64(mont_mul((v), (v), &M), c, n)
    do {
        x = y;
        for (u64 i = 0; i < r; i++) y = SQF_F(y);
        it += r;
        u64 k = 0;
        do {
            ys = y;
            u64 lim = (mblk < r - k) ? mblk : r - k;
            for (u64 i = 0; i < lim; i++) { y = SQF_F(y); q = mont_mul(q, x > y ? x - y : y - x, &M); }
            g = gcd64(q, n);
            k += lim; it += lim;
        } while (k < r && g == 1);
        r <<= 1;
    } while (g == 1 && it < budget);
    if (g == n) {
        for (u64 i = 0; i <= mblk; i++) {
            ys = SQF_F(ys);
            g = gcd64(x > ys ? x - ys : ys - x, n);
            if (g != 1) break;
        }
    }
#undef SQF_F
    sqf_st.rho_iters += it;
    if (g == 1 || g == n) return 0;
    return g;
}

static int sqf_rough(u64 c, u32 pidx);

// c: composite, not a perfect square, c >= P^3, every prime factor >= P = sqf_primes[pidx]
static int sqf_hard(u64 c, u32 pidx) {
    sqf_st.hard_calls++;
    if (is_cube64(c)) { sqf_st.cube_no++; return SF_NO; }
    for (int a = 0; a < sqf_rho_tries; a++) {
        sqf_st.rho_calls++;
        u64 f = sqf_rho(c, 1 + 2 * (u64)a, sqf_rho_budget);
        if (f > 1 && f < c) {
            u64 g = c / f;
            if ((u128)f * g != c) { sqf_st.unk++; return SF_UNK; }
            sqf_st.rho_found++;
            if (gcd64(f, g) != 1) return SF_NO;
            int r1 = sqf_rough(f, pidx);
            if (r1 != SF_YES) return r1;
            return sqf_rough(g, pidx);
        }
    }
    // deterministic fallback: trial division over [P, icbrt(c)]
    sqf_st.td2_calls++;
    u64 L = icbrt64(c);
    for (u32 i = pidx; i < sqf_nprimes && sqf_primes[i] <= L; i++) {
        u64 q = c * sqf_div[i].inv;
        if (q <= sqf_div[i].lim) {
            sqf_st.td2_found++;
            if (q * sqf_div[i].inv <= sqf_div[i].lim) return SF_NO;
            return sqf_rough(q, i + 1);
        }
    }
    if (L >= sqf_primes[sqf_nprimes - 1]) { sqf_st.unk++; return SF_UNK; }  // table too short (impossible)
    sqf_st.td2_none_yes++;
    return SF_YES;
}

// c >= 1, every prime factor of c >= P = sqf_primes[pidx]
static int sqf_rough(u64 c, u32 pidx) {
    if (c == 1) { sqf_st.cof_small_yes++; return SF_YES; }
    if (pidx >= sqf_nprimes) { sqf_st.unk++; return SF_UNK; }
    u128 P = sqf_primes[pidx], P2 = P * P;
    if ((u128)c < P2) { sqf_st.cof_small_yes++; return SF_YES; }
    if (is_square64(c)) { sqf_st.sq_no++; return SF_NO; }
    if ((u128)c < P2 * P) { sqf_st.lt_p3_yes++; return SF_YES; }
    sqf_st.mr_calls++;
    if (sqf_is_prime(c)) { sqf_st.mr_prime++; return SF_YES; }
    return sqf_hard(c, pidx);
}

static int sqf_cert(u64 m) {
    sqf_st.calls++;
    if ((m & 3) == 0) { sqf_st.td_no++; return SF_NO; }  // includes m == 0
    if (!(m & 1)) m >>= 1;
    const sqf_divtab_t *D = sqf_div;
    const u32 nB1 = sqf_nB1;
    for (u32 i = 0; i < nB1; i++) {
        u64 q = m * D[i].inv;
        if (__builtin_expect(q <= D[i].lim, 0)) {
            m = q;
            if (m * D[i].inv <= D[i].lim) { sqf_st.td_no++; return SF_NO; }
            if (m < sqf_psq_next[i]) { sqf_st.td_early_yes++; return SF_YES; }  // m is 1 or prime
        }
    }
    return sqf_rough(m, nB1);
}

// Quick partial test: returns 1 if p^2 | m for some prime p <= 47 (=> definitely NOT squarefree).
// Returns 0 if no such p (says nothing about squarefreeness).
static inline int sqf_quick_no(u64 m) {
    sqf_st.quick_calls++;
    if ((m & 3) == 0) { sqf_st.quick_no++; return 1; }
    for (int i = 0; i < SQF_NQ; i++) if (m * sqf_qinv[i] <= sqf_qlim[i]) { sqf_st.quick_no++; return 1; }
    return 0;
}

#endif

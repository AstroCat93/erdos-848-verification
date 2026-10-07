// test_v.c -- cross-checks of the vertex set V_N = {n <= N : n^2+1 not squarefree}:
//  (1) main sieve (vsieve.h: sqrt(-1) mod q^2 by Hensel lifting) vs independent sieve (vsieve2.h:
//      Tonelli-Shanks mod p + direct test of p^2 | n^2+1) over ALL n <= N;
//  (2) direct reference test ref_sf(n^2+1) for every n <= D;
//  (3) direct reference test for R random n <= N and R random members of V.
// usage: test_v N [D=1e7] [R=200000] [threads=10]
#include "common.h"
#include "ref.h"
#include "vsieve.h"
#include "vsieve2.h"
static u64 *V; static u64 N;
typedef struct { u64 lo, hi, bad, cnt; } dj_t;
static void *dworker(void *a) {
    dj_t *J = a;
    for (u64 n = J->lo; n < J->hi; n++) { int inv = !ref_sf(n * n + 1); J->cnt += inv; if (inv != vbit(V, n)) { if (!J->bad++) printf("MISMATCH direct n=%llu\n", (unsigned long long)n); } }
    return NULL;
}
typedef struct { u64 seed, R, bad, inV; } rj_t;
static void *rworker(void *a) {
    rj_t *J = a; u64 s = J->seed;
    for (u64 i = 0; i < J->R; i++) {
        u64 n = 1 + rng_below(&s, N);
        if (i & 1) { while (!vbit(V, n)) { n++; if (n > N) n = 1; } }  // odd i: a random member of V
        int inv = !ref_sf(n * n + 1); J->inV += inv;
        if (inv != vbit(V, n)) { if (!J->bad++) printf("MISMATCH random n=%llu\n", (unsigned long long)n); }
    }
    return NULL;
}
int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    N = argc > 1 ? (u64)strtod(argv[1], 0) : 100000000ULL;
    u64 D = argc > 2 ? (u64)strtod(argv[2], 0) : 10000000ULL; if (D > N) D = N;
    u64 R = argc > 3 ? (u64)strtod(argv[3], 0) : 200000;
    int NT = argc > 4 ? atoi(argv[4]) : 10;
    ref_init();
    double t0 = now_sec(); u64 nq;
    V = vsieve_build(N, NT, &nq);
    double t1 = now_sec(); u64 np2;
    u64 *V2 = vsieve2_build(N, NT, &np2);
    double t2 = now_sec();
    u64 words = (N >> 6) + 1, mism = 0, c1 = 0, c2 = 0;
    for (u64 w = 0; w < words; w++) { c1 += __builtin_popcountll(V[w]); c2 += __builtin_popcountll(V2[w]); if (V[w] != V2[w]) { if (!mism) printf("first word mismatch at n~%llu\n", (unsigned long long)(w * 64)); mism++; } }
    printf("[V sieve vs independent sieve] N=%llu: |V|=%llu vs %llu, primes %llu vs %llu, mismatching words=%llu  (%.1fs + %.1fs)\n",
           (unsigned long long)N, (unsigned long long)c1, (unsigned long long)c2, (unsigned long long)nq, (unsigned long long)np2, (unsigned long long)mism, t1 - t0, t2 - t1);
    pthread_t th[64]; dj_t J[64]; u64 st = D / NT + 1, bad = 0, cnt = 0;
    for (int t = 0; t < NT; t++) { J[t].lo = 1 + t * st; J[t].hi = 1 + (t + 1) * st; if (J[t].hi > D + 1) J[t].hi = D + 1; if (J[t].lo > D + 1) J[t].lo = D + 1; J[t].bad = J[t].cnt = 0; pthread_create(&th[t], NULL, dworker, &J[t]); }
    for (int t = 0; t < NT; t++) { pthread_join(th[t], NULL); bad += J[t].bad; cnt += J[t].cnt; }
    printf("[direct ref_sf(n^2+1)] all n <= %llu: %llu in V (sieve count %llu), mismatches=%llu  (%.1fs)\n", (unsigned long long)D, (unsigned long long)cnt, (unsigned long long)vcount_upto(V, D), (unsigned long long)bad, now_sec() - t2);
    double t3 = now_sec(); rj_t RJ[64]; u64 rbad = 0, rin = 0;
    for (int t = 0; t < NT; t++) { RJ[t].seed = 777 + t; RJ[t].R = R / NT; RJ[t].bad = RJ[t].inV = 0; pthread_create(&th[t], NULL, rworker, &RJ[t]); }
    for (int t = 0; t < NT; t++) { pthread_join(th[t], NULL); rbad += RJ[t].bad; rin += RJ[t].inV; }
    printf("[random direct] %llu samples (half uniform in [1,N], half random members of V; %llu in V): mismatches=%llu  (%.1fs)\n", (unsigned long long)(R / NT * NT), (unsigned long long)rin, (unsigned long long)rbad, now_sec() - t3);
    int ok = !mism && !bad && !rbad && c1 == c2;
    printf("TEST_V %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

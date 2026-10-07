// bench_sqf.c -- throughput of sqf_cert on u*v+1 for various B1 / rho budgets
#include "common.h"
#include "sqf.h"
int main(int argc, char **argv) {
    u64 lo = argc > 1 ? strtoull(argv[1], 0, 10) : 500000000ULL, hi = argc > 2 ? strtoull(argv[2], 0, 10) : 1000000000ULL;
    long N = argc > 3 ? atol(argv[3]) : 2000000;
    u64 B1s[] = {4096, 16384, 65536};
    setvbuf(stdout, NULL, _IOLBF, 0);
    u64 *ms = xmalloc(8 * N); u64 s = 99;
    for (long i = 0; i < N; i++) { u64 u = lo + rng_below(&s, hi - lo), v = lo + rng_below(&s, hi - lo); ms[i] = u * v + 1; }
    // keep only those passing the quick filter (as in the colouring code)
    sqf_init(4096);
    long M = 0; for (long i = 0; i < N; i++) if (!sqf_quick_no(ms[i])) ms[M++] = ms[i];
    printf("u,v in [%llu,%llu): %ld of %ld pass quick filter\n", (unsigned long long)lo, (unsigned long long)hi, M, N);
    for (int bi = 0; bi < 3; bi++) {
        // re-init tables with new B1
        sqf_nB1 = 0; while (sqf_primes[sqf_nB1] <= B1s[bi]) sqf_nB1++;
        u64 budgets[] = {0, 64, 128, 256};
        for (int rb = 0; rb < 4; rb++) {
            sqf_rho_budget = budgets[rb]; sqf_rho_tries = budgets[rb] ? 1 : 0;
            memset(&sqf_st, 0, sizeof sqf_st);
            double t0 = now_sec(); long y = 0;
            for (long i = 0; i < M; i++) y += sqf_cert(ms[i]) == SF_YES;
            double dt = now_sec() - t0;
            printf("B1=%6llu rho=%5llu: %.3f us/cert  yes=%.4f  mr=%.3f prime=%.3f hard=%.4f rho_found=%.4f td2=%.5f rho_it/cert=%.1f ltP3=%.3f\n",
                   (unsigned long long)B1s[bi], (unsigned long long)budgets[rb], dt / M * 1e6, (double)y / M, (double)sqf_st.mr_calls / M, (double)sqf_st.mr_prime / M,
                   (double)sqf_st.hard_calls / M, (double)sqf_st.rho_found / M, (double)sqf_st.td2_calls / M, (double)sqf_st.rho_iters / M, (double)sqf_st.lt_p3_yes / M);
        }
    }
    return 0;
}

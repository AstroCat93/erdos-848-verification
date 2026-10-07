// checker.c -- independent verification of colour's output.  Shares no code with colour.c, sqf.h, vsieve.h:
// vertices come from vsieve2.h (Tonelli-Shanks + direct p^2 | n^2+1 tests), squarefreeness from ref.h
// (trial division to the cube root with hardware '%', plus binary-search integer square root).
//
//   checker replay <logfile> <N> [threads]
//       Replays the move log written by `colour --log`.  Records are (u32 v, u32 k): "v is now in class k".
//       The record whose v is the next vertex of V (in increasing order) ends the step for that vertex;
//       other records are moves of already-placed vertices.  After every step it checks:
//         - every vertex of V up to the current vertex was placed, in order, none skipped;
//         - a C7 vertex (7 mod 25) opens exactly class t(v) = T+1; C7 vertices never move;
//         - other vertices go to classes 1..T (T = t(v));
//         - every vertex placed/moved during the step is squarefree-compatible (ref_sf(uv+1)) with every
//           other member of its final class (these pair checks are batched and run in parallel).
//       By induction the state after every step is a proper colouring of G_n with t(n) classes.
//
//   checker dump <dumpfile> [threads]
//       Verifies a colouring dump (colour --dump-at): T = t(N), class k starts with 25k-18, every member is
//       in V and <= N, every vertex of V_N appears exactly once, every same-class pair is ref-squarefree.
#include "common.h"
#include "ref.h"
#include "vsieve2.h"

static int NTHR = 10;
static inline int inV(const u64 *V, u64 n) { return (V[n >> 6] >> (n & 63)) & 1; }
static inline u64 tN(u64 n) { return (n + 18) / 25; }

// ---------------- parallel pair verification
typedef struct { u32 a, b; } pair_t;
typedef struct { pair_t *P; u64 lo, hi, bad; pair_t first; } pj_t;
static void *pair_worker(void *arg) {
    pj_t *J = arg;
    for (u64 i = J->lo; i < J->hi; i++) {
        u64 a = J->P[i].a, b = J->P[i].b;
        if (!ref_sf(a * b + 1)) { if (!J->bad++) J->first = J->P[i]; }
    }
    return NULL;
}
static u64 total_pairs, total_bad;
static void verify_pairs(pair_t *P, u64 np) {
    if (!np) return;
    pthread_t th[64]; pj_t J[64]; int nt = np < 10000 ? 1 : NTHR;
    u64 step = np / nt + 1;
    for (int t = 0; t < nt; t++) { J[t].P = P; J[t].lo = t * step < np ? t * step : np; J[t].hi = (t + 1) * step < np ? (t + 1) * step : np; J[t].bad = 0; }
    if (nt == 1) pair_worker(&J[0]);
    else { for (int t = 0; t < nt; t++) pthread_create(&th[t], NULL, pair_worker, &J[t]); for (int t = 0; t < nt; t++) pthread_join(th[t], NULL); }
    for (int t = 0; t < nt; t++) if (J[t].bad) {
        printf("CHECK FAIL: pair (%u,%u): %u*%u+1 is NOT squarefree (reference test)\n", J[t].first.a, J[t].first.b, J[t].first.a, J[t].first.b);
        total_bad += J[t].bad;
    }
    total_pairs += np;
}

// ---------------- class storage for replay (independent of colour.c)
#define CCAP 10   /* colour never exceeds 8 members (CAP 7 + c_k); a larger class is reported as an error */
typedef struct { u32 cnt; u32 m[CCAP]; } ccls_t;
// vertex -> class map indexed by the rank of the vertex in V (memory |V|*4 bytes instead of N*4)
static u32 *rank_blk;  // rank_blk[w] = #{V elements < 64*w}
static inline u64 vrank(const u64 *V, u64 n) { return rank_blk[n >> 6] + __builtin_popcountll(V[n >> 6] & ((1ULL << (n & 63)) - 1)); }

static int replay(const char *logp, u64 N) {
    double t0 = now_sec(); u64 np;
    u64 *V = vsieve2_build(N, NTHR, &np);
    u64 nV = 0; for (u64 w = 0; w <= (N >> 6); w++) nV += __builtin_popcountll(V[w]);
    printf("independent V: %llu primes p=1 mod 4, |V_N|=%llu (%.1fs)\n", (unsigned long long)np, (unsigned long long)nV, now_sec() - t0);
    FILE *f = fopen(logp, "rb"); if (!f) die("cannot open %s", logp);
    setvbuf(f, NULL, _IOFBF, 1 << 24);
    u64 TM = tN(N);
    ccls_t *C = xcalloc(TM + 2, sizeof(ccls_t));
    u64 nw = (N >> 6) + 2; rank_blk = xmalloc(sizeof(u32) * nw);
    { u64 acc = 0; for (u64 w = 0; w < nw; w++) { rank_blk[w] = (u32)acc; acc += __builtin_popcountll(V[w]); } }
    u32 *whr = xcalloc(nV + 1, sizeof(u32));  // class of each placed vertex (0 = not placed), by rank
#define where_(n) whr[vrank(V, (n))]
    u64 PB = 1 << 22; pair_t *P = xmalloc(sizeof(pair_t) * (PB + 64)); u64 np2 = 0;
    u32 step_v[256]; int nstep = 0;
    u64 x = 0;  // last fully processed vertex
    u64 next = 1; while (next <= N && !inV(V, next)) next++;
    u64 T = 0, recs = 0, errors = 0, moves = 0, steps = 0;
    u32 r[2];
#define ERR(...) do { if (errors++ < 20) { printf("CHECK FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)
    while (fread(r, 4, 2, f) == 2) {
        recs++;
        u64 v = r[0], k = r[1];
        if (v == next) {
            // placement of the new vertex; ends the step
            if (v % 25 == 7) {
                if (k != tN(v) || k != T + 1) ERR("C7 vertex %llu opened class %llu (T=%llu)", (unsigned long long)v, (unsigned long long)k, (unsigned long long)T);
                T = k;
            } else if (k < 1 || k > T || T != tN(v)) ERR("vertex %llu placed in class %llu (T=%llu, t(v)=%llu)", (unsigned long long)v, (unsigned long long)k, (unsigned long long)T, (unsigned long long)tN(v));
            if (where_(v)) ERR("vertex %llu placed twice", (unsigned long long)v);
            if (k >= 1 && k <= TM) {
                if (C[k].cnt >= CCAP) ERR("class %llu overflow", (unsigned long long)k); else C[k].m[C[k].cnt++] = (u32)v;
                where_(v) = (u32)k;
            }
            if (nstep < 256) step_v[nstep++] = (u32)v; else ERR("step too long at %llu", (unsigned long long)v);
            // verify every vertex placed/moved in this step against its final class
            for (int i = 0; i < nstep; i++) {
                u64 u = step_v[i], kk = where_(u);
                int dup = 0; for (int j = 0; j < i; j++) if (step_v[j] == u) dup = 1;
                if (dup) continue;
                for (int j = 0; j < C[kk].cnt; j++) { u32 w = C[kk].m[j]; if (w != u) P[np2++] = (pair_t){(u32)u, w}; }
                if (C[kk].cnt == 0 || C[kk].m[0] != 25 * kk - 18) ERR("class %llu does not start with its C7 element", (unsigned long long)kk);
            }
            if (np2 >= PB) { verify_pairs(P, np2); np2 = 0; }
            nstep = 0; x = v; steps++;
            next = v + 1; while (next <= N && !inV(V, next)) next++;
        } else if (v < next && v >= 1 && inV(V, v) && where_(v)) {
            // move of an existing vertex
            if (v % 25 == 7) ERR("C7 vertex %llu moved", (unsigned long long)v);
            if (k < 1 || k > T) ERR("move of %llu to invalid class %llu", (unsigned long long)v, (unsigned long long)k);
            else {
                ccls_t *A = &C[where_(v)]; int found = 0;
                for (int j = 0; j < A->cnt; j++) if (A->m[j] == v) { A->m[j] = A->m[--A->cnt]; found = 1; break; }
                if (!found) ERR("internal: %llu not in its class", (unsigned long long)v);
                if (C[k].cnt >= CCAP) ERR("class %llu overflow", (unsigned long long)k); else C[k].m[C[k].cnt++] = (u32)v;
                where_(v) = (u32)k; moves++;
                if (nstep < 256) step_v[nstep++] = (u32)v; else ERR("step too long");
            }
        } else {
            ERR("unexpected record (%llu,%llu): next V vertex is %llu", (unsigned long long)v, (unsigned long long)k, (unsigned long long)next);
            break;
        }
        if (errors > 20) break;
    }
    fclose(f);
    verify_pairs(P, np2);
    int complete = (next > N) && nstep == 0;
    if (!complete) ERR("log ends after vertex %llu; next V vertex %llu <= N not placed", (unsigned long long)x, (unsigned long long)next);
    if (complete && T != tN(N)) ERR("final T=%llu != t(N)=%llu", (unsigned long long)T, (unsigned long long)tN(N));
    printf("REPLAY %s: N=%llu records=%llu steps(vertices)=%llu moves=%llu pairs_checked=%llu bad_pairs=%llu errors=%llu  => %s  (%.1fs)\n",
           logp, (unsigned long long)N, (unsigned long long)recs, (unsigned long long)steps, (unsigned long long)moves,
           (unsigned long long)total_pairs, (unsigned long long)total_bad, (unsigned long long)errors,
           (errors == 0 && total_bad == 0 && complete && steps == nV) ? "PASS (proper t(n)-colouring of G_n for every n <= N)" : "FAIL",
           now_sec() - t0);
    return (errors == 0 && total_bad == 0 && complete && steps == nV) ? 0 : 1;
}

static int check_dump(const char *path) {
    double t0 = now_sec();
    FILE *f = fopen(path, "rb"); if (!f) die("cannot open %s", path);
    setvbuf(f, NULL, _IOFBF, 1 << 24);
    char mg[8]; u64 N = 0, T = 0;
    if (fread(mg, 1, 8, f) != 8 || memcmp(mg, "E848DUMP", 8) || fread(&N, 8, 1, f) != 1 || fread(&T, 8, 1, f) != 1) die("bad dump");
    u64 np; u64 *V = vsieve2_build(N, NTHR, &np);
    u64 nV = 0; for (u64 w = 0; w <= (N >> 6); w++) nV += __builtin_popcountll(V[w]);
    printf("dump N=%llu T=%llu; independent |V_N|=%llu (%.1fs)\n", (unsigned long long)N, (unsigned long long)T, (unsigned long long)nV, now_sec() - t0);
    u64 errors = 0;
    if (T != tN(N)) ERR("T=%llu != t(N)=%llu", (unsigned long long)T, (unsigned long long)tN(N));
    u64 *seen = xcalloc((N >> 6) + 2, 8);
    u64 PB = 1 << 24; pair_t *P = xmalloc(sizeof(pair_t) * (PB + 1024)); u64 npairs = 0, nmem = 0, maxsz = 0;
    for (u64 k = 1; k <= T; k++) {
        u32 sz; if (fread(&sz, 4, 1, f) != 1 || sz < 1 || sz > 1000) die("dump truncated/corrupt at class %llu", (unsigned long long)k);
        u32 m[1000]; if (fread(m, 4, sz, f) != sz) die("dump truncated");
        if (m[0] != 25 * k - 18) ERR("class %llu first member %u != %llu", (unsigned long long)k, m[0], (unsigned long long)(25 * k - 18));
        if (sz > maxsz) maxsz = sz;
        for (u32 i = 0; i < sz; i++) {
            u64 w = m[i];
            if (w < 1 || w > N || !inV(V, w)) ERR("class %llu member %llu not in V_N", (unsigned long long)k, (unsigned long long)w);
            else if ((seen[w >> 6] >> (w & 63)) & 1) ERR("vertex %llu appears twice", (unsigned long long)w);
            else seen[w >> 6] |= 1ULL << (w & 63);
            nmem++;
            for (u32 j = i + 1; j < sz; j++) P[npairs++] = (pair_t){m[i], m[j]};
        }
        if (npairs >= PB) { verify_pairs(P, npairs); npairs = 0; }
        if (errors > 20) break;
    }
    verify_pairs(P, npairs);
    fclose(f);
    if (nmem != nV) ERR("%llu members but |V_N|=%llu", (unsigned long long)nmem, (unsigned long long)nV);
    for (u64 w = 0; w <= (N >> 6); w++) if (seen[w] != V[w]) { ERR("coverage mismatch near %llu", (unsigned long long)(w * 64)); break; }
    int ok = errors == 0 && total_bad == 0;
    printf("DUMP CHECK %s: N=%llu classes=%llu members=%llu max_class=%llu pairs_checked=%llu bad=%llu errors=%llu => %s (%.1fs)\n",
           path, (unsigned long long)N, (unsigned long long)T, (unsigned long long)nmem, (unsigned long long)maxsz, (unsigned long long)total_pairs,
           (unsigned long long)total_bad, (unsigned long long)errors, ok ? "PASS" : "FAIL", now_sec() - t0);
    return ok ? 0 : 1;
}

int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    ref_init();
    if (argc >= 4 && !strcmp(argv[1], "replay")) { if (argc > 4) NTHR = atoi(argv[4]); return replay(argv[2], strtoull(argv[3], 0, 10)); }
    if (argc >= 3 && !strcmp(argv[1], "dump")) { if (argc > 3) NTHR = atoi(argv[3]); return check_dump(argv[2]); }
    fprintf(stderr, "usage: checker replay <log> <N> [threads] | checker dump <dumpfile> [threads]\n");
    return 2;
}

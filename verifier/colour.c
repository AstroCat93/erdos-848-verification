// colour.c -- online proper colouring of G_N with exactly t(N) = floor((N+18)/25) classes, for N = 1..NMAX.
//
// G_N: vertices V_N = {n <= N : n^2+1 not squarefree}; edge a~b (a != b) iff ab+1 not squarefree.
// Invariant maintained after processing every n (i.e. for every N):
//   (i)   every vertex of V_N lies in exactly one class;
//   (ii)  there are exactly t(N) classes, class k contains c_k = 25k-18 (the k-th element of C7 = 7 mod 25);
//   (iii) every pair u != v in a common class has uv+1 CERTIFIED squarefree (sqf_cert == SF_YES).
// Hence omega(G_N) <= chi(G_N) <= t(N) for every N <= NMAX.
//
// Processing n in V: n = 7 mod 25 opens class t(n) = {n}; otherwise n is inserted into an existing class,
// possibly after eviction chains (each evicted vertex is re-homed elsewhere; C7 elements never move).
// Every move of a vertex into a class is (re)certified against every member of the final destination class.
//
// Performance: vertices are processed in chunks.  Phase A (parallel, read-only): for each vertex of the
// chunk the SAME deterministic planner is run "dry" on the pre-chunk state, and every certificate it
// computes is recorded in the vertex's slot.  Phase B (sequential): the real planner runs on the real
// state; certificate lookups hit the recorded results (pure function of the pair), misses are computed.
#include "common.h"
#include "sqf.h"
#include "ref.h"
#include "vsieve.h"
#include <pthread.h>
#include <stdatomic.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>

#define CAP 7      /* max members per class besides c_k (a full class is simply not used as a destination) */
#define MAXE 2     /* max evictions per class */
#define DMAX 4     /* max eviction-chain depth for the randomized search */
#define MAXREC 256 /* certificate records per vertex slot */
#define MAXPLAN 64
#define MAXLOCK 64

typedef struct { u32 cnt; u32 m[CAP]; } cls_t;  // 32 bytes; c_k is implicit

// ------------------------------------------------------------------ global state
static cls_t *CL;           // CL[1..TMAX]
static u64 TMAX, T, freeC18, nC18;
static u64 *V; static u64 NMAX;
static int NTHR = 10;

static inline u64 tN(u64 n) { return (n + 18) / 25; }
static inline u64 c7(u64 k) { return 25 * k - 18; }
static inline int isC18(u64 n) { return n % 25 == 18; }
static inline int isC7(u64 n) { return n % 25 == 7; }

// ------------------------------------------------------------------ statistics (checkpointed)
#define NB 10
typedef struct {
    u64 nV, nC7, nC18, nOut;
    u64 depth_hist[3][DMAX + 3];     // [type 0=C18,1=outsider][depth at which plan found; DMAX+1 = exhaustive]
    u64 plan_size_hist[3][8];        // number of moves (1 = direct)
    u64 tries_hist[3][NB];           // classes examined in Phase B: 1,2,3,4,5-8,9-16,17-32,33-64,65-128,>128
    u64 repairs[3];                  // plans with >= 1 eviction
    u64 c18_moves, out_moves;        // existing vertices moved
    u64 exhaustive, failures;
    u64 max_class, max_class_n;      // largest class size seen (incl. c_k) and where
    u64 cache_hit, cache_miss;       // Phase-B certificate lookups
    u64 commit_checks;
    u64 cert_yes, cert_no, cert_unk; // results of sqf_cert in phase B misses + phase A
    sqf_stats_t sq;                  // accumulated certificate path statistics (all threads)
    double cpu_phaseA, cpu_phaseB;
    u64 audits, audit_pairs;
} stats_t;
static stats_t ST;

// ------------------------------------------------------------------ certificate slots / planner context
typedef struct { u32 n, nrec; u64 key[MAXREC]; uint8_t res[MAXREC]; } pslot_t;
typedef struct { u32 v, from, to; } mv_t;
typedef struct {
    int dry; pslot_t *slot; u64 Tn, free, n; u32 counter;
    mv_t plan[MAXPLAN]; int nplan;
    u32 locked[MAXLOCK]; int nlocked;
    u64 tries; int exhaustive; int level;
    u64 yes, no, unk;
} pctx_t;

static int pair_cert(pctx_t *X, u64 a, u64 b) {
    u64 key = a < b ? (a << 32 | b) : (b << 32 | a);
    pslot_t *S = X->slot;
    if (S) for (u32 i = 0; i < S->nrec; i++) if (S->key[i] == key) { if (!X->dry) ST.cache_hit++; return S->res[i]; }
    int r = sqf_cert(a * b + 1);   // a,b < 2^32 so ab+1 < 2^64
    if (r == SF_YES) X->yes++; else if (r == SF_NO) X->no++; else X->unk++;
    if (S && S->nrec < MAXREC) { S->key[S->nrec] = key; S->res[S->nrec] = (uint8_t)r; S->nrec++; }
    if (!X->dry) ST.cache_miss++;
    return r;
}

static inline int is_locked(const pctx_t *X, u64 k) { for (int i = 0; i < X->nlocked; i++) if (X->locked[i] == k) return 1; return 0; }

// Conflict set of y with class k: members w with y*w+1 not certified squarefree.
// Returns -1 if c_k conflicts or more than maxE conflicts; else the number of conflicts (stored in E).
static int conflicts(pctx_t *X, u64 y, u64 k, int maxE, u32 *E) {
    const cls_t *c = &CL[k];
    u32 cnt = c->cnt; u32 mem[CAP]; for (u32 i = 0; i < cnt; i++) mem[i] = c->m[i];
    u64 ck = c7(k);
    int nE = 0; uint8_t dec[CAP] = {0};
    if (sqf_quick_no(y * ck + 1)) return -1;
    for (u32 i = 0; i < cnt; i++)
        if (sqf_quick_no(y * mem[i] + 1)) { if (nE >= maxE) return -1; E[nE++] = mem[i]; dec[i] = 1; }
    if (pair_cert(X, y, ck) != SF_YES) return -1;
    for (u32 i = 0; i < cnt; i++)
        if (!dec[i] && pair_cert(X, y, mem[i]) != SF_YES) { if (nE >= maxE) return -1; E[nE++] = mem[i]; }
    return nE;
}

static const int KD[DMAX + 1] = {24, 32, 16, 8, 6};  // random candidate classes per depth
static inline u64 seed_of(u64 a, u64 b, u64 c) { u64 s = a * 0x9E3779B97F4A7C15ULL ^ (b + 1) * 0xC2B2AE3D27D4EB4FULL ^ (c + 7) * 0x165667B19E3779F9ULL; splitmix64(&s); return s; }

// Find a home for homeless vertex y (currently in class `from`, 0 = new) using eviction depth <= d.
// On success appends moves to X->plan (and locks the destination classes) and returns 1.
static int plan_home(pctx_t *X, u64 y, u64 from, int d, u64 seed) {
    u64 s = seed;
    int c18 = isC18(y);
    int top = (X->level == 0);
    u64 K = (top && X->exhaustive) ? X->Tn : (u64)KD[d];
    X->level++;
    for (int64_t ci = -1; ; ci++) {
        u64 k;
        if (ci == -1) {  // C18 vertices: the C18-free class first
            if (!c18 || !X->free || is_locked(X, X->free)) continue;
            k = X->free;
        } else if ((u64)ci == K) {  // other vertices: the C18-free class only as a last resort
            if (c18 || !X->free || is_locked(X, X->free)) break;
            k = X->free;
        } else {
            if ((u64)ci > K) break;
            if (c18 && d == 0) break;  // every non-free class already holds a C18 element (adjacent to y)
            k = (top && X->exhaustive) ? (u64)ci + 1 : 1 + rng_below(&s, X->Tn);
            if (k == X->free || is_locked(X, k)) continue;  // free class: reserved for C18 / tried at ci=-1
        }
        X->tries++;
        u32 E[MAXE]; int maxE = d > 0 ? MAXE : 0;
        int nE = conflicts(X, y, k, maxE, E);
        if (nE < 0) continue;
        if (CL[k].cnt + 1 - nE > CAP) continue;
        int mark = X->nplan, lmark = X->nlocked;
        if (X->nplan >= MAXPLAN || X->nlocked >= MAXLOCK) continue;
        X->plan[X->nplan++] = (mv_t){(u32)y, (u32)from, (u32)k};
        X->locked[X->nlocked++] = (u32)k;
        if (nE == 0) { X->level--; return 1; }
        int ok = 1;
        for (int e = 0; e < nE && ok; e++)
            if (!plan_home(X, E[e], k, d - 1, seed_of(E[e], d - 1, X->n * 64 + (++X->counter)))) ok = 0;
        if (ok) { X->level--; return 1; }
        X->nplan = mark; X->nlocked = lmark;
    }
    X->level--;
    return 0;
}

static inline u64 free_pred(u64 n) { u64 nc18 = (n + 6) / 25; return nc18 < tN(n) ? tN(n) : 0; }  // #C18 < n vs t(n)

// ------------------------------------------------------------------ replay log
static FILE *LOGF; static u64 log_records;
static inline void log_rec(u32 v, u32 k) { u32 r[2] = {v, k}; if (LOGF) { fwrite(r, 4, 2, LOGF); log_records++; } }

// ------------------------------------------------------------------ class operations
static inline int hasC18(u64 k) { for (u32 i = 0; i < CL[k].cnt; i++) if (CL[k].m[i] % 25 == 18) return 1; return 0; }
static void cls_remove(u64 k, u32 v) {
    cls_t *c = &CL[k];
    for (u32 i = 0; i < c->cnt; i++) if (c->m[i] == v) { c->m[i] = c->m[--c->cnt]; return; }
    die("cls_remove: vertex %u not in class %llu", v, (unsigned long long)k);
}
static void cls_add(u64 k, u32 v, u64 n) {
    cls_t *c = &CL[k];
    if (c->cnt >= CAP) die("cls_add: class %llu full", (unsigned long long)k);
    c->m[c->cnt++] = v;
    if (c->cnt + 1 > ST.max_class) { ST.max_class = c->cnt + 1; ST.max_class_n = n; }
}

static void commit(pctx_t *X, u64 n) {
    for (int i = 0; i < X->nplan; i++) if (X->plan[i].from) cls_remove(X->plan[i].from, X->plan[i].v);
    for (int i = 0; i < X->nplan; i++) cls_add(X->plan[i].to, X->plan[i].v, n);
    // re-certify every moved vertex against every member of its final destination class
    for (int i = 0; i < X->nplan; i++) {
        u64 v = X->plan[i].v, k = X->plan[i].to;
        ST.commit_checks++;
        if (pair_cert(X, v, c7(k)) != SF_YES) die("COMMIT CHECK FAILED n=%llu v=%llu class=%llu vs c_k", (unsigned long long)n, (unsigned long long)v, (unsigned long long)k);
        for (u32 j = 0; j < CL[k].cnt; j++) {
            u64 w = CL[k].m[j]; if (w == v) continue;
            ST.commit_checks++;
            if (pair_cert(X, v, w) != SF_YES) die("COMMIT CHECK FAILED n=%llu v=%llu class=%llu vs %llu", (unsigned long long)n, (unsigned long long)v, (unsigned long long)k, (unsigned long long)w);
        }
    }
    // replay log: moves of existing vertices first, then the new vertex
    for (int i = 0; i < X->nplan; i++) if (X->plan[i].from) {
        log_rec(X->plan[i].v, X->plan[i].to);
        if (isC18(X->plan[i].v)) ST.c18_moves++; else ST.out_moves++;
    }
    for (int i = 0; i < X->nplan; i++) if (!X->plan[i].from) log_rec(X->plan[i].v, X->plan[i].to);
    // C18 bookkeeping: exactly T - nC18 (0 or 1) classes lack a C18 element
    if (isC18(n)) nC18++;
    u64 cand[2 * MAXPLAN + 1]; int nc = 0;
    if (freeC18) cand[nc++] = freeC18;
    for (int i = 0; i < X->nplan; i++) { cand[nc++] = X->plan[i].to; if (X->plan[i].from) cand[nc++] = X->plan[i].from; }
    u64 nf = 0, f = 0;
    for (int i = 0; i < nc; i++) {
        int dup = 0; for (int j = 0; j < i; j++) if (cand[j] == cand[i]) dup = 1;
        if (dup) continue;
        if (!hasC18(cand[i])) { nf++; f = cand[i]; }
    }
    if (nf != T - nC18) die("C18 bookkeeping mismatch at n=%llu: found %llu free, expected %llu", (unsigned long long)n, (unsigned long long)nf, (unsigned long long)(T - nC18));
    freeC18 = nf ? f : 0;
}

static inline int tries_bucket(u64 t) {
    if (t <= 4) return t ? (int)t - 1 : 0;
    if (t <= 8) return 4; if (t <= 16) return 5; if (t <= 32) return 6; if (t <= 64) return 7; if (t <= 128) return 8; return 9;
}

static void dump_failure(u64 n);

// Phase B: process one vertex for real
static void process_vertex(u64 n, pslot_t *slot) {
    ST.nV++;
    if (isC7(n)) {
        ST.nC7++;
        u64 k = tN(n);
        if (k != T + 1) die("C7 vertex %llu: expected class %llu, T=%llu", (unsigned long long)n, (unsigned long long)k, (unsigned long long)T);
        if (T != nC18) die("C7 vertex %llu arrives while a C18-free class exists", (unsigned long long)n);
        if (k > TMAX) die("TMAX exceeded");
        T = k; CL[k].cnt = 0; freeC18 = k;
        log_rec((u32)n, (u32)k);
        return;
    }
    int type = isC18(n) ? 0 : 1;
    if (type == 0) ST.nC18++; else ST.nOut++;
    static pctx_t X;
    X.dry = 0; X.slot = slot; X.Tn = T; X.free = freeC18; X.n = n; X.tries = 0; X.exhaustive = 0; X.yes = X.no = X.unk = 0;
    if (T != tN(n)) die("T mismatch at n=%llu", (unsigned long long)n);
    int found = -1;
    for (int d = 0; d <= DMAX && found < 0; d++) {
        X.counter = 0; X.nplan = X.nlocked = 0; X.level = 0;
        if (plan_home(&X, n, 0, d, seed_of(n, d, 0))) found = d;
    }
    if (found < 0) {  // emergency: exhaustive top-level search at depth 1, then 2
        ST.exhaustive++;
        fprintf(stdout, "NOTE: n=%llu needed exhaustive repair search (tries so far %llu)\n", (unsigned long long)n, (unsigned long long)X.tries);
        fflush(stdout);
        for (int d = 1; d <= 2 && found < 0; d++) {
            X.counter = 0; X.nplan = X.nlocked = 0; X.level = 0; X.exhaustive = 1;
            if (plan_home(&X, n, 0, d, seed_of(n, d, 99))) found = DMAX + 1;
        }
        X.exhaustive = 0;
    }
    ST.cert_yes += X.yes; ST.cert_no += X.no; ST.cert_unk += X.unk;
    if (found < 0) {
        ST.failures++;
        dump_failure(n);
        die("PLACEMENT FAILURE at n=%llu (no proper placement found even after exhaustive search)", (unsigned long long)n);
    }
    ST.depth_hist[type][found]++;
    ST.plan_size_hist[type][X.nplan < 8 ? X.nplan : 7]++;
    ST.tries_hist[type][tries_bucket(X.tries)]++;
    if (X.nplan > 1) ST.repairs[type]++;
    commit(&X, n);
}

// ------------------------------------------------------------------ Phase A worker pool
static u32 *chunkV; static u64 nchunk; static pslot_t *slots; static u64 slots_cap;
static atomic_ullong phaseA_next; static atomic_int phaseA_done;
static pthread_mutex_t pool_mx = PTHREAD_MUTEX_INITIALIZER; static pthread_cond_t pool_cv = PTHREAD_COND_INITIALIZER, done_cv = PTHREAD_COND_INITIALIZER;
static u64 pool_gen; static int pool_quit;
static sqf_stats_t *thr_sq; static u64 *thr_yes, *thr_no, *thr_unk;
#define DPRE 2

static void phaseA_vertex(u64 idx, pctx_t *X) {
    u64 n = chunkV[idx];
    pslot_t *S = &slots[idx]; S->n = (u32)n; S->nrec = 0;
    if (isC7(n)) return;
    X->dry = 1; X->slot = S; X->Tn = tN(n); X->free = free_pred(n); X->n = n; X->tries = 0; X->exhaustive = 0;
    for (int d = 0; d <= DPRE; d++) {
        X->counter = 0; X->nplan = X->nlocked = 0; X->level = 0;
        if (plan_home(X, n, 0, d, seed_of(n, d, 0))) break;
    }
}
static void *pool_worker(void *arg) {
    int id = (int)(intptr_t)arg; u64 seen = 0;
    pctx_t *X = xcalloc(1, sizeof(pctx_t));
    for (;;) {
        pthread_mutex_lock(&pool_mx);
        while (pool_gen == seen && !pool_quit) pthread_cond_wait(&pool_cv, &pool_mx);
        if (pool_quit) { pthread_mutex_unlock(&pool_mx); break; }
        seen = pool_gen;
        pthread_mutex_unlock(&pool_mx);
        X->yes = X->no = X->unk = 0;
        for (;;) {
            u64 i0 = atomic_fetch_add(&phaseA_next, 8);
            if (i0 >= nchunk) break;
            u64 i1 = i0 + 8 < nchunk ? i0 + 8 : nchunk;
            for (u64 i = i0; i < i1; i++) phaseA_vertex(i, X);
        }
        thr_sq[id] = sqf_st; thr_yes[id] += X->yes; thr_no[id] += X->no; thr_unk[id] += X->unk;
        pthread_mutex_lock(&pool_mx);
        if (atomic_fetch_add(&phaseA_done, 1) + 1 == NTHR) pthread_cond_signal(&done_cv);
        pthread_mutex_unlock(&pool_mx);
    }
    free(X);
    return NULL;
}
static pthread_t *pool_th;
static void pool_start(void) {
    pool_th = xmalloc(sizeof(pthread_t) * NTHR);
    thr_sq = xcalloc(NTHR, sizeof(sqf_stats_t)); thr_yes = xcalloc(NTHR, 8); thr_no = xcalloc(NTHR, 8); thr_unk = xcalloc(NTHR, 8);
    for (int t = 0; t < NTHR; t++) pthread_create(&pool_th[t], NULL, pool_worker, (void *)(intptr_t)t);
}
static void pool_run_phaseA(void) {
    atomic_store(&phaseA_next, 0); atomic_store(&phaseA_done, 0);
    pthread_mutex_lock(&pool_mx);
    pool_gen++;
    pthread_cond_broadcast(&pool_cv);
    while (atomic_load(&phaseA_done) < NTHR) pthread_cond_wait(&done_cv, &pool_mx);
    pthread_mutex_unlock(&pool_mx);
}
static void pool_stop(void) {
    pthread_mutex_lock(&pool_mx); pool_quit = 1; pthread_cond_broadcast(&pool_cv); pthread_mutex_unlock(&pool_mx);
    for (int t = 0; t < NTHR; t++) pthread_join(pool_th[t], NULL);
}
// total certificate statistics (main thread + workers)
static void sum_sq(sqf_stats_t *out, u64 *yes, u64 *no, u64 *unk) {
    *out = ST.sq;
    u64 *a = (u64 *)out; const u64 *b = (const u64 *)&sqf_st;
    for (size_t k = 0; k < sizeof(sqf_stats_t) / 8; k++) a[k] += b[k];
    *yes = ST.cert_yes; *no = ST.cert_no; *unk = ST.cert_unk;
    for (int t = 0; t < NTHR && thr_sq; t++) {
        b = (const u64 *)&thr_sq[t]; for (size_t k = 0; k < sizeof(sqf_stats_t) / 8; k++) a[k] += b[k];
        *yes += thr_yes[t]; *no += thr_no[t]; *unk += thr_unk[t];
    }
}

// ------------------------------------------------------------------ audit (independent reference squarefree test)
typedef struct { u64 k0, k1; u64 pairs, bad; u64 first_bad_a, first_bad_b; } aud_t;
static double audit_sample = 1.0;
static void *audit_worker(void *arg) {
    aud_t *A = arg;
    u64 s = A->k0 * 0x9E3779B97F4A7C15ULL + 12345;
    for (u64 k = A->k0; k < A->k1; k++) {
        if (audit_sample < 1.0 && (double)(splitmix64(&s) >> 11) / 9007199254740992.0 >= audit_sample) continue;
        u64 mem[CAP + 1]; u32 cnt = CL[k].cnt; mem[0] = c7(k);
        for (u32 i = 0; i < cnt; i++) mem[i + 1] = CL[k].m[i];
        for (u32 i = 0; i <= cnt; i++) for (u32 j = i + 1; j <= cnt; j++) {
            A->pairs++;
            if (!ref_sf(mem[i] * mem[j] + 1)) { if (!A->bad++) { A->first_bad_a = mem[i]; A->first_bad_b = mem[j]; } }
        }
    }
    return NULL;
}
// full structural audit + pair audit; n = last processed integer
static int audit(u64 n, int verbose) {
    double t0 = now_sec();
    int ok = 1;
    if (T != tN(n)) { printf("AUDIT FAIL: T=%llu != t(%llu)=%llu\n", (unsigned long long)T, (unsigned long long)n, (unsigned long long)tN(n)); ok = 0; }
    u64 words = (n >> 6) + 1;
    u64 *seen = xcalloc(words + 1, 8);
    u64 nmem = 0, nfree = 0, maxc = 0, hist[CAP + 2] = {0};
    for (u64 k = 1; k <= T; k++) {
        u64 ck = c7(k);
        if (!vbit(V, ck) || ck > n) { printf("AUDIT FAIL: c_k=%llu invalid\n", (unsigned long long)ck); ok = 0; }
        if ((seen[ck >> 6] >> (ck & 63)) & 1) { printf("AUDIT FAIL: duplicate %llu\n", (unsigned long long)ck); ok = 0; }
        seen[ck >> 6] |= 1ULL << (ck & 63); nmem++;
        int c18 = 0;
        if (CL[k].cnt > CAP) { printf("AUDIT FAIL: class %llu cnt %u\n", (unsigned long long)k, CL[k].cnt); ok = 0; continue; }
        for (u32 i = 0; i < CL[k].cnt; i++) {
            u64 w = CL[k].m[i];
            if (w > n || w == 0 || !vbit(V, w) || isC7(w)) { printf("AUDIT FAIL: bad member %llu in class %llu\n", (unsigned long long)w, (unsigned long long)k); ok = 0; continue; }
            if ((seen[w >> 6] >> (w & 63)) & 1) { printf("AUDIT FAIL: duplicate %llu\n", (unsigned long long)w); ok = 0; }
            seen[w >> 6] |= 1ULL << (w & 63); nmem++;
            c18 += isC18(w);
        }
        if (c18 > 1) { printf("AUDIT FAIL: class %llu has %d C18 elements\n", (unsigned long long)k, c18); ok = 0; }
        if (!c18) nfree++;
        if (CL[k].cnt + 1 > maxc) maxc = CL[k].cnt + 1;
        hist[CL[k].cnt + 1]++;
    }
    // seen must equal V restricted to [1,n]
    u64 mism = 0;
    for (u64 w = 0; w < words; w++) {
        u64 vw = V[w];
        if (w == words - 1) { u64 rem = (n & 63) + 1; if (rem < 64) vw &= (1ULL << rem) - 1; }
        if (w == 0) vw &= ~1ULL;
        if (vw != seen[w]) mism++;
    }
    if (mism) { printf("AUDIT FAIL: coverage mismatch in %llu words\n", (unsigned long long)mism); ok = 0; }
    u64 nv = vcount_upto(V, n);
    if (nv != nmem) { printf("AUDIT FAIL: |V cap [1,n]|=%llu but %llu members\n", (unsigned long long)nv, (unsigned long long)nmem); ok = 0; }
    u64 c18cnt = (n + 7) / 25;
    if (nfree != T - c18cnt) { printf("AUDIT FAIL: %llu C18-free classes, expected %llu\n", (unsigned long long)nfree, (unsigned long long)(T - c18cnt)); ok = 0; }
    free(seen);
    // pairs with the independent reference test
    int nt = (T > 1000) ? NTHR : 1;
    pthread_t th[64]; aud_t A[64];
    u64 step = T / nt + 1;
    for (int t = 0; t < nt; t++) { A[t] = (aud_t){1 + t * step, 1 + (t + 1) * step > T + 1 ? T + 1 : 1 + (t + 1) * step, 0, 0, 0, 0}; if (A[t].k0 > T + 1) A[t].k0 = T + 1; }
    if (nt == 1) audit_worker(&A[0]);
    else { for (int t = 0; t < nt; t++) pthread_create(&th[t], NULL, audit_worker, &A[t]); for (int t = 0; t < nt; t++) pthread_join(th[t], NULL); }
    u64 pairs = 0, bad = 0;
    for (int t = 0; t < nt; t++) { pairs += A[t].pairs; if (A[t].bad && !bad) printf("AUDIT FAIL: pair %llu,%llu not squarefree by reference test\n", (unsigned long long)A[t].first_bad_a, (unsigned long long)A[t].first_bad_b); bad += A[t].bad; }
    if (bad) ok = 0;
    ST.audits++; ST.audit_pairs += pairs;
    if (verbose || !ok)
        printf("AUDIT N=%llu: %s  classes=%llu (=t(N)) members=%llu (=|V_N|) maxclass=%llu C18free=%llu pairs_checked=%llu%s bad=%llu (%.1fs)\n",
               (unsigned long long)n, ok ? "PASS" : "FAIL", (unsigned long long)T, (unsigned long long)nmem, (unsigned long long)maxc,
               (unsigned long long)nfree, (unsigned long long)pairs, audit_sample < 1.0 ? " (sampled)" : "", (unsigned long long)bad, now_sec() - t0);
    if (verbose || !ok) {
        printf("   class size histogram (incl. c_k):");
        for (int s = 1; s <= CAP + 1; s++) printf(" %d:%llu", s, (unsigned long long)hist[s]);
        printf("\n");
    }
    fflush(stdout);
    return ok;
}

// dump colouring: header "E848DUMP", u64 N, u64 T; then per class: u32 size, size x u32 members (c_k first)
static void dump_colouring(const char *path, u64 n) {
    FILE *f = fopen(path, "wb"); if (!f) die("cannot open %s", path);
    setvbuf(f, NULL, _IOFBF, 1 << 22);
    fwrite("E848DUMP", 1, 8, f); fwrite(&n, 8, 1, f); fwrite(&T, 8, 1, f);
    for (u64 k = 1; k <= T; k++) {
        u32 sz = CL[k].cnt + 1, ck = (u32)c7(k);
        fwrite(&sz, 4, 1, f); fwrite(&ck, 4, 1, f); fwrite(CL[k].m, 4, CL[k].cnt, f);
    }
    if (fclose(f)) die("write error %s", path);
    printf("DUMP N=%llu written to %s\n", (unsigned long long)n, path); fflush(stdout);
}
static char dump_dir[1024] = ".";
static void dump_failure(u64 n) {
    char p[1200]; snprintf(p, sizeof p, "%s/FAILURE_state_%llu.bin", dump_dir, (unsigned long long)n);
    dump_colouring(p, n - 1);
}

// ------------------------------------------------------------------ checkpoint
#define CKPT_VER 2
typedef struct { char magic[8]; u64 ver, n_next, T, freeC18, nC18, log_off, log_records, cap; stats_t st; } ckhdr_t;
static void save_ckpt(const char *path, u64 n_next) {
    char tmp[1200]; snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "wb"); if (!f) die("cannot open %s", tmp);
    ckhdr_t H; memset(&H, 0, sizeof H); memcpy(H.magic, "E848CKPT", 8);
    H.ver = CKPT_VER; H.n_next = n_next; H.T = T; H.freeC18 = freeC18; H.nC18 = nC18; H.cap = CAP;
    if (LOGF) {
        fflush(LOGF); off_t o = ftello(LOGF);
        if (o < 0 || (u64)o != 8 * log_records) die("log offset %lld inconsistent with %llu records", (long long)o, (unsigned long long)log_records);
        H.log_off = (u64)o;
    }
    H.log_records = log_records;
    sum_sq(&H.st.sq, &H.st.cert_yes, &H.st.cert_no, &H.st.cert_unk);
    stats_t s2 = ST; s2.sq = H.st.sq; s2.cert_yes = H.st.cert_yes; s2.cert_no = H.st.cert_no; s2.cert_unk = H.st.cert_unk; H.st = s2;
    fwrite(&H, sizeof H, 1, f);
    u64 off = 1;
    while (off <= T) { u64 c = T + 1 - off; if (c > (1 << 22)) c = 1 << 22; if (fwrite(&CL[off], sizeof(cls_t), c, f) != c) die("ckpt write"); off += c; }
    if (fflush(f) || fsync(fileno(f)) || fclose(f)) die("ckpt write error");
    if (rename(tmp, path)) die("rename ckpt");
    printf("CHECKPOINT saved: %s (next n=%llu, T=%llu)\n", path, (unsigned long long)n_next, (unsigned long long)T); fflush(stdout);
}

static volatile sig_atomic_t stop_flag = 0;
static void on_signal(int s) { (void)s; stop_flag = 1; }

// ------------------------------------------------------------------ progress
static double T0, wall_prev; static u64 n_prev;
static void progress(u64 n) {
    double t = now_sec();
    sqf_stats_t S; u64 y, no, un; sum_sq(&S, &y, &no, &un);
    u64 tot18 = ST.depth_hist[0][0] + ST.depth_hist[0][1] + ST.depth_hist[0][2] + ST.depth_hist[0][3] + ST.depth_hist[0][4] + ST.depth_hist[0][5];
    (void)tot18;
    printf("N=%llu verified | t(N)=%llu |V_N|=%llu (C7 %llu, C18 %llu, other %llu) | %.1fs (+%.1fs, %.2f M/s) | repairs C18=%llu out=%llu | depth C18[%llu %llu %llu %llu %llu|ex %llu] out[%llu %llu %llu|ex %llu] | maxclass=%llu | certs yes=%llu no=%llu UNCERTAIN=%llu td2=%llu | phaseB hit=%llu miss=%llu\n",
           (unsigned long long)n, (unsigned long long)T, (unsigned long long)ST.nV, (unsigned long long)ST.nC7, (unsigned long long)ST.nC18, (unsigned long long)ST.nOut,
           t - T0, t - wall_prev, (n - n_prev) / (t - wall_prev) / 1e6,
           (unsigned long long)ST.repairs[0], (unsigned long long)ST.repairs[1],
           (unsigned long long)ST.depth_hist[0][0], (unsigned long long)ST.depth_hist[0][1], (unsigned long long)ST.depth_hist[0][2], (unsigned long long)ST.depth_hist[0][3], (unsigned long long)ST.depth_hist[0][4], (unsigned long long)ST.depth_hist[0][5],
           (unsigned long long)ST.depth_hist[1][0], (unsigned long long)ST.depth_hist[1][1], (unsigned long long)ST.depth_hist[1][2], (unsigned long long)ST.depth_hist[1][5],
           (unsigned long long)ST.max_class, (unsigned long long)y, (unsigned long long)no, (unsigned long long)un, (unsigned long long)S.td2_calls,
           (unsigned long long)ST.cache_hit, (unsigned long long)ST.cache_miss);
    fflush(stdout);
    wall_prev = t; n_prev = n;
}
static void final_report(u64 n) {
    sqf_stats_t S; u64 y, no, un; sum_sq(&S, &y, &no, &un);
    const char *tn[2] = {"C18", "other"};
    printf("\n==== FINAL REPORT: all N <= %llu verified (online invariant held after every n) ====\n", (unsigned long long)n);
    printf("t(N)=%llu  |V_N|=%llu  C7=%llu C18=%llu other=%llu\n", (unsigned long long)T, (unsigned long long)ST.nV, (unsigned long long)ST.nC7, (unsigned long long)ST.nC18, (unsigned long long)ST.nOut);
    for (int ty = 0; ty < 2; ty++) {
        u64 tot = ty ? ST.nOut : ST.nC18;
        printf("[%s] placements=%llu repairs(plans with evictions)=%llu (%.3f%%)\n", tn[ty], (unsigned long long)tot, (unsigned long long)ST.repairs[ty], tot ? 100.0 * ST.repairs[ty] / tot : 0);
        printf("   depth at which plan found:"); for (int d = 0; d <= DMAX; d++) printf(" d%d=%llu", d, (unsigned long long)ST.depth_hist[ty][d]); printf(" exhaustive=%llu\n", (unsigned long long)ST.depth_hist[ty][DMAX + 1]);
        printf("   plan size (#moves):"); for (int s = 1; s < 8; s++) printf(" %d%s=%llu", s, s == 7 ? "+" : "", (unsigned long long)ST.plan_size_hist[ty][s]); printf("\n");
        const char *bn[NB] = {"1", "2", "3", "4", "5-8", "9-16", "17-32", "33-64", "65-128", ">128"};
        printf("   classes examined:"); for (int b = 0; b < NB; b++) printf(" %s:%llu", bn[b], (unsigned long long)ST.tries_hist[ty][b]); printf("\n");
        u64 gt[NB]; u64 acc = 0; for (int b = NB - 1; b >= 0; b--) { acc += ST.tries_hist[ty][b]; gt[b] = acc; }
        if (tot) printf("   fraction needing >1: %.4f  >2: %.4f  >4: %.5f  >8: %.6f  >16: %.7f  >32: %.8f\n", (double)gt[1] / tot, (double)gt[2] / tot, (double)gt[4] / tot, (double)gt[5] / tot, (double)gt[6] / tot, (double)gt[7] / tot);
    }
    printf("existing vertices moved: C18=%llu other=%llu ; exhaustive searches=%llu ; placement failures=%llu\n", (unsigned long long)ST.c18_moves, (unsigned long long)ST.out_moves, (unsigned long long)ST.exhaustive, (unsigned long long)ST.failures);
    printf("max class size (incl. c_k)=%llu (first reached at n=%llu)\n", (unsigned long long)ST.max_class, (unsigned long long)ST.max_class_n);
    printf("certificates: computed=%llu YES=%llu NO=%llu UNCERTAIN=%llu ; phase-B lookups hit=%llu miss=%llu ; commit re-checks=%llu\n",
           (unsigned long long)(y + no + un), (unsigned long long)y, (unsigned long long)no, (unsigned long long)un, (unsigned long long)ST.cache_hit, (unsigned long long)ST.cache_miss, (unsigned long long)ST.commit_checks);
    printf("cert paths: calls=%llu td_no=%llu td_early_yes=%llu cof_small_yes=%llu sq_no=%llu <P^3_yes=%llu MR=%llu(prime %llu) cube_no=%llu hard=%llu rho=%llu(found %llu, iters %llu) td2=%llu(found %llu, none->yes %llu) unk=%llu quick=%llu(no %llu)\n",
           (unsigned long long)S.calls, (unsigned long long)S.td_no, (unsigned long long)S.td_early_yes, (unsigned long long)S.cof_small_yes, (unsigned long long)S.sq_no,
           (unsigned long long)S.lt_p3_yes, (unsigned long long)S.mr_calls, (unsigned long long)S.mr_prime, (unsigned long long)S.cube_no, (unsigned long long)S.hard_calls,
           (unsigned long long)S.rho_calls, (unsigned long long)S.rho_found, (unsigned long long)S.rho_iters, (unsigned long long)S.td2_calls, (unsigned long long)S.td2_found,
           (unsigned long long)S.td2_none_yes, (unsigned long long)S.unk, (unsigned long long)S.quick_calls, (unsigned long long)S.quick_no);
    printf("audits=%llu (pairs re-checked with reference test: %llu)\n", (unsigned long long)ST.audits, (unsigned long long)ST.audit_pairs);
    fflush(stdout);
}

// ------------------------------------------------------------------ main
static u64 parse_u64(const char *s) {  // accepts 1e9, 4294967295, 1000000
    char *e; double d = strtod(s, &e);
    if (strchr(s, 'e') || strchr(s, 'E')) return (u64)(d + 0.5);
    return strtoull(s, NULL, 10);
}
static u64 *parse_list(const char *s, int *cnt) {
    u64 *a = xmalloc(sizeof(u64) * 1024); *cnt = 0;
    char *dup = strdup(s), *tok = strtok(dup, ",");
    while (tok && *cnt < 1024) { a[(*cnt)++] = parse_u64(tok); tok = strtok(NULL, ","); }
    free(dup); return a;
}

int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    u64 prog = 10000000, ckpt_every = 0, audit_every = 0, audit_vertex_upto = 0, span = 32768, B1 = 4096;
    const char *ckpt = NULL, *resume = NULL, *logp = NULL;
    u64 *audit_at = NULL, *dump_at = NULL; int n_audit_at = 0, n_dump_at = 0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i]; const char *v = i + 1 < argc ? argv[i + 1] : NULL;
#define ARG(name) (!strcmp(a, name) && v && (i++, 1))
        if (ARG("--nmax")) NMAX = parse_u64(v);
        else if (ARG("--threads")) NTHR = atoi(v);
        else if (ARG("--progress")) prog = parse_u64(v);
        else if (ARG("--ckpt")) ckpt = v;
        else if (ARG("--ckpt-every")) ckpt_every = parse_u64(v);
        else if (ARG("--resume")) resume = v;
        else if (ARG("--log")) logp = v;
        else if (ARG("--audit-every")) audit_every = parse_u64(v);
        else if (ARG("--audit-at")) audit_at = parse_list(v, &n_audit_at);
        else if (ARG("--audit-sample")) audit_sample = atof(v);
        else if (ARG("--audit-every-vertex-upto")) audit_vertex_upto = parse_u64(v);
        else if (ARG("--dump-at")) dump_at = parse_list(v, &n_dump_at);
        else if (ARG("--dump-dir")) snprintf(dump_dir, sizeof dump_dir, "%s", v);
        else if (ARG("--span")) span = parse_u64(v);
        else if (ARG("--B1")) B1 = parse_u64(v);
        else if (ARG("--rho-budget")) sqf_rho_budget = parse_u64(v);
        else if (ARG("--rho-tries")) sqf_rho_tries = atoi(v);
        else die("unknown/incomplete argument: %s", a);
    }
    if (NMAX < 1 || NMAX >= (1ULL << 32)) die("--nmax must be in [1, 2^32-1]");
    if (NTHR < 1) NTHR = 1;
    if (NTHR > 64) NTHR = 64;
    if (span < 1) span = 1;
    signal(SIGTERM, on_signal); signal(SIGINT, on_signal);
    T0 = now_sec(); wall_prev = T0;
    printf("colour: NMAX=%llu threads=%d B1=%llu rho=%llux%d span=%llu CAP=%d MAXE=%d DMAX=%d KD={%d,%d,%d,%d,%d} DPRE=%d\n",
           (unsigned long long)NMAX, NTHR, (unsigned long long)B1, (unsigned long long)sqf_rho_budget, sqf_rho_tries, (unsigned long long)span, CAP, MAXE, DMAX, KD[0], KD[1], KD[2], KD[3], KD[4], DPRE);
    sqf_init(B1); ref_init();
    u64 nq = 0; double ts = now_sec();
    V = vsieve_build(NMAX, NTHR, &nq);
    printf("V sieve: %llu primes q=1 mod 4 used, |V_NMAX|=%llu, %.1fs\n", (unsigned long long)nq, (unsigned long long)vcount_upto(V, NMAX), now_sec() - ts);
    TMAX = tN(NMAX);
    CL = (cls_t *)xcalloc(TMAX + 2, sizeof(cls_t));
    u64 n_start = 1;
    if (resume) {
        FILE *f = fopen(resume, "rb"); if (!f) die("cannot open checkpoint %s", resume);
        ckhdr_t H; if (fread(&H, sizeof H, 1, f) != 1 || memcmp(H.magic, "E848CKPT", 8) || H.ver != CKPT_VER || H.cap != CAP) die("bad checkpoint");
        if (H.T > TMAX) die("checkpoint has more classes than t(NMAX)");
        T = H.T; freeC18 = H.freeC18; nC18 = H.nC18; ST = H.st; n_start = H.n_next;
        if (fread(&CL[1], sizeof(cls_t), T, f) != T) die("checkpoint truncated");
        fclose(f);
        if (logp) {
            struct stat sb;
            if (H.log_off != 8 * H.log_records) die("checkpoint log offset %llu inconsistent with %llu records", (unsigned long long)H.log_off, (unsigned long long)H.log_records);
            if (stat(logp, &sb) || (u64)sb.st_size < H.log_off) die("log %s missing or shorter than checkpoint offset %llu", logp, (unsigned long long)H.log_off);
            if (truncate(logp, (off_t)H.log_off)) die("cannot truncate log %s", logp);
            LOGF = fopen(logp, "ab"); log_records = H.log_records;
        }
        printf("RESUMED from %s at n=%llu (T=%llu)\n", resume, (unsigned long long)n_start, (unsigned long long)T);
        if (!audit(n_start - 1, 1)) die("audit of resumed state failed");
    } else if (logp) LOGF = fopen(logp, "wb");
    if (logp && !LOGF) die("cannot open log %s", logp);
    if (LOGF) setvbuf(LOGF, NULL, _IOFBF, 1 << 24);
    n_prev = n_start - 1;
    slots_cap = span + 64;
    slots = (pslot_t *)xmalloc(sizeof(pslot_t) * slots_cap);
    chunkV = (u32 *)xmalloc(sizeof(u32) * slots_cap);
    pool_start();

    u64 n = n_start;
    while (n <= NMAX && !stop_flag) {
        // next boundary at which some action happens (after processing it)
        u64 hi = n + span - 1; if (hi > NMAX) hi = NMAX;
#define BOUND(every) do { if (every) { u64 b_ = ((n - 1) / (every) + 1) * (every); if (b_ < hi) hi = b_; } } while (0)
        BOUND(prog);
        if (ckpt) BOUND(ckpt_every);
        BOUND(audit_every);
        for (int i = 0; i < n_audit_at; i++) if (audit_at[i] >= n && audit_at[i] < hi) hi = audit_at[i];
        for (int i = 0; i < n_dump_at; i++) if (dump_at[i] >= n && dump_at[i] < hi) hi = dump_at[i];
        if (audit_vertex_upto && n <= audit_vertex_upto) { if (hi > n) hi = n; }
        // collect vertices in [n, hi]
        nchunk = 0;
        for (u64 x = n; x <= hi;) {
            u64 w = V[x >> 6] >> (x & 63);
            if (!w) { x = (x | 63) + 1; continue; }
            x += __builtin_ctzll(w);
            if (x > hi) break;
            chunkV[nchunk++] = (u32)x; x++;
        }
        double ta = now_sec();
        if (nchunk > 4 && NTHR > 1) pool_run_phaseA();
        else for (u64 i = 0; i < nchunk; i++) { slots[i].n = chunkV[i]; slots[i].nrec = 0; }
        double tb = now_sec();
        for (u64 i = 0; i < nchunk; i++) process_vertex(chunkV[i], &slots[i]);
        ST.cpu_phaseA += tb - ta; ST.cpu_phaseB += now_sec() - tb;
        // actions at hi
        if (audit_vertex_upto && hi <= audit_vertex_upto && nchunk) { if (!audit(hi, 0)) die("AUDIT FAILED at N=%llu", (unsigned long long)hi); }
        int do_audit = audit_every && hi % audit_every == 0;
        for (int i = 0; i < n_audit_at; i++) if (audit_at[i] == hi) do_audit = 1;
        for (int i = 0; i < n_dump_at; i++) if (dump_at[i] == hi) { char p[1200]; snprintf(p, sizeof p, "%s/colouring_%llu.bin", dump_dir, (unsigned long long)hi); dump_colouring(p, hi); }
        if (do_audit && !audit(hi, 1)) die("AUDIT FAILED at N=%llu", (unsigned long long)hi);
        if (prog && hi % prog == 0) progress(hi);
        if (ckpt && ckpt_every && hi % ckpt_every == 0) save_ckpt(ckpt, hi + 1);
        n = hi + 1;
    }
    u64 last = n - 1;
    if (stop_flag) { printf("STOP requested: stopping after N=%llu\n", (unsigned long long)last); if (ckpt) save_ckpt(ckpt, n); }
    pool_stop();
    if (last % (prog ? prog : 1) != 0 || !prog) progress(last);
    if (!stop_flag && ckpt) save_ckpt(ckpt, n);   // must precede fclose(LOGF): the checkpoint records the log offset
    if (LOGF) { fflush(LOGF); fclose(LOGF); LOGF = NULL; printf("replay log: %llu records\n", (unsigned long long)log_records); }
    final_report(last);
    printf("phase times: A(parallel prefetch)=%.1fs B(sequential)=%.1fs total=%.1fs\n", ST.cpu_phaseA, ST.cpu_phaseB, now_sec() - T0);
    return 0;
}

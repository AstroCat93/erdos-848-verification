// common.h -- small shared utilities (timing, RNG, allocation, fatal errors).
#ifndef COMMON_H
#define COMMON_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef unsigned __int128 u128;

static inline double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void die(const char *fmt, ...) {
    va_list ap; va_start(ap, fmt);
    fprintf(stderr, "FATAL: "); vfprintf(stderr, fmt, ap); fprintf(stderr, "\n");
    va_end(ap); fflush(stderr); exit(3);
}

static void *xmalloc(size_t n) { void *p = malloc(n ? n : 1); if (!p) die("malloc(%zu) failed", n); return p; }
static void *xcalloc(size_t n, size_t s) { void *p = calloc(n ? n : 1, s ? s : 1); if (!p) die("calloc(%zu,%zu) failed", n, s); return p; }

// splitmix64: tiny, fast, state = one u64 (easy to checkpoint)
static inline u64 splitmix64(u64 *s) {
    u64 z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
// uniform in [0, n) (n < 2^63), bias < 2^-30 for our n -- only used for heuristics
static inline u64 rng_below(u64 *s, u64 n) { return (u64)(((u128)splitmix64(s) * n) >> 64); }

#endif

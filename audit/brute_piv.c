// Audit brute force (independent of the authors' code).
// For N given, computes V_N, O, and for every pivot b in V_N (b<=N) and every x<=N the smallest and
// largest prime p with p^2 | xb+1.  Prints per pivot, per x-set, the counts needed to test Lemma 1/2/5.
// usage: brute_piv N  > out.txt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
typedef unsigned long long ull;

static ll modinv(ll a, ll m) { // a, m coprime, returns a^{-1} mod m in [0,m)
    ll g = m, x = 0, x1 = 1, a1 = a % m; if (a1 < 0) a1 += m;
    ll b1 = a1;
    // ext euclid on (m, a1)
    ll r0 = m, r1 = b1, s0 = 0, s1 = 1;
    while (r1) { ll q = r0 / r1, t = r0 - q * r1; r0 = r1; r1 = t; t = s0 - q * s1; s0 = s1; s1 = t; }
    (void)g; (void)x; (void)x1;
    if (r0 != 1) return -1;
    s0 %= m; if (s0 < 0) s0 += m; return s0;
}
static ll powmod(ll a, ll e, ll m) { __int128 r = 1, b = a % m; while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; } return (ll)r; }

#define NSETS 12
static const char *setname[NSETS] = {"C7","C18","C7o","C18o","r7_100","r43_100","r57_100","r93_100","Oany","Oodd","all","Oeven"};
static int P1s[] = {13, 30, 60, 100, 200, 1000};
#define NP1 6
static ll Ys[] = {100, 300, 1000, 3000, 10000, 30000};
#define NY 6

int main(int argc, char **argv) {
    ll N = atoll(argv[1]);
    ll PM = 2 * N + 10; // primes needed: p^2 <= N*N+1  => p <= N+1
    char *isc = calloc(PM + 1, 1); int *pr = malloc(sizeof(int) * (PM / 2 + 10)); int np = 0;
    for (ll i = 2; i <= PM; i++) { if (!isc[i]) { pr[np++] = (int)i; for (ll j = i * i; j <= PM; j += i) isc[j] = 1; } }
    // V: x^2+1 non-squarefree
    char *isV = calloc(N + 1, 1);
    for (int k = 0; k < np; k++) {
        ll p = pr[k]; if ((ull)p * p > (ull)N * N + 1) break;
        if (p == 2 || p % 4 == 3) continue;
        ll R;
        if (p == 5) R = 7;
        else {
            ll c = 2; while (powmod(c, (p - 1) / 2, p) != p - 1) c++;
            ll r = powmod(c, (p - 1) / 4, p);
            ll t = (r * r + 1) / p; ll inv2r = modinv((2 * r) % p, p);
            ll kk = ((-(t % p) * inv2r) % p + p) % p; R = r + kk * p;
        }
        ll p2 = p * p;
        if ((ll)(((__int128)R * R + 1) % p2) != 0) { fprintf(stderr, "hensel fail %lld\n", p); return 1; }
        for (ll x = R % p2; x <= N; x += p2) if (x) isV[x] = 1;
        for (ll x = (p2 - R % p2) % p2; x <= N; x += p2) if (x) isV[x] = 1;
    }
    // sanity: brute check isV for x <= 3000
    for (ll x = 1; x <= N && x <= 3000; x++) {
        ll v = x * x + 1; int nsf = 0; for (ll d = 2; d * d <= v; d++) if (v % (d * d) == 0) { nsf = 1; break; }
        if (nsf != isV[x]) { fprintf(stderr, "isV mismatch %lld\n", x); return 1; }
    }
    unsigned short *setmask = calloc(N + 1, sizeof(unsigned short));
    for (ll x = 1; x <= N; x++) {
        unsigned m = 0; ll r25 = x % 25, r100 = x % 100, r50 = x % 50;
        if (r25 == 7) m |= 1 << 0;
        if (r25 == 18) m |= 1 << 1;
        if (r50 == 7) m |= 1 << 2;
        if (r50 == 43) m |= 1 << 3;
        if (r100 == 7) m |= 1 << 4;
        if (r100 == 43) m |= 1 << 5;
        if (r100 == 57) m |= 1 << 6;
        if (r100 == 93) m |= 1 << 7;
        int inO = isV[x] && r25 != 7 && r25 != 18;
        if (inO) m |= 1 << 8;
        if (inO && (x & 1)) m |= 1 << 9;
        m |= 1 << 10;
        if (inO && !(x & 1)) m |= 1 << 11;
        setmask[x] = (unsigned short)m;
    }
    int *sp = malloc(sizeof(int) * (N + 1)), *lp = malloc(sizeof(int) * (N + 1));
    printf("# N=%lld sets:", N); for (int s = 0; s < NSETS; s++) printf(" %s", setname[s]); printf("\n");
    printf("# per line: b  then per set: tot  small[P1=13,30,60,100,200,1000]  large[same]; then Ycounts[100,300,1000,3000,10000,30000]\n");
    for (ll b = 1; b <= N; b++) {
        if (!isV[b]) continue;
        memset(sp, 0, sizeof(int) * (N + 1)); memset(lp, 0, sizeof(int) * (N + 1));
        ull lim = (ull)N * b + 1;
        for (int k = 0; k < np; k++) {
            ll p = pr[k]; ull p2 = (ull)p * p; if (p2 > lim) break;
            if (b % p == 0) continue;
            ll inv = modinv(b % (ll)p2, (ll)p2);
            ll x0 = (ll)p2 - inv; if (x0 == (ll)p2) x0 = 0; // x0*b = -1 mod p2
            if (x0 == 0) x0 = p2;
            for (ll x = x0; x <= N; x += p2) { if (!sp[x]) sp[x] = (int)p; lp[x] = (int)p; }
        }
        // verify a few
        printf("%lld", b);
        ll tot[NSETS] = {0}, sm[NSETS][NP1], lg[NSETS][NP1]; memset(sm, 0, sizeof sm); memset(lg, 0, sizeof lg);
        ll ycnt[NY] = {0};
        for (ll x = 1; x <= N; x++) {
            if (!sp[x]) continue;
            for (int j = 0; j < NY; j++) if (lp[x] > Ys[j]) ycnt[j]++;
            unsigned m = setmask[x];
            for (int s = 0; s < NSETS; s++) if (m >> s & 1) {
                tot[s]++;
                for (int i = 0; i < NP1; i++) { if (sp[x] <= P1s[i]) sm[s][i]++; if (lp[x] > P1s[i]) lg[s][i]++; }
            }
        }
        for (int s = 0; s < NSETS; s++) {
            printf(" %lld", tot[s]); for (int i = 0; i < NP1; i++) printf(" %lld", sm[s][i]); for (int i = 0; i < NP1; i++) printf(" %lld", lg[s][i]);
        }
        for (int j = 0; j < NY; j++) printf(" %lld", ycnt[j]);
        printf("\n");
        // spot-check squarefreeness for a few x by trial division
        if (b % 997 == 0) {
            for (ll x = 1; x <= N; x += 4999) {
                ull v = (ull)x * b + 1; int nsf = 0;
                for (ull d = 2; d * d <= v; d++) if (v % (d * d) == 0) { nsf = 1; break; }
                if (nsf != (sp[x] != 0)) { fprintf(stderr, "sqf mismatch b=%lld x=%lld\n", b, x); return 1; }
            }
        }
    }
    return 0;
}

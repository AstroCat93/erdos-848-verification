"""
Independent re-implementation (audit) of the Part-A bounds, written from notes/proof.md, NOT from rigorous.py.
Differences by design:
  * statuses for p=2,5,q are derived by brute-force enumeration over the allowed pivot residues (mod 4, mod 25),
    not hand-coded;
  * the matched-q gain is the exact integer knapsack over sets M with prod q <= sqrt(N^2+1) (resp. /5),
    maximised jointly with the pivot's residues mod 4 and mod 25;
  * own prime sieve, own Y / Yq search, Pell sum evaluated term by term (float64, with 1e-9 relative slack).
"""
import math, sys
from fractions import Fraction as Fr
from functools import lru_cache
import numpy as np

LIM = 130_000_000
_s = np.ones(LIM + 1, dtype=np.bool_); _s[:2] = False
for i in range(2, int(LIM ** 0.5) + 1):
    if _s[i]: _s[i * i::i] = False
PRIMES = np.nonzero(_s)[0].astype(np.int64)
del _s
def pi(x):
    assert x <= LIM
    return int(np.searchsorted(PRIMES, x, side='right'))
SP = [int(p) for p in PRIMES[:2000]]

P0, P1, QS, QSPLIT = 13, 1000, 400, 60
QLIST = [q for q in SP if q % 4 == 1 and 13 <= q <= QS]

def cdiv(a, b): return -(-a // b)

def sqrt_m1(q):
    """the two roots of r^2 = -1 mod q^2"""
    q2 = q * q
    # brute force mod q then Hensel
    r = next(r for r in range(1, q) if (r * r + 1) % q == 0)
    t = (r * r + 1) // q
    k = (-t * pow(2 * r, -1, q)) % q
    R = r + k * q
    assert (R * R + 1) % q2 == 0
    return R % q2, (-R) % q2

def frac_hits(p, L, r, bres, bmod):
    """for x = r (mod L), pivot b = bres (mod bmod) where p^2 | bmod or p^2 | gcd: fraction of the class with
    p^2 | x b + 1.  Requires bmod multiple of p^2 (so b mod p^2 is determined)."""
    p2 = p * p
    assert bmod % p2 == 0
    b = bres % p2
    if b % p == 0: return Fr(0)
    d = math.gcd(L, p2)
    xs = [x for x in range(p2) if x % d == r % d]
    hit = sum(1 for x in xs if (x * b + 1) % p2 == 0)
    return Fr(hit, len(xs))

def class_bound(N, L, fr):
    """Lemma 1.  fr: dict p -> fraction in [0,1] (0 never, 1 always, else sub g) for p | special; other p<=P1 -> 1/p^2."""
    full = cdiv(N, L)
    small, mid = [], []
    for p in SP:
        if p > P1: break
        g = fr.get(p)
        if g is None:
            assert L % p != 0, (p, L)
            g = Fr(1, p * p)
        if g == 1: return Fr(full)
        if g == 0: continue
        (small if p <= P0 else mid).append(g)
    prod = Fr(1)
    for g in small: prod *= 1 - g
    b = Fr(N, L) * (1 - prod) + (2 ** len(small) - 1)
    for g in mid: b += cdiv(N * g.numerator, L * g.denominator)
    return min(Fr(full), b)

# pivot residue sets
def piv_res(kind, parity):
    r25 = {'C7': [7], 'C18': [18], 'O': [s for s in range(25) if s not in (7, 18)]}[kind]
    r4 = [0, 2] if parity == 'even' else [1, 3]
    return r25, r4

def C_class(N, r, L, kind, parity, r4only=None):
    """worst case over pivots of the given kind/parity of Lemma-1 bound for class r mod L (L | 100)."""
    r25s, r4s = piv_res(kind, parity)
    if r4only is not None: r4s = [r4only]
    best = Fr(0)
    for b25 in r25s:
        for b4 in r4s:
            fr = {2: frac_hits(2, L, r, b4, 4), 5: frac_hits(5, L, r, b25, 25)}
            best = max(best, class_bound(N, L, fr))
    return best

def rho_max(N, parity, avoid5):
    def kmax(X):
        k, prod = 0, 1
        for p in SP[1:]:
            if avoid5 and p == 5: continue
            if prod * p > X: return k
            prod *= p; k += 1
    if parity == 'odd': return 2 ** kmax(N)
    return max(c * 2 ** kmax(N // 2 ** v) for v, c in ((1, 1), (2, 2), (3, 4)) if 2 ** v <= N)

def glob(N, parity, avoid5):
    rho = rho_max(N, parity, avoid5)
    sq = math.isqrt(N)
    best = None
    ps_all = PRIMES[PRIMES > P1]
    for e in np.linspace(0.5, 0.97, 300):
        Y = int(N ** e)
        if Y <= max(P1, sq) or Y > LIM: continue
        ps = ps_all[ps_all <= sq]
        s = int(((N + ps * ps - 1) // (ps * ps)).sum())
        s += pi(Y) - pi(max(sq, P1))
        v = s + rho * (Fr(2 * (N + 1), Y) + Fr(N * N + 1, Y * Y))
        if best is None or v < best[0]: best = (v, Y)
    return best

def pell(N, Yq):
    M = (N * N + 1) // (Yq * Yq)
    if M < 2: return 0.0
    m = np.arange(M + 1)
    adm = np.ones(M + 1, dtype=bool); adm[:2] = False; adm[0::4] = False
    for p in PRIMES[(PRIMES <= M) & (PRIMES % 4 == 3)]:
        adm[p::p] = False
    ms = m[adm].astype(np.float64)
    t = 0.5 + math.log(2 * N + 1) / (2 * np.log1p(np.sqrt(ms)))
    return float(t.sum()) * (1 + 1e-9) + 1e-6

def O_tail(N, div):
    best = None
    sq = math.isqrt(N)
    qs_all = PRIMES[(PRIMES > QS) & (PRIMES % 4 == 1)]
    for e in np.linspace(0.55, 0.95, 120):
        Yq = int(N ** e)
        if Yq <= QS or Yq > LIM: continue
        M = (N * N + 1) // (Yq * Yq)
        if M > 3 * 10 ** 7: continue
        qs = qs_all[qs_all <= Yq]
        s = int((2 * ((N + div * qs * qs - 1) // (div * qs * qs))).sum())
        v = s + pell(N, Yq)
        if best is None or v < best: best = v
    return best

def knap_exact(items, cap):
    """items: list of (q, gain); maximize sum gain over subsets with prod q <= cap (exact DFS)."""
    items = sorted(items)
    best = [Fr(0)]
    def dfs(i, prod, val):
        if val > best[0]: best[0] = val
        for j in range(i, len(items)):
            q, g = items[j]
            if prod * q > cap: break  # sorted ascending
            dfs(j + 1, prod * q, val + g)
    dfs(0, 1, Fr(0))
    return best[0]

def G(N, kind, parity, xpar):
    """Lemma 5, worst case over pivots (jointly over b mod 25, b mod 4), excluding the global term."""
    div = 1 if xpar == 'any' else 2
    r25s, r4s = piv_res(kind, parity)
    cap = math.isqrt(N * N + 1) if kind == 'O' else math.isqrt((N * N + 1) // 25)
    roots = {q: sqrt_m1(q) for q in QLIST}
    best = None
    # x-class residues mod 4 (relevant for p=2): parity class
    for b25 in r25s:
        for b4 in r4s:
            base, items = Fr(0), []
            for q in QLIST:
                Lq = q * q * div
                per_eps = []
                for eps_r in roots[q]:
                    # class x = eps_r (q^2), x of parity xpar
                    if div == 2:
                        xr = eps_r if (eps_r % 2 == 1) else eps_r + q * q   # odd representative mod 2q^2
                        assert xr % 2 == 1
                    else:
                        xr = eps_r
                    if q <= QSPLIT:
                        one = Fr(0)
                        for s in range(25):
                            if s in (7, 18): continue
                            L = Lq * 25
                            # representative of x = xr (Lq), x = s (25)
                            xx = next(x for x in range(xr, L, Lq) if x % 25 == s)
                            fr = {2: frac_hits(2, L, xx, b4, 4), 5: frac_hits(5, L, xx, b25, 25), q: Fr(0)}
                            one += class_bound(N, L, fr)
                        full = 23 * cdiv(N, Lq * 25)
                    else:
                        fr = {2: frac_hits(2, Lq, xr, b4, 4), 5: Fr(1, 25), q: Fr(0)}
                        # p=5 with L coprime to 5: worst case 1/25 (5 not | b)
                        one = class_bound(N, Lq, fr)
                        full = cdiv(N, Lq)
                    per_eps.append((one, full))
                base += per_eps[0][0] + per_eps[1][0]
                gain = max(per_eps[0][1] - per_eps[0][0], per_eps[1][1] - per_eps[1][0])
                items.append((q, gain))
            v = base + knap_exact(items, cap)
            if best is None or v > best: best = v
    return best + Fr(O_tail(N, div)).limit_denominator(10 ** 6) + 1

def O_crude(N, div):
    return sum(2 * 23 * cdiv(N, 25 * div * q * q) for q in QLIST) + Fr(O_tail(N, div)).limit_denominator(10 ** 6) + 1

def cases(N):
    gOe = glob(N, 'even', False)[0]; gOo = glob(N, 'odd', False)[0]
    gCe = glob(N, 'even', True)[0]; gCo = glob(N, 'odd', True)[0]
    res = {}
    # II.a
    res['II.a'] = C_class(N, 7, 25, 'O', 'even') + C_class(N, 18, 25, 'O', 'even') + G(N, 'O', 'even', 'any') + gOe
    # II.b.2b
    Oo = O_crude(N, 2)
    GC7 = G(N, 'C7', 'even', 'odd'); GOo = G(N, 'O', 'odd', 'odd')
    res['II.b.2b'] = C_class(N, 43, 50, 'C7', 'even') + C_class(N, 7, 25, 'O', 'odd') + gCe + gOo + min(Oo, GC7, GOo)
    parts = dict(gOe=gOe, gOo=gOo, gCe=gCe, gCo=gCo, Oo=Oo, GC7=GC7, GOo=GOo,
                 C18o_by_C7e=C_class(N, 43, 50, 'C7', 'even'), C7_by_Oo=C_class(N, 7, 25, 'O', 'odd'),
                 C7_by_Oe=C_class(N, 7, 25, 'O', 'even'), GOe_any=G(N, 'O', 'even', 'any'))
    return res, parts

if __name__ == '__main__':
    for a in sys.argv[1:]:
        N = int(float(a))
        res, parts = cases(N)
        tgt = Fr(N - 6, 25) + 1
        print(f"N={N}")
        for k, v in parts.items(): print(f"   {k:14s} /N = {float(v)/N:.6f}")
        for k, v in res.items(): print(f"   {k:8s} UB/N = {float(v)/N:.6f}  target/N = {float(tgt)/N:.6f}  {'ok' if v < tgt else 'FAIL'}")

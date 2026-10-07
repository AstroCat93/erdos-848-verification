"""
Rigorous explicit bound for Erdos #848 for all N >= N1.

Everything returned as an "upper bound" is a true upper bound:
 * exact integers / Fractions for all counting main terms and ceilings,
 * prime counts exact (sieve) up to PLIM, Rosser-Schoenfeld above,
 * sums over primes beyond the table bounded by integral/tail bounds,
 * logarithms only in the Pell tail, with explicit outward rounding slack.

See ../notes/proof.md for the mathematics (lemma numbers referenced below).
"""
from fractions import Fraction as Fr
import math, sys, itertools
import numpy as np

PLIM = 2 * 10**8   # exact prime table bound

def _sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    s[4::2] = False
    for i in range(3, int(n**0.5) + 1, 2):
        if s[i]: s[i*i::2*i] = False
    return np.nonzero(s)[0].astype(np.int64)

PR = _sieve(PLIM)
SMALLP = [int(p) for p in PR[:20000]]
Q1 = [q for q in SMALLP if q % 4 == 1 and q >= 13]

def primepi_ub(x):
    """upper bound for pi(x); exact if x <= PLIM."""
    x = int(x)
    if x < 2: return 0
    if x <= PLIM:
        return int(np.searchsorted(PR, x, side='right'))
    # Rosser & Schoenfeld (1962), (3.6): pi(x) < 1.25506 x / ln x for x > 1.
    lx = math.log(x) * (1 - 1e-12)
    return math.ceil(1.25506 * x / lx) + 1

# ---- sum_{p > P} 1/p^2, rigorous upper bound -----------------------------------
_inv_sq_cache = {}
def inv_sq_tail(P, residue4=None):
    """upper bound on sum over primes p > P (optionally p = 1 mod 4) of 1/p^2."""
    key = (P, residue4)
    if key in _inv_sq_cache: return _inv_sq_cache[key]
    arr = PR[PR > P]
    if residue4 is not None: arr = arr[arr % 4 == residue4]
    a = arr.astype(np.float64)
    s = math.fsum((1.0 / (a * a)).tolist())   # correctly rounded sum of rounded terms
    s = s * (1 + 1e-12)                      # covers per-term rounding (rel. err <= 2^-52 each)
    # tail beyond PLIM: sum over odd n > PLIM of 1/n^2 <= 1/(2*(PLIM-1))
    s += 1.0 / (2 * (PLIM - 1))
    r = Fr(s).limit_denominator(10**18) + Fr(1, 10**17)
    _inv_sq_cache[key] = r
    return r

# ---- rho* -----------------------------------------------------------------------
def _max_omega_odd(X, avoid5=False):
    k, prod = 0, 1
    for p in SMALLP[1:]:
        if avoid5 and p == 5: continue
        if prod * p > X: break
        prod *= p; k += 1
    return k

def rho_star(N, parity, avoid5=False):
    """max over b<=N with given parity (and 5 not | b if avoid5) of 2^{omega_odd(b)} * c2(v2(b)),
    c2 = 1 (v2<=1), 2 (v2=2), 4 (v2>=3).  This bounds #{x mod b : m x^2 = 1 (mod b)} for every m (Lemma 2)."""
    if parity == 'odd':
        return 2 ** _max_omega_odd(N, avoid5)
    best = 0
    for v, c2 in ((1, 1), (2, 2), (3, 4)):
        if 2**v <= N:
            best = max(best, c2 * 2 ** _max_omega_odd(N // 2**v, avoid5))
    return best

def E_switch(N, Y, rho):
    """Lemma 2: #{a<=N : exists prime p>Y, p^2 | ab+1} <= rho*(2(N+1)/Y + (N^2+1)/Y^2)."""
    return rho * (Fr(2 * (N + 1), Y) + Fr(N * N + 1, Y * Y))

def sum_ceil_N_over_p2(N, lo, hi):
    """upper bound on sum_{lo < p <= hi, p prime} ceil(N/p^2)."""
    if hi <= lo: return 0
    sq = math.isqrt(N)
    if min(sq, hi) <= PLIM:
        ps = PR[(PR > lo) & (PR <= min(sq, hi))]
        s = 0
        if len(ps):
            ps_obj = ps.astype(object)
            s = sum((N + p * p - 1) // (p * p) for p in ps_obj.tolist())
        lo2 = max(lo, sq)
        if hi > lo2:
            s += primepi_ub(hi) - (int(np.searchsorted(PR, lo2, side='right')) if lo2 <= PLIM else 0)
        return s
    # large N: ceil(x) <= x + 1
    return N * inv_sq_tail(lo) + primepi_ub(hi)

def glob(N, P1, rho, Y=None):
    """Global part for one sieving element: sum_{P1<p<=Y} ceil(N/p^2) + E(N,Y)  (bounds #{x<=N: exists p>P1, p^2|xb+1})."""
    if Y is not None:
        return sum_ceil_N_over_p2(N, P1, Y) + E_switch(N, Y, rho)
    best = None
    lo = max(P1 + 1, math.isqrt(N) + 1)
    for e in [0.5 + 0.01 * k for k in range(46)]:
        Yc = max(lo, int(round(N ** e)))
        v = sum_ceil_N_over_p2(N, P1, Yc) + E_switch(N, Yc, rho)
        if best is None or v < best[0]: best = (v, Yc)
    # local refinement
    v0, Y0 = best
    for f in [0.8, 0.9, 0.95, 1.05, 1.1, 1.25]:
        Yc = max(lo, int(Y0 * f))
        v = sum_ceil_N_over_p2(N, P1, Yc) + E_switch(N, Yc, rho)
        if v < best[0]: best = (v, Yc)
    return best

# ---- class bound (Lemma 1) ------------------------------------------------------
def class_bound(N, L, stat, P0, P1):
    """Upper bound for #{x in [N]: x = r (mod L), exists prime p <= P1 with p^2 | xb+1},
    returning the full class size ceil(N/L) if some prime is 'always'.
    stat[p] in {('always',), ('never',), ('sub', g)}; unspecified p -> ('sub', 1/p^2)."""
    full = -(-N // L)
    small, mid = [], []
    for p in SMALLP:
        if p > P1: break
        st = stat.get(p, ('sub', Fr(1, p * p)))
        if st[0] == 'always': return full
        if st[0] == 'never': continue
        (small if p <= P0 else mid).append(st[1])
    prod = Fr(1)
    for g in small: prod *= (1 - g)
    main = Fr(N, L) * (1 - prod) + (2 ** len(small) - 1)
    for g in mid:
        # #{x<=N in a single class mod L/g} <= ceil(N g / L)
        main += -(-(N * g.numerator) // (L * g.denominator))
    return min(Fr(full), main)

def stat_C(r, L, elt_parity, elt_mod4=None):
    """statuses for a class of C7 / C18 (r mod L, L in {25,50,100}) sieved by an element that is
    NOT in the same mod-25 class family as forces 25 | xb+1 (caller's responsibility): p=5 never."""
    st = {5: ('never',)}
    if elt_parity == 'even':
        st[2] = ('never',)
    else:
        if L == 100:
            st[2] = ('always',) if ((r % 4) * elt_mod4 + 1) % 4 == 0 else ('never',)
        elif L == 50:
            st[2] = ('sub', Fr(1, 2)) if r % 2 == 1 else ('never',)
        else:
            st[2] = ('sub', Fr(1, 4))
    return st

def SB(N, classes, elt_parity, P0, P1, elt_mod4=None):
    return sum(class_bound(N, L, stat_C(r, L, elt_parity, elt_mod4), P0, P1) for (r, L) in classes)

# ---- the set O (Lemma 3, Lemma 4) ----------------------------------------------
def pell_tail(N, Yq):
    """#{x<=N : q^2 | x^2+1 for some q > Yq} <= sum_{2<=m<=M, m admissible} (1/2 + ln(2N+1)/(2 ln(1+sqrt m))),
    M = floor((N^2+1)/Yq^2).  Bounded rigorously by dyadic blocks."""
    M = (N * N + 1) // (Yq * Yq)
    if M < 2: return Fr(0)
    if M > 4 * 10**7:
        # analytic: sum_{2<=m<=M} (1/2 + lnN2/(2 ln(1+sqrt m))) with no admissibility restriction;
        # split at m <= K = isqrt(M): terms <= 1/2 + lnN2/(2 ln(1+sqrt 2)); m > K: <= 1/2 + lnN2/(2 ln(1+K^(1/2)))
        lnN = math.log(2 * N + 1) * (1 + 1e-12)
        K = math.isqrt(M)
        t1 = K * (0.5 + lnN / (2 * math.log1p(math.sqrt(2)) * (1 - 1e-12)))
        t2 = M * (0.5 + lnN / (2 * math.log1p(math.sqrt(K)) * (1 - 1e-12)))
        return Fr(t1 + t2).limit_denominator(10**6) + 1
    adm = np.ones(M + 1, dtype=bool); adm[:2] = False
    adm[4::4] = False
    for p in PR[PR <= M]:
        if p % 4 == 3: adm[p::p] = False
    lnN = math.log(2 * N + 1) * (1 + 1e-12)
    total = Fr(0)
    lo = 2
    while lo <= M:
        hi = min(M, 2 * lo - 1)
        cnt = int(adm[lo:hi + 1].sum())
        if cnt:
            denom = math.log1p(math.sqrt(lo)) * (1 - 1e-12)   # lower bound for ln(1+sqrt m), m >= lo
            per = 0.5 + lnN / (2 * denom)
            total += Fr(cnt) * (Fr(per).limit_denominator(10**12) + Fr(1, 10**9))
        lo = hi + 1
    return total

def O_tail(N, Qs, parity, Yq=None):
    """upper bound #{x<=N of the given parity: exists prime q = 1 mod 4, q > Qs, q^2 | x^2+1}."""
    div = 1 if parity == 'any' else 2
    if Yq is None:
        best = None
        for e in [0.55 + 0.01 * k for k in range(40)]:
            Yc = max(Qs + 1, int(round(N ** e)))
            v = O_tail(N, Qs, parity, Yc)
            if v is None: continue
            if best is None or v < best: best = v
        return best
    pt = pell_tail(N, Yq)
    if pt is None: return None
    sq = math.isqrt(N)
    s = Fr(0)
    if min(sq, Yq) <= PLIM:
        qs = PR[(PR > Qs) & (PR <= min(sq, Yq)) & (PR % 4 == 1)].astype(object).tolist()
        s += sum(2 * (-(-N // (div * q * q))) for q in qs)
        lo2 = max(Qs, sq)
        if Yq > lo2:
            s += 2 * (primepi_ub(Yq) - (int(np.searchsorted(PR, lo2, side='right')) if lo2 <= PLIM else 0))
    else:
        s += Fr(2 * N, div) * inv_sq_tail(Qs, 1) + 2 * primepi_ub(Yq)
    return s + pt

def O_crude(N, parity, Qs):
    """|O cap [N] cap parity| <= sum_{q<=Qs} 2*23*ceil(N/(25*div*q^2)) + O_tail."""
    div = 1 if parity == 'any' else 2
    s = sum(2 * 23 * (-(-N // (25 * div * q * q))) for q in Q1 if q <= Qs)
    return s + O_tail(N, Qs, parity)

def fractional_knapsack_ub(items, cap):
    """items: list of (gain>=0 Fraction, weight>0 float); returns an upper bound on max sum of gains of a subset with
    total weight <= cap.  Uses LP duality: for ANY lam >= 0,  OPT <= lam*cap + sum_i max(0, g_i - lam*w_i).
    lam is chosen heuristically (critical ratio of the greedy order) but the bound is valid for every lam.
    All arithmetic exact (weights/cap converted to Fractions; they are already rounded outward by the caller)."""
    its = [(Fr(g), Fr(w)) for g, w in items if g > 0]
    if not its: return Fr(0)
    capF = Fr(cap)
    its.sort(key=lambda z: z[0] / z[1], reverse=True)
    w, lam = Fr(0), Fr(0)
    for g, wt in its:
        if w + wt <= capF:
            w += wt
        else:
            lam = g / wt
            break
    return lam * capF + sum(max(Fr(0), g - lam * wt) for g, wt in its)

def G_bound(N, elt_kind, elt_parity, x_parity, P0, P1, Qs, Qsplit):
    """Upper bound (WITHOUT the element's global term) for #{x in O cap [N], x of x_parity, xb+1 non-squarefree},
    worst case over sieving elements b of the given kind ('O','C7','C18') and parity.  (Lemma 5)"""
    div = 1 if x_parity == 'any' else 2
    base, items = Fr(0), []
    for q in Q1:
        if q > Qs: break
        Lq = q * q * div
        st = {q: ('never',)}
        if x_parity == 'even' or elt_parity == 'even':
            st[2] = ('never',)
        else:
            st[2] = ('sub', Fr(1, 2) if x_parity == 'odd' else Fr(1, 4))
        if q <= Qsplit:
            one = Fr(0)
            for s25 in range(25):
                if s25 in (7, 18): continue
                st5 = dict(st)
                # 25 | xb+1 for all x in the class iff s25 = -b^{-1} (mod 25).  For b in O this is one class of R25
                # (worst case 5 not | b); for b in C7 it is s25=7, for b in C18 it is s25=18: excluded.
                st5[5] = ('always',) if (elt_kind == 'O' and s25 == 1) else ('never',)   # representative class s0=-b^{-1} (any unit class; bound is independent of r)
                one += class_bound(N, Lq * 25, st5, P0, P1)
            full = 23 * (-(-N // (25 * Lq)))
        else:
            st5 = dict(st); st5[5] = ('sub', Fr(1, 25))
            one = class_bound(N, Lq, st5, P0, P1)
            full = -(-N // Lq)
        base += 2 * one
        items.append((Fr(full) - one, math.log(q) * (1 - 1e-12)))
    # matched primes: b = eps r_q (mod q^2) makes that class fully counted; prod of matched q^2 divides b^2+1 <= N^2+1
    # (and 25 | b^2+1 additionally if b in C7 or C18).
    cap = math.log(N * N + 1) / 2 * (1 + 1e-12)
    if elt_kind in ('C7', 'C18'): cap -= math.log(5) * (1 - 1e-12)
    gain = fractional_knapsack_ub(items, cap)
    tail = O_tail(N, Qs, x_parity)
    return base + gain + tail

# ---- the cases ------------------------------------------------------------------
C7_25 = [(7, 25)]; C18_25 = [(18, 25)]
C7_o50 = [(7, 50)]; C18_o50 = [(43, 50)]

def all_cases(N, P0=13, P1=1000, Qs=400, Qsplit=60, verbose=False):
    t_plus_1_lb = Fr(N - 6, 25) + 1          # |A| >= t(N)+1 > (N-6)/25 + ... ; need UB < (N-6)/25 + 1
    gl = {}
    for kind, par, a5 in (('O', 'even', False), ('O', 'odd', False), ('C', 'even', True), ('C', 'odd', True)):
        gl[(kind, par)] = glob(N, P1, rho_star(N, par, a5))[0]
    c100 = -(-N // 100)
    R = {}
    # pair bound for (r1, r2) odd classes mod 100 of C7 / C18 that are mutually incompatible mod 4 and 25
    def pairbound(r1, r2):
        m4_1, m4_2 = r1 % 4, r2 % 4
        sieved = SB(N, [(r1, 100)], 'odd', P0, P1, m4_2) + SB(N, [(r2, 100)], 'odd', P0, P1, m4_1) \
                 + gl[('C', 'odd')] * 2
        return max(Fr(c100), sieved)
    # Case I (A* empty)
    R['I.1'] = SB(N, C7_25, 'even', P0, P1) + SB(N, C18_25, 'even', P0, P1) + 2 * gl[('C', 'even')]
    R['I.2'] = SB(N, C18_o50, 'even', P0, P1) + SB(N, C7_25, 'odd', P0, P1) + gl[('C', 'even')] + gl[('C', 'odd')]
    R['I.3'] = SB(N, C7_o50, 'even', P0, P1) + SB(N, C18_25, 'odd', P0, P1) + gl[('C', 'even')] + gl[('C', 'odd')]
    R['I.4'] = pairbound(7, 43) + pairbound(57, 93)
    # Case II.a (A* has an even element b)
    R['II.a'] = SB(N, C7_25, 'even', P0, P1) + SB(N, C18_25, 'even', P0, P1) \
                + G_bound(N, 'O', 'even', 'any', P0, P1, Qs, Qsplit) + gl[('O', 'even')]
    # Case II.b (A* nonempty, all odd; b in A* odd)
    Oodd = O_crude(N, 'odd', Qs)
    G_O_odd = G_bound(N, 'O', 'odd', 'odd', P0, P1, Qs, Qsplit)
    G_C7_even = G_bound(N, 'C7', 'even', 'odd', P0, P1, Qs, Qsplit)
    G_C18_even = G_bound(N, 'C18', 'even', 'odd', P0, P1, Qs, Qsplit)
    worst = Fr(0)
    for m4 in (1, 3):
        if m4 == 1: free, siev = (7, 43), [(57, 100), (93, 100)]
        else:       free, siev = (57, 93), [(7, 100), (43, 100)]
        v = pairbound(*free) + SB(N, siev, 'odd', P0, P1, m4) + gl[('O', 'odd')] + min(Oodd, G_O_odd)
        worst = max(worst, v)
    R['II.b.1'] = worst
    R['II.b.2a'] = SB(N, C18_25, 'even', P0, P1) + SB(N, C7_25, 'even', P0, P1) + 2 * gl[('C', 'even')] \
                   + min(Oodd, G_C7_even, G_C18_even)
    R['II.b.2b'] = SB(N, C18_o50, 'even', P0, P1) + SB(N, C7_25, 'odd', P0, P1) + gl[('C', 'even')] + gl[('O', 'odd')] \
                   + min(Oodd, G_C7_even, G_O_odd)
    R['II.b.3b'] = SB(N, C7_o50, 'even', P0, P1) + SB(N, C18_25, 'odd', P0, P1) + gl[('C', 'even')] + gl[('O', 'odd')] \
                   + min(Oodd, G_C18_even, G_O_odd)
    margins = {k: t_plus_1_lb - v for k, v in R.items()}
    if verbose:
        print(f"N = {N:.4e}:  globals/N: " + ", ".join(f"{k[0]}{k[1][0]}={float(v)/N:.5f}" for k, v in gl.items()))
        for k in R:
            print(f"   {k:8s} UB/N = {float(R[k])/N:.6f}   margin/N = {float(margins[k])/N:+.6f}  {'ok' if margins[k] > 0 else 'FAIL'}")
    return R, margins

if __name__ == '__main__':
    Ns = [int(float(a)) for a in sys.argv[1:]] or [10**8, 2 * 10**8, 10**9]
    for N in Ns:
        all_cases(N, verbose=True)

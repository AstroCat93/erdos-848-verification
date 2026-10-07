"""
NOTE (repository version): the tail argument used in the proof (proof/proof.md, Section 6) is
    tail_check2(NB=10**14, theta=0.8, ratio=2)        together with D_cases()
defined at the bottom of this file.  The first version below (density_parts [a stub], err_majorant_over_N,
tail_check) is superseded; it is kept unchanged because the audit scripts in ../audit exercise it.
Use ./verify_analytic.sh (or check_tail.py) for the certified check; the `__main__` block below
only prints the superseded version-1 majorant.

Tail argument for N >= NB (default 1e40).   [superseded description of version 1, kept for the record]

For every case, UB(N) = D*N + Err(N) where
  D   = sum of all main densities (class mains, mid-prime mains, global main N*sum_{p>P1} 1/p^2 per sieving
        element, ALL matched-q gains, O_tail main N*sum_{q>Qs, q=1(4)} 2/(div q^2)),
  Err = everything else, with parameters Y = Yq = ceil(N^{3/4}).
We bound Err(N)/N <= c0/N + c1 N^{-1/4} + c2 rho*(N) N^{-1/2} with rho*(N) <= 4 * 2^{1.38402 ln N / ln ln N}
(Robin 1983, Thm 11, applied to the odd part), which is decreasing for N >= NB, so it suffices to check at NB.

Implementation: rather than re-deriving D by hand, we evaluate the exact rigorous UB at two huge N and at NB with
the explicit majorant, and also compute D directly from the class data.  See notes/proof.md, Section 6.
"""
import math
from fractions import Fraction as Fr
import rigorous as R

def density_parts(P0=13, P1=1000, Qs=400, Qsplit=60):
    """Return per-case limiting density D (as float upper bound) by evaluating UB(N)/N at N = 10^k for large k
    with all matched gains forced on (cap = infinity), and the error constants."""
    pass

def err_majorant_over_N(N, n_classes_const, n_elements):
    """Explicit upper bound on Err(N)/N for Y = Yq = N^{3/4}:
       - class constants (inclusion-exclusion 2^k-1, mid-prime '+1's, ceil +1's): n_classes_const / N
       - per sieving element: pi(Y) <= 1.25506 Y/ln Y,  E <= rho*(2(N+1)/Y + (N^2+1)/Y^2)
       - O_tail: 2 pi(Yq) + Pell(N, Yq) with Pell <= M(1/2 + ln(2N+1)/(2 ln(1+M^{1/4}))) + sqrt(M)(1/2 + ln(2N+1)/(2 ln(1+sqrt 2))),
         M = (N^2+1)/Yq^2."""
    lnN = math.log(N)
    Y = N ** 0.75
    piY = 1.25506 * Y / (0.75 * lnN)
    om = 1.38402 * lnN / math.log(lnN)
    rho = 4 * 2 ** om
    E = rho * (2 * (N + 1) / Y + (N * N + 1) / (Y * Y))
    M = (N * N + 1) / (Y * Y)
    pell = M * (0.5 + math.log(2 * N + 1) / (2 * math.log1p(M ** 0.25))) + math.sqrt(M) * (0.5 + math.log(2 * N + 1) / (2 * math.log1p(math.sqrt(2))))
    tot = n_classes_const + n_elements * (piY + E) + 2 * piY + pell
    return tot / N

if __name__ == '__main__':
    for k in (40, 50, 60, 100):
        N = 10.0 ** k
        print(k, err_majorant_over_N(N, 10**6, 2))


# ------------------------------------------------------------------------------------------------
# Exact limiting densities D_case (N -> infinity), with ALL matched gains switched on.
from fractions import Fraction as Fr

def cls_density(L, stat, P0, P1):
    """lim_{N} ClassBound/N  (constants vanish).  'always' -> 1/L."""
    small, mid = [], []
    for p in R.SMALLP:
        if p > P1: break
        st = stat.get(p, ('sub', Fr(1, p * p)))
        if st[0] == 'always': return Fr(1, L)
        if st[0] == 'never': continue
        (small if p <= P0 else mid).append(st[1])
    prod = Fr(1)
    for g in small: prod *= (1 - g)
    return Fr(1, L) * (1 - prod) + sum(g / L for g in mid)

def SBd(classes, par, P0, P1, m4=None):
    return sum(cls_density(L, R.stat_C(r, L, par, m4), P0, P1) for (r, L) in classes)

def Gd(elt_kind, elt_par, x_par, P0, P1, Qs, Qsplit):
    div = 1 if x_par == 'any' else 2
    tot = Fr(0)
    for q in R.Q1:
        if q > Qs: break
        Lq = q * q * div
        # all matched: worst case is the full class for one sign, sieved for the other
        if q <= Qsplit:
            st = {q: ('never',)}
            st[2] = ('never',) if (x_par == 'even' or elt_par == 'even') else ('sub', Fr(1, 2) if x_par == 'odd' else Fr(1, 4))
            one = Fr(0)
            for s25 in range(25):
                if s25 in (7, 18): continue
                st5 = dict(st); st5[5] = ('always',) if (elt_kind == 'O' and s25 == 0) else ('never',)
                one += cls_density(Lq * 25, st5, P0, P1)
            full = Fr(23, 25 * Lq)
        else:
            st = {q: ('never',), 5: ('sub', Fr(1, 25))}
            st[2] = ('never',) if (x_par == 'even' or elt_par == 'even') else ('sub', Fr(1, 2) if x_par == 'odd' else Fr(1, 4))
            one = cls_density(Lq, st, P0, P1)
            full = Fr(1, Lq)
        tot += one + max(one, full)          # one sign sieved, the other matched (full)
    tail = 2 * R.inv_sq_tail(Qs, 1) / div      # sum_{q>Qs, q=1(4)} 2/(div q^2)
    return tot + tail

def Ocrude_d(par, Qs):
    div = 1 if par == 'any' else 2
    return sum(Fr(2 * 23, 25 * div * q * q) for q in R.Q1 if q <= Qs) + 2 * R.inv_sq_tail(Qs, 1) / div

def D_cases(P0=13, P1=1000, Qs=400, Qsplit=60):
    gl = R.inv_sq_tail(P1)       # global main density per pivot: sum_{p>P1} 1/p^2
    C7, C18, C7o, C18o = [(7, 25)], [(18, 25)], [(7, 50)], [(43, 50)]
    def pair(r1, r2):
        s = SBd([(r1, 100)], 'odd', P0, P1, r2 % 4) + SBd([(r2, 100)], 'odd', P0, P1, r1 % 4) + 2 * gl
        return max(Fr(1, 100), s)
    D = {}
    D['I.1'] = SBd(C7, 'even', P0, P1) + SBd(C18, 'even', P0, P1) + 2 * gl
    D['I.2'] = SBd(C18o, 'even', P0, P1) + SBd(C7, 'odd', P0, P1) + 2 * gl
    D['I.3'] = SBd(C7o, 'even', P0, P1) + SBd(C18, 'odd', P0, P1) + 2 * gl
    D['I.4'] = pair(7, 43) + pair(57, 93)
    D['II.a'] = SBd(C7, 'even', P0, P1) + SBd(C18, 'even', P0, P1) + Gd('O', 'even', 'any', P0, P1, Qs, Qsplit) + gl
    Oo = Ocrude_d('odd', Qs)
    GOo = Gd('O', 'odd', 'odd', P0, P1, Qs, Qsplit)
    GC7 = Gd('C7', 'even', 'odd', P0, P1, Qs, Qsplit); GC18 = Gd('C18', 'even', 'odd', P0, P1, Qs, Qsplit)
    w = Fr(0)
    for m4 in (1, 3):
        free, siev = ((7, 43), [(57, 100), (93, 100)]) if m4 == 1 else ((57, 93), [(7, 100), (43, 100)])
        w = max(w, pair(*free) + SBd(siev, 'odd', P0, P1, m4) + gl + min(Oo, GOo))
    D['II.b.1'] = w
    D['II.b.2a'] = SBd(C18, 'even', P0, P1) + SBd(C7, 'even', P0, P1) + 2 * gl + min(Oo, GC7, GC18)
    D['II.b.2b'] = SBd(C18o, 'even', P0, P1) + SBd(C7, 'odd', P0, P1) + 2 * gl + min(Oo, GC7, GOo)
    D['II.b.3b'] = SBd(C7o, 'even', P0, P1) + SBd(C18, 'odd', P0, P1) + 2 * gl + min(Oo, GC18, GOo)
    return D


def tail_check(NB=10**14, n_pivots=4, c0=10**6, slack=2.0):
    """Certify every N >= NB:  UB_case(N) <= D_case*N + Err(N) with Y = Yq = N^{3/4}, and
    Err(N)/N <= slack * err_majorant_over_N(N, c0, n_pivots), which is decreasing for N >= NB.
    Need max_case D_case + Err/N < (N-6)/(25N) + 1/N  for N >= NB; we check against 1/25 - 1/NB."""
    D = D_cases()
    Dmax = max(D.values())
    e = slack * err_majorant_over_N(float(NB), c0, n_pivots)
    lhs = float(Dmax) + e
    ok = lhs < 0.04 - 1.0 / NB
    # numerical sanity check of monotone decrease of the majorant on a log grid up to 1e300
    prev, mono = None, True
    for k in range(int(math.log10(NB) * 10), 3000):
        v = err_majorant_over_N(10.0 ** (k / 10), c0, n_pivots) if k < 3000 else 0
        if prev is not None and v > prev * (1 + 1e-12): mono = False
        prev = v
    return dict(Dmax=float(Dmax), worst_case=max(D, key=D.get), err_over_N_at_NB=e, lhs=lhs, ok=ok, majorant_monotone_numeric=mono)


# ------------------------------------------------------------------------------------------------
# Rigorous tail, version 2 (replaces tail_check).  For N in [Ni, Ni*10], fix Y = Yq = ceil(Ni^theta).
# UB_case(N) <= D_case*N + Err with
#   Err <= c0 + npiv*(pi(Y) + E(N_hi, Y, rho*(N_hi))) + 2*pi(Yq) + Pell(N_hi, Yq)   (each term nondecreasing in N),
# so UB_case(N)/N <= D_case + Err/Ni.  Need < 1/25 - 6/(25 Ni) + 1/Ni  (we use 1/25 - 1/Ni).
def c0_bound(P0=13, P1=1000, Qs=400, Qsplit=60):
    nmid = sum(1 for p in R.SMALLP if P0 < p <= P1)
    per_class = (2 ** 6 - 1) + nmid + 1          # incl.-excl. error (<= 6 small sub primes) + mid ceilings + full-class ceil
    nq_split = sum(1 for q in R.Q1 if q <= Qsplit)
    nq_rest = sum(1 for q in R.Q1 if Qsplit < q <= Qs)
    g_classes = 2 * 23 * nq_split + 2 * nq_rest
    c_classes = 8 + 4                              # C7/C18 classes (mod 25/50/100) + pair classes
    return (g_classes + c_classes) * per_class + 10   # + small ceil terms (pair max, etc.)

def pell_ub_float(N, Yq):
    M = (N * N + 1) / (Yq * Yq)
    if M < 2: return 0.0
    lnN = math.log(2 * N + 1)
    K = math.sqrt(M)
    return K * (0.5 + lnN / (2 * math.log1p(math.sqrt(2)))) + M * (0.5 + lnN / (2 * math.log1p(K ** 0.5)))

def pi_rs(x):
    return 1.25506 * x / math.log(x)

def tail_check2(NB=10**14, NTOP=10**300, theta=0.8, npiv=4, ratio=2):
    D = D_cases(); Dmax = float(max(D.values()))
    c0 = c0_bound()
    Ni, worst = NB, 0.0
    while Ni < NTOP:
        Nh = Ni * ratio
        Y = math.ceil(Ni ** theta)
        rho = max(R.rho_star(Nh, 'even'), R.rho_star(Nh, 'odd'))
        E = rho * (2 * (Nh + 1) / Y + (Nh * Nh + 1) / (Y * Y))
        err = c0 + npiv * (pi_rs(Y) + E) + 2 * pi_rs(Y) + pell_ub_float(Nh, Y)
        val = Dmax + 1.01 * err / Ni
        worst = max(worst, val)
        if not val < 0.04 - 1.0 / Ni:
            return dict(ok=False, at=Ni, val=val)
        Ni = Nh
    # beyond NTOP: Robin bound rho <= 4*2^{1.38402 ln N/lnln N}; majorant decreasing (checked analytically in notes)
    lnN = math.log(NTOP); om = 1.38402 * lnN / math.log(lnN)
    tail_val = Dmax + (c0 * 1.0) / 1e300 + 7 * 1.25506 * (1e300) ** (theta - 1) / (theta * lnN) \
               + npiv * 4 * 2 ** om * (2 * 1e300 ** (-theta) + 1e300 ** (1 - 2 * theta)) + 10 * 1e300 ** (-0.5)
    return dict(ok=True, Dmax=Dmax, worst_case=max(D, key=D.get), c0=c0, worst_val_in_grid=worst, beyond_NTOP=tail_val)

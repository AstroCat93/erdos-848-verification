import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analytic'))
import rigorous as R
from sympy import factorint
from math import gcd
def rho_brute(b):
    best = 0
    sq = {}
    for x in range(b):
        sq.setdefault(x * x % b, 0); sq[x * x % b] += 1
    # m x^2 = 1  <=> x^2 = m^{-1} for units m
    for u, c in sq.items():
        if gcd(u, b) == 1: best = max(best, c)
    return best if b > 1 else 1
def bound(b):
    f = factorint(b); v2 = f.get(2, 0); w = sum(1 for p in f if p != 2)
    return 2**w * (1 if v2 <= 1 else (2 if v2 == 2 else 4))
bad = [b for b in range(1, 3001) if rho_brute(b) > bound(b)]
eq = sum(1 for b in range(1, 3001) if rho_brute(b) == bound(b))
print('rho*(b) <= 2^w c2 violations for b<=3000:', bad[:10], ' equality in', eq, 'cases')
# max over b <= N by parity / avoid5 vs rho_star
for N in [10, 100, 1000, 5000, 30030, 100000, 510510*2, 10**7]:
    for par in ('odd', 'even'):
        for a5 in (False, True):
            if N <= 100000:
                cands = [b for b in range(1, N + 1) if (b % 2 == (1 if par == 'odd' else 0)) and (not a5 or b % 5)]
                m = max(bound(b) for b in cands)
            else:
                m = None
            r = R.rho_star(N, par, a5)
            if m is not None and m != r: print('MISMATCH', N, par, a5, m, r)
    print(N, {(p, a): R.rho_star(N, p, a) for p in ('odd', 'even') for a in (False, True)})

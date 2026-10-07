import sys, random
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analytic'))
import numpy as np
import rigorous as R
import indep as I
N = int(float(sys.argv[1])); Qs = int(sys.argv[2]); Qsplit = int(sys.argv[3])
P0 = 13
QL = [q for q in R.Q1 if q <= Qs]
roots = {q: I.sqrt_m1(q) for q in QL}
SPR = [int(p) for p in R.SMALLP if p <= 1000]
def classes(xpar):
    out = []
    for q in QL:
        q2 = q * q
        for r in roots[q]:
            xs = np.arange(r if r else q2, N + 1, q2, dtype=np.int64)
            if xpar == 'odd': xs = xs[xs % 2 == 1]
            if q <= Qsplit:
                m = xs % 25; xs = xs[(m != 7) & (m != 18)]
            out.append(xs)
    return np.concatenate(out)
X = {'odd': classes('odd'), 'any': classes('any')}
# V membership for pivots (b^2+1 non-squarefree) via sieve over q
isV = np.zeros(N + 1, dtype=bool)
for q in [5] + [q for q in R.SMALLP if q % 4 == 1 and q > 5 and q * q <= N * N + 1 and q <= N]:
    if q == 5: rs = (7, 18)
    elif q in roots: rs = roots[q]
    else: rs = I.sqrt_m1(q)
    for r in rs:
        isV[r if r else q * q::q * q] = True
def nmatched(b): return sum(1 for q in QL if (b * b + 1) % (q * q) == 0)
def T(b, xpar, P1):
    xs = X[xpar]; v = xs * b + 1
    hit = np.zeros(len(xs), dtype=bool)
    for p in SPR:
        if p > P1: break
        hit |= (v % (p * p) == 0)
    return int(hit.sum())
random.seed(1)
for kind, par, xpar in (('O', 'odd', 'odd'), ('O', 'even', 'any'), ('C7', 'even', 'odd'), ('C18', 'even', 'odd')):
    r25 = {'C7': (7,), 'C18': (18,)}.get(kind)
    allb = [b for b in range(1, N + 1) if isV[b] and (b % 2 == (1 if par == 'odd' else 0))
            and ((b % 25 in r25) if r25 else (b % 25 not in (7, 18)))]
    multi = [b for b in allb if nmatched(b) >= 2]
    m13 = [b for b in allb if (b * b + 1) % 169 == 0]
    piv = sorted(set(multi + random.sample(m13, min(400, len(m13))) + random.sample(allb, min(400, len(allb)))))
    for P1 in (13, 50, 200):
        ub = R.G_bound(N, kind, par, xpar, P0, P1, Qs, Qsplit) - R.O_tail(N, Qs, xpar)
        vals = [(T(b, xpar, P1), b) for b in piv]
        t, bw = max(vals)
        print(f"{kind:3s} {par:4s} x={xpar:3s} P1={P1:4d}: pivots={len(piv):5d} (multi-matched {len(multi)})  max truth={t} at b={bw} (matched {nmatched(bw)})  bound={float(ub):.1f}  ratio={t/float(ub):.3f}  {'OK' if t <= ub else 'VIOLATION'}")

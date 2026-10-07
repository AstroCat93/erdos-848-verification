import sys, math
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analytic'))
import rigorous as R
from fractions import Fraction as Fr
from sympy import factorint

SETS = "C7 C18 C7o C18o r7_100 r43_100 r57_100 r93_100 Oany Oodd all Oeven".split()
P1s = [13, 30, 60, 100, 200, 1000]
Ys = [100, 300, 1000, 3000, 10000, 30000]

def load(fn):
    rows = {}
    for line in open(fn):
        if line.startswith('#'): continue
        v = list(map(int, line.split()))
        b = v[0]; k = 1; d = {}
        for s in SETS:
            tot = v[k]; sm = v[k+1:k+7]; lg = v[k+7:k+13]; k += 13
            d[s] = (tot, sm, lg)
        d['Y'] = v[k:k+6]
        rows[b] = d
    return rows

def kind(b):
    r = b % 25
    return 'C7' if r == 7 else ('C18' if r == 18 else 'O')

def rho_b(b):
    f = factorint(b); v2 = f.get(2, 0); w = sum(1 for p in f if p != 2)
    return 2**w * (1 if v2 <= 1 else (2 if v2 == 2 else 4))

def check(N, fn, Qs=400, Qsplit=60):
    rows = load(fn)
    P0 = 13
    viol = []
    # map: (pivot predicate, set, bound fn(P1))
    tests = []
    def add(name, pred, setn, bfun, pkind, ppar, a5):
        tests.append((name, pred, setn, bfun, pkind, ppar, a5))
    ev = lambda b: b % 2 == 0
    od = lambda b: b % 2 == 1
    # C classes sieved by pivots (Lemma 1 part)
    add('C7 by O even', lambda b: kind(b)=='O' and ev(b), 'C7', lambda P1: R.SB(N, R.C7_25, 'even', P0, P1), 'O','even',False)
    add('C18 by O even', lambda b: kind(b)=='O' and ev(b), 'C18', lambda P1: R.SB(N, R.C18_25, 'even', P0, P1), 'O','even',False)
    add('C7 by O odd', lambda b: kind(b)=='O' and od(b), 'C7', lambda P1: R.SB(N, R.C7_25, 'odd', P0, P1), 'O','odd',False)
    add('C18 by O odd', lambda b: kind(b)=='O' and od(b), 'C18', lambda P1: R.SB(N, R.C18_25, 'odd', P0, P1), 'O','odd',False)
    add('C7 by C18 even', lambda b: kind(b)=='C18' and ev(b), 'C7', lambda P1: R.SB(N, R.C7_25, 'even', P0, P1), 'C','even',True)
    add('C18 by C7 even', lambda b: kind(b)=='C7' and ev(b), 'C18', lambda P1: R.SB(N, R.C18_25, 'even', P0, P1), 'C','even',True)
    add('C18o by C7 even', lambda b: kind(b)=='C7' and ev(b), 'C18o', lambda P1: R.SB(N, R.C18_o50, 'even', P0, P1), 'C','even',True)
    add('C7o by C18 even', lambda b: kind(b)=='C18' and ev(b), 'C7o', lambda P1: R.SB(N, R.C7_o50, 'even', P0, P1), 'C','even',True)
    add('C7 by C18 odd', lambda b: kind(b)=='C18' and od(b), 'C7', lambda P1: R.SB(N, R.C7_25, 'odd', P0, P1), 'C','odd',True)
    add('C18 by C7 odd', lambda b: kind(b)=='C7' and od(b), 'C18', lambda P1: R.SB(N, R.C18_25, 'odd', P0, P1), 'C','odd',True)
    # mod-100 pair and II.b.1 classes
    add('r7 by r43', lambda b: b%100==43, 'r7_100', lambda P1: R.SB(N, [(7,100)], 'odd', P0, P1, 3), 'C','odd',True)
    add('r43 by r7', lambda b: b%100==7, 'r43_100', lambda P1: R.SB(N, [(43,100)], 'odd', P0, P1, 3), 'C','odd',True)
    add('r57 by r93', lambda b: b%100==93, 'r57_100', lambda P1: R.SB(N, [(57,100)], 'odd', P0, P1, 1), 'C','odd',True)
    add('r93 by r57', lambda b: b%100==57, 'r93_100', lambda P1: R.SB(N, [(93,100)], 'odd', P0, P1, 1), 'C','odd',True)
    add('r57 by O 1mod4', lambda b: kind(b)=='O' and b%4==1, 'r57_100', lambda P1: R.SB(N, [(57,100)], 'odd', P0, P1, 1), 'O','odd',False)
    add('r93 by O 1mod4', lambda b: kind(b)=='O' and b%4==1, 'r93_100', lambda P1: R.SB(N, [(93,100)], 'odd', P0, P1, 1), 'O','odd',False)
    add('r7 by O 3mod4', lambda b: kind(b)=='O' and b%4==3, 'r7_100', lambda P1: R.SB(N, [(7,100)], 'odd', P0, P1, 3), 'O','odd',False)
    add('r43 by O 3mod4', lambda b: kind(b)=='O' and b%4==3, 'r43_100', lambda P1: R.SB(N, [(43,100)], 'odd', P0, P1, 3), 'O','odd',False)
    # O sieved (Lemma 5)
    add('O any by O even', lambda b: kind(b)=='O' and ev(b), 'Oany', lambda P1: R.G_bound(N, 'O', 'even', 'any', P0, P1, Qs, Qsplit), 'O','even',False)
    add('O odd by O odd', lambda b: kind(b)=='O' and od(b), 'Oodd', lambda P1: R.G_bound(N, 'O', 'odd', 'odd', P0, P1, Qs, Qsplit), 'O','odd',False)
    add('O odd by C7 even', lambda b: kind(b)=='C7' and ev(b), 'Oodd', lambda P1: R.G_bound(N, 'C7', 'even', 'odd', P0, P1, Qs, Qsplit), 'C','even',True)
    add('O odd by C18 even', lambda b: kind(b)=='C18' and ev(b), 'Oodd', lambda P1: R.G_bound(N, 'C18', 'even', 'odd', P0, P1, Qs, Qsplit), 'C','even',True)
    out = []
    for (name, pred, setn, bfun, pk, pp, a5) in tests:
        piv = [b for b in rows if pred(b)]
        if not piv: continue
        rho = R.rho_star(N, pp, a5)
        line = [f"{name:20s} n={len(piv):5d}"]
        for i, P1 in enumerate(P1s):
            tsm = max(rows[b][setn][1][i] for b in piv)
            ub = bfun(P1)
            g = R.glob(N, P1, rho)[0]
            tlg = max(rows[b]['all'][2][i] for b in piv)
            ttot = max(rows[b][setn][0] for b in piv)
            ok1 = tsm <= ub; ok2 = tlg <= g; ok3 = ttot <= ub + g
            if not (ok1 and ok2 and ok3): viol.append((name, P1, tsm, float(ub), tlg, float(g), ttot))
            line.append(f"P1={P1}: {tsm}/{float(ub):.0f} gl {tlg}/{float(g):.0f}")
        out.append("  ".join(line))
    for l in out: print(l)
    # Lemma 2 direct: #{x: lp > Y} <= E(N,Y,rho(b))
    bad = 0; worst = 0
    for b, d in rows.items():
        rb = rho_b(b)
        for j, Y in enumerate(Ys):
            if Y < 2: continue
            e = R.E_switch(N, Y, rb)
            if d['Y'][j] > e: bad += 1; print('Lemma2 VIOL', b, Y, d['Y'][j], float(e))
            if e > 0: worst = max(worst, d['Y'][j] / float(e))
    print(f"Lemma 2 per-pivot check: violations={bad}, max truth/bound={worst:.4f}")
    print("VIOLATIONS:", viol)

if __name__ == '__main__':
    N = int(sys.argv[1]); fn = sys.argv[2]
    Qs = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    Qsplit = int(sys.argv[4]) if len(sys.argv) > 4 else 60
    check(N, fn, Qs, Qsplit)

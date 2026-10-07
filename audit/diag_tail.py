import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analytic'))
import rigorous as R
import indep as I
from fractions import Fraction as Fr
N = int(float(sys.argv[1]))
for par, div in (('odd', 2), ('any', 1)):
    tR = R.O_tail(N, 400, par); tI = I.O_tail(N, div)
    print(par, 'O_tail rig/N', float(tR)/N, ' indep/N', tI/N)
    # best Yq for rigorous & pieces
    best=None
    for e in [0.55 + 0.01*k for k in range(40)]:
        Yc = max(401, int(round(N**e)))
        pt = R.pell_tail(N, Yc); tot = R.O_tail(N, 400, par, Yc)
        if best is None or tot < best[0]: best=(tot, Yc, pt)
    print('   rig best Yq', best[1], 'pell/N', float(best[2])/N, 'pell-part M', (N*N+1)//best[1]**2)
for kind, par, xp, div in (('O','odd','odd',2), ('C7','even','odd',2), ('O','even','any',1)):
    gR = R.G_bound(N, kind, par, xp, 13, 1000, 400, 60) - R.O_tail(N, 400, xp)
    gI = I.G(N, kind, par, xp) - (Fr(I.O_tail(N, div)).limit_denominator(10**6) + 1)
    print(kind, par, xp, 'base+gain rig/N', float(gR)/N, ' indep(exact knapsack)/N', float(gI)/N)

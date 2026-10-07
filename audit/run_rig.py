import sys, time
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analytic'))
t0=time.time()
import rigorous as R
print('import', time.time()-t0)
from fractions import Fraction as Fr
for N in [42840000, 10**9]:
    t0=time.time()
    Rb, m = R.all_cases(N, verbose=True)
    print('time', time.time()-t0)
    # components
    P0,P1,Qs,Qsplit=13,1000,400,60
    gl = {}
    for kind, par, a5 in (('O','even',False),('O','odd',False),('C','even',True),('C','odd',True)):
        v,Y = R.glob(N, P1, R.rho_star(N, par, a5)); gl[(kind,par)]=(v,Y,R.rho_star(N,par,a5))
        print(kind,par,'rho*',R.rho_star(N,par,a5),'Y',Y,'gl/N',float(v)/N)
    print('SB C7 even /N', float(R.SB(N,R.C7_25,'even',P0,P1))/N)
    print('SB C7 odd /N', float(R.SB(N,R.C7_25,'odd',P0,P1))/N)
    print('SB C18o50 even /N', float(R.SB(N,R.C18_o50,'even',P0,P1))/N)
    for args in [('O','even','any'),('O','odd','odd'),('C7','even','odd'),('C18','even','odd')]:
        print('G',args, float(R.G_bound(N,*args,P0,P1,Qs,Qsplit))/N)
    print('Oodd crude /N', float(R.O_crude(N,'odd',Qs))/N)
    print('Otail odd /N', float(R.O_tail(N,Qs,'odd'))/N, 'any', float(R.O_tail(N,Qs,'any'))/N)

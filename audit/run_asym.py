import sys, time
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analytic'))
import asymptotic as A
import rigorous as R
print('density_parts returns:', A.density_parts())
D = A.D_cases()
for k, v in D.items(): print(f"  D[{k}] = {float(v):.6f}")
print('max D', max(float(v) for v in D.values()))
for k in (40, 50, 60, 100):
    N = 10.0**k
    print(k, 'err/N (as called in __main__)', A.err_majorant_over_N(N, 10**6, 2), ' with 4 elements', A.err_majorant_over_N(N, 10**6, 4))
for N in (10**12, 10**16, 10**20, 10**40, 10**100):
    t=time.time(); Rb,_ = R.all_cases(N); w=max(Rb,key=lambda c:Rb[c]); print(f"N=1e{len(str(N))-1}: worst {w} UB/N={float(Rb[w])/N:.6f}  ({time.time()-t:.0f}s)", {k: round(float(v)/N,5) for k,v in Rb.items()})

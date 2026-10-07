import sys, math
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analytic'))
import asymptotic as A
D = A.D_cases()
nel = {'I.1': 2, 'I.2': 2, 'I.3': 2, 'I.4': 4, 'II.a': 1, 'II.b.1': 3, 'II.b.2a': 2, 'II.b.2b': 2, 'II.b.3b': 2}
for k in (12, 13, 14, 15, 16, 17, 20, 40):
    N = 10.0 ** k
    worst = max(float(D[c]) + A.err_majorant_over_N(N, 10**6, nel[c]) for c in D)
    wc = max(D, key=lambda c: float(D[c]) + A.err_majorant_over_N(N, 10**6, nel[c]))
    print(f"N_B=1e{k}: max_case D + Err/N = {worst:.5f} (case {wc})  need < 0.04")

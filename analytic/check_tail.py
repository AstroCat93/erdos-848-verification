"""
Re-run the certified tail check of proof/proof.md, Section 6 (all N >= 10^14).

  * D_cases(): exact limiting densities of every case (all matched gains switched on);
  * tail_check2(NB=1e14, theta=0.8, ratio=2): doubling grid N_i -> 2 N_i up to 10^300, checking
    D_max + 1.01 Err/N_i < 1/25 - 1/N_i on each [N_i, 2 N_i];
  * the closed-form majorant beyond 10^300 (Robin's bound for omega).

Exit status 0 iff every check passes.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import asymptotic as A  # noqa: E402

t0 = time.time()
D = A.D_cases()
print("Limiting densities D_case (exact, printed as floats):")
for k in D:
    print(f"   {k:8s} {float(D[k]):.7f}")
wc = max(D, key=D.get)
print(f"max D_case = {float(D[wc]):.7f}  (case {wc})")

res = A.tail_check2(NB=10**14, theta=0.8, ratio=2)
print("tail_check2(NB=1e14, theta=0.8, ratio=2):", res)
ok = bool(res.get("ok")) and res["worst_val_in_grid"] < 0.04 and res["beyond_NTOP"] < 0.04
print(f"worst value of D_max + 1.01*Err/N_i on [1e14, 1e300]: {res.get('worst_val_in_grid', float('nan')):.5f}  (need < 0.04 - 1/N_i)")
print(f"value beyond 1e300: {res.get('beyond_NTOP', float('nan')):.7f}  (need < 0.04)")
print(f"TAIL {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.1f}s)")
sys.exit(0 if ok else 1)

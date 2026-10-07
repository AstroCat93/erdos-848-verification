"""
Interval certification (proof/proof.md, Section 5, "Interval grid").

For N_lo <= N <= N_hi: a valid A in [N] is also valid in [N_hi]; if A is not contained in C7 or C18, the case
analysis run AT N_hi gives |A| <= UB(N_hi).  Hence  max_case UB(N_hi) < (N_lo - 6)/25 + 1 <= t(N) + 1
certifies |A| <= t(N) for every N in [N_lo, N_hi].  No monotonicity of UB in N is used.
(An earlier version of this docstring appealed to monotonicity; that argument was replaced after the audit.
The code below is unchanged: it only ever evaluates UB at the right endpoint N_hi.)

Usage: python3 grid.py N1 NEND logfile      (appends to logfile)
"""
import sys, time, math
from fractions import Fraction as Fr
import rigorous as R

def main():
    N1 = int(float(sys.argv[1])); NEND = int(float(sys.argv[2])); log = open(sys.argv[3], 'a')
    lo = N1; ratio = 1.02; k = 0
    t0 = time.time()
    log.write(f"# grid start N1={N1} NEND={NEND:.3e}\n"); log.flush()
    while lo < NEND:
        hi = int(lo * ratio)
        Rb, _ = R.all_cases(hi)
        worst_case = max(Rb, key=lambda c: Rb[c])
        ub = Rb[worst_case]
        target = Fr(lo - 6, 25) + 1
        if ub < target:
            k += 1
            slack = float((target - ub) / lo)
            log.write(f"OK [{lo:.6e}, {hi:.6e}] worst={worst_case} UB/N_hi={float(ub)/hi:.6f} slack/N_lo={slack:.6f} t={time.time()-t0:.0f}s\n"); log.flush()
            lo = hi
            # adapt: aim to use ~half the slack
            ratio = max(1.002, min(1.5, 1 + 0.5 * slack / (float(ub) / hi + 1e-9)))
        else:
            log.write(f"retry [{lo:.6e}, {hi:.6e}] worst={worst_case} UB={float(ub):.6e} target={float(target):.6e}\n"); log.flush()
            ratio = 1 + (ratio - 1) / 2
            if ratio < 1.0005:
                log.write(f"FAIL at lo={lo}\n"); log.flush(); return
    log.write(f"# DONE: certified all N in [{N1}, {lo}] in {k} intervals, {time.time()-t0:.0f}s\n"); log.flush()

if __name__ == '__main__':
    main()

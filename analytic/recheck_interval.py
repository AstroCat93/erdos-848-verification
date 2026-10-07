"""
Re-certify grid intervals independently of grid.py's chaining.

usage: python3 recheck_interval.py LO1 HI1 [LO2 HI2 ...]

For each [LO, HI] it evaluates every case bound at HI with rigorous.all_cases (exact Fractions / integers)
and checks   max_case UB(HI) < (LO - 6)/25 + 1,   which certifies |A| <= t(N) for all N in [LO, HI]
(proof/proof.md, Section 5, "Interval grid").  The endpoints printed in the grid logs are rounded to
7 significant digits, so callers should round LO down and HI up (a superset interval is harder, not easier).

Exit status 0 iff every interval is certified.
"""
import os
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rigorous as R  # noqa: E402


def main(argv):
    vals = [int(float(a)) if ("e" in a or "E" in a) else int(a) for a in argv]
    if len(vals) < 2 or len(vals) % 2:
        print(__doc__)
        return 2
    all_ok = True
    for lo, hi in zip(vals[::2], vals[1::2]):
        t0 = time.time()
        Rb, _ = R.all_cases(hi)
        worst = max(Rb, key=lambda c: Rb[c])
        ub = Rb[worst]
        target = Fr(lo - 6, 25) + 1
        ok = ub < target
        all_ok &= ok
        print(f"[{lo}, {hi}]  worst={worst}  UB/N_hi={float(ub) / hi:.6f}  slack/N_lo={float((target - ub) / lo):.6f}  "
              f"{'CERTIFIED' if ok else 'NOT CERTIFIED'}  ({time.time() - t0:.0f}s)", flush=True)
    print("INTERVALS", "PASS" if all_ok else "FAIL")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

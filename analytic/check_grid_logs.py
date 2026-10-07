"""
Consistency check of the certified grid logs grid2/log_{A,B,C}.txt (proof/proof.md, Section 5).

Checks, for each segment: every interval line is 'OK' (no 'retry'/'FAIL'), consecutive intervals are
contiguous (printed hi of one == printed lo of the next), the '# DONE' line is present; and across segments:
the first interval starts at <= 4.2e7, each segment starts no later than the previous one ends, and the last
one ends at >= 1e14 (where the tail of Section 6 takes over).

This is a bookkeeping check only.  To re-certify intervals, use recheck_interval.py (a sample) or
`verify_analytic.sh --full` (re-runs grid.py over all three segments and compares with these logs).
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OK_RE = re.compile(r"^OK \[([0-9.e+]+), ([0-9.e+]+)\] worst=(\S+) UB/N_hi=([0-9.]+) slack/N_lo=([0-9.]+)")
DONE_RE = re.compile(r"^# DONE: certified all N in \[(\d+), (\d+)\] in (\d+) intervals")


def check(path):
    n_ok, prev_hi, done, errors, min_slack = 0, None, None, [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("# grid start"):
                continue
            m = OK_RE.match(line)
            if m:
                lo, hi = m.group(1), m.group(2)
                if prev_hi is not None and lo != prev_hi:
                    errors.append(f"gap: previous hi {prev_hi} != lo {lo}")
                prev_hi = hi
                n_ok += 1
                s = float(m.group(5))
                min_slack = s if min_slack is None else min(min_slack, s)
                continue
            m = DONE_RE.match(line)
            if m:
                done = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
                continue
            errors.append(f"unexpected line: {line}")
    if done is None:
        errors.append("no DONE line")
    elif done[2] != n_ok:
        errors.append(f"DONE reports {done[2]} intervals, found {n_ok}")
    return n_ok, done, min_slack, errors


def main():
    segs = [os.path.join(HERE, "grid2", f"log_{s}.txt") for s in "ABC"]
    all_ok, prev_end = True, None
    for p in segs:
        n, done, min_slack, errs = check(p)
        name = os.path.relpath(p, HERE)
        if done:
            print(f"{name}: {n} intervals, certified [{done[0]}, {done[1]}], min slack/N_lo = {min_slack}, errors = {len(errs)}")
        for e in errs:
            print("   ERROR", e)
        all_ok &= not errs and done is not None
        if done:
            if prev_end is None and done[0] > 42_000_000:
                print("   ERROR first segment starts above 4.2e7"); all_ok = False
            if prev_end is not None and done[0] > prev_end:
                print(f"   ERROR segment starts at {done[0]} > previous end {prev_end}"); all_ok = False
            prev_end = done[1]
    if prev_end is None or prev_end < 10**14:
        print("ERROR grid does not reach 1e14"); all_ok = False
    print("GRID LOGS", "PASS" if all_ok else "FAIL", f"(coverage [42000000, {prev_end}])")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())

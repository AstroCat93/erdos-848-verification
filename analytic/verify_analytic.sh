#!/usr/bin/env bash
# Quick re-verification of Part A (analytic, N >= 4.2e7).  See proof/proof.md, Sections 5-6.
#
#   ./verify_analytic.sh          quick mode (a few minutes):
#                                   1. bookkeeping check of the certified grid logs grid2/log_{A,B,C}.txt
#                                   2. exact limiting densities + certified tail for N >= 1e14 (tail_check2)
#                                   3. independent re-certification of 5 sample grid intervals, including
#                                      the tightest one [4.2e7, 4.284e7] and the last one ending at 1.14e14
#   ./verify_analytic.sh --full   additionally re-runs grid.py over all of [4.2e7, 1e14] (three segments in
#                                 parallel, ~25 min wall on an M4 Pro) and diffs the result against grid2/.
#
# Requires python3 (>= 3.9) with numpy.  Exit status 0 iff everything passes.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}

echo "== 1. grid log consistency"
$PY check_grid_logs.py

echo "== 2. tail N >= 1e14"
$PY check_tail.py

echo "== 3. sample grid intervals (endpoints rounded outward from the 7-digit log values)"
# first interval (tightest, case II.b.2b); last of segment A; first of B; last of B; last of C
$PY recheck_interval.py \
    42000000 42840000 \
    8876111000 10036750979 \
    10000000000 10200000000 \
    877991000000 1004273536264 \
    99373500000000 114170693749961

if [ "${1:-}" = "--full" ]; then
    echo "== 4. full grid re-run (writes to regrid/)"
    rm -rf regrid && mkdir -p regrid
    $PY grid.py 42000000 1e10 regrid/log_A.txt &
    $PY grid.py 1e10 1e12 regrid/log_B.txt &
    $PY grid.py 1e12 1e14 regrid/log_C.txt &
    wait
    status=0
    for s in A B C; do
        # drop the timing fields before comparing
        if diff <(sed -E 's/ t=[0-9]+s$//; s/ intervals, [0-9]+s$/ intervals/' grid2/log_$s.txt) \
                <(sed -E 's/ t=[0-9]+s$//; s/ intervals, [0-9]+s$/ intervals/' regrid/log_$s.txt) >/dev/null; then
            echo "segment $s: re-run identical to grid2/log_$s.txt (timings ignored)"
        else
            echo "segment $s: re-run DIFFERS from grid2/log_$s.txt"; status=1
        fi
    done
    [ $status -eq 0 ] || { echo "FULL GRID FAIL"; exit 1; }
    echo "FULL GRID PASS"
fi
echo "ANALYTIC CHECKS PASS"

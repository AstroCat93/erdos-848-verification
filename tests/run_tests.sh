#!/usr/bin/env bash
# Independent small-N checks and certifier unit tests (~1 min on an M4 Pro).  Exit status 0 iff all pass.
#
#  1. clique_check.c : exact max-clique induction, |A| <= t(N) for every N <= 1e5 (no colouring involved)
#  2. brute.py       : pure-Python Bron-Kerbosch, omega(G_N) == t(N) at a few N
#  3. verifier unit tests: test_sqf (squarefree certifier vs reference), test_v (two sieves for V agree)
#  4. checker replay of the included runs/log_1e5.bin, and brute_cmp (edge-sieve graph, no squarefree test)
#  5. negative controls: a truncated log and a log with a planted non-squarefree pair must both be REJECTED
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
VER="$HERE/../verifier"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
CC=${CC:-clang}
PY=${PYTHON:-python3}
expect() { grep -q -- "$1" "$2" || { echo "FAILED: expected '$1' in $2"; cat "$2"; exit 1; }; }

echo "== 1. clique_check (all N <= 100000)"
$CC -O2 -o "$TMP/clique_check" "$HERE/clique_check.c"
"$TMP/clique_check" 100000 2>/dev/null | tee "$TMP/cc.txt"
expect "violations=0" "$TMP/cc.txt"

echo "== 2. brute.py (omega(G_N) vs t(N))"
for n in 100 500 1000 2000 3000; do
    out=$($PY "$HERE/brute.py" $n)
    echo "   $out" | cut -c1-60
    om=$(echo "$out" | awk '{print $6}'); t=$(echo "$out" | awk '{print $8}')
    [ "$om" = "$t" ] || { echo "FAILED: omega $om != t(N) $t at N=$n"; exit 1; }
done

echo "== 3. verifier unit tests"
make -s -C "$VER"
"$VER/test_sqf" 100000 10 > "$TMP/sqf.txt"; tail -1 "$TMP/sqf.txt"
expect "TOTAL FAILURES: 0" "$TMP/sqf.txt"
"$VER/test_v" 1e8 1e6 200000 10 > "$TMP/v.txt"; tail -1 "$TMP/v.txt"
expect "TEST_V PASS" "$TMP/v.txt"

echo "== 4. replay of runs/log_1e5.bin"
"$VER/checker" replay "$VER/runs/log_1e5.bin" 100000 10 | tee "$TMP/r.txt"
expect "=> PASS" "$TMP/r.txt"
"$VER/brute_cmp" 1e5 "$VER/runs/log_1e5.bin" clique 2>/dev/null | grep -v progress > "$TMP/b.txt"; tail -1 "$TMP/b.txt"
expect "BRUTE_CMP PASS" "$TMP/b.txt"

echo "== 5. negative controls (must FAIL)"
head -c 96488 "$VER/runs/log_1e5.bin" > "$TMP/trunc.bin"
if "$VER/checker" replay "$TMP/trunc.bin" 100000 2 > "$TMP/n1.txt"; then echo "FAILED: truncated log accepted"; exit 1; fi
grep "CHECK FAIL" "$TMP/n1.txt" | head -1
# record #4 of the log is (41 -> class 2); put 41 into class 1 instead, next to 7: 41*7+1 = 288 = 2^5 * 3^2
$PY - "$VER/runs/log_1e5.bin" "$TMP/planted.bin" <<'EOF'
import struct, sys
b = bytearray(open(sys.argv[1], 'rb').read())
assert struct.unpack_from('<II', b, 32) == (41, 2)
struct.pack_into('<II', b, 32, 41, 1)
open(sys.argv[2], 'wb').write(b)
EOF
if "$VER/checker" replay "$TMP/planted.bin" 100000 2 > "$TMP/n2.txt"; then echo "FAILED: planted bad pair accepted"; exit 1; fi
grep "CHECK FAIL" "$TMP/n2.txt" | head -1

echo "ALL TESTS PASS"

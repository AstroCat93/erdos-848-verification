#!/usr/bin/env bash
# Reproduce Part B (all N <= NMAX) from scratch and verify it independently.  See proof/proof.md, Section 7.
#
#   ./reproduce.sh [quick|1e8|1e9|full] [threads]        (default: quick, 10 threads)
#
#   level   NMAX            colour     checker replay + dump    disk (log + dump)   (M4 Pro, 10 threads)
#   quick   10^7            ~3 s       ~2 s                     ~10 MB
#   1e8     10^8            ~15 s      ~11 s                    ~155 MB
#   1e9     10^9            ~5 min     ~5 min                   ~1.6 GB
#   full    2^32-1          ~25 min    ~70 min                  ~6.7 GB  (class table needs ~5.5 GB RAM)
#
# Steps: build; run `colour` (online colouring, writes the move log and a dump of the final colouring);
# replay the log with the independent `checker` (re-certifies every same-class pair with code that shares
# nothing with colour.c); check the dump; compare the SHA-256 of the move log with the published value
# (colour is deterministic and independent of the thread count, so the log is a byte-exact prefix of the
# published log_2p32_full.bin; see runs/SHA256SUMS.txt).
set -euo pipefail
cd "$(dirname "$0")"
LEVEL=${1:-quick}
THR=${2:-10}
case "$LEVEL" in
    quick) N=10000000;   EXP=5e2a72780d804b2bf92649d9eba2ecec9b61fa0efc665235f3b40c31c5293714 ;;
    1e8)   N=100000000;  EXP=210e198ebee43b647e0057ae1e0d8f7c0a0db041dcd11072644eafaaac5f5274 ;;
    1e9)   N=1000000000; EXP=65b214d314e28103c608a072620ecdef5394cc714a3db41d42e6099fc3b18618 ;;
    full)  N=4294967295; EXP=7d933b60f12c49a78b2c934cc3cb0da192acde34e18ea6bd4e92c3b7979cb7c3 ;;
    *) echo "usage: $0 [quick|1e8|1e9|full] [threads]"; exit 2 ;;
esac
sha256() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }

make -s
OUT="out_$LEVEL"
mkdir -p "$OUT"
echo "== colour --nmax $N --threads $THR"
start=$(date +%s)
./colour --nmax "$N" --threads "$THR" --log "$OUT/log.bin" --dump-at "$N" --dump-dir "$OUT" > "$OUT/colour.log" 2>&1
echo "   $(grep 'FINAL REPORT' "$OUT/colour.log")  [$(( $(date +%s) - start ))s]"
echo "== checker replay"
./checker replay "$OUT/log.bin" "$N" "$THR" 2>&1 | tee "$OUT/checker_replay.log"
echo "== checker dump"
./checker dump "$OUT/colouring_$N.bin" "$THR" 2>&1 | tee "$OUT/checker_dump.log"
grep -q '=> PASS' "$OUT/checker_replay.log" && grep -q '=> PASS' "$OUT/checker_dump.log" || { echo "PART B FAIL"; exit 1; }
got=$(sha256 "$OUT/log.bin")
echo "== move log SHA-256: $got"
if [ "$got" = "$EXP" ]; then echo "   matches the published value"; else echo "   DIFFERS from the published value $EXP"; exit 1; fi
echo "PART B PASS for all N <= $N  [total $(( $(date +%s) - start ))s]"

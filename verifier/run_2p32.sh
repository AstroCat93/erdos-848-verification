#!/bin/bash
# HISTORICAL SCRIPT, kept for the record: it resumes from runs/ckpt_1e9.bin and runs/log_1e9.bin, which are
# NOT in this repository (see runs/SHA256SUMS.txt).  For a from-scratch run use ./reproduce.sh full.
#
# Stretch run: extend the verified colouring from N=1e9 to N=2^32-1 (resume from 1e9 checkpoint),
# then independently re-verify with checker (full replay of the move log + final dump check).
# NOTE (2026-10-07): the first execution used a colour binary with a bug (final checkpoint written after the
# log was closed -> log offset 0 in ckpt_1e9.bin), so runs/log_2p32.bin only holds records for n > 1e9 and the
# "checker replay" step above failed (correctly).  The complete log is runs/log_2p32_full.bin =
# log_1e9.bin ++ log_2p32.bin, replayed independently from n=1 (runs/checker_2p32_replay_full.log).
# colour.c is now fixed (checkpoint before fclose, offset consistency checks) and the ckpt headers repaired.
cd "$(dirname "$0")"
set -o pipefail
echo "START $(date)" > runs/STATUS_2p32.txt
cp runs/ckpt_1e9.bin runs/ckpt_2p32.bin && cp runs/log_1e9.bin runs/log_2p32.bin || { echo "COPY FAILED" >> runs/STATUS_2p32.txt; exit 1; }
./colour --nmax 4294967295 --threads 10 --progress 1e7 --resume runs/ckpt_2p32.bin --ckpt runs/ckpt_2p32.bin --ckpt-every 2.5e8 \
         --log runs/log_2p32.bin --dump-at 4294967295 --dump-dir runs > runs/colour_2p32.log 2>&1
echo "colour exit $? $(date)" >> runs/STATUS_2p32.txt
./checker replay runs/log_2p32.bin 4294967295 10 > runs/checker_2p32_replay.log 2>&1
echo "checker replay exit $? $(date)" >> runs/STATUS_2p32.txt
./checker dump runs/colouring_4294967295.bin 10 > runs/checker_2p32_dump.log 2>&1
echo "checker dump exit $? $(date)" >> runs/STATUS_2p32.txt
echo "DONE $(date)" >> runs/STATUS_2p32.txt

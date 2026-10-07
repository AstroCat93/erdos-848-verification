# Part B verifier.  `make` builds every program; tested with Apple clang 21 on arm64 (macOS).
# Elsewhere: `make CC=gcc` (needs a 64-bit compiler with unsigned __int128 and pthreads).
CC = clang
CFLAGS = -O3 -march=native -Wall -Wno-unused-function -Wno-unused-but-set-variable
HDR = common.h sqf.h ref.h vsieve.h vsieve2.h
PROGS = colour checker test_sqf test_v brute_cmp bench_sqf

.PHONY: all clean

all: $(PROGS)

%: %.c $(HDR)
	$(CC) $(CFLAGS) -o $@ $< -lpthread

clean:
	rm -f $(PROGS)

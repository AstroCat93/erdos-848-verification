# Erdős Problem #848 for all N: an explicit analytic bound plus a colouring certificate

This repository contains a computer-assisted proof that answers [Erdős Problem #848](https://www.erdosproblems.com/848)
for **every** N, together with all code and logs needed to check it.

> **Status (2026-10-07).** The proof is complete but computer-assisted, and it has **not been peer reviewed**.
> It was produced with AI assistance; see the [disclosure](#ai-assistance-disclosure).
> It is **not the first resolution** of the problem. Two earlier, also unrefereed, full-resolution claims
> exist (A. C. Li; I. Pitchford / OpenAI Codex), and see [Prior work](#prior-work-and-credit) for the rest.
> This repository gives an **independent route** to the same answer, which may be easier to check.

## The problem

Let A ⊆ {1, …, N} be a set such that ab + 1 is **not squarefree** for all a, b ∈ A, including a = b.
Erdős and Sárközy asked whether |A| can exceed the number of n ≤ N with n ≡ 7 (mod 25), which is

  t(N) = #{n ≤ N : n ≡ 7 (mod 25)} = ⌊(N + 18)/25⌋.

The set {n ≤ N : n ≡ 7 (mod 25)} attains t(N): 7·7 + 1 = 50, so 25 | ab + 1 whenever a ≡ b ≡ 7 (mod 25).

## Result

**Theorem.** For every N ≥ 1, every such A satisfies |A| ≤ t(N) = ⌊(N + 18)/25⌋.

The proof (`proof/proof.md`) has two parts that overlap by a factor of about 100.

| part | range of N | method | where |
|---|---|---|---|
| **A** (analytic) | N ≥ N₁ = 4.2·10⁷ | Sawhney's case split, made fully explicit (see below). Every inequality is checked in exact rational arithmetic on a finite grid of intervals [4.2·10⁷, 1.14·10¹⁴], and a certified tail covers N ≥ 10¹⁴. | proof §1–6, `analytic/` |
| **B** (computational) | N ≤ 2³² − 1 ≈ 4.29·10⁹ | An online proper colouring of the "conflict graph" G_N with exactly t(N) colour classes, maintained for every N. Every same-class pair u, v has uv + 1 *certified* squarefree, so no valid A can have more than t(N) elements. Every N is replayed by an independent checker. | proof §7, `verifier/` |

The analytic part starts at N₁ = 4.2·10⁷: exact-rational interval checks cover [4.2·10⁷, 1.14·10¹⁴] and a
uniform argument covers N ≥ 10¹⁴ (explicit error bounds on a doubling grid up to 10³⁰⁰, evaluated in floating
point with a 1% safety factor, and a closed form beyond). (For comparison, Li's interval-wise checks start at
5·10⁶ and his uniform argument at 5·10¹¹; see the table in the paper.) Because N₁ is small, the finite part
fits on a laptop. All N ≤ 10⁸, which already exceeds N₁, take under half a minute including the independent
replay. All N ≤ 2³² − 1 take about 30 min of colouring plus about 70 min of independent checking.

Part A follows Sawhney's case split. The tools it adds are standard, and several appear in earlier explicit
work (proof §0):

1. **Switching for large primes.** Large prime squares are counted through the cofactor m in ab + 1 = mp², which
   costs about 2^ω(b)·N^{2/3} instead of π(N). The device is classical (Estermann, Hooley, Heath-Brown) and is
   used in explicit form by Sothanaphan, Li and Pitchford.
2. **Self-sieving.** The elements of A outside the classes 7, 18 (mod 25) are sieved by the same pivot elements.
3. **A mod-100 parity analysis.** It uses the fact that x ≡ 7 and y ≡ 43 (mod 100) give xy + 1 ≡ 2 (mod 4)
   and xy + 1 ≡ 2 (mod 25).

The elements outside 7, 18 (mod 25) with a large prime-square divisor of x² + 1 are bounded with the negative
Pell equation, as in Sothanaphan and Pitchford.

## Repository layout

```
proof/proof.md              the proof (Parts A and B), with lemma numbers used in the code comments
analytic/                   Part A
  rigorous.py                 all case bounds UB_case(N), with exact Fractions/integers (needs numpy)
  grid.py                     chains certified intervals [N_lo, N_hi] (the run behind grid2/)
  asymptotic.py               D_cases() (exact limiting densities) and tail_check2 (N >= 1e14)
  grid2/log_{A,B,C}.txt       the certified grid: [4.2e7,1.0037e10], [1e10,1.0043e12], [1e12,1.1417e14]
  grid_log*.txt               earlier, superseded runs (greedy knapsack; not used by the proof)
  verify_analytic.sh          quick re-verification (about 80 s); --full re-runs the whole grid
  check_grid_logs.py, check_tail.py, recheck_interval.py   helpers used by verify_analytic.sh
verifier/                   Part B (C with unsigned __int128, C11 atomics, pthreads)
  colour.c                    online colouring with t(N) classes and certified squarefreeness
  sqf.h                       squarefree certifier (trial division, Miller-Rabin, Pollard-Brent, cube-root fallback)
  vsieve.h                    sieve for V_N = {n <= N : n^2+1 not squarefree}
  checker.c                   independent replay/dump checker; shares no code with colour.c/sqf.h/vsieve.h
  ref.h, vsieve2.h            the checker's squarefree test and V sieve (ref.h is also compiled into colour.c,
                                but only for its optional internal audits, never for its certificates)
  common.h                    small shared utilities (timing, allocation, RNG); no squarefree or sieve logic
  test_sqf.c, test_v.c        unit tests of the certifier and of the two V sieves
  brute_cmp.c                 compares against an edge-sieve graph (no squarefree test) + exact clique induction
  bench_sqf.c                 certifier benchmark
  reproduce.sh                build, colour, replay, dump-check and hash-compare at 1e7 / 1e8 / 1e9 / 2^32-1
  run_2p32.sh                 the historical script of the 2^32-1 run (resumed from 1e9; see "Caveats")
  runs/                       logs of the original runs, log_1e5.bin, SHA256SUMS.txt for the large artifacts
tests/
  run_tests.sh                small-N checks, unit tests, negative controls (about 50 s)
  clique_check.c              exact max-clique induction for all N <= 1e5 (independent of the colouring)
  brute.py                    pure-Python Bron-Kerbosch max clique for one N
audit/                      adversarial audit of Part A: independent re-implementation of the two binding cases (II.a, II.b.2b), brute-force lemma tests (see audit/README.md)
notes/literature.md         literature and forum notes (as of 2026-10-07)
```

## Reproducing

**Hardware and software used** for every number below: Apple M4 Pro (12 cores: 8 performance + 4 efficiency),
24 GB RAM, macOS 27.0.1, Apple clang 21.0.0, Python 3.9.6, numpy 2.0.2 (and sympy 1.14.0 for `audit/`).
Part B used 10 threads.

**Requirements:** a C compiler with `unsigned __int128` and pthreads, `make`, and Python ≥ 3.9 with numpy
(`pip install -r requirements.txt`). Only Apple clang on arm64 was tested. The Makefile defaults to `clang`;
elsewhere set `CC` in the environment, e.g. `CC=gcc verifier/reproduce.sh quick` and `CC=gcc tests/run_tests.sh`
(both scripts pass it on to `make`). `-march=native` is set in `verifier/Makefile`.

Every command below was run on 2026-10-07 from a **fresh copy of this repository** (no binaries, no prior
outputs), except `reproduce.sh 1e9`, `reproduce.sh full` and `verify_analytic.sh --full`. The runtimes are
wall-clock times from that test; for the three exceptions they come from the logs of the original runs.

### Part A: analytic, N ≥ 4.2·10⁷

```sh
analytic/verify_analytic.sh          # about 80 s
```

This runs three checks: (1) the bookkeeping of the certified grid logs; (2) the tail N ≥ 10¹⁴; (3) an
independent exact re-certification of five grid intervals. Those five include the tightest interval
[4.2·10⁷, 4.284·10⁷] and the last interval, which ends at 1.14·10¹⁴. Expected output (abridged):

```
grid2/log_A.txt: 83 intervals, certified [42000000, 10036750979], min slack/N_lo = 0.000891, errors = 0
grid2/log_B.txt: 36 intervals, certified [10000000000, 1004273536264], min slack/N_lo = 0.006039, errors = 0
grid2/log_C.txt: 35 intervals, certified [1000000000000, 114170693749961], min slack/N_lo = 0.006363, errors = 0
GRID LOGS PASS (coverage [42000000, 114170693749961])
max D_case = 0.0281946  (case II.a)
tail_check2(NB=1e14, theta=0.8, ratio=2): {'ok': True, ..., 'c0': 78658, 'worst_val_in_grid': 0.029189943166037777, ...}
TAIL PASS  (2.3s)
[42000000, 42840000]  worst=II.b.2b  UB/N_hi=0.038342  slack/N_lo=0.000891  CERTIFIED  (7s)
[8876111000, 10036750979]  worst=II.a  UB/N_hi=0.028690  slack/N_lo=0.007558  CERTIFIED  (11s)
[10000000000, 10200000000]  worst=II.a  UB/N_hi=0.028683  slack/N_lo=0.010744  CERTIFIED  (12s)
[877991000000, 1004273536264]  worst=II.a  UB/N_hi=0.027883  slack/N_lo=0.008107  CERTIFIED  (11s)
[99373500000000, 114170693749961]  worst=II.a  UB/N_hi=0.027647  slack/N_lo=0.008236  CERTIFIED  (35s)
ANALYTIC CHECKS PASS
```

Every case bound must stay below 1/25 = 0.04 per unit of N. The binding case is II.b.2b near N₁, with
slack 0.00089·N. For large N the worst case is II.a, with limiting density 0.02819.

`analytic/verify_analytic.sh --full` also re-runs `grid.py` over all of [4.2·10⁷, 10¹⁴] and diffs the
result with `grid2/`. The three segments run in parallel. The original segment runs took 1547 s, 766 s and
1387 s, so expect about 26 min of wall time. The full mode was **not** re-run in the fresh-copy test; the logs
in `grid2/` come from the original run. A single interval can be certified with
`python3 analytic/recheck_interval.py LO HI`.

### Part B: computational, N ≤ 2³² − 1

```sh
verifier/reproduce.sh quick      # all N <= 1e7:   2.6 s
verifier/reproduce.sh 1e8        # all N <= 1e8:   24 s  (already covers N_1 = 4.2e7)
verifier/reproduce.sh 1e9        # all N <= 1e9:    about 10 min (estimate from the original logs)
verifier/reproduce.sh full       # all N <= 2^32-1: about 30 min colour + 70 min checker (estimate);
                                 #   about 6.7 GB of disk; the class table alone needs 5.5 GB of RAM
```

Each level does the following:

1. Builds the programs with `make`.
2. Runs `colour`, which writes the move log and a dump of the final colouring.
3. Replays the whole log with the independent `checker`. For every n, the checker verifies that the state is
   a proper colouring of G_n with t(n) classes. It re-certifies every pair that becomes same-class with its
   own squarefree test.
4. Checks the dump.
5. Compares the SHA-256 of the move log with the published value.

Expected output of `quick`:

```
   ==== FINAL REPORT: all N <= 10000000 verified (online invariant held after every n) ====
REPLAY out_quick/log.bin: N=10000000 records=1204766 steps(vertices)=1051583 moves=153183 pairs_checked=1249886 bad_pairs=0 errors=0  => PASS (proper t(n)-colouring of G_n for every n <= N)
DUMP CHECK out_quick/colouring_10000000.bin: N=10000000 classes=400000 members=1051583 max_class=8 pairs_checked=1010633 bad=0 errors=0 => PASS
== move log SHA-256: 5e2a72780d804b2bf92649d9eba2ecec9b61fa0efc665235f3b40c31c5293714
   matches the published value
PART B PASS for all N <= 10000000
```

and of `1e8`:

```
REPLAY out_1e8/log.bin: N=100000000 records=12046324 steps(vertices)=10515898 moves=1530426 pairs_checked=12493370 bad_pairs=0 errors=0  => PASS (...)
DUMP CHECK out_1e8/colouring_100000000.bin: N=100000000 classes=4000000 members=10515898 max_class=8 pairs_checked=10103256 bad=0 errors=0 => PASS
== move log SHA-256: 210e198ebee43b647e0057ae1e0d8f7c0a0db041dcd11072644eafaaac5f5274
   matches the published value
```

**Determinism.** `colour` is deterministic: its randomised search is seeded by the vertex index. Its output
does not depend on the thread count, and a run resumed from a checkpoint gives the same output as a fresh
run. So the move log for any NMAX is a byte-exact prefix of the published `log_2p32_full.bin`. These
properties were tested directly:

- The separate runs to 10⁵, 10⁷, 10⁸ and 10⁹ are prefixes of each other.
- At 10⁷, runs with 10 and with 3 threads give identical logs.
- At 3·10⁷, a run resumed from a 10⁷ checkpoint matches a fresh run.
- A fresh 10⁸ run reproduces the published hashes of `log_1e8.bin` and `colouring_100000000.bin`.

The hashes that `reproduce.sh 1e9` and `full` compare against are those of the original runs. A from-scratch
run to 2³² − 1 should reproduce `log_2p32_full.bin`, but it was not repeated in the fresh-copy test.

### Small-N checks, unit tests and negative controls

```sh
tests/run_tests.sh               # about 50 s; ends with "ALL TESTS PASS"
```

This script runs the following checks:

- `clique_check` checks every N ≤ 10⁵ by exact max-clique induction (`violations=0`), with no colouring involved.
- `brute.py` finds ω(G_N) = t(N) for N = 100, 500, 1000, 2000, 3000.
- `test_sqf` reports 0 failures.
- `test_v` passes.
- `checker` replays the included `runs/log_1e5.bin`.
- `brute_cmp` builds G_N from an edge sieve that uses no squarefree test, and finds 0 mismatches and 0 clique
  violations for N ≤ 10⁵.
- Two negative controls must be rejected, and are: a truncated log ("log ends after vertex 99982 ...") and a
  log with a planted bad pair ("pair (41,7): 41*7+1 is NOT squarefree").

## Results of the original runs

All of the following are from the logs in `verifier/runs/` (original runs, 2026-10-07, 10 threads).

| NMAX | \|V_N\| | t(N) | log records | moves of placed vertices | pairs re-certified by `checker` | `colour` wall | `checker replay` wall |
|---|---|---|---|---|---|---|---|
| 10⁸ | 10,515,898 | 4,000,000 | 12,046,324 | 1,530,426 | 12,493,370 | 15.3 s | 5.4 s |
| 10⁹ | 105,158,686 | 40,000,000 | 120,457,382 | 15,298,696 | 124,904,274 | 450 s (incl. 2 internal audits) | 156 s |
| 2³² − 1 | 451,653,474 | 171,798,692 | 517,358,696 | 65,705,222 | 536,466,170 | 1,327 s (resumed at 10⁹) | 2,399 s |

The checker's dump checks passed at N = 4.2·10⁷, 10⁸, 10⁹ and 2³² − 1. The check at 2³² − 1 took 1,813 s
and covered 433,903,423 same-class pairs.

Over the whole run, `colour` computed 578,428,609 squarefree certificates: 576,183,405 squarefree,
2,245,204 not squarefree, and **0 uncertain**. A C18 vertex needed a repair, a single eviction, 34.5% of the
time. No placement ever failed, and the maximum class size was 8.

Further checks recorded in the logs:

- The two V sieves agree up to 2³² − 1 (`test_v_2p32.log`).
- Direct tests of n² + 1 agree for n ≤ 10⁸.
- The certifier unit tests passed on about 18.5 M random and adversarial inputs plus exhaustive ranges, with
  0 failures (`test_sqf.log`).
- `brute_cmp` passed up to 3·10⁵ (`brute_cmp_3e5.log`).

### Large artifacts not included

The move logs and colouring dumps of the original runs (the seven files in `verifier/runs/SHA256SUMS.txt`
other than `log_1e5.bin` and `log_1e7.bin`) total about 8.4 GB and are not in the repository.
`reproduce.sh` regenerates them, except the dump at N₁, which needs one extra command:

```sh
cd verifier && make -s && mkdir -p out_n1 && \
  ./colour --nmax 42000000 --log out_n1/log.bin --dump-at 42000000 --dump-dir out_n1 > out_n1/colour.log && \
  ./checker dump out_n1/colouring_42000000.bin && shasum -a 256 out_n1/colouring_42000000.bin
```

Their SHA-256 (full list in `verifier/runs/SHA256SUMS.txt`):

| file | bytes | SHA-256 |
|---|---|---|
| `colouring_42000000.bin` (dump at N₁) | 24,386,684 | `16da05ba07522f9007fff433a8867940295ac27d946015d831e127c99b99d93a` |
| `log_1e8.bin` | 96,370,592 | `210e198ebee43b647e0057ae1e0d8f7c0a0db041dcd11072644eafaaac5f5274` |
| `log_1e9.bin` | 963,659,056 | `65b214d314e28103c608a072620ecdef5394cc714a3db41d42e6099fc3b18618` |
| `log_2p32_full.bin` | 4,138,869,568 | `7d933b60f12c49a78b2c934cc3cb0da192acde34e18ea6bd4e92c3b7979cb7c3` |

## Caveats and known issues

- **The first 2³² − 1 replay failed, correctly.** `colour` originally wrote its final checkpoint after closing
  the log, so the 10⁹ checkpoint recorded log offset 0. The resumed run's `log_2p32.bin` therefore holds only
  the records for n > 10⁹, and the checker rejected it (`runs/checker_2p32_replay.log`, kept on purpose). The
  complete log `log_2p32_full.bin` = `log_1e9.bin` ++ `log_2p32.bin` was then replayed from n = 1 and passed
  (`runs/checker_2p32_replay_full.log`). `colour.c` has since been fixed: it now checkpoints before `fclose`
  and checks offset consistency on resume.
- **Floating point in Part A.** All counting main terms and ceilings in `rigorous.py` are exact
  (`Fraction`/integers). Floating point is used only in these places, each rounded outward with explicit
  slack: sums of 1/p² over large primes, the logarithms in the Pell tail and in the knapsack weights, and the
  Rosser–Schoenfeld bound for π(x) above 2·10⁸.
  The tail check `tail_check2` (N ≥ 10¹⁴) is evaluated in double precision with a 1.01 safety factor. Its
  worst value, 0.02919, is far below the target 0.04.
- **Superseded code is kept.** `asymptotic.py` still contains the superseded v1 tail (`tail_check`, the
  `density_parts` stub), because the audit scripts use it. The certified tail is `tail_check2`. The
  docstring of `grid.py` was corrected to the interval argument actually used; its code is unchanged.
- **Untested platforms.** Part B was built and run only with Apple clang on arm64.
- **No formal verification.** This work is not formalised (Li's Lean development is a separate proof).

## Prior work and credit

- **Problem source.** P. Erdős and A. Sárközy; see [Er92b, p. 239] in the erdosproblems.com bibliography,
  and [erdosproblems.com/848](https://www.erdosproblems.com/848) with its
  [forum thread](https://www.erdosproblems.com/forum/thread/848) (from August 2025).
- **First density bounds** (as recorded on the [problem page](https://www.erdosproblems.com/848)).
  - **W. van Doorn**'s argument gives |A| ≤ (0.108… + o(1))N.
  - **D. Weisenberg** noted that the same proof gives about 0.105·N.
  - **S. Cambie** is credited in Sawhney's note, together with van Doorn and Weisenberg, with substantial
    progress towards the problem.
- **M. Sawhney (2025; acknowledging M. Sellke, with assistance from ChatGPT 5 Pro)** proved the statement for
  all sufficiently large N: [Problem_848.pdf](https://www.math.columbia.edu/~msawhney/Problem_848.pdf). Part A
  here follows his case split.
- **N. Sothanaphan, with GPT-5.4 (March 2026),** made Sawhney's argument explicit, with N₀ = 2.64·10¹⁷.
  The notes were posted in the erdosproblems.com forum thread.
- **Two full-resolution claims precede this work.** Both are unrefereed.
  - **Alex Chengyu Li** (Jul–Aug 2026), using OpenAI ChatGPT 5.5 and ChatGPT 5.6 and a personally developed
    AI-assisted proof-engineering system (as stated in the claim), has a Lean 4 formalisation. It is proof
    claim #204 on erdosproblems.com, the code is at [github.com/crabsatellite/erdos-848-squarefree-product](https://github.com/crabsatellite/erdos-848-squarefree-product),
    and the paper is SSRN 7230480. Its interval-wise checks start at 5·10⁶ and its uniform argument at
    5·10¹¹; below 5·10⁶ it uses colouring certificates.
  - **Ian Pitchford / OpenAI Codex** (Jul 2026) has a certificate chain up to 2.64·10¹⁷ and then
    Sothanaphan's bound. The code is at [github.com/ipitchford/erdos-848-all-n](https://github.com/ipitchford/erdos-848-all-n)
    and the archive is Zenodo [10.5281/zenodo.21647629](https://doi.org/10.5281/zenodo.21647629).
- **Other computations.** Forum participants ran independent exact checks for small N. For example, the
  user ScottHughes checked N ≤ 100,006 by max clique. See `notes/literature.md`.
- **Explicit number theory used.** Rosser–Schoenfeld (1962): π(x) < 1.25506 x / ln x. Robin (1983):
  ω(n) ≤ 1.38402 ln n / ln ln n.

This argument does not use the results, code or certificates of Li or Pitchford. Both claims were found during
the literature search for this work, and their public write-ups were read for comparison and attribution. Its
aim is a short, checkable proof whose finite part runs on a laptop. It is **not** a priority claim.

## AI-assistance disclosure

The mathematics, code, computations, audit, and this write-up were produced by **Claude**, an AI model by
Anthropic (model Claude Opus 5.5, running in Claude Code), at the request of [AUTHOR NAME].
[EDIT: state accurately what [AUTHOR NAME] did, e.g. posed the task, which parts they have read and checked.]
The adversarial audit in `audit/` was also AI-generated. The proof has not been peer reviewed. Please treat it as an
unrefereed, AI-generated, computer-assisted proof, and report any error you find (issues are welcome).

## How to cite

[AUTHOR NAME] (2026). *Erdős Problem #848 for all N: an explicit analytic bound plus a colouring certificate.*
[REPOSITORY URL] (AI-assisted; Claude Opus 5.5).

## License

The code is under the MIT License. The documentation, proof, notes, logs and data are under CC BY 4.0.
See `LICENSE`.

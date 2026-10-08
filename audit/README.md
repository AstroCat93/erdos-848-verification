# Adversarial audit of Part A

This directory holds the code and the small data files of the adversarial audit of the analytic part
(proof Sections 1–6). The audit's findings are summarised in `proof/proof.md` (Section 6, paragraph "Audit")
and below. Like the rest of this repository, the audit was carried out by an AI agent (a separate Claude
instance told to look for errors). It does not replace human review.

## Findings

- No fatal issue was found. The audit re-derived Lemmas 1–5 and the case split from the text of the proof.
  It also re-implemented the bounds of the two binding cases (II.a, II.b.2b) independently (`indep.py`). The
  proof reports that the two implementations agree to within about 2% and that `rigorous.py` always gives the
  larger bound. For example, at N = 42,840,000
  `indep.py` gives 0.037765·N for case II.b.2b and `rigorous.py` gives 0.038342·N; the target is 0.04·N.
- Brute-force tests found no bound exceeded (`brute_piv.c` + `cmp_brute.py`, `lemma5_test.py`, `rho_check.py`).
- These minor issues were raised and fixed in the current proof and code:
  1. **Monotonicity in N.** The original grid argument assumed the case bounds are monotone in N.
     `mono_check.py` probes this for a single class. The assumption was replaced by the "a valid A in [N]
     is valid in [N_hi]" argument (proof Section 5), which needs no monotonicity.
  2. **Matched-prime knapsack.** It now uses the exact LP-dual upper bound instead of a greedy sort.
  3. **Class label.** One mislabelled class was corrected.
  4. **Tail.** The tail argument was incomplete. It is now completed by `tail_check2`, which uses a doubling
     grid from 1e14 to 1e300 and a closed-form majorant beyond that. The v1 `tail_check` is superseded.

  These fixes were made after the audit and were re-checked by later separate review agents (also Claude
  instances), not by the original audit.

## Files

| file | what it does | usage / runtime (M4 Pro) |
|---|---|---|
| `indep.py` | Independent re-implementation of the Part A case bounds, written from `proof.md` and not from `rigorous.py`. It derives the statuses of p = 2, 5, q by enumerating pivot residues, uses an exact integer knapsack for matched primes, and has its own sieve and Pell sum. | `python3 indep.py 42840000` (27 s) |
| `brute_piv.c` | For every pivot b ≤ N and every x ≤ N, finds the smallest and largest p with p² \| xb+1, and prints per-pivot counts for the sets used in Lemmas 1, 2 and 5. | `cc -O2 -o brute_piv brute_piv.c && ./brute_piv 3000 > data/bp3000.txt` |
| `data/bp3000.txt`, `data/bp20000.txt` | Output of `brute_piv` for N = 3000 and 20000. The N = 1e5 and 4e5 files (6 MB and 27 MB) are left out; regenerate them with `brute_piv`. | |
| `cmp_brute.py` | Checks every brute-force count against the corresponding `rigorous.py` bound and prints `VIOLATIONS: [...]`. | `python3 cmp_brute.py 3000 data/bp3000.txt` prints `VIOLATIONS: []` (42 s) |
| `lemma5_test.py` | Lemma 5 against the truth for random pivots and for every pivot with ≥ 2 matched primes. | `python3 lemma5_test.py N Qs Qsplit` |
| `rho_check.py` | Brute-force check of ρ*(b) ≤ 2^ω_odd(b)·c₂(v₂(b)) for b ≤ 3000, and of the maximisation in `rigorous.rho_star`. | `python3 rho_check.py` (4 s) |
| `diag_tail.py` | Compares the O-tail and G bounds of `rigorous.py` and `indep.py` at one N. | `python3 diag_tail.py N` |
| `mono_check.py` | Counts how often a single class's gain decreases between consecutive N, for N near 4.2e7 (finding 1). | `python3 mono_check.py` |
| `run_rig.py` | Component breakdown of `rigorous.all_cases` at N = 42,840,000 and 1e9. | `python3 run_rig.py` |
| `run_asym.py`, `tail_nb.py` | Exercise the superseded v1 tail in `asymptotic.py` (finding 4). | |

All scripts locate `../analytic` relative to their own path, so they can be run from any directory. They need
numpy, and `cmp_brute.py` and `rho_check.py` also need sympy.

# Erdős #848 literature notes (compiled 2026-10-07)

## 0. Forum state (fetched 2026-10-07)
- https://www.erdosproblems.com/forum/thread/848 : 48 comments, newest = sproutseeds 05 Apr 2026 (post 5255). "The ability to post new comments has been suspended."
- Problem page "last edited 06 December 2025" (still says solved for sufficiently large N by Sawhney).
- **Proof claims tab (1)**: https://www.erdosproblems.com/forum/thread/848/proof-claims
  - Claim #204 by Alex Chengyu Li (alexchengyuli), submitted 2026-08-16 03:37:45, "using OpenAI ChatGPT 5.5; OpenAI ChatGPT 5.6; and a personally developed AI-assisted research and proof-engineering system."
  - Claims max = floor((N+18)/25) for every N >= 1.
  - Summary: Hall-equivalence reformulation; N <= 5e6 by certified prefix colourings; above that, square-divisor estimates + matching reduce Hall defect to structured residual partitioned by five 2-adic types and residues mod 9; three-/four-pivot arguments + bounds on solutions of a transformed quadratic equation; "exact rational inequalities over finitely many exhaustive structural cases, followed by a uniform argument for N >= 5*10^11".
  - Proof: https://papers.ssrn.com/abstract=7230480
  - Formalisation: https://github.com/crabsatellite/erdos-848-squarefree-product (tag v1.0.5-kernel; Lean 4; 75 OLean cache shards; ~32 GiB RAM, 200 GB disk)
  - Claim comments (https://www.erdosproblems.com/forum/proof-claims/204/comments):
    - ScottHughes 19 Aug 2026: independent recheck; f(N)=floor((N+18)/25) for all N<=100,006 (own max-clique solver); "every constant in Sothanaphan's 2.64e17 note re-derives exactly (its true crossover is ~2.636e17)"; tightest middle-range gate (E2 case) clears by ~1e-5; mentions a SECOND independent all-N determination by ipitchford / "OpenAI Codex": paper https://doi.org/10.5281/zenodo.21647629 , repo https://github.com/ipitchford/erdos-848-all-n ("certificate chain plus envelope arguments to 2.64e17 and then Sothanaphan"). Caveat: OLean cache is Windows-x86_64 only.
    - MalekZ 24 Aug 2026: "author-approved Significance record" https://hjyuh.github.io/significance/2026-alexchengyuli-erdos-848/

## 1. Other sources checked (Q1)
- teorth/erdosproblems (cloned 2026-10-07, HEAD 98e8a47): data/problems.yaml #848 status `decidable`, last_update 2025-10-19; formalized yes (2026-02-01). Not updated after the Jul/Aug 2026 all-N claims.
- Issues: #151 (Tao, 2025-11-03) lists 848 as "solved asymptotically, but perhaps some combination of explicit analytic number theory and numerical methods can be used to convert this to a complete solution." #196 wiki suggestion. No PRs mentioning 848.
- Wiki (teorth/erdosproblems.wiki, last commit 2026-06-30): AI-contributions table rows for 848: Sothanaphan explicit bound (5-23 Mar 2026), MalekZ "Incorrect proof found", Sawhney/Sellke GPT-5 "Full solution" (asymptotic), erdosbanger Lean (28 Jan 2026). No row for Li or ipitchford.
- arXiv API searches (ab+1 squarefree; Erdos Sarkozy squarefree; non-squarefree; Erdos problem 848; au:Sothanaphan; etc.): NO arXiv paper on #848 in 2026. (2511.16072 = Bubeck et al., "Early science acceleration experiments with GPT-5", contains Sawhney–Sellke 848 story.)
- SproutSeeds/erdos-problems (npm workspace): pack packs/number-theory/problems/848; last 848 commit 2026-04-25; PROGRESS.md (2026-04-19) "Global Decision: not_decided"; raw exact coverage N<=10,000 (local 40,500); uses imported N0=2.64e17; repo-audited threshold only exp(1420). No new threshold. Audited MalekZ repo as blocked.
- ipitchford/erdos-848-all-n (created 2026-07-28; Zenodo 10.5281/zenodo.21647629 published 2026-07-28): "unrefereed candidate computer-assisted result", f(N)=floor((N+18)/25) all N. Ranges: 1..100000006 coloring certificates; 1e8..1e9 structural split; 1e9..1e12 exact-rational short-shift envelopes; 1e12..2.64e17 exact-rational rank envelopes; >=2.64e17 Sothanaphan's PDF (SHA-256 8162113a...796f, which matches the Drive file I downloaded). Attribution: OpenAI Codex workflow; Ian Pitchford. No external reproduction/peer review/formalization claimed.
- crabsatellite/erdos-848-squarefree-product (Alex Chengyu Li; created 2026-07-24; releases v1.0.0-kernel 2026-07-30 ... v1.0.5-kernel 2026-08-04; last commit 2026-09-05 "Correct diagonal selector descriptions"). SSRN 7230480. Lean 4 endpoint `Erdos848.PaperGeneratedCertificateProvider.all_N`, axioms only propext/Classical.choice/Quot.sound, 30,638 modules. Zenodo concept DOI 10.5281/zenodo.21750213.

### Li paper structure (paper/erdos_848_kernel_checked_exact_extremal_bound.tex, "August 2026")
- Hall reformulation: C_N = 7 mod 25 class; for compatible outside B, need |B|+|T_N(B)| <= |C_N|; target |H_N(B)|/N <= 1/25 - 7/(25N).
- N<=5e6: prefix colouring certificate.
- Pure opposite-base matching (D_N=18 mod 25 into C_N, N>=5e6): primes split p<=47 / 47<p<=N/26 / p>N/26.
  - Lemma (high prime, SWITCHING TRICK): "Write bx+1=p^2 m. Since p>N/26 ... bx+1 <= N^2+1 < 677(floor(N/26)+1)^2. Hence m<=676." x=r+50t => p^2 m=(br+1)+50bt. For q|b prime (q not 2,5): p^2 m ≡ 1 (mod q) so m is a nonzero QR mod q. For fixed m, root classes s mod 50b with s^2 m ≡ br+1 (mod 50b). Spacing: two solutions p1<p2 in one class => p2-p1>=50b, and 50b(t2-t1)=m(p2-p1)(p2+p1) => t2-t1 >= 2m(N/26). Coefficient 11264064/3125 (count, then /N, worst at N=5e6).
  - Medium primes: (1/25)Σ_{47<p<=N/26} 1/p^2 < 3887e-6/25 plus endpoint π(N/26)/N <= 59/10^4 (cites Rosser–Schoenfeld 1962, Dusart arXiv:1002.0442).
- Uniform mixed degree (N>=5e6): residual pivot b has > N/50+N/525 squarefree neighbours in each of C_N, D_N; split p in {2,3,7} / 7<p<=N/55 / p>N/55 (m<=3025; "normal" and "mod-5-twisted" transformed equations); tau_7, rho_7 constants. => |R_N(B)| > 2N/525 - 1.
- Three-pivot "half decomposition" inequality: 1 <= (1/2)1_{k=3} + (1/2)1_{k>=2} + (1/2) l(x).
- Ten-case decomposition by 2-adic class (E1,E2,E3,O1,O3) and mod 9 cells; budget D_beta + F_beta + Q_beta < 1/25 - 7/(25L).
- Intervals: [5e6,1e7), [1e7,2e7), [2e7,4e7), [4e7,2e8) (6 blocks), [2e8,2e9) (splits s=75,85,100,125), [2e9,5e11) (rows 2,2.3,2.8,4.2,6,12,500 billion), and N>=5e11 uniform: split 120 (p>N/120, m<=14400), "840^4 <= 5e11", factor (N^{1/4}+1)/N <= 841/840^4, sieve primes up to 47, common transformed-root bound 8062340e-9, diagonal bounds 1e-9(26005710,...,917036), counting moduli (1821,2892,2627,2295,3643,4590,3643,3643).
- Displayed example margin [2e7,4e7) E2-generic: slack = 75490557093924693317991/7596668444022826249000000000 ≈ 9.94e-6 (in units of N).
- Li Lean statement (lean4/Erdos848/ProblemCore.lean), read 2026-10-07:
  `def NonSquarefreeProductProp (A : Finset ℕ) : Prop := ∀ a ∈ A, ∀ b ∈ A, ¬ Squarefree (a * b + 1)`
  `def OriginalProblem848Statement (N : ℕ) : Prop := ∀ A : Finset ℕ, A ⊆ Finset.Icc 1 N → NonSquarefreeProductProp A → A.card ≤ ((Finset.Icc 1 N).filter fun n => n % 25 = 7).card`
  `theorem all_N : ∀ N, OriginalProblem848Statement N` (PaperGeneratedCertificateProvider.lean). Statement looks faithful (diagonal included). I did NOT build it.
- Significance record (hjyuh.github.io, MalekZ): "Significance has not independently rebuilt the repository, audited the finite certificates, or assessed the mathematical argument." 0 reviews.
- No GitHub issues/PRs on either all-N repo. teorth wiki/yaml not updated for them.
- Computation coverage claims (all sources): MalekZ 1e7 (withdrawn/likely incorrect); SproutSeeds 10,000 (local 40,500); erdosbanger Lean N=50,100; ScottHughes max-clique N<=100,006; Li prefix colouring N<=5e6 (Lean); ipitchford coloring certs N<=100,000,006.
- Sub-2.64e17 explicit thresholds: Li's analytic "unbounded tail" N>=5e11 (different method: Hall/pivots, needs his structural reductions); ipitchford rank envelopes 1e12..2.64e17 (exact-rational, per-interval, NOT a uniform threshold). old-bielefelder (22 Mar 2026) AI ping-pong "at least 3.0e17, perhaps 2.94e17" (unpublished, superseded). ScottHughes: Sothanaphan's inequalities' true crossover ≈ 2.636e17.

## 2. Sothanaphan explicit-threshold notes (Drive), read in full
Files (downloaded via drive uc?export=download): 2.64e17 = 1ujhm4_WYpgRV_rd1rJXIfHyvx16COEKe (dated Mar 24 2026, 12pp, SHA-256 8162113a...796f); 3.3e17 = 1SR1U_9Rp4B7nQwljjhCDvdCRZ7ln9Jtb (Mar 23); 7e17 = 1-8Xo-CNtRA9LnLaSZslJjVMgMkOQIbPa (Mar 22); exp(1420) = 1lSEzb9l-ThEc4ygZUg6KQGIIYvedd6gx (Mar 7); exp(1958) = 1nYBzwaBcMVWWW_Ncn38TERGalpuPz2hW; variant (distinct a,b) = 15yo2FQFm5sZkgu5Ny1WiAWN3F5XMER2p; althofer.de/exp1185_or_exp1420.pdf (Mar 5).

### 2.64e17 note structure
- Lemma 1: ρ(q) := max_{(u,q)=1} #{y mod q: y^2≡u} ≤ λ q^{1/4}, λ = 128/(8·3·5·7·11·13)^{1/4} < 6.876.
- Prop 1 (explicit squarefree count in AP; Prachar/Hooley style): for q ≤ x, (a,q)=1, D ≥ 1:
  |S2(x;q,a) − (x/q)∏_{p∤q}(1−1/p^2)| ≤ D + x/(qD) + λ q^{1/4}(2x/(qD) + x/D^2).
  Tail Σ_{d>D} A_d handled by SWITCHING: n = m d^2, m ≤ x/D^2, d^2 ≡ a m̄ (mod q) has ≤ ρ(q) classes, d ≤ sqrt(x/m) => ≤ ρ(q) Σ_m (sqrt(x/m)/q + 1). (All d > D, not only primes.)
- Cor 1 (replaces Sawhney L2.2): q0 ∈ {25,100}, b ≤ N, p^2 ∤ bt+1 for p | q0, N ≥ 2q0:
  #{a ≤ N: a ≡ t (q0), µ(ab+1)=0} ≤ (N/q0)(1 − ∏_{(p,q0 b)=1}(1−1/p^2)) + κ_{q0} N^{3/4}, κ25 = 4.7, κ100 = 5.3 (D = ⌈α N^{3/4}⌉, α25 = 25/8, α100 = 7/2; E25(2.64e17) < 4.700, E100 < 5.275). Modulus fed to Prop 1 is d1 = q0 b/g ≤ q0 N, so the crude λ d1^{1/4} costs N^{1/4}.
- Lemma 2: Cquad := 1 − ∏_{p≡1(4), p≥13}(1−2/p^2) < 0.02734669.
- Prop 2 (replaces Sawhney L2.1), q ∈ {25,50}, t ≢ 7,18 (mod 25), R ≥ 2:
  #{x ≤ N: x ≡ t (q), µ(x^2+1)=0} ≤ X·Cquad + R(1+log R) + 2X(log R + 2)/R + ((N^2+1)/R^2)(1 + log_{2+√3} N).
  Large r > R: x^2+1 = d r^2, d ≤ (N^2+1)/R^2, negative Pell x^2 − d u^2 = −1; consecutive solutions grow by > 2+√3 => ≤ 1 + log_{2+√3} N per d. (Pell tail counts ALL x ≤ N, not just the class.)
- Cor 2: R = ⌈(259/200) N^{2/3}⌉; per-class error < 8.577e-5 at 2.64e17; ×23 classes => 0.001973. My recomputation: per class t1 = R(1+log R)/N = 5.653e-5, Pell tail = 2.924e-5, middle term 4e-12.
- Constants: C_{2,5} = 1 − (6/π^2)/((3/4)(24/25)) = 0.1556568029; C5 = 1 − 25/(4π^2) = 0.3667426022.
- Case margins (|A|/N < 0.04 needed), at N = 2.64e17:
  Case 1 (A* has even elt): 23/25 Cquad + 2/25 C25 + 0.001973 + 9.4 N^{-1/4} => 0.0376115 + 0.001973 + 0.000415 = 0.0399995; slack ≈ 5e-7. BINDING. (main-term slack 0.0023885)
  Case 2a: 0.0356927 + 0.001973 + 0.000468 (10.6 N^{-1/4}) + 1e-17 = 0.0381337; slack 0.00187.
  Case 2b: 0.0334755 + 0.001973 + 0.000415 = 0.0358635; slack 0.00414.
  Case 3: 2/25 C5 + 9.4N^{-1/4} = 0.0293395 + 0.000415; slack 0.01025.
- Cost ranking (Case 1 budget 0.0023885): x^2+1 lemma (Prop2/Cor2) 0.001973 (~83%: Möbius-truncation endpoint 23·R(1+logR)/N ≈ 0.00130; Pell tail 23·2.92e-5 ≈ 0.00067); Cor 1 (ab+1 in AP) 0.000415.
  [Inference, not in the note: the Pell tail T_R is a bound over all x ≤ N, so it need not be paid 23 times; also R(1+log R) uses Σ_{r≤R} d(r) though only r composed of primes ≡1 (4), ≥13 contribute.]
- Evolution: 7e17 (Mar 22): ρ(q) ≤ 12 q^{1/4}, κ25=7, κ100=8, R=⌈(8/5)N^{2/3}⌉, Pell count log_2(2N), Case 1 0.0376119+0.0019+0.000485. 3.3e17 (Mar 23): λ<11.564, κ25=5.6, κ100=6.3, α25=15/4, α100=21/5, R=⌈1.3N^{2/3}⌉, log_{2+√3}. 2.64e17 (Mar 24): λ<6.876 (includes the 2-part 8), κ 4.7/5.3, R=⌈1.295N^{2/3}⌉.

### exp(1420) note (Mar 7) + Althofer PDF (Mar 5)
- No switching: Prop 2 (CRT-periodic): small p ≤ T exact CRT density with boundary ≤ q∏_{p≤T}p^2; large p > T union bound N/(qp^2) + 1 each => + π(N). Lemma 2(a) for A*: (23/25)C1 + 46π(√N)/N + 2π(N;4,1)/N (each p > √N gives ≤ 2 x ≤ N).
- C1 = Σ_{p≥13, p≡1(4)} 2/p^2 < 0.02762792. M1 = 23/25 C1 + 2/25 C2 < 0.037870231, δ = 1/25 − M1 > 0.002129769.
- E1(N) at L=1420, T=697: 1.5159e-5 + 1.40945e-3 (2π(N)/N, Dusart) + 7.04729e-4 (2π(N;4,1)/N, BMOR) + ... < 0.00212934 vs δ 0.00212977. The cost is the "+1 per prime" boundary term π(N)/N ~ 1/log N; exp(1185) fails (Ẽ = 0.00255273).
- Prime bounds used: Dusart [Ramanujan J. 45 (2018) 227–251, Eq. (5.4)] as quoted: π(x) ≤ x/log x (1 + 1/log x + 2.53816/log^2 x) "for all x ≥ 88789" (I could not verify the constant/range independently); BMOR 2018 Thm 1.3 with constant 1/160 for x ≥ 8e9 (note: for 3 ≤ q ≤ 10^4 the paper's c0(q) is actually 1/840; 1/160 is the q > 10^4 value — weaker but valid).

## 3. Citable explicit inputs
(a) ab+1 / bn+1 in AP with p^2 | bn+1, p > Y, uniform in b:
- No peer-reviewed explicit version found. Non-explicit: Prachar, Monatsh. Math. 62 (1958) 173–176; Hooley, Bull. LMS 7 (1975) 133–138; Nunes, Mathematika 63 (2017) 483–498 (arXiv:1605.03347).
- Explicit (AI-generated/unrefereed) forms:
  * Sothanaphan 2.64e17 note Prop 1 + Cor 1 (above).
  * ipitchford proofs/CORRECTED_LEMMAS.md "Lemma A": for q a perfect square, p^2 ∤ bt+1 for p | q, D ≥ 2,
    |S − (N/q)∏_{p∤qb}(1 − 1/p^2)| ≤ D + N/(qD) + C_b(2(bN+1)/(bD) + (bN+1)/D^2), C_b = 2^{ω(b)+1} (safe) or ρ(b) = 2^{ω(b_odd)} r(v2(b)), r(0)=r(1)=1, r(2)=2, r(e≥3)=4. (Switching mod b instead of mod q0 b.)
  * Li paper: large primes p > N/s, write bx+1 = p^2 m with m ≤ s^2 (s=26: N^2+1 < 677(⌊N/26⌋+1)^2 => m ≤ 676); m must be a nonzero QR mod every prime q | b (q ≠ 2,5); for fixed m, p lies in root classes s^2 m ≡ br+1 (mod 50b); spacing p2−p1 ≥ 50b and 50b(t2−t1) = m(p2−p1)(p2+p1) => t2−t1 ≥ 2m(N/26). [Inference: the QR restriction on m offsets the 2^{ω(b)} growth in root classes, which is why the bound is uniform in b; the paper handles this by a finite enumeration indexed by the prime-divisor set S(b).] Medium primes 47 < p ≤ N/26: union bound with endpoint π(N/26)/N ≤ 59/10^4 at N ≥ 5e6.
(b) n ≤ N with q^2 | n^2+1 for some q > Y:
- Published non-explicit: Estermann, Math. Ann. 105 (1931) 653–662 (O(x^{2/3} log x)); Heath-Brown, Acta Arith. 155 (2012) 1–13, arXiv:1010.6217 (O(x^{7/12+ε})). No explicit published version found.
- Explicit (unrefereed): Sothanaphan Prop 2 (Pell, ratio 2+√3). ipitchford "Lemma B": T(N,P) ≤ C_8(M)·K, q = P+1, M = ⌊(N^2+1)/q^2⌋, H = N + sqrt(N^2+1), L = sqrt(2q^2−1) + q√2, R = 3+2√2, K = 0 if H < L else 1 + ⌊log(H/L)/log R⌋; C_8(M) = 3a + 1_{r≥1} + 1_{r≥2} + 1_{r≥5} − ⌈s/2⌉ for M = 8a + r, s = ⌊√M⌋ (counts m ≤ M with m ≡ 1,2,5 mod 8 minus squares). Basis: all positive solutions of n^2 − m p^2 = −1 are odd powers of one fundamental solution, ratio ε^2 ≥ (1+√2)^2 = 3+2√2.
  Also ipitchford's union-bound version: #{n ≤ N: n mod Q ∈ T, ∃p ≥ 13, p ≡ 1 (4), p^2 | n^2+1} ≤ (rN/Q)S + 2rπ(P) + ⌊(N^2+1)/P^2⌋(1 + log(2N+1)/log(2+√3)), S = 2Σ_{p≥13, p≡1(4)} 1/p^2 < 59/2000.
  Structure of negative Pell solutions: standard textbook fact (e.g. Niven–Zuckerman–Montgomery, Intro. Theory of Numbers, 5th ed., §7.8 — check theorem number).
(c) ω(n), d(n), π(x):
- Robin, Acta Arith. 42 (1983) 367–389: Thm 11 ω(n) ≤ 1.38402 log n/log log n (n ≥ 3) [often quoted 1.3841]; Thm 12 ω(n) ≤ log n/log log n + 1.45743 log n/(log log n)^2 (n ≥ 3); Thm 13 ω(n) ≤ log n/(log log n − 1.1714) (n ≥ 26); Thm 16 ω(n) ≤ log n/loglog n + log n/(loglog n)^2 + 2.89726 log n/(loglog n)^3 (n ≥ 3). (Verified via survey arXiv:2404.17165, Thms 30–33.)
- Nicolas–Robin 1983 (Canad. Math. Bull. 26): log d(n)/log 2 ≤ 1.5379 log n/log log n (n ≥ 3).
- Rosser–Schoenfeld, Illinois J. Math. 6 (1962) 64–94: π(x) < 1.25506 x/log x (x > 1); π(x) > x/log x (x ≥ 17) [Cor. 1, (3.5)–(3.6)].
- Dusart arXiv:1002.0442 Thm 6.9: x/ln x (1+1/ln x) ≤ π(x) (x ≥ 599) and π(x) ≤ x/ln x (1 + 1.2762/ln x); x/(ln x − 1) ≤ π(x) (x ≥ 5393); π(x) ≤ x/(ln x − 1.1) (x ≥ 60184); x/ln x(1 + 1/ln x + 2/ln^2 x) ≤ π(x) (x ≥ 88783); π(x) ≤ x/ln x(1 + 1/ln x + 2.334/ln^2 x) (x ≥ 2,953,652,287).
- Bennett–Martin–O'Bryant–Rechnitzer, Illinois J. Math. 62 (2018) 427–532, arXiv:1802.00085: Thm 1.3 |π(x;q,a) − Li(x)/φ(q)| < c_π(q) x/(log x)^2 for x ≥ x_π(q), c_π(q) ≤ 1/840 (3 ≤ q ≤ 10^4), x_π ≤ 8·10^9 (q ≤ 10^5). Cor 1.6: for 1 ≤ q ≤ 1200, x ≥ 50q^2: x/(φ(q) log x) < π(x;q,a) < x/(φ(q) log x)(1 + 5/(2 log x)).
- Squarefree counts: Cohen–Dress–El Marraki (Funct. Approx. 37 (2007)): |Σ_{n≤x} µ^2(n) − 6x/π^2| ≤ 0.02767 √x for x ≥ 438653; Cohen–Dress 1988: |R(x+y) − R(x)| ≤ 1.6749√y + 0.6212 x/y; Ramaré 2019: Σ µ^2(n) = 6x/π^2 + O*(1.06 √(x/log x)) (x > 1). (Quoted from TME-EMT database, tme-emt-1f5ce3.pages.math.cnrs.fr/Art10.html.)

## 4. Tao–van Doorn arXiv:2512.01087 (v2, 7 Dec 2025) Appendix A
- Theorem 16: A ⊆ [M+1, M+N] avoiding ω(p^2) < p^2 classes mod p^2 for each p; for any Q ≥ 1: |A| ≤ (N + Q^4)/Σ_{q≤Q} h(q), h(q) = µ^2(q)∏_{p|q} ω(p^2)/(p^2 − ω(p^2)). Proof via Montgomery–Vaughan arithmetic large sieve (fractions a/q' with q | q' | q^2 are 1/Q^4-separated). Fully explicit.
- Usability (inference): it is an UPPER bound for a SIFTED set (one avoiding classes). For #848 one needs upper bounds on the set of a with p^2 | ab+1 for SOME p (union of classes), i.e. a lower bound for survivors — the large sieve doesn't give that. Also Q^4 error forces Q ≤ N^{1/4}, and Σ h(q) converges, so it can't beat the trivial density anyway. Matches Sothanaphan/GPT's "couldn't make use of this large sieve".

# Erdős Problem #848 for all N — an independent verification (explicit analytic bound plus colouring certificate)

*AI disclosure: the mathematics, code and text were produced by Claude (Anthropic, Claude Opus 5.5), with adversarial review by separate Claude instances and computer verification. Nathan Durfee (composer and prompter) conceived and directed the project; he has not personally verified the mathematics. See README, "AI-assistance disclosure".*

*Draft, 2026-10-07. Status: an independent verification, by a complete computer-assisted proof for all N, of the resolution previously claimed by A. C. Li and by I. Pitchford / OpenAI Codex (Section 8); analytic part audited by a separate AI agent (no fatal issues; minor issues fixed); finite part replayed by an independently written checker. Not peer reviewed.*

## 0. Statement

For $N\ge1$ let $t(N)=\#\{n\le N:\ n\equiv 7\pmod{25}\}=\lfloor (N+18)/25\rfloor$.

**Theorem.** If $A\subseteq[N]$ and $ab+1$ is not squarefree for all $a,b\in A$ (including $a=b$), then $|A|\le t(N)$.
Equality holds for $A=\{n\le N: n\equiv7 \pmod{25}\}$, since $7\cdot 7+1=50$ and so $25\mid ab+1$ whenever $a\equiv b\equiv 7\pmod{25}$.

**Structure of the proof.**
- Part A (analytic, Sections 1–6): the theorem holds for every $N\ge N_1:=4.2\cdot10^{7}$. The numerical inequalities are checked by exact rational arithmetic on a finite grid of intervals covering $[4.2\cdot10^7,1.14\cdot10^{14}]$ (Section 5), plus a uniform tail argument for $N\ge10^{14}$ (Section 6: explicit error bounds evaluated on a doubling grid up to $10^{300}$, in floating point with a 1% safety factor, and a closed form beyond).
- Part B (computational, Section 7): the theorem holds for every $N\le 2^{32}-1$ (in particular every $N\le N_1$), certified by an online proper colouring.

Sawhney (2025; acknowledging M. Sellke, with assistance from ChatGPT 5 Pro) proved the theorem for all sufficiently large $N$. Part A follows his case split. The tools it adds are standard, and several appear in earlier explicit work:
1. **Switching.** Large prime squares are counted through the cofactor ($ab+1=mp^2$ with $m$ small), not one prime at a time. This replaces an error of size $\pi(N)$ by about $2^{\omega(b)}N^{2/3}$. The device is classical (Estermann, Hooley, Heath-Brown) and is used in explicit form by Sothanaphan, Li and Pitchford.
2. **Self-sieving of $A^*$.** The elements of $A$ outside $7,18\pmod{25}$ are sieved by the same "pivot" elements that sieve the residue classes.
3. **A mod-100 parity analysis.** This exploits that $x\equiv7$ and $y\equiv43\pmod{100}$ have $xy+1\equiv 2\pmod 4$ and $xy+1\equiv 2\pmod{25}$.

The elements outside $7,18\pmod{25}$ with a large prime-square divisor of $x^2+1$ are bounded with the negative Pell equation (Lemma 4), as in Sothanaphan and Pitchford.

Notation:
- $V_N=\{n\le N:\ \mu(n^2+1)=0\}$. Every element of a valid $A$ lies in $V_N$.
- $C_7=\{n\equiv7\ (25)\}$, $C_{18}=\{n\equiv18\ (25)\}$, $O=V\setminus(C_7\cup C_{18})$.
- For a valid $A$: $A_7=A\cap C_7$, $A_{18}=A\cap C_{18}$, $A^*=A\cap O$.
- $R_{25}=\mathbb Z/25\setminus\{7,18\}$.

## 1. Sieving one residue class by one element (Lemma 1)

Fix $b\in[N]$ (the *pivot*), a modulus $L$, a class $r\bmod L$, and parameters $P_0\le P_1$. For a prime $p$ let
$$B_p=\{x\in[N]:\ x\equiv r\ (L),\ p^2\mid xb+1\}.$$

- If $p\mid b$ then $B_p=\emptyset$.
- Otherwise $p^2\mid xb+1\iff x\equiv -b^{-1}\pmod{p^2}$. Together with $x\equiv r\ (L)$, this is either:
  - inconsistent ("never"),
  - automatic ("always": $p^2\mid\gcd(L,p^2)$ and $p^2\mid rb+1$),
  - or a single class modulo $\operatorname{lcm}(L,p^2)$ ("sub"). In that case the class has relative density $g_p=\gcd(L,p^2)/p^2$ inside the class $r\bmod L$.

**Lemma 1.** If some $p\le P_1$ is "always", then $|\bigcup_{p\le P_1}B_p|\le\lceil N/L\rceil$. Otherwise, writing $S=\{p\le P_0 \text{ sub}\}$ and $k=|S|$,
$$\Big|\bigcup_{p\le P_1}B_p\Big|\ \le\ \frac NL\Big(1-\prod_{p\in S}(1-g_p)\Big)+(2^k-1)+\sum_{\substack{P_0<p\le P_1\\ \text{sub}}}\Big\lceil \frac{N g_p}{L}\Big\rceil .$$

*Proof.*
- **Primes $p\in S$.** For each $\emptyset\ne T\subseteq S$, the congruences $x\equiv r\ (L)$ and $x\equiv -b^{-1}\ (p^2)$ for $p\in T$ are pairwise compatible. The $p^2$ are pairwise coprime, and each is compatible with $L$ by the "sub" hypothesis. So by the Chinese remainder theorem they define one class modulo $L\prod_{p\in T}g_p^{-1}$.
- Its number of elements in $[N]$ is $\frac NL\prod_{T}g_p+\theta_T$ with $|\theta_T|<1$.
- Inclusion–exclusion is exact. Summing the main terms gives $\frac NL(1-\prod_S(1-g_p))$, and the errors total less than $2^k-1$.
- **Primes $P_0<p\le P_1$.** The union bound applies, with $\#\{x\le N \text{ in one class mod } M\}\le\lceil N/M\rceil$. $\square$

## 2. Large primes: the switching lemma (Lemma 2)

For $b\ge1$ let
$$\rho^*(b)=\max_{m}\#\{x\bmod b:\ mx^2\equiv1\ (b)\}.$$
Since $mx^2\equiv1$ forces $\gcd(m,b)=1$ and then $x^2\equiv m^{-1}$, the count is a number of square roots of a unit. Hence $\rho^*(b)\le 2^{\omega_{\rm odd}(b)}c_2(v_2(b))$, with $c_2(0)=c_2(1)=1$, $c_2(2)=2$, $c_2(\ge3)=4$.

**Lemma 2.** For $b\in[N]$ and $Y\ge1$,
$$E(b;N,Y):=\#\{a\in[N]:\ \exists\text{ prime }p>Y,\ p^2\mid ab+1\}\ \le\ \rho^*(b)\Big(\frac{2(N+1)}{Y}+\frac{N^2+1}{Y^2}\Big).$$

*Proof.*
- Map $a\mapsto (m,p)$ with $ab+1=mp^2$. This is injective, and the image satisfies $mp^2\equiv1\ (b)$ and $Y<p\le X_m:=\sqrt{(Nb+1)/m}$, so $m<M:=(Nb+1)/Y^2$.
- For fixed $m$, $p$ lies in at most $\rho^*(b)$ classes modulo $b$, so there are at most $\rho^*(b)(X_m/b+1)$ choices.
- Summing over $m<M$ and using $\sum_{m\le M}m^{-1/2}\le2\sqrt M$, the total is $\le\rho^*(b)\,(2(Nb+1)/(Yb)+M)$.
- Finally use $(Nb+1)/b\le N+1$ and $Nb+1\le N^2+1$. $\square$

**Corollary (global part of one pivot).** For $P_1<Y$:
$$\#\{x\in[N]:\ \exists p>P_1,\ p^2\mid xb+1\}\le \sum_{P_1<p\le Y}\lceil N/p^2\rceil+E(b;N,Y)=:\mathrm{Gl}(b).$$

**Combining.** Let $\mathcal C$ be a family of pairwise disjoint residue classes. Every $x$ with $\mu(xb+1)=0$ has a prime $p$ with $p^2\mid xb+1$, and either $p\le P_1$ or $p>P_1$. Therefore
$$\#\{x\in\textstyle\bigcup\mathcal C\cap[N]:\ \mu(xb+1)=0\}\le\sum_{C\in\mathcal C}\mathrm{Lemma1}(C)+\mathrm{Gl}(b).$$
In particular, **the global part is paid once per pivot**, however many classes that pivot sieves.

We need $\max\rho^*(b)$ over $b\le N$ of a prescribed parity (and $5\nmid b$ for pivots in $C_7\cup C_{18}$). It is computed exactly: maximise $2^{\omega_{\rm odd}}c_2$ over $b=2^v\cdot(\text{odd part})\le N$.

## 3. The set $O$ (Lemmas 3–4)

**Lemma 3.** If $x\in O$, then $x\bmod 25\in R_{25}$, and there is a prime $q\equiv1\ (4)$ with $q\ge13$ and $x\equiv\pm r_q\pmod{q^2}$, where $\pm r_q$ are the two roots of $r^2\equiv-1\pmod{q^2}$.

*Proof.*
- $p^2\nmid x^2+1$ for $p=2$, since $x^2+1\equiv1,2\pmod 4$.
- $p^2\nmid x^2+1$ for $p\equiv3\ (4)$, since $-1$ is not a square modulo such $p$.
- $25\mid x^2+1\iff x\equiv\pm7\pmod{25}$, which is excluded for $x\in O$. $\square$

**Lemma 4 (tail).** Let $\mathrm{div}=2$ if a parity is prescribed and $1$ otherwise. For $Q_s<Y_q$:
$$\#\{x\le N \text{ of the parity}:\ \exists q>Q_s,\ q\equiv1\ (4),\ q^2\mid x^2+1\}\le\sum_{\substack{Q_s<q\le Y_q\\ q\equiv1(4)}}2\Big\lceil\frac{N}{\mathrm{div}\,q^2}\Big\rceil+\mathrm{Pell}(N,Y_q).$$
Here
$$\mathrm{Pell}(N,Y_q)=\sum_{\substack{2\le m\le (N^2+1)/Y_q^2\\ 4\nmid m,\ p\mid m\Rightarrow p\not\equiv3(4)}}\Big(\tfrac12+\tfrac{\ln(2N+1)}{2\ln(1+\sqrt m)}\Big).$$

*Proof.*
- For $q\le Y_q$, the condition defines 2 classes modulo $\mathrm{div}\cdot q^2$.
- For $q>Y_q$, write $x^2+1=mq^2$, so $m\le(N^2+1)/Y_q^2$. Here $m$ divides $x^2+1$, so $4\nmid m$ and $m$ has no prime factor $\equiv 3\ (4)$.
- The solutions of $x^2-mq^2=-1$ with $x,q>0$ are the norm $-1$ units $\eta^{2k+1}$, $k\ge0$, of the order $\mathbb Z[\sqrt m]$, where $\eta\ge1+\sqrt m$.
- From $x=(u-u^{-1})/2\le N$ we get $u\le2N+1$. Hence the number of $k$ is at most $\frac12+\frac{\ln(2N+1)}{2\ln\eta}$.
- For $m$ a square there are no solutions. $\square$

## 4. Sieving $O$ by a pivot (Lemma 5)

Fix a pivot $b$, a parity for $x$, and parameters $Q_{\rm split}\le Q_s$. For each prime $q\equiv1\ (4)$ with $13\le q\le Q_s$ and each sign $\varepsilon=\pm1$, consider the classes $x\equiv\varepsilon r_q\ (q^2)$ of the given parity.
- If $q\le Q_{\rm split}$, refine further into the 23 classes $x\bmod 25\in R_{25}$.
- If $q> Q_{\rm split}$, drop the $R_{25}$ restriction. This only overcounts.

In each class apply Lemma 1, with these statuses:
- **$p=q$:** "never", unless $b\equiv \varepsilon r_q\ (q^2)$. In that case the whole $\varepsilon$-class is compatible with $b$ ("matched").
- **$p=5$:**
  - In split classes, "always" for exactly one $s\in R_{25}$ if $b\in O$ (namely $s=-b^{-1}$, which lies in $R_{25}$ because $b\not\equiv\pm7$). If instead $b\in C_7$ or $C_{18}$, then $s=-b^{-1}\in\{7,18\}\notin R_{25}$, so "never".
  - In unsplit classes, "sub" with $g=1/25$.
- **$p=2$:** "never" if $b$ or $x$ is even; otherwise $g=1/2$ (odd $x$) or $g=1/4$ (any $x$).
- **Other $p$:** $g=1/p^2$, the worst case $p\nmid b$.

**Matched primes.** Each matched $q$ contributes (full class) − (sieved bound) extra, and matches only one sign. If $q$ is matched then $q^2\mid b^2+1$. Hence $\prod_{q\in M}q^2\le N^2+1$, and also $25\prod_{q\in M}q^2\le N^2+1$ if $b\in C_7\cup C_{18}$. The worst case over admissible $M$ is bounded above by the fractional-knapsack (LP) relaxation, with weights $\ln q$ and capacity $\frac12\ln(N^2+1)$ (minus $\ln5$).

The remaining $q>Q_s$ are handled by Lemma 4. Call the resulting bound $G(b)$. It excludes $\mathrm{Gl}(b)$, which is added once per pivot.

## 5. The case analysis (Theorem A)

Let $A$ be valid with $|A|\ge t(N)+1$. Then $A\not\subseteq C_7$, and $A\not\subseteq C_{18}$ because $|C_{18}\cap[N]|\le t(N)$. Every bound below is a sum of Lemma 1/5 class bounds plus one $\mathrm{Gl}$ per pivot used; each must be shown $<\frac{N-6}{25}+1\le t(N)+1$.

Odd classes mod 100:
- $C_7^{\rm odd}=\{7,57\}$, $C_{18}^{\rm odd}=\{43,93\}$.
- $7\equiv43\equiv3$ and $57\equiv93\equiv1\pmod 4$.

**Key observation.** If $x\equiv7,\ y\equiv43\pmod{100}$, or $x\equiv57,\ y\equiv93\pmod{100}$, then $xy+1\equiv2\pmod4$ and $xy+1\equiv 2\pmod{25}$. So such a pair needs a prime $p\ne2,5$.

Consequently the two classes of such a pair contribute at most
$$\mathrm{pair}=\max\Big(\lceil N/100\rceil,\ \text{(each class sieved by an element of the other)}+2\,\mathrm{Gl}\Big).$$

**Case I: $A^*=\emptyset$** (so $A_7\neq\emptyset\neq A_{18}$).

| Case | Situation | Bound used |
|---|---|---|
| I.1 | even $c\in A_7$, even $d\in A_{18}$ | $C_7$ sieved by $d$, $C_{18}$ sieved by $c$; $p\ne2,5$ |
| I.2 | even $c\in A_7$, $A_{18}$ all odd | $A_{18}\subseteq 43\ (50)$ sieved by $c$; $C_7$ sieved by odd $d\in A_{18}$ |
| I.3 | $A_7$ all odd, even $d\in A_{18}$ | symmetric to I.2 |
| I.4 | everything odd | $\mathrm{pair}(7,43)+\mathrm{pair}(57,93)$ |

**Case II.a: $A^*$ has an even element $b$.** $C_7$ and $C_{18}$ are sieved by $b$ with $p\ne 2,5$, and $A^*$ is sieved by $b$ via Lemma 5 (any parity). One pivot.

**Case II.b: $A^*\ne\emptyset$ is all odd; fix $b\in A^*$.**

| Case | Situation | Bound used |
|---|---|---|
| II.b.1 | $A_7,A_{18}$ all odd | If $b\equiv1\ (4)$, classes $\equiv3\ (4)$ are always compatible with $b$, so the pair $(7,43)$ uses $\mathrm{pair}$, while $57,93$ (where $xb+1\equiv2\ (4)$) are sieved by $b$ with $p\ne2,5$. Symmetric if $b\equiv3\ (4)$. $A^*$ uses $\min(\text{crude }\lvert O^{\rm odd}\rvert, G(b))$. |
| II.b.2a | even $c\in A_7$, even $d\in A_{18}$ | $C_{18}$ sieved by $c$, $C_7$ sieved by $d$, $A^*$ uses $\min(\lvert O^{\rm odd}\rvert,G(c),G(d))$ |
| II.b.2b | even $c\in A_7$, $A_{18}$ all odd | $A_{18}\subseteq43\ (50)$ sieved by $c$, $C_7$ sieved by $b$, $A^*$ uses $\min(\lvert O^{\rm odd}\rvert,G(c),G(b))$ |
| II.b.3b | $A_7$ all odd, even $d\in A_{18}$ | symmetric to II.b.2b |

**Checks on the case split.**
- The cases are exhaustive.
- **The 5-statuses are correct.** $25\mid xb+1$ with $x\equiv7$ needs $b\equiv7$, and with $x\equiv18$ needs $b\equiv18$. Neither happens for $b\in O$. For the $C_7$–$C_{18}$ cross pairs, $7\cdot18+1=127\equiv2\pmod{25}$.

**Parameters.**
- $P_0=13$, $P_1=1000$, $Q_s=400$, $Q_{\rm split}=60$.
- $Y$ is optimised per pivot type, and $Y_q$ per tail.
- $\pi(x)$ is exact up to $2\cdot10^8$; above that, Rosser–Schoenfeld $\pi(x)<1.25506\,x/\ln x$.
- Code: `analytic/rigorous.py`, using exact `Fraction` arithmetic for all main terms and integers for all ceilings.

**Interval grid (no monotonicity needed).** Let $N_{lo}\le N\le N_{hi}$ and let $A\subseteq[N]$ be valid. Then $A\subseteq[N_{hi}]$ is also valid. If $A\subseteq C_7$ or $A\subseteq C_{18}$, then $|A|\le t(N)$ directly. Otherwise the case analysis *at $N_{hi}$* gives $|A|\le \mathrm{UB}(N_{hi})$. So $\mathrm{UB}(N_{hi})<(N_{lo}-6)/25+1\le t(N)+1$ certifies every $N\in[N_{lo},N_{hi}]$. `analytic/grid.py` chains such intervals over $[4.2\cdot10^7,10^{14}]$, in three segments with logs `analytic/grid2/log_{A,B,C}.txt`: $[4.2\cdot10^7,1.0037\cdot10^{10}]$ (83 intervals), $[10^{10},1.0043\cdot10^{12}]$ (36 intervals), $[10^{12},1.1417\cdot10^{14}]$ (35 intervals), all certified with zero retries; minimum slack $0.00089\,N_{lo}$ (first interval, case II.b.2b).
(An earlier single run, logged in `analytic/grid_log.txt` and using a greedy knapsack sort, covered $[4.2\cdot10^7,4.7\cdot10^{14}]$. The certified rerun uses the exact dual LP bound $\mathrm{OPT}\le\lambda\,\mathrm{cap}+\sum_q\max(0,g_q-\lambda w_q)$, valid for every $\lambda\ge0$.)

Sample margins, in units of $N$ (the worst case is II.b.2b at small $N$ and II.a at large $N$):

| $N$ | worst UB$/N$ | needed |
|---|---|---|
| $4.2\cdot10^7$ | 0.03845 | $<0.04$ |
| $10^8$ | 0.03450 | |
| $10^9$ | 0.03052 | |
| $10^{12}$ | 0.02788 | |
| $10^{40}$ | 0.02806 | |
| $10^{100}$ | 0.02819 | |

## 6. Tail $N\ge10^{14}$ (`analytic/asymptotic.py: tail_check2`)

**Limiting densities.** Let $D_{\rm case}$ be the sum of all main densities: the class mains, the mid-prime mains, $\sum_{p>P_1}p^{-2}$ per pivot, the $O$-tail density $\sum_{q>Q_s,\,q\equiv1(4)}2/(\mathrm{div}\,q^2)$, and **all** matched gains switched on. These are computed in exact arithmetic by `D_cases()`. Then
$$\max_{\rm case}D_{\rm case}=0.0281945\ldots \quad\text{(case II.a)}.$$
Every case bound satisfies $\mathrm{UB}_{\rm case}(N)\le D_{\rm case}N+\mathrm{Err}$.

**Error terms.** With $Y=Y_q=\lceil N_i^{0.8}\rceil$ held fixed on the interval $[N_i,2N_i]$,
$$\mathrm{Err}\le c_0+4\big(\pi(Y)+E(2N_i,Y,\rho^*(2N_i))\big)+2\pi(Y_q)+\mathrm{Pell}(2N_i,Y_q).$$
- **The constant $c_0$.** It counts the inclusion–exclusion errors, the mid-prime ceilings and the full-class ceilings: 226 per class across at most 348 classes, plus 10 for the pair bounds, so $c_0\le 78{,}658$.
- **Pivots.** At most 4 pivots are used in any case (I.4 uses 4; II.b.1 uses 3).
- **Prime counts.** $\pi$ is bounded by Rosser–Schoenfeld: $\pi(x)<1.25506\,x/\ln x$.
- **$\rho^*$.** It is computed exactly.

**The check.** Doubling from $N_i=10^{14}$ to $10^{300}$, we verify $D_{\max}+1.01\,\mathrm{Err}/N_i<1/25-1/N_i$. The worst value attained is $0.02919$.

**Beyond $10^{300}$.**
- Use Robin's bound $\omega(n)\le1.38402\ln n/\ln\ln n$ for $n\ge3$, so that $\rho^*(N)\le4\cdot2^{1.38402\ln N/\ln\ln N}$.
- The error majorant then has the form
$$\frac{c_0}{N}+7\cdot1.25506\,\frac{N^{-0.2}}{0.8\ln N}+17\cdot2^{1.38402\ln N/\ln\ln N}\,\big(2N^{-0.8}+N^{-0.6}\big)+10N^{-1/2}.$$
- **Each term decreases in $N$.** For the third term, $\frac{d}{du}\big[0.9594\,u/\ln u-0.6u\big]<0$ for $u=\ln N\ge 690$.
- The majorant is $<10^{-61}$ at $10^{300}$.

**Combined with the grid.** The grid covers $[4.2\cdot10^7,10^{14}]$, so Part A holds for all $N\ge4.2\cdot10^{7}$.

**Audit.** An adversarial audit by a separate AI agent (another Claude instance; audit/) re-derived Lemmas 1–5 and the case split, and re-implemented the two binding cases (II.a, II.b.2b) from this text. Its values for these cases lie within about 2% of ours, and rigorous.py is always the larger (for example $0.037765N$ against $0.038342N$ for case II.b.2b at $N=42{,}840{,}000$). Brute-force tests found no bound exceeded. It found no fatal issue. Its minor points are addressed:
- monotonicity has been replaced by the argument above,
- the knapsack now uses the dual LP bound,
- a class label was fixed,
- the tail is now complete.

The fixes it prompted (interval principle, dual LP bound, tail_check2) were made after the audit and were re-checked by later separate review agents, not by the original audit.

## 7. Part B: all $N\le 2^{32}-1$ by online colouring (Theorem B)

Let $G_N$ be the graph on $V_N$ with $a\sim b$ iff $\mu(ab+1)=0$; then $\omega(G_N)\le\chi(G_N)$. We increase $n$ one at a time and, at every $N$, maintain a proper colouring of $G_N$ satisfying:
- the classes are indexed by $C_7\cap[N]$, so there are exactly $t(N)$ of them;
- class $k$ contains $c_k=25k-18$;
- every same-class pair $u\ne v$ has $uv+1$ *certified* squarefree.

**Updates.**
- A new $n\in C_7$ opens a class.
- Any other new vertex is inserted into an existing class, possibly after a repair. At most two earlier vertices move per insertion: typically a $C_{18}$ element swaps with another class's $C_{18}$ element (34.5% of $C_{18}$ insertions); in about 6.4M cases one further vertex also moves. Every repair was found at search depth 1.
- No static colouring can work, because the nested property already fails at $n=43$. Hence the colouring must be updated online.

**Certification.** Squarefreeness of $m=uv+1<2^{64}$ is certified as follows:
1. Trial division by small primes.
2. For the cofactor $c$: deterministic Miller–Rabin primality testing (Jaeschke and Sinclair bases), an exact perfect-square test when $c<B^3$, and Pollard–Brent factorisation with a deterministic cube-root trial-division fallback.

Uncertain cases are treated as edges, which is sound. No uncertain case occurred.

**Results** (`verifier/`, 10 threads):
- The invariant held at every $N\le 2^{32}-1=4{,}294{,}967{,}295$, with $|V|=451{,}653{,}474$; `colour` computed $578{,}428{,}609$ squarefreeness certificates (0 uncertain).
- The maximum repair depth was 1.
- `colour` took 450 s to reach $10^9$ and a further 1,327 s (resumed from a checkpoint) to reach $2^{32}-1$, about 30 min in all; the independent replay and dump check took about 70 min.

**Independent checks.**
- A separate program (`checker.c`) replayed the complete operation log from $n=1$. It re-certified every new same-class pair (536,466,170 pair checks) using a squarefree routine that shares no code with the main program, and dump checks were run at $4.2\cdot10^7$, $10^8$, $10^9$ and $2^{32}-1$.
- Two independently written sieves for $V$ agree up to $2^{32}-1$, and direct tests of $n^2+1$ for $n\le10^8$ also agree.
- An exact clique-number induction confirms $\omega(G_N)\le t(N)$ for all $N\le3\cdot10^5$ without any colouring, and a colouring check against an edge-sieve graph (no squarefree test) agrees for $N\le3\cdot10^5$.
- Certifier unit tests: about $18.5$M random and adversarial inputs plus exhaustive ranges, with 0 failures.
- Negative controls: a planted bad pair and a truncated log are both rejected.

**Overlap.** Part B covers $N\le4.29\cdot10^9$, about $100\times$ beyond the analytic threshold $N_1=4.2\cdot10^7$, so the two parts overlap with a wide margin.

## 8. Relation to other work

- Density bounds (per the problem page erdosproblems.com/848): van Doorn's argument gives $|A|\le(0.108\ldots+o(1))N$, and Weisenberg noted that the same proof gives about $0.105N$. Sawhney's note also credits Cambie, with van Doorn and Weisenberg, with substantial progress.
- Sawhney (2025; acknowledging M. Sellke, with assistance from ChatGPT 5 Pro): proof for all large $N$.
- Sothanaphan + GPT‑5.4 (Mar 2026): explicit $N_0=2.64\cdot10^{17}$.
- Two unrefereed full-resolution claims (Jul/Aug 2026):
  - A. C. Li (with OpenAI ChatGPT 5.5 and ChatGPT 5.6 and a personally developed AI-assisted research and proof-engineering system): Lean 4 proof; colourings for $N\le5\cdot10^6$, interval-wise checks from $5\cdot10^6$, uniform argument from $5\cdot10^{11}$.
  - I. Pitchford / OpenAI Codex: certificate chain up to $2.64\cdot10^{17}$.

For comparison, our interval checks cover $[4.2\cdot10^7,1.14\cdot10^{14}]$ and our uniform argument covers $N\ge10^{14}$, so Li's thresholds are lower on both counts. The aim here is an independent, checkable verification whose finite part runs on a laptop. We became aware of both claims during this work, after the approach of Part A had been designed; we then read their public write-ups (including their explicit switching and Pell lemmas) for comparison and attribution. All lemmas here are proved in full, and no result, code or certificate of theirs is used. This is not a priority claim.

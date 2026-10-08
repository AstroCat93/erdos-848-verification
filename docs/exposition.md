
# Erdős Problem #848 for every $N$: an independent verification

*Nathan Durfee (composer and prompter), with Claude (Anthropic)*

> **Authorship and AI disclosure.** Written by Claude (Anthropic, model Claude Opus 5.5) at the direction of Nathan Durfee, who composed and prompted the project. Claude produced the mathematics, the code and the text; a separate Claude instance carried out an adversarial audit. Below, "we" means these Claude instances. No affiliation. The work has not been refereed, and no human, including Nathan Durfee, has checked it line by line. Nathan Durfee is responsible for releasing it. Code, logs and an issue tracker for corrections: https://github.com/AstroCat93/erdos-848-verification.

## 0. Summary

Erdős and Sárközy asked whether a set $A\subseteq\{1,\dots,N\}$ in which $ab+1$ is never squarefree can be larger than the residue class $7 \bmod 25$. We give a computer-assisted proof that it cannot, for every $N$: $|A|\le t(N)=\lfloor (N+18)/25\rfloor$.

This is an **independent verification** of a result already claimed: A. C. Li (with a Lean 4 formalisation) and I. Pitchford with OpenAI Codex posted full-resolution claims in July–August 2026; neither has been refereed. We make no priority claim, and our thresholds are not lower than Li's. We became aware of both claims during this work, after the approach of Part A had been designed; we then read their public write-ups (including their explicit switching and Pell lemmas) for comparison and attribution. All lemmas here are proved in full, and no result, code or certificate of theirs is used.

The proof has two parts.

- **Part A** follows the case split of Sawhney's large-$N$ argument and makes it explicit for $N\ge 4.2\cdot 10^7$, with some additional devices. If $A$ is contained in neither class $7$ nor class $18 \bmod 25$, some $b\in A$ (a *pivot*) lies outside the class being sieved, so $25\nmid xb+1$ for every $x$ in it; an even pivot, or in the all-odd cases a pairing of classes mod $100$, also rules out $4$. Then $\mu(xb+1)=0$ needs some $p^2$ with $p\ne2,5$, which happens on at most about $15.6\%$ of the class. Primes $p>1000$ cost one switching bound per pivot, across all classes. A nine-case analysis gives a worst-case limiting density of $0.02819\ldots$ against the $1/25$ required; the inequalities are certified exactly on $154$ intervals up to $1.14\cdot 10^{14}$, and by explicit error terms beyond.
- **Part B** covers $N\le 2^{32}-1$. It uses a proper colouring of the conflict graph with exactly $t(N)$ colours, maintained online as $N$ grows and replayed by an independent checker.

What may be new (§8): the pivots also sieve the elements outside the two classes, with a knapsack bound on "matched" primes; the mod-$100$ pairing; and the online move-log colouring with an independent replay.

## 1. The problem, the extremal example, and why $25$ and $\pm7$

Call $A\subseteq[N]=\{1,\dots,N\}$ *admissible* if $\mu(ab+1)=0$ for all $a,b\in A$, **including $a=b$**. This is the natural reading of the problem page and the one in Li's formalisation; we do not treat the variant with $a\ne b$. Let
$$t(N)=\#\{n\le N:\ n\equiv 7 \ (\mathrm{mod}\ 25)\}=\Big\lfloor\frac{N+18}{25}\Big\rfloor .$$

**Theorem.** For every $N\ge 1$, every admissible $A\subseteq[N]$ satisfies $|A|\le t(N)$.

**Equality.** Since $7\cdot 7+1=50$, we have $25\mid ab+1$ whenever $a\equiv b\equiv 7 \pmod{25}$. Similarly $18^2+1=325=13\cdot 25$. The class $18 \bmod 25$ has $\lfloor (N+7)/25\rfloor$ elements, so it is also extremal whenever this count equals $t(N)$. The two classes cannot be combined: for $x\equiv7$, $y\equiv18\pmod{25}$ we have $xy+1\equiv2\pmod{25}$, so $25$ never helps, and indeed $7\cdot18+1=127$ is squarefree.

**Why these residues.** Taking $a=b$ shows that every element of $A$ lies in
$$V_N=\{n\le N:\ \mu(n^2+1)=0\}.$$
Suppose $p^2\mid n^2+1$. Then $p\ne 2$, because $n^2+1\equiv 1,2 \pmod 4$, and $p\not\equiv 3\pmod 4$, because $-1$ is not a square modulo such $p$. For each prime $p\equiv 1\pmod 4$ the congruence $r^2\equiv-1 \pmod{p^2}$ has exactly two roots $\pm r_p$. Each *root class* $x\equiv r_p\pmod{p^2}$ is admissible on its own, since $xy+1\equiv r_p^2+1\equiv 0\pmod{p^2}$.

So $V$ is a union of root classes. Its density is
$$1-\prod_{p\equiv 1(4)}\bigl(1-2p^{-2}\bigr)=0.10516\ldots.$$
Since $A\subseteq V_N$, this gives $|A|\le(0.10516\ldots+o(1))N$, matching the $\approx0.105N$ noted by Weisenberg; the problem page attributes $|A|\le(0.108\ldots+o(1))N$ to van Doorn. The densest root classes are those for $p=5$, namely $7$ and $18\pmod{25}$, each of density $1/25$. The theorem says that no admissible set, whether or not it is a union of root classes, does better.

**Decomposition.** We split
$$V_N=C_7\sqcup C_{18}\sqcup O,$$
where $C_r=\{n\le N: n\equiv r\ (25)\}$ and $O$ is the rest of $V_N$. For an admissible $A$ we write $A_7=A\cap C_7$, $A_{18}=A\cap C_{18}$ and $A^*=A\cap O$. The densities are $1/25$ each for $C_7$ and $C_{18}$, and
$$\tfrac{23}{25}\Bigl(1-\prod_{q\equiv1(4),\,q\ge13}(1-2q^{-2})\Bigr)\approx 0.0252$$
for $O$.

**Earlier work.**

- The problem is due to Erdős and Sárközy [Er92b, p. 239] and is erdosproblems.com/848.
- Sawhney (2025) proved the theorem for all sufficiently large $N$, together with a stability version. His note credits van Doorn, Weisenberg and Cambie with substantial progress and thanks M. Sellke; the argument was obtained with assistance from ChatGPT 5 Pro.
- N. Sothanaphan, working with GPT-5.4, made the argument explicit for $N\ge N_0=2.64\cdot10^{17}$ (March 2026).
- The two full-resolution claims, by A. C. Li and by I. Pitchford (with OpenAI Codex), neither refereed, are compared with ours in §8.

## 2. Strategy

Suppose $|A|\ge t(N)+1$. Since $|C_{18}|\le t(N)$, the set $A$ is contained in neither $C_7$ nor $C_{18}$. Together with admissibility, this is all that §§3–5 use about $A$.

**Pivots.** Every $b\in A$ can act as a *pivot*: each $x\in A$ has some $p^2\mid xb+1$. For $x\in C_7$ and $b\notin C_7$ we have $25\nmid xb+1$, because $25\mid 7b+1$ forces $b\equiv 7\pmod{25}$; if $b$ is also even, $4\nmid xb+1$. The remaining primes $p\ne 2,5$ catch at most a fraction
$$1-\prod_{p\ne2,5}(1-p^{-2})\approx 0.156$$
of $C_7$ (a prime $p\mid b$ catches nothing). So an even pivot outside $C_7\cup C_{18}$ shrinks $C_7\cup C_{18}$ from density $0.08$ to at most about $0.0125$ ($0.0126$ with the explicit bound of Lemma 1).

**Ingredients.** Small primes $p\le 1000$ inside a residue class are handled by exact CRT inclusion–exclusion and a union bound (Lemma 1). Large primes $p>1000$ are handled by a *global* switching bound that does not depend on the class, so each pivot pays it once however many classes it sieves (Lemma 2). The pivots also sieve $A^*$, which is covered by the root classes of primes $q\ge 13$ (Lemma 3); the root classes of very large $q$ are handled by a negative-Pell bound (Lemma 4), and "matched" primes by a knapsack bound (Lemma 5). Depending on the parities and residues present in $A$, a case analysis (§5) chooses one to four pivots and assigns to each the sets it sieves. Parameters: $P_0=13$, $P_1=1000$, $Q_s=400$, $Q_{\rm sp}=60$.

## 3. Counting lemmas

### 3.1 The class sieve

Fix a pivot $b$ and a class $K=\{x\le N: x\equiv r\ (L)\}$. For a prime $p$ let
$$B_p=\{x\in K:\ p^2\mid xb+1\}.$$
If $p\mid b$ then $B_p=\emptyset$. Otherwise $p^2\mid xb+1$ is equivalent to $x\equiv -b^{-1}\pmod{p^2}$. Combined with $x\equiv r\ (L)$, the CRT says this condition is one of three kinds:

- **never:** it has no solutions in $K$;
- **always:** it holds on all of $K$, which happens exactly when $p^2\mid L$ and $p^2\mid rb+1$;
- **sub:** it defines a single class modulo $\mathrm{lcm}(L,p^2)$, of relative density $g_p=\gcd(L,p^2)/p^2$ in $K$.

**Lemma 1 (class sieve).** Suppose no $p\le P_1$ is *always* (otherwise we use the trivial bound $\lceil N/L\rceil$), and let $S$ be the set of *sub* primes $p\le P_0$. Then
$$\Big|\bigcup_{p\le P_1}B_p\Big|\ \le\ \frac NL\Big(1-\prod_{p\in S}(1-g_p)\Big)+\bigl(2^{|S|}-1\bigr)+\sum_{\substack{P_0<p\le P_1\\ p\ \text{sub}}}\Big\lceil \frac{Ng_p}{L}\Big\rceil .$$

*Proof.* Let $\emptyset\ne T\subseteq S$. The congruences for $p\in T$ are pairwise compatible, so $\bigcap_{p\in T}B_p$ is one class modulo $L\prod_T g_p^{-1}$. Its size is $\frac NL\prod_T g_p+\theta_T$ with $|\theta_T|<1$. Inclusion–exclusion is exact: the main terms sum to $\frac NL\bigl(1-\prod_S(1-g_p)\bigr)$, and the errors total less than $2^{|S|}-1$. For the primes in $(P_0,P_1]$, use the union bound. $\square$

The right-hand side only grows when a *never* prime is declared *sub*. So every prime $p\nmid L$ may be given its worst case $g_p=p^{-2}$, uniformly over all pivots, and only the statuses of $2$, $5$ (and of $q$ in Lemma 5) need thought. For $L\in\{25,50,100\}$ write $\Lambda_L=\Lambda_L(N)$ for the right-hand side of Lemma 1 on a class of modulus $L$ with $2$ and $5$ both *never* and $g_p=p^{-2}$ for every other $p\le P_1$, and $\Lambda'_{25}$ for the variant in which $2$ is *sub* with $g_2=1/4$ (an odd pivot). For example, $\Lambda_{25}/N\to 0.0063191$, which is about $0.158/25$ (slightly above $0.156/25$ because the primes in $(13,1000]$ are union-bounded).

### 3.2 Large primes: switching

Let
$$\rho^*(b)=\max_m\#\{x \bmod b:\ mx^2\equiv 1\ (\mathrm{mod}\ b)\}$$
be the largest number of square roots of a unit modulo $b$. Then $\rho^*(b)\le 2^{\omega_{\rm odd}(b)}c_2(v_2(b))$, where $c_2=1,1,2,4$ for $v_2(b)=0,1,2,\ge 3$.

**Lemma 2 (switching).** For $b\le N$ and $Y\ge 1$,
$$E(b;N,Y):=\#\{a\le N:\ \exists\ \text{prime } p>Y,\ p^2\mid ab+1\}\ \le\ \rho^*(b)\Big(\frac{2(N+1)}{Y}+\frac{N^2+1}{Y^2}\Big).$$

*Proof.* Write $ab+1=mp^2$ with $p>Y$.

- The map $a\mapsto(m,p)$ is injective, and $m<M:=(Nb+1)/Y^2$.
- For fixed $m$ we have $p\le X_m=\sqrt{(Nb+1)/m}$ and $mp^2\equiv1\pmod b$. So $p$ lies in at most $\rho^*(b)$ classes mod $b$, which allows at most $\rho^*(b)(X_m/b+1)$ values.
- Sum over $m<M$, using $\sum_{m<M} m^{-1/2}\le 2\sqrt M$. This gives $\rho^*(b)\bigl(2(Nb+1)/(bY)+(Nb+1)/Y^2\bigr)$.
- Finally use $(Nb+1)/b\le N+1$ and $Nb+1\le N^2+1$. $\square$

Counting by the cofactor $m$ instead of by $p$ avoids the "$+1$ per prime" of a union bound over $p\le N$. That union bound costs about $\pi(N)\approx N/\log N$, which fits inside the final slack of about $0.0118N$ only once $\log N\gtrsim 85$, that is $N\gtrsim10^{37}$. Taking $Y\approx N^{2/3}$ already makes the error $O(\rho^*(b)\,N^{2/3})$. Switching is classical (Estermann, Hooley, Heath-Brown), and explicit versions appear in Sothanaphan, Li and Pitchford.

**One global term per pivot.** For any $Y>P_1$, let
$$\mathrm{Gl}(b)=\sum_{P_1<p\le Y}\lceil N/p^2\rceil+E(b;N,Y).$$
Every $x$ with $\mu(xb+1)=0$ has some $p^2\mid xb+1$, either with $p\le P_1$ (counted by Lemma 1 in whichever class contains $x$) or with $p>P_1$ (counted by $\mathrm{Gl}(b)$). Hence, for any classes $K_1,\dots,K_s$ that need not be disjoint,
$$\#\{x\in K_1\cup\dots\cup K_s:\ \mu(xb+1)=0\}\ \le\ \sum_i(\text{Lemma 1 bound for }K_i)+\mathrm{Gl}(b).$$
We take the worst case over all $b\le N$ of a given kind and parity (with $5\nmid b$ for pivots in $C_7\cup C_{18}$), maximising the bound $2^{\omega_{\rm odd}(b)}c_2(v_2(b))$ for $\rho^*(b)$ exactly; for $b\le 4.2\cdot10^7$ it is at most $512=2^7\cdot4$, attained by $b=8\cdot3\cdot5\cdot7\cdot11\cdot13\cdot17\cdot19$. $Y$ is optimised numerically. The main density of $\mathrm{Gl}$ is only $\sum_{p>1000}p^{-2}\approx 1.27\cdot10^{-4}$. Its error term is what makes small $N$ expensive: for the two pivots of the binding case at $N=4.2\cdot10^7$, $\mathrm{Gl}$ is $0.0064N$ and $0.0051N$, whereas at $N=10^{14}$ every pivot type costs about $2\cdot10^{-4}N$.

### 3.3 The structure of $O$ and the negative-Pell tail

**Lemma 3 (structure of $O$).** If $x\in O$, then $x\bmod 25\in R_{25}:=\mathbb Z/25\setminus\{7,18\}$, and $q^2\mid x^2+1$, that is $x\equiv\pm r_q\pmod{q^2}$, for some prime $q\equiv 1\pmod 4$ with $q\ge 13$.

The prime $5$ is excluded because $25\mid x^2+1$ forces $x\equiv\pm7\pmod{25}$. So $O$ is covered by the $46$ classes
$$\{x\equiv\pm r_q\ (q^2),\ x\equiv s\ (25)\},\qquad s\in R_{25},$$
for each prime $q\le Q_s=400$, together with a tail of $x$ for which some $q>400$ works.

**Lemma 4 (tail).** Let $\delta=2$ if a parity of $x$ is prescribed and $\delta=1$ otherwise, and let $Q_s<Y_q$. Then
$$\#\{x\le N\ \text{(of the prescribed parity)}:\ \exists\, q>Q_s,\ q\equiv1\ (4),\ q^2\mid x^2+1\}\ \le \sum_{\substack{Q_s<q\le Y_q\\ q\equiv1(4)}}2\Big\lceil\frac{N}{\delta q^2}\Big\rceil+\sum_{m}\Big(\frac12+\frac{\ln(2N+1)}{2\ln(1+\sqrt m)}\Big),$$
where $m$ runs over $2\le m\le (N^2+1)/Y_q^2$ with $4\nmid m$ and no prime factor $\equiv 3\pmod 4$.

*Proof.* For $q\le Y_q$ there are two classes modulo $\delta q^2$. For $q>Y_q$, write $x^2+1=mq^2$, so that $m$ is in the stated range. If $m$ is a square there are no solutions. Otherwise let $\varepsilon=x_1+y_1\sqrt m$ be the least solution with $x_1,y_1\ge1$ of $x^2-my^2=\pm1$ (the fundamental unit of the order $\mathbb Z[\sqrt m]$), so $\varepsilon\ge1+\sqrt m$. If its norm is $-1$, the positive solutions of $x^2-my^2=-1$ are exactly the $\varepsilon^{2j+1}$; otherwise there are none. The number $u=x+q\sqrt m$ satisfies $2x=u-u^{-1}$, so $x\le N$ forces $u<2N+1$. That leaves at most $\frac12+\ln(2N+1)/(2\ln\varepsilon)$ solutions for each $m$. $\square$

After optimising $Y_q$, the error beyond the main term is $O(N^{2/3})$. Non-explicit bounds of this kind go back to Estermann; see also Heath-Brown. The explicit Pell argument also appears in Sothanaphan and Pitchford.

## 4. Sieving $A^*$ by a pivot: matched primes

This step is the main source of slack. Fix a pivot $b$ and a parity type $\tau$: $\tau$ is "all $x$" for an even pivot $b\in O$ (Case II.a), and "odd $x$" otherwise (Case II.b, where $A^*$ is all odd). Write $O^\tau$ for the elements of $O$ of type $\tau$. Apply Lemma 1 to each $O$-class $K_{q,\varepsilon,s}=\{x\equiv\varepsilon r_q\ (q^2),\ x\equiv s\ (25)\}$ of type $\tau$. For $q>Q_{\rm sp}=60$ we drop the condition on $s$, which only enlarges the class. The statuses are as follows.

- **$p=q$.** Since $q^2$ divides the modulus, $q$ is *always* or *never*. On the class, $xb+1\equiv \varepsilon r_q b+1\pmod{q^2}$. Since $(\varepsilon r_q)^{-1}=-\varepsilon r_q$, this vanishes exactly when $b\equiv\varepsilon r_q\pmod{q^2}$. So $q$ is *always* exactly when the pivot lies in the same root class as the class being sieved. We then call $q$ **matched**. Otherwise $q$ is *never*.
- **$p=5$.** In a split class ($q\le Q_{\rm sp}$) and with $5\nmid b$, we have $25\mid xb+1$ exactly when $s\equiv-b^{-1}\pmod{25}$. Since $-7^{-1}\equiv 7$ and $-18^{-1}\equiv 18$, if $b\in O$ this singles out exactly one $s\in R_{25}$, and that class is counted in full. If $b\in C_7\cup C_{18}$, then $-b^{-1}\in\{7,18\}$, so $5$ is *never*. In an unsplit class ($q>Q_{\rm sp}$), $25$ does not divide the modulus, so $5$ is *sub* with $g_5=1/25$. If $5\mid b$, then $5$ is *never*.
- **$p=2$.** It is *never* if $b$ is even. If $b$ is odd (so $\tau$ is "odd $x$"), it is *sub* with $g_2=1/2$.
- **Other $p$.** They are given the worst case $g_p=p^{-2}$.

**Why matching happens, and why it is rare.** A matched prime is the $O$-analogue of the extremal example. The root class $x\equiv r_q\pmod{q^2}$ is itself admissible, just as $C_7$ is, and a pivot inside it is compatible with all of it. So a matched $q$ gives back
$$\mathrm{gain}_q=(\text{full size of that class})-(\text{its sieved bound}),$$
and for at most one sign, since $r_q\not\equiv-r_q\pmod{q^2}$. How many $q$ can be matched at once? A matched $q$ has $b\equiv\pm r_q\pmod{q^2}$, so $q^2\mid b^2+1$. These moduli are coprime, so
$$\prod_{q\ \text{matched}}q^2\ \Big|\ b^2+1\le N^2+1,\qquad\text{hence}\qquad \sum_{q\ \text{matched}}\ln q\ \le\ \mathrm{cap}:=\tfrac12\ln(N^2+1)\ \bigl(-\ln 5\ \text{if } b\in C_7\cup C_{18}\bigr),$$
where the extra $\ln 5$ appears because $25\mid b^2+1$ for $b\in C_7\cup C_{18}$. So the matched primes form a feasible set for a knapsack with weights $\ln q$ and capacity about $\ln N$. Since $\sum_{\rm matched}\ln q\le\mathrm{cap}$, for every $\nu\ge0$ (this is LP weak duality, though it is immediate)
$$\sum_{q\ \text{matched}}\mathrm{gain}_q\ \le\ \nu\cdot\mathrm{cap}+\sum_{13\le q\le Q_s}\max\bigl(0,\ \mathrm{gain}_q-\nu\ln q\bigr).$$
At $N=4.2\cdot10^7$ the capacity is $\approx17.55$, so at most five primes can be matched at once ($13\cdot17\cdot29\cdot37\cdot41\approx 9.7\cdot10^6$, while a sixth prime would push the product past $5\cdot10^8>N$), and at most four for a pivot in $C_7\cup C_{18}$. As $N\to\infty$ all of them can be, and our limiting densities count every gain.

**Lemma 5 (sieving $O$ by a pivot).** Let $G^\tau(b)$ be the sum of the Lemma 1 bounds for the classes $K_{q,\varepsilon,s}$ with $q\le Q_s$, the knapsack term above and the tail of Lemma 4. Then
$$\#\{x\in O^\tau:\ \mu(xb+1)=0\}\le G^\tau(b)+\mathrm{Gl}(b).$$

*Proof.* By Lemma 3, every $x\in O^\tau$ lies in some $K_{q,\varepsilon,s}$ with $q\le Q_s$ or in the set of Lemma 4. Apply Lemma 1 with the statuses above to each class, the knapsack bound to the matched ones, and §3.2 for $p>P_1$. $\square$

For an even pivot in $O$ sieving all of $O$, $G/N\to 0.0154$, whereas the crude cover of $O$ by its $46$ classes per $q$ has density $0.0254$. In the worst case (II.a below), bounding $A^*$ by $|O|$ alone would give a limiting density of about $0.038$ and leave a margin of only about $0.002N$. Self-sieving gives $0.0282$, a margin of $\approx0.0118N$. It is also what makes $4.2\cdot10^7$ reachable: in the binding case there (§5) it saves about $0.005N$, more than three times the remaining slack.

## 5. The case analysis

**Two elementary facts** decide the statuses of $2$ and $5$:

- For $x\in C_7$, we have $25\mid xb+1$ only if $b\in C_7$, and likewise for $C_{18}$.
- If $x$ and $b$ are odd, then $4\mid xb+1$ if and only if $x\equiv-b\pmod 4$.

**Odd elements mod 100.** The odd elements of $C_7$ are those $\equiv 7 \pmod{50}$, that is, $\equiv 7$ or $57 \pmod{100}$. The odd elements of $C_{18}$ are those $\equiv 43\pmod{50}$, that is, $\equiv 43$ or $93\pmod{100}$. Modulo $4$, $7\equiv43\equiv3$ and $57\equiv93\equiv1$.

**The mod-100 observation.**
$$7\cdot43+1=302,\qquad 57\cdot93+1=5302,$$
and both are $\equiv 2\pmod{100}$, so $\equiv 2 \pmod 4$ and $\equiv 2\pmod{25}$. Hence a pair $x\equiv 7$, $y\equiv 43\pmod{100}$ (or $57$, $93$) can only be compatible through a prime square other than $4$ and $25$.

For such a pair, either one class misses $A$, and the pair contributes at most $\lceil N/100\rceil$, or each class is sieved by an element of the other with $p\notin\{2,5\}$, leaving about $2\cdot 0.158\cdot N/100$ plus two global terms. So the pair contributes at most
$$\Pi=\max\bigl(\lceil N/100\rceil,\ 2\Lambda_{100}+2\,\mathrm{Gl}\bigr)\sim N/100.$$
The cross pairs are different: $7\cdot93+1=652=4\cdot163$, so $x\equiv 7$ and $y\equiv 93$ can be compatible through $4$. This is why the classes are paired as above.

Likewise, an odd pivot $b\in O$ with $b\equiv1\pmod 4$ makes $xb+1\equiv2\pmod 4$ for $x\equiv57,93\pmod{100}$. It therefore sieves those two classes without $2$ or $5$, while $7$ and $43$ form a pair. If $b\equiv 3\pmod 4$, the roles are reversed.

**The cases.** Below, $c\in A_7$, $d\in A_{18}$ and $b\in A^*$. "All odd" allows the set to be empty unless stated otherwise. The cases are exhaustive:

- If $A^*=\emptyset$, both $A_7$ and $A_{18}$ are nonempty. The cases are then split by whether each contains an even element.
- Otherwise, either $A^*$ has an even element, or it is all odd and the same split applies.

Each bound is a sum of Lemma 1 and Lemma 5 bounds plus one $\mathrm{Gl}$ per pivot, with worst-case statuses. We need $\mathrm{UB}<(N+19)/25$, that is, $\mathrm{UB}/N$ essentially below $0.04$.

| Case | Situation | Who sieves what | $\mathrm{UB}/N$ at $4.2\cdot10^7$ | limit $D$ |
|---|---|---|---|---|
| I.1 | $A^*=\emptyset$; even $c$ and even $d$ | $d$ sieves $C_7$; $c$ sieves $C_{18}$ | 0.025397 | 0.012893 |
| I.2, I.3 | $A^*=\emptyset$; even $c$; $A_{18}\ne\emptyset$ all odd (and mirror) | $c$ sieves $A_{18}\subseteq 43\ (50)$; an odd $d$ sieves $C_7$ | 0.029552 | 0.018318 |
| I.4 | $A^*=\emptyset$; $A_7\cup A_{18}$ all odd | pairs $(7,43)$ and $(57,93)$ mod $100$, bounded by $\Pi$ | 0.026753 | 0.020000 |
| II.a | $A^*$ has an even element $b$ | $b$ sieves $C_7$, $C_{18}$ and $O$ | 0.036817 | **0.028195** |
| II.b.1 | $A^*\ne\emptyset$ all odd; $A_7\cup A_{18}$ all odd | $b$ sieves two mod-$100$ classes and $O$; the other two form a pair | 0.033701 | 0.023608 |
| II.b.2a | $A^*\ne\emptyset$ all odd; even $c$, even $d$ | $d$ sieves $C_7$; $c$ sieves $C_{18}$; $c$ or $d$ sieves $O$ | 0.034297 | 0.020400 |
| II.b.2b, 3b | $A^*\ne\emptyset$ all odd; even $c$; $A_{18}$ all odd (and mirror) | $c$ sieves $A_{18}\subseteq43\ (50)$; $b$ sieves $C_7$; $c$ or $b$ sieves $O$ | **0.038451** | 0.025825 |

Wherever $O$ can be sieved by either of two pivots, or bounded crudely, we take the minimum of the available bounds. The maximum number of pivots in any case is four (I.4); II.b.1 uses three, and the other cases at most two.

**Where the budget goes.**

- *Limit.* The worst limiting case is II.a:
$$D_{\max}=2(0.0063191)+0.0154294+0.0001270=0.0281946,$$
made up of $C_7$, $C_{18}$, $O$ and one global main term.
- *At $4.2\cdot10^7$.* The binding case is II.b.2b, which has two pivots and an expensive odd-pivot class. Its components, in units of $N$, are:
  - $\Lambda_{50}=0.003162$;
  - $\Lambda'_{25}=0.014906$ (here the odd pivot $b$ can still use $4$);
  - two global terms, $0.006377+0.005106$;
  - for $A^*$, the minimum of $0.013954$ (crude bound on the odd part of $O$), $0.008899$ ($c$ sieves $O$) and $0.012054$ ($b$ sieves $O$).

  The total is $0.038451$.

## 6. Certification for $4.2\cdot10^7\le N\le 1.14\cdot10^{14}$, and the tail

**Interval principle.** Let $N_{lo}\le N\le N_{hi}$, and let $A\subseteq[N]$ be admissible.

- If $A\subseteq C_7$ or $A\subseteq C_{18}$, then $|A|\le t(N)$ directly.
- Otherwise $A$ is also an admissible subset of $[N_{hi}]$ and is contained in neither class there, so the case analysis *at $N_{hi}$* gives $|A|\le \mathrm{UB}(N_{hi})$.

Since $t(N)\ge (N-6)/25$, the single exact inequality
$$\mathrm{UB}(N_{hi})<\frac{N_{lo}+19}{25}\le t(N)+1$$
certifies every $N$ in the interval, and no monotonicity of $\mathrm{UB}$ in $N$ is needed.

**The grid.** Main terms are exact rationals and all ceilings are computed in integers. The few real numbers that enter are rounded outward and then converted to rationals: logarithms in the Pell sum and the knapsack weights, and Rosser–Schoenfeld values of $\pi(x)$ above $2\cdot10^8$ (below that, $\pi$ is exact, and every prime sum on the grid is evaluated exactly). $154$ intervals, chained with an adaptive ratio and no retries, cover $[4.2\cdot10^7,\,1.1417\cdot10^{14}]$. The minimum slack is $8.91\cdot10^{-4}N_{lo}$, which is at least $8.7\cdot10^{-4}N$, on the first interval $[4.2\cdot10^7,\,4.284\cdot10^7]$ (case II.b.2b). The binding case is II.b.2b up to about $8.7\cdot10^8$ and II.a beyond.

**The tail $N\ge10^{14}$.**

- **Limiting densities.** For each case, $D_{\rm case}$ is the coefficient of $N$ in the case bound once ceilings, inclusion–exclusion errors and the non-main parts of the global and $O$-tail terms are dropped and all matched gains are switched on. These are computed in exact rational arithmetic, with the two prime sums that enter, $\sum_{p>1000}p^{-2}$ and $2\sum_{q>400,\,q\equiv1(4)}q^{-2}$, replaced by rigorous rational upper bounds. Then $\max D=0.0281946$ (case II.a).
- **Error term.** In every case
$$\mathrm{UB}_{\rm case}(N)\le D_{\rm case}N+c_0+4\bigl(\pi(Y)+E^*(N,Y)\bigr)+2\pi(Y_q)+\mathrm{Pell}(N,Y_q),$$
where $E^*(N,Y)$ is the bound of Lemma 2 with $\rho^*$ maximised over $b\le N$, and $c_0=78{,}658$. This constant collects the inclusion–exclusion errors and ceilings over at most $348$ classes, and the factor $4$ reflects that no case uses more than $4$ pivots.
- **Doubling grid.** Set $N_i=2^i\cdot10^{14}$ and fix $Y=Y_q=\lceil N_i^{0.8}\rceil$ on $[N_i,2N_i]$. The exponent $0.8$ is a convenient, non-optimised choice: the bound on $\rho^*$ is already $8192$ at $N=10^{14}$, so the switching term $\rho^*N^2/Y^2$ needs $Y$ well above $N^{2/3}$; at $N=10^{14}$ it and the $\pi(Y)$ terms are each about $5\cdot10^{-4}N$. Every error term is nondecreasing in $N$ at fixed $Y$, so with $\mathrm{Err}_i$ the error term evaluated at $2N_i$ it suffices to check
$$D_{\max}+1.01\cdot\mathrm{Err}_i/N_i<\tfrac1{25}-\tfrac1{N_i}$$
for each $i$, which we do up to $10^{300}$. This check is in floating point. The factor $1.01$ covers rounding in $\mathrm{Err}_i$, whose relative error in double precision is many orders of magnitude below $1\%$, and $D_{\max}$ itself is computed exactly. The worst value is $0.02919$, at $i=0$.
- **Beyond $10^{300}$.** Robin's bound $\omega(n)\le1.38402\ln n/\ln\ln n$ and Rosser–Schoenfeld give a closed-form majorant for $\mathrm{Err}(N)/N$. Each of its terms decreases, and the majorant is below $10^{-61}$ at $10^{300}$.

## 7. Part B: $N\le 2^{32}-1$ by an online colouring

**The graph.** Let $G_N$ have vertex set $V_N$, with $a\sim b$ when $\mu(ab+1)=0$. Admissible sets are exactly the cliques, so $|A|\le\omega(G_N)\le\chi(G_N)$. Since $C_7\cap[N]$ is itself a clique of size $t(N)$, a colouring with $t(N)$ colours is optimal, and each of its classes contains exactly one element of $C_7$. We therefore index the classes by the *anchors* $c_k=25k-18$. Also, $t$ increases only at $n\equiv7\pmod{25}$, and all such $n$ are vertices. So it suffices to colour $G_v$ with $t(v)$ colours for every vertex $v\le 2^{32}-1$.

**Why the colouring must be allowed to change.** We need a colouring with $t(N)$ colours for *every* $N$, not just one large $N$. A colouring in which vertices never change colour already fails at $n=43$. The vertices up to $43$ are $7,18,32,38,41,43$, with $t(18)=1$ and $t(32)=t(43)=2$, so $7$ and $18$ are forced into class $1$ and $32$ opens class $2$. But $43\cdot18+1=775=5^2\cdot31$ and $43\cdot32+1=1377=3^4\cdot17$, so $43$ fits in neither class, wherever $38$ and $41$ were placed, unless an earlier vertex moves. In the computed run, $38$ joined class $1$ and $41$ joined class $2$ ($41\cdot7+1=288=2^5\cdot3^2$ rules out class $1$). When $43$ arrives, $18$ moves to class $2$ ($18\cdot32+1=577$ and $18\cdot41+1=739$ are prime) and $43$ joins class $1$ ($43\cdot7+1=302=2\cdot151$, $43\cdot38+1=1635=3\cdot5\cdot109$). So the certificate is a **log of moves**, not a fixed colouring.

**Tight moments.** $C_{18}$ is also a clique, so right after $25k-7$ arrives each of the $k=t$ classes holds exactly one element of $C_{18}$. These are the hard moments.

**Statistics.** About a third ($34.5\%$) of $C_{18}$ insertions needed a local repair, almost always a swap; no repair searched deeper than one level, and no class ever had more than $8$ members.

**Certification.** The program that builds the colouring certifies each same-class pair squarefree ($uv+1<2^{64}$) with a fast routine: trial division, deterministic Miller–Rabin, an exact square test and Pollard–Brent. It need not be trusted, since the checker re-proves every claim. Uncertain cases are treated as edges, which is sound; none occurred among the $578{,}428{,}609$ certificates computed. (This exceeds the $536{,}466{,}170$ replay re-checks below because it includes the $2{,}245{,}204$ conflicts found and pairs from placements tried and abandoned; the replay re-checks only pairs created by the moves made.)

**Independent replay.** A separate checker, sharing no code with the program that builds the colouring, builds $V_N$ with its own sieve and tests squarefreeness by trial division up to the cube root of the cofactor followed by an exact square test. It replays the $517{,}358{,}696$-record log from $n=1$ and checks, after each vertex, that vertices arrive in order with none skipped, that there are exactly $t(n)$ classes, that anchors never move, and that every same-class pair created in the step is squarefree. That is $536{,}466{,}170$ pair re-checks, and the replay passes. Soundness follows by induction on the steps: the state after each vertex's step is a proper $t(n)$-colouring of $G_n$. The final colourings at $4.2\cdot10^7$, $10^8$, $10^9$ and $2^{32}-1$ were also checked in full.

**Further checks.** Two independently written sieves agree on $V_N$ for $N=2^{32}-1$ ($|V_N|=451{,}653{,}474$). An exact clique-number induction, using $\omega(G_n)\le\max(\omega(G_{n-1}),\,1+\omega(\text{earlier neighbours of } n))$, confirms $\omega(G_n)\le t(n)$ for all $n\le3\cdot10^5$ with no colouring involved. For $N=10^5$ and $3\cdot10^5$ the colouring was also checked against an edge-sieved graph that uses no squarefree test, and two negative controls (a truncated log and a planted bad pair) are rejected.

Part B reaches $4.29\cdot10^9$, about $100$ times the analytic threshold.

## 8. What is new, how it compares, and limitations

**What is standard.** The case split by whether $A$ meets $O$ goes back to Sawhney. Switching is classical, with explicit forms in Sothanaphan, Li and Pitchford. The negative-Pell tail appears in Sothanaphan and Pitchford, and colouring certificates for an initial range are used by Li and Pitchford.

**What may be new.** The three devices listed at the end of §0. We have not checked whether equivalent devices appear in the work of Sawhney, Sothanaphan, Li or Pitchford.

**Comparison.** The entries for other work are taken from its public documentation; we have not audited it.

| Work | Finite computation | Interval-wise explicit | Uniform explicit |
|---|---|---|---|
| Sothanaphan (GPT-5.4) | — | — | $N\ge2.64\cdot10^{17}$ |
| Pitchford / OpenAI Codex | colourings, $N\le100{,}000{,}006$ | $[10^8,\ 2.64\cdot10^{17}]$ (certificate chain) | via Sothanaphan |
| Li (Lean 4) | colouring, $N\le5\cdot10^6$ | $[5\cdot10^6,\ 5\cdot10^{11}]$ | $N\ge5\cdot10^{11}$ |
| This work | colouring, $N\le2^{32}-1$ | $[4.2\cdot10^7,\ 1.14\cdot10^{14}]$ | $N\ge10^{14}$ |

- Li's claim comes with a **Lean 4 formalisation**; if it builds as documented, that is a much stronger guarantee than anything here. Li's analytic ranges also start lower on both counts ($5\cdot10^6$ and $5\cdot10^{11}$, against our $4.2\cdot10^7$ and $10^{14}$). Li used ChatGPT 5.5/5.6 and an AI-assisted proof-engineering system (erdosproblems.com proof claim #204, August 2026).

Any merits of the present route are modest: it is independent, it is short (five lemmas and nine cases), its margins are large (slack at least $8.7\cdot10^{-4}N$, tending to $\approx0.0118N$), and its finite part has a small trusted checker and reruns on one machine.

**Limitations.**

- The proof is computer-assisted and not refereed, and it has no formalisation.
- The audit was done by another instance of the same model, so its errors may be correlated with ours. It re-derived the lemmas and the case split, re-implemented the two binding cases (II.a, II.b.2b) and brute-force tested the lemma bounds; it found no fatal issue, and its four minor issues have been fixed.
- In Part A, floating point without outward rounding is used only in the tail check for $N\ge10^{14}$ (§6), where even $\mathrm{Err}_0/N_0\approx10^{-3}$ is under a tenth of the available slack.
- The trusted base is about $550$ lines of Python for Part A, and about $350$ lines of C for the Part B checker, plus the compiler and hardware.

**Reproducing.** Part A's quick re-check takes under two minutes, and Part B up to $10^8$ (already past the analytic threshold) under a minute. Commands, full-run timings and run logs, including one correctly rejected replay caused by a since-fixed checkpoint bug, are in the repository README.

Corrections, especially to the attributions above, are welcome.

### References

- [Er92b] P. Erdős, as cited at p. 239 on erdosproblems.com/848; T. F. Bloom, *Erdős Problem #848*, erdosproblems.com/848, and its forum thread.
- van Doorn, $|A|\le(0.108\ldots+o(1))N$, and Weisenberg, $\approx0.105N$: density bounds as reported on erdosproblems.com/848.
- M. Sawhney, note on Problem 848 (2025).
- N. Sothanaphan (with GPT-5.4), explicit threshold $N_0=2.64\cdot10^{17}$ (March 2026).
- A. C. Li, Lean 4 formalisation and manuscript; erdosproblems.com proof claim #204 (August 2026).
- I. Pitchford and OpenAI Codex, all-$N$ certificate chain, Zenodo doi:10.5281/zenodo.21647629 (July 2026).
- T. Estermann, Math. Ann. 105 (1931); C. Hooley, Bull. LMS 7 (1975); D. R. Heath-Brown, Acta Arith. 155 (2012).
- J. B. Rosser and L. Schoenfeld, Illinois J. Math. 6 (1962); G. Robin, Acta Arith. 42 (1983).

# Exact reference sizes and success probabilities

In the [one-use model](model.md), $`N`$ sensor qubits contain an unknown $`k`$-element subset subjected to a common $`Z`$ rotation. Assume $`1\le k\le N/2`$. The angle varies in a nonempty open interval, the receiver may abstain, and every conclusive subset must be correct throughout that interval. The preparation and readout conserve total excitation number. Only the $`r`$ bypass qubits can carry charge outside the sensors during the interrogation.

For a receiver, define $`\eta_S=\langle Z_S\psi,E_SZ_S\psi\rangle`$ on a pure input $`\psi`$, where $`E_S`$ is its conclusive effect and $`Z_S=\prod_{i\in S}Z_i`$. The [moment isolation argument](model.md) gives

```math
\Pr\{\text{output }S\mid S,\theta\}
=\sin^{2k}(\theta/2)\,\eta_S.
```

For a mixed input the coefficient is averaged over its pure components. Let $`\eta^*_{N,k}(r)`$ be the largest worst-subset coefficient. The same optimal preparation and receiver attain the corresponding pointwise success for every angle in the interval.

```math
\eta^*_{N,k}(r)=0\quad(r<k).
```

```math
\eta^*_{N,k}(k)=\frac{1}{\binom Nk}.
```

```math
\eta^*_{N,1}(r)=\min\{1,r^2/N\}\quad(r\geq0).
```

Thus $`k`$ bypass qubits are necessary and sufficient for a positive probability at every nonidentity angle. For one fault, $`\lceil\sqrt N\rceil`$ bypass qubits suffice to attain the unrestricted success coefficient of one.

## A bound from total sensor charge

Write $`q=\sum_{i=1}^N(1-Z_i)/2`$ for sensor charge. The operator $`\sum_{|A|=j}Z_A`$ is the polynomial $`K_j(q)`$, where

```math
K_j(q)=\sum_{\ell=0}^{j}(-2)^\ell
\binom{N-\ell}{j-\ell}\binom q\ell.
```

This follows by expanding each product of $`j`$ factors $`1-2x_i`$ on a bit string $`x`$, then counting the extensions of each $`\ell`$-subset. The polynomials $`K_0,\ldots,K_{k-1}`$ span all polynomials of degree at most $`k-1`$. The leading coefficient of $`K_k`$ is $`(-2)^k/k!`$.

Moment isolation says that $`E_S`$ annihilates $`Z_A\psi`$ for every $`|A|\le k`$ except $`A=S`$. For every polynomial $`g`$ of degree at most $`k-1`$, the vector $`f=(K_k(q)-g(q))\psi`$ therefore satisfies $`E_S^{1/2}f=E_S^{1/2}Z_S\psi`$ for each $`S`$. Since the conclusive effects sum to at most the identity,

```math
\sum_{|S|=k}\eta_S\le \langle\psi,\lvert K_k(q)-g(q)\rvert^2\psi\rangle.
```

In each total-charge sector the sensor charges lie in a consecutive interval of width at most $`r`$. For $`r<k`$, interpolate $`K_k`$ on the sensor charges in that interval by a polynomial $`g`$ of degree at most $`r`$. The bound is then zero. For $`r=k`$, choose a consecutive interval $`a,\ldots,a+k`$ containing the sensor charge support. The interpolation polynomial $`F`$ with values $`F(a+j)=(-1)^j`$ has $`k`$th forward difference $`(-2)^k`$. Thus $`K_k-F`$ has degree at most $`k-1`$, and taking it as $`g`$ gives $`\sum_S\eta_S\le 1`$.

These sector bounds also apply to mixed inputs and classical randomizations. For number-conserving effects and interrogation, the probabilities depend only on the input blocks of fixed total charge. A zero wrong-label probability is a sum of nonnegative probabilities, so it is zero on every pure component of each occupied block. The bound on the mean coefficient over all subsets holds componentwise and hence for their mixture. The smallest subset coefficient is at most that mean.

## A receiver with exactly $`k`$ bypass qubits

For a sensor string with occupied set $`X\subseteq[N]`$ and $`|X|\le k`$, let $`|b_X\rangle`$ be that string together with a known bypass string of weight $`k-|X|`$. These vectors are orthonormal and all have total charge $`k`$. Prepare

```math
w_j=2^{-k}\binom{k}{j},\qquad
p_X=\frac{w_{|X|}}{\binom N{|X|}},\qquad
|\psi\rangle=\sum_{|X|\le k}\sqrt{p_X}\,|b_X\rangle.
```

For each $`|S|=k`$, define the reciprocal vector

```math
|u_S\rangle=(-2)^{-k}\sum_{X\subseteq S}
\frac{(-1)^{k-|X|}}{\sqrt{p_X}}\,|b_X\rangle.
```

The signed sum in $`\langle u_S,Z_A\psi\rangle`$ is a finite difference on the Boolean cube of subsets of $`S`$. Using $`Z_A|b_X\rangle=(-1)^{|A\cap X|}|b_X\rangle`$, it factors over the $`k`$ positions in $`S`$. It vanishes if one of those positions is outside $`A`$, and equals $`(-2)^k`$ before the prefactor if all are in $`A`$. Hence $`\langle u_S,Z_A\psi\rangle=\delta_{S,A}`$ for all $`|A|\le k`$.

The Gram matrix of the $`u_S`$ has nonnegative entries. Its row sum, obtained by counting the $`k`$-subsets containing each $`j`$-subset of $`S`$, is

```math
4^{-k}\sum_{j=0}^k\binom{k}{j}\binom{N-j}{k-j}
\frac{\binom Nj}{w_j}
=4^{-k}\binom Nk\sum_{j=0}^k\frac{\binom{k}{j}^2}{w_j}
=\binom Nk.
```

The constant row sum is the largest eigenvalue, since it is both an eigenvalue and an upper bound from the maximum absolute row sum. It is also the norm of $`\sum_S|u_S\rangle\langle u_S|`$. The measurement takes

```math
E_S=\frac{|u_S\rangle\langle u_S|}{\binom Nk},\qquad
E_{?}=I-\sum_{|S|=k}E_S.
```

Every effect preserves total charge. The reciprocal identities isolate the term with the correct subset in the rotation expansion. Each subset has coefficient $`1/\binom Nk`$, which attains the preceding bound.

## Exact tradeoff for one fault

For $`k=1`$, $`K_1=N-2q`$. On each pure component in a fixed total-charge sector, choose $`g=\langle K_1\rangle`$. The charge bound gives $`\sum_i\eta_i\le 4\mathrm{Var}(q)\le r^2`$, using the variance bound for a variable in an interval of width at most $`r`$. Each $`\eta_i`$ is also at most one. The preceding mixture argument gives the claimed upper bound.

For the construction, let $`q`$ follow a chosen distribution on consecutive allowed charges. Within each charge use the uniform superposition of sensor strings, paired with one known bypass string of complementary charge. The resulting pure preparation is permutation symmetric on the sensors. Define

```math
v_i=(I-|\psi\rangle\langle\psi|)Z_i|\psi\rangle,\qquad
G_{ij}=\langle v_i,v_j\rangle.
```

The matrix $`G`$ is the covariance of the exchangeable signs $`1-2x_i`$ under the preparation's sensor-string probabilities. It is constant on the diagonal and separately constant off the diagonal. Its eigenvalue along the all-ones vector and its eigenvalue on the orthogonal complement are

```math
\lambda_{\parallel}=\frac{4\mathrm{Var}(q)}{N},\qquad
\lambda_{\perp}=\frac{4\,\mathbb E[q(N-q)]}{N(N-1)}.
```

The first identity follows from summing the signs. The second follows by evaluating the covariance on the difference of two coordinate vectors and counting oppositely occupied pairs at each $`q`$. When the smaller eigenvalue $`\eta`$ is positive, take the reciprocal vectors $`u_i`$ in the span of the $`v_i`$, so that $`\langle u_i,v_j\rangle=\delta_{ij}`$. Their frame operator has norm $`1/\eta`$. The conclusive effects $`\eta\,|u_i\rangle\langle u_i|`$, with the remaining effect for abstention, attain coefficient $`\eta`$ at every site. They are orthogonal to $`\psi`$ and confined to the chosen total-charge sector.

### When $`1\le r\le\sqrt N`$

Choose integer endpoints $`a,a+r\in[0,N]`$ with midpoint $`N/2+\delta`$ and $`|\delta|\le 1/2`$. Put probability one half on each. Then

```math
\lambda_{\parallel}=\frac{r^2}{N},\qquad
\lambda_{\perp}=\frac{N^2-r^2-4\delta^2}{N(N-1)}.
```

The transverse eigenvalue is at least $`r^2/N`$ exactly when $`N(N-r^2)\ge 4\delta^2`$. If $`r^2<N`$, the integer difference is at least one. If $`r^2=N`$, $`r`$ and $`N`$ have the same parity, so the midpoint can be exactly $`N/2`$. The charge sector with total charge $`a+r`$ supports both endpoints and yields the desired gain. For $`r=0`$, always abstaining attains the zero optimum.

### When $`r\ge\sqrt N`$

Put $`s=\min(r,N)`$ and choose an interval $`[a,a+s]`$ centered to within $`1/2`$ of $`\mu=N/2`$. Its endpoint distribution with mean $`\mu`$ has variance

```math
(\mu-a)(a+s-\mu)=s^2/4-\delta^2\ge N/4.
```

Here $`\delta`$ is the displacement of the interval midpoint from $`\mu`$. The inequality again follows from integrality if $`s^2>N`$ and from matching parity if $`s^2=N`$. The distribution on the nearest central integers with mean $`\mu`$ has variance at most $`1/4`$. Both distributions fit in the same charge interval. Mixing their probability weights gives a distribution with mean $`\mu`$ and variance exactly $`N/4`$, used coherently in the preparation at total charge $`a+s`$. It makes both eigenvalues equal to one, so the reciprocal measurement becomes projective on the span of the $`v_i`$.

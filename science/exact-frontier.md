# Optimal receivers at finite reference width

The [model](model.md) fixes one simultaneous interrogation of $`N`$ sensors, $`r`$ idle reference qubits and $`1\le k\le N/2`$ unknown affected sites. Its moment-isolation lemma gives a pointwise success $`\eta_S\sin^{2k}(\theta/2)`$. The conclusive effect for $`S`$ annihilates every lower-degree moment and each degree-$`k`$ moment except its own.

Sensor permutations preserve the resources and are transitive on the true labels. A uniform random permutation converts an average coefficient into the same coefficient for each label. For a fixed strategy, random thinning of larger conclusive coefficients can also equalize them at their minimum. This section uses a common coefficient $`\eta`$ to characterize the optimum.

## A finite semidefinite program and its receiver

Index a Gram matrix by all sets of size at most $`k`$. Define

```math
G_{A,B}=\langle\psi|Z_{A\mathbin{\triangle}B}|\psi\rangle,
\qquad
P_{A,B}=\mathbf 1_{A=B,\ |A|=k}.
```

For a common success coefficient $`\eta\ge0`$, feasibility is equivalent to
$`G-\eta P\succeq0`$. To see necessity, use the contraction
$`|v\rangle\mapsto\bigoplus_{|S|=k}E_S^{1/2}|v\rangle`$. The squared norm of the
image of $`\sum_A x_A Z_A|\psi\rangle`$ is $`\eta\sum_{|S|=k}|x_S|^2`$ by the
annihilation conditions. The squared norm before applying the contraction is
$`x^*Gx`$.

For sufficiency, the inequality makes the map $`L`$ below well-defined on the span
of the moment vectors even if they are dependent.

```math
L\left(\sum_A x_A Z_A|\psi\rangle\right)=\sqrt{\eta}\sum_{|S|=k}x_S|S\rangle.
```

It is a contraction. Extend it by zero to the orthogonal complement in sector $`Q`$.
Take $`E_S=L^*|S\rangle\langle S|L`$ and put the remaining weight in the
inconclusive effect. The pulse expansion gives exact no-error readout and the
stated coefficient. An isometry with components $`L`$ and
$`(I-L^*L)^{1/2}`$ realizes the receiver with a neutral outcome register within
sector $`Q`$. For the conclusive component, leave the original charged system in
one fixed vector of sector $`Q`$ and put the label in the neutral register. Use
a distinct neutral flag for the inconclusive component. This construction
specifies the measurement as a finite linear map.

The Gram matrix only depends on the distribution of the sensor bit string in the
computational basis. Averaging that distribution over permutations preserves
$`G-\eta P\succeq0`$ because the inequality is linear in $`G`$ and $`P`$ is invariant.
The averaged distribution is specified by probabilities $`p_w`$ for sensor charge
$`w`$. Given $`w`$, all bit strings of weight $`w`$ have equal probability. In sector
$`Q`$, the charge interval is

```math
I_Q=[\max(0,Q-r),\min(N,Q)]\cap\mathbb Z.
```

Let $`|D_{N,w}\rangle`$ be the normalized equal-amplitude state on sensor bit
strings of weight $`w`$. Choose one normalized idle state $`|R_u\rangle`$ of charge
$`u`$ for each $`0\le u\le r`$. Then the single pure state

```math
|\psi(p,Q)\rangle=\sum_{w\in I_Q}\sqrt{p_w}\,|D_{N,w}\rangle|R_{Q-w}\rangle
```

realizes the averaged Gram matrix in the same total sector and uses no additional
charged registers. This removes any need for a coherent phase reference external to
the probe and idle qubits.

Define the polynomials $`K_\ell`$ by

```math
K_\ell(w)=\sum_{h=0}^{\ell}(-1)^h\binom wh\binom{N-w}{\ell-h}.
```

Use the convention that a binomial coefficient with an impossible lower index
is zero. For each fixed set of size $`\ell`$, the expected eigenvalue of its $`Z`$
product on the weight-$`w`$ slice is $`K_\ell(w)/\binom N\ell`$. This follows
by counting how many sets of size $`\ell`$ intersect a fixed weight-$`w`$ string in
$`h`$ ones. For each $`Q=0,\ldots,N+r`$, solve

```math
\begin{aligned}
\text{maximize}\quad & \eta\\
\text{subject to}\quad & p_w\ge0,\quad\sum_{w\in I_Q}p_w=1,\quad \eta\ge0,\\
&G(p)-\eta P\succeq0,\\
&G(p)_{A,B}=\sum_{w\in I_Q}p_w\,
 \frac{K_{|A\mathbin{\triangle}B|}(w)}{\binom N{|A\mathbin{\triangle}B|}}.
\end{aligned}
```

Take the largest sector optimum to obtain $`\eta^*(N,k,r)`$. Classical mixtures
over $`Q`$ cannot improve this result. Their mean coefficient over a uniform true
label is a weighted sum of the separate sector mean coefficients. Their worst-label
coefficient is at most that mean. Optimizing over $`Q`$ corresponds to choosing the
charge sector during preparation. With a prescribed sector, use just its program.

## Reducing the Gram matrix by sensor permutations

Let $`V_s`$ be the space of vectors indexed by $`s`$-subsets. Define
$`U_s=B_{s+1,s}:V_s\to V_{s+1}`$ and $`D_s=U_{s-1}^*:V_s\to V_{s-1}`$, with
$`D_0=0`$. Counting the common supersets and common subsets of two $`s`$-sets gives

```math
D_{s+1}U_s-U_{s-1}D_s=(N-2s)I_{V_s}.
```

On the diagonal the two counts are $`N-s`$ and $`s`$. Off the diagonal they agree,
with both equal to one when the sets differ in one element and zero otherwise.
In particular, $`U_s`$ is injective for $`s<N/2`$ by taking the quadratic form
of this identity. Since $`B_{a,b}=U_{a-1}\cdots U_b/(a-b)!`$, this proves the full
column rank of $`B_{k,j}`$.

Define $`H_j=\ker D_j`$. For $`0\le j\le s\le k`$, the images $`B_{s,j}H_j`$
give orthogonal subspaces that fill $`V_s`$. The commutator gives the norm of the
lift and the lowering identities used to verify this decomposition. For $`h\in H_j`$,
suppress the level index on repeated applications of $`U`$ and $`D`$. Inductively,
$`DU^t h=t(N-2j-t+1)U^{t-1}h`$. With the factorial normalization of $`B`$, this
implies

```math
\begin{aligned}
\|B_{s,j}h\|^2&=\binom{N-2j}{s-j}\|h\|^2,\\
B_{s,i}^*B_{s,j}h
&=\begin{cases}
0 & i<j,\\
\binom{N-i-j}{s-i}B_{i,j}h & j\le i\le s.
\end{cases}
\end{aligned}
```

The zero in this identity proves orthogonality. Injectivity of the raising map
shows that $`\dim H_j=\binom Nj-\binom N{j-1}`$, with $`\binom N{-1}=0`$.
These dimensions sum to $`\dim V_s`$ as $`j`$ runs from 0 to $`s`$.

The $`(a,b)`$ block of $`G(p)`$ has entry $`f_{a,b}(t)`$ at intersection size
$`t=|A\cap B|`$. Here $`f_{a,b}(t)`$ is the Gram expression above evaluated at
$`\ell=a+b-2t`$. The possible $`t`$ range is 0 through $`\min(a,b)`$ because
$`a+b\le N`$. Newton's difference formula gives

```math
G(p)_{a,b}=\sum_{i=0}^{\min(a,b)}c_{a,b,i}\,B_{a,i}B_{b,i}^*,
\qquad c_{a,b,i}=(\Delta^i f_{a,b})(0),\quad
(\Delta f)(t)=f(t+1)-f(t).
```

Indeed, $`B_{a,i}B_{b,i}^*`$ has entry $`\binom ti`$, and these binomial
polynomials form the Newton basis. Also
$`B_{a,i}B_{i,j}=\binom{a-j}{i-j}B_{a,j}`$ by counting the intermediate sets.
Combining these facts with the lowering identity shows that the normalized lifts
of one unit vector in $`H_j`$ are acted on through the matrix

```math
M_j[a,b]=
\sqrt{\frac{\binom{N-2j}{a-j}}{\binom{N-2j}{b-j}}}
\sum_{i=j}^{\min(a,b)}c_{a,b,i}\binom{N-i-j}{b-i}\binom{a-j}{i-j},
\qquad a,b=j,\ldots,k.
```

The same matrix occurs for each orthonormal direction in $`H_j`$. The operator
$`P`$ restricts to the projection onto the level-$`k`$ entry of each such block.
Thus the full Gram constraint is equivalent to $`k+1`$ simultaneous positive
semidefinite constraints $`M_j-\eta e_k e_k^*\succeq0`$, of sizes
$`k-j+1`$, where $`e_k`$ is the last coordinate of the block indexed by
$`j,\ldots,k`$.

## Sharp scaling at fixed trajectory size

The unit coefficient is an upper bound because $`G_{S,S}=1`$. The [radial inequality](model.md#the-radial-bound) bounds the mean coefficient by the squared norm of any polynomial of the sensor charge with leading coefficient $`(-2)^k/k!`$, divided by $`\binom Nk`$. The [exact endpoints](exact-endpoints.md) apply this inequality below and at the minimum useful width.

For $`r>0`$, choose a real interval of width $`r`$ with center $`c`$ containing
the allowed sensor charges. Take $`q`$ so that $`K_k(w)-q(w)`$ is a multiple of the
Chebyshev polynomial $`T_k(2(w-c)/r)`$ with the same leading coefficient as $`K_k`$.
Since $`T_k`$ has leading coefficient $`2^{k-1}`$ and modulus at most 1 on
$`[-1,1]`$, the total bound is

```math
\eta^*(N,k,r)\le
\min\!\left\{1,\frac{r^{2k}}{2^{2k-2}(k!)^2\binom Nk}\right\}.
```

### A modular receiver using the same pulse

The [one-location receiver](exact-endpoints.md#exact-tradeoff-for-one-fault) achieves $`\min(1,b^2/n)`$ with $`n`$ sensors and $`b`$ idle qubits. Its conclusive effects annihilate the unrotated probe, so an empty group always abstains.

For $`r\ge k`$, partition the $`N`$ sensors uniformly into $`k`$ labeled groups of
balanced sizes $`n_1,\ldots,n_k`$. Divide the $`r`$ idle qubits among the groups as
evenly as possible, with sizes $`r_i`$. Run the single-location receiver in each
group, with all groups exposed to the same one use of the physical pulse. The
preparations have fixed total charge within each group, and their sum is fixed
independently of the randomly chosen partition. Report the union of the group's
labels only when all groups return a location. If the true trajectory does not split
one per group, at least one group is empty and cannot report a location. If it
does split, the group outcomes are independent conditional on $`S,\theta`$ and the
partition. This gives a common guaranteed coefficient

```math
\eta_{\rm part}=\frac{\prod_i n_i}{\binom Nk}\prod_{i=1}^k\min(1,r_i^2/n_i).
```

The first factor is the probability a fixed $`k`$-set splits one per group. This is
the fraction of split sets in any fixed balanced partition, by permutation symmetry.
Here $`N/(2k)\le n_i\le 3N/(2k)`$ and $`r_i\ge r/(2k)`$. Together with
$`\binom Nk\le N^k/k!`$, these bound the first factor below by $`k!/(2k)^k`$
and each factor in the second product below by $`\min(1,r^2/(6kN))`$. The
converse gives the same order of upper bound since
$`\binom Nk\ge (N/2)^k/k!`$. Thus, uniformly over $`N\ge2k`$ and integers $`r\ge k`$,

```math
\eta^*(N,k,r)=\Theta_k\!\left(\min\{1,(r^2/N)^k\}\right).
```

## Recovering the unrestricted coefficient exactly

The preceding estimate guarantees a constant fraction of the unrestricted receiver
at width of order $`\sqrt N`$ for fixed $`k`$. The same order of width can also
recover the unrestricted coefficient exactly. If the first $`2k`$ moments of $`w`$
agree with those of $`\operatorname{Bin}(N,1/2)`$, the Gram matrix is the identity.
Indeed, $`K_\ell(w)`$ has degree $`\ell`$ and for the binomial distribution the bit
string is uniform. Its sign products for all nonempty sets have zero expectation.
The resulting charge-tagged probe therefore has orthonormal moment vectors up to
degree $`k`$. Direct projection onto those of degree $`k`$ attains $`\eta=1`$.

To see that the moments can be matched with integer support in a band of width
$`O_k(\sqrt N)`$, define $`y=(2w-N)/\sqrt N`$ and
$`v(y)=(y,y^2,\ldots,y^{2k})`$. Let $`m`$ be the mean of $`v(y)`$ for the standard
Gaussian. There exists $`L_k<\infty`$ such that $`m`$ is in the interior of the convex
hull of $`\{v(y):|y|\le L_k\}`$.

For this claim, take the monic degree-$`k+1`$ polynomial orthogonal to all polynomials
of degree at most $`k`$ in the Gaussian inner product. It has $`k+1`$ distinct real
roots. Otherwise the product of its sign-changing linear factors would have degree
at most $`k`$ and its inner product with the original polynomial would have a fixed,
nonzero sign. Division by the orthogonal polynomial followed by interpolation at its
roots yields an exact Gaussian quadrature for degree at most $`2k+1`$. Each quadrature
weight is positive because it is the Gaussian integral of the square of the
corresponding Lagrange basis polynomial, to which exactness applies. Choose $`L_k`$
to put all the roots in the interior of the interval. The quadrature represents $`m`$
as a positive mixture of the associated points $`v(y)`$. The convex hull has full
dimension because a polynomial of degree at most $`2k`$ constant on an interval
has all its nonconstant coefficients zero. If $`m`$ were on its boundary, a nonzero
supporting linear functional would give a nonnegative polynomial on the interval
vanishing at every one of the $`k+1`$ roots. Those are interior roots, each of even
multiplicity. This contradicts the polynomial's degree bound of $`2k`$.

Now restrict the support to the grid $`y=(2w-N)/\sqrt N`$ for integer $`0\le w\le N`$
with $`|y|\le L_k`$. For fixed $`k`$ the grid fills the interval as $`N`$ grows.
Continuity of $`v`$ makes the support functions of the grid convex hull converge
uniformly on unit directions to those of the full convex hull. A fixed smaller ball
around $`m`$ is therefore contained in every sufficiently fine grid convex hull,
by the separating-hyperplane criterion for convex containment. The moments of the
normalized binomial converge to $`m`$ up to degree $`2k`$. For instance, write its
normalization as a sum of $`N`$ independent symmetric signs divided by $`\sqrt N`$.
In an even moment the leading terms pair all indices, yielding the Gaussian moment,
and every term with fewer distinct indices is of smaller order. Odd moments vanish.
The normalized binomial moment vector eventually lies in the smaller ball, hence in
the grid hull. This gives the required nonnegative weights. They are obtained by a
finite linear feasibility problem imposing the binomial moments on the grid.

The grid spans sensor charges of width at most $`L_k\sqrt N`$. Use its maximum
charge as $`Q`$ and give each occupied sensor weight the corresponding idle charge.
For large $`N`$ this constructs the exact optimal coefficient one with
$`r\le \lceil L_k\sqrt N\rceil`$ idle qubits. Finitely many smaller arrays can use
the full binomial distribution with $`r=N`$, so the implied constant can be increased
to cover all $`N\ge2k`$. The converse gives an order-$`\sqrt N`$ lower bound on the
width needed for this coefficient at fixed $`k`$.

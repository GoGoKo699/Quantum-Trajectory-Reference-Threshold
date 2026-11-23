# Charge-conserving circuits and resource comparison

Each sensing round uses $`N`$ queried qubits. A bypass qubit passes around the queried devices and can participate coherently in preparation and readout. A neutral register has charge operator $`0`$ on its entire state space. Its coherent pointer states, classical records, and work bits supply degeneracy within the same total-charge sector. Preparation work bits can be reset before interrogation. Readout pointers can be introduced after interrogation.

The native charged operations used here are single-wire occupation phases and two-wire exchanges. An exchange acts as a two-level rotation on $`|10\rangle,|01\rangle`$ and leaves the zero- and two-excitation states alone. A charged wire may control a neutral bit through its occupation, and a neutral bit may control an exchange. Each operation commutes with total excitation. Local gates on neutral bits can compute and uncompute several controls. If controlled exchange is unavailable as a primitive, write $`X_{ij}=|10\rangle\langle01|+|01\rangle\langle10|`$ and let $`D`$ be a neutral-controlled occupation sign on wire $`i`$. The product $`e^{-i\alpha X_{ij}/2} D e^{i\alpha X_{ij}/2} D`$ is the controlled $`e^{-i\alpha X_{ij}}`$. Occupation phases change the exchange axis in its two-dimensional sector. These gates presume calibrated internal control of exchange phases. They require no operations that create coherence between different total excitation numbers.

## Local Bell benchmark

Prepare an independent pair $`(|10\rangle+|01\rangle)/\sqrt2`$ at each site. The first member goes through the queried device and the second bypasses it. The balanced pair exchange prepares this state from a single excitation. After interrogation the inverse exchange converts the symmetric and antisymmetric states to different occupation readings. A hit has antisymmetric probability $`\sin^2(\theta/2)`$ and an idle site has antisymmetric probability zero. Report the subset exactly when $`k`$ sites register the antisymmetric outcome. This implements the local anomaly receiver of [Llorens et al.](https://arxiv.org/abs/2605.04981) for this rotation family. It uses $`N`$ bypass qubits, $`N`$ pair exchanges at each of preparation and readout, and ideal parallel depth one at each end when all pair couplings are available. The conclusive probability is $`\sin^{2k}(\theta/2)`$.

## A linear-size receiver with one bypass qubit

Consider first a single hit among $`n\ge 2`$ sites, with one bypass wire $`b`$. In the one-excitation sector let $`|j\rangle`$ put the excitation on sensor $`j`$, let $`|b\rangle`$ put it on the bypass, and define

```math
 |s\rangle=\frac{1}{\sqrt n}\sum_{j=1}^n|j\rangle,
 \qquad |\psi\rangle=\frac{|b\rangle+|s\rangle}{\sqrt2},
 \qquad |d\rangle=\frac{|b\rangle-|s\rangle}{\sqrt2}.
```

Choose a product $`D`$ of $`n-1`$ sensor exchanges such that $`D|s\rangle=|1\rangle`$. Such a product follows by Givens elimination of the $`n-1`$ amplitudes on the other sensor rails. It works with a star centered on the first sensor, or with a sweep on adjacent sensors. Prepare $`|\psi\rangle`$ by a balanced exchange between $`b`$ and sensor 1 followed by $`D^\dagger`$.

At readout apply $`D`$ and the pair transform which maps $`(|b\rangle+|1\rangle)/\sqrt2`$ to $`|b\rangle`$ and $`(|b\rangle-|1\rangle)/\sqrt2`$ to $`|1\rangle`$. Set the acceptance amplitude of a neutral pointer to zero when the excitation is on $`b`$, to one on sensor 1, and to $`-1/\sqrt2`$ on any other sensor. At most $`n`$ occupation-controlled real pointer rotations implement these amplitudes, since only one wire is occupied. Apply $`D^\dagger`$ and measure the sensor occupations. The accepted occupied rail is the label. Reject the other pointer result.

To check this circuit, let $`P_{\rm ang}`$ project onto the sensor space perpendicular to $`|s\rangle`$. The accepted map $`K`$ kills $`|\psi\rangle`$, maps $`|d\rangle`$ to $`|s\rangle`$, and acts as $`-1/\sqrt2`$ on the angular subspace. Since $`Z_j|\psi\rangle=|\psi\rangle-\sqrt{2/n}|j\rangle`$, it follows that

```math
 K Z_j|\psi\rangle = \frac{|s\rangle}{n}+\frac{P_{\rm ang}|j\rangle}{\sqrt n}
 =\frac{|j\rangle}{\sqrt n}.
```

This gives exactly $`\sin^2(\theta/2)/n`$, with no wrong labels and a single neutral pointer bit. At $`n=2`$, the whole probe uses two exchanges and the readout uses three exchanges and at most two controlled pointer rotations. The sensor pair can use a common beam splitter for $`D`$ and $`D^\dagger`$. The one-excitation design deliberately uses the $`0,1`$ sensor-charge pair, which suffices for the exact one-bypass law even when $`n`$ is large.

## Compiling the one-hit optimum at larger bypass width

The same one-hit optimum, $`\eta=\min(1,r^2/n)`$, is available with a structured, though larger, circuit at general $`r\ge1`$. Use the optimal symmetric probe and its sensor-charge distribution from the one-hit result. It has at most four occupied sensor weights, with a specified normalized bypass basis word at each compensating charge. Only $`r_{\rm use}=\min(r,\lceil\sqrt n\rceil)`$ bypass wires need take part. In total charge $`Q`$ write the probe as

```math
 |\psi\rangle=\sum_w\sqrt{p_w}\,|D_{n,w}\rangle|b_{Q-w}\rangle,
 \qquad z=\langle\psi|Z_i|\psi\rangle,
 \qquad u_i=(Z_i-z)|\psi\rangle,
```

where $`|D_{n,w}\rangle`$ is the equal superposition of sensor strings of weight $`w`$. The Gram matrix of the $`u_i`$ has eigenvalue $`\lambda_0`$ on the uniform label vector and $`\lambda_1`$ on its orthogonal complement. Choose a real orthonormal basis $`c^{(a)}`$ of the label space. Take $`c^{(0)}`$ uniform and, for $`1\le a\le n-1`$, take $`c^{(a)}`$ with its first $`a`$ entries equal to $`1/\sqrt{a(a+1)}`$, its next entry $`-a/\sqrt{a(a+1)}`$, and all others zero. Set $`\lambda_a=\lambda_1`$ for $`a>0`$ and $`v_a=\sum_i c_i^{(a)}u_i/\sqrt{\lambda_a}`$. The $`v_a`$ are orthonormal and orthogonal to the probe. The accepted map to neutral label space is

```math
 K=\sum_{a=0}^{n-1}\sqrt{\eta/\lambda_a}\,|c^{(a)}\rangle\langle v_a|,
 \qquad K|\psi\rangle=0,\qquad Ku_i=\sqrt\eta\,|i\rangle.
```

To prepare each $`v_a`$ from a known charge-$`Q`$ word, group strings by their weight $`w`$ and, for $`a>0`$, by the weight $`u`$ in the first $`a`$ positions and occupation $`b`$ in position $`a+1`$. The remaining sensor positions have weight $`w-u-b`$. The coefficient multiplying the probe amplitude on each string depends only on these counts. It is $`(n-2w-nz)/\sqrt n`$ for the uniform label vector and $`2(ab-u)/\sqrt{a(a+1)}`$ for $`a>0`$. The multiplicity for the latter counts is $`\binom a u\binom{n-a-1}{w-u-b}`$. Prepare a neutral register over the feasible counts with signed amplitudes obtained by multiplying the probe amplitude by this coefficient and the square root of the multiplicity, then dividing by $`\sqrt{\lambda_a}`$. There are $`O(n)`$ possible count tuples because the sensor-weight support is bounded. Controlled swaps place the charge-$`Q`$ anchor into a sorted word for each tuple. Blockwise Dicke preparation fills the sensor blocks uniformly at their prescribed weights. Erase the neutral tuple register by computing the same counts from the occupations. Probe preparation uses the same steps with just the weight label and its amplitudes $`\sqrt{p_w}`$.

An explicit Dicke preparation step on a block of length $`m`$ uses the identity

```math
 |D_{m,w}\rangle=
 \sqrt{w/m}\,|1\rangle|D_{m-1,w-1}\rangle+
 \sqrt{1-w/m}\,|0\rangle|D_{m-1,w}\rangle .
```

On sorted words, the two representatives for $`0<w<m`$ differ only by exchanging wires 1 and $`w+1`$. Coherently count the current suffix weight, perform the corresponding exchange, and uncompute the count, which this exchange preserves. Sweep through progressively shorter suffixes. The two extreme weights require no exchange.

For readout let $`V_a`$ prepare $`v_a`$ from the common anchor word. Reserve a neutral pointer flag for refusal and one for each $`a`$. For each $`a`$, apply $`V_a^\dagger`$ on the refusal flag only, then rotate that flag into flag $`a`$ on the anchor with amplitude $`\sqrt{\eta/\lambda_a}`$. Apply $`V_a`$ on the refusal flag only. Accepted branches retain the anchor and skip later steps. Finally rotate the accepted flag basis to the vectors $`c^{(a)}`$ and measure the resulting label. This realizes $`K`$ and never borrows a charged ancilla.

For a coarse gate bound, put $`L=n+r_{\rm use}`$ and assume arbitrary wire-pair contacts, native arbitrary-angle rotations, and elementary neutral reversible gates. For each $`V_a`$, preparing the $`O(n)`$-entry sparse neutral register by two-level rotations, controlled sorting into charge representatives, and the suffix construction together take $`O(L^2\log L)`$ elementary gates with straightforward reversible count and pattern controls. The resulting readout bound is $`O(nL^2\log L)`$ gates and $`O(\log L)`$ neutral work and pointer bits. Contact restrictions, physical coherence time, and precision synthesis affect the physical cost.

## Minimal width and modular receivers

At bypass width $`r=k`$, the explicit joint reciprocal-vector POVM can be implemented wholly in the total-charge-$`k`$ sector. The occupied logical basis consists of strings $`|x\rangle|b_{k-|x|}\rangle`$, $`|x|\le k`$. For each sensor weight, a bypass word with the corresponding charge is sufficient. For any $`\eta`$ delivered by the moment-Gram contraction, write $`K`$ for the accepted map into a neutral label register. With $`|q_0\rangle`$ a known charged word in the same sector, an implementation is the isometry

```math
 |\phi\rangle|\varnothing\rangle\longmapsto
 \big(\sqrt{I-K^\dagger K}|\phi\rangle\big)|\varnothing\rangle
 +|q_0\rangle K|\phi\rangle .
```

The neutral pointer outcomes $`\varnothing`$ and $`S`$ realize abstention and the conclusive effects. If $`D`$ is the number of supported logical basis vectors and $`M=\binom Nk`$, the displayed map extends to a unitary on the span of the $`D`$ input words with the refusal flag and $`M`$ anchor words with accepted flags. A two-level decomposition needs $`O((D+M)^2)`$ arbitrary rotations in this basis. To implement a two-level rotation between two fixed-charge words, match the positions where the first word has an extra excitation with those where the second lacks one. Occupation-pattern-controlled exchanges and neutral flag flips give a path between the words. Conjugating a final rotation by the preceding controlled exchanges or flips restores the other basis words. With $`H=N+r+\lceil\log_2(M+1)\rceil`$ charged wires and pointer bits combined, the path has length $`O(H)`$, and each pattern control can count mismatches in $`O(H\log H)`$ elementary gates using $`O(\log H)`$ neutral work bits. This compiler works in any chosen total sector. At $`r=k`$, $`D=\sum_{j=0}^k\binom Nj`$.

Choose a uniform random partition into groups of prescribed sizes $`n_1,\ldots,n_k\ge2`$, and allocate $`r_j\ge1`$ bypass qubits to group $`j`$, with $`\sum_j r_j\le r`$. In each group run the one-hit receiver. Report a subset only if every group reports a label. If some group contains multiple hits, another has no hits and rejects with certainty. For any fixed size-$`k`$ hidden subset, the probability of one hit in each group is $`\prod_j n_j/\binom Nk`$. The resulting normalized conclusive probability is

```math
 \frac{1}{\binom Nk}\prod_{j=1}^k\min(n_j,r_j^2).
```

In particular, take $`r_j=1`$ at $`r=k`$. The single-excitation circuits above implement the exact minimum-width law with $`O(N)`$ exchanges and pointer rotations in total, plus the classical randomness and ability to route or schedule the chosen groups. One neutral pointer per active group suffices. Groups can be read in parallel when their couplings are available. A fixed physical wiring restricted to one particular partition serves only its promise of one hit per group. The randomized claim for an arbitrary size-$`k`$ hidden subset requires configurable grouping.

## Costs charged by the comparison

All comparators use the same $`N`$ device interrogations, once per device, in the same parallel round. The bypass-width saving counts every coherent spare wire able to hold compensating excitation. The local Bell benchmark has a short independent readout and costs one spare wire per site. The modular circuit exchanges bypass width for abstention probability and uses configurable routing. The joint compiler uses additional coherent control and charge-neutral work and outcome storage. The model tracks these as distinct resources. The bypass register may be consumed. Each bound is an ideal-gate statement, and exact abstention relies on the stated unitary fault model and exact control.

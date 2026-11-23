# Finite receiver checks

Run from the repository root with Python 3.10 or later and NumPy.

```bash
python3 verification/check_receiver.py
```

The check uses deterministic instances, including the identity pulse and four nonzero angles. A successful run prints these totals.

```text
Minimal-width receiver and balanced-split checks: 21 (N, k) cases
Below-threshold radial identities: 240 charge sectors
One-hit covariance tradeoff: 4272 (N, r) cases
One-hit full measurement and Gram certificates: 60 (N, r) cases
Single-excitation circuit checks: 31 arrays, five angles each
```

| Check | Evidence |
| --- | --- |
| [Minimum-width receiver](../science/exact-endpoints.md) | The Möbius vectors are reciprocal to the required moments. Their normalized conclusive effects leave a positive-semidefinite abstention effect. Direct pulse evolution gives the predicted diagonal matrix of label probabilities. |
| [Below-threshold obstruction](../science/exact-endpoints.md) | For $`2\le N\le 8`$ and each tested $`r\lt k`$, the summed degree-$`k`$ moment is in the span of the lower moments on every allowed charge sector. |
| [One-hit tradeoff](../science/exact-endpoints.md) | The two covariance eigenvalues give the claimed optimum for $`2\le N\le 90`$ and all tested widths from zero through $`N+1`$. For $`N\le9`$, the script constructs the reciprocal measurement explicitly for positive gain and compares label probabilities after pulse evolution. |
| [Moment-Gram criterion](../science/exact-frontier.md) | At the explicit minimum-width and one-hit preparations, the script checks the positive-semidefinite feasibility matrix and independently evaluates the limiting gain after projecting away lower moments. No optimization package is required. |
| [Charge-conserving implementation](../science/implementation.md) | The script composes the one-excitation preparation, exchanges and accepted pointer map on $`N+1`$ rails. It checks the correct label and abstention behavior for $`2\le N\le 32`$. It also enumerates the one-hit-per-group event in balanced partitions and compares the resulting minimum-width gain with the joint receiver. |

These computations check finite instances of the ideal model with a numerical tolerance of $`10^{-10}`$. The linked derivations establish the all-size and all-angle claims. The circuit check works in the one-excitation sector and uses the ideal accepted pointer amplitudes.

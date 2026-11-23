# Quantum Trajectory Reference Threshold

**How many idle qubits does a sensor array need to identify an unknown disturbance with zero wrong locations?** One pulse rotates an unknown set of $`k`$ of $`N`$ sensors about a known axis, by an unknown angle $`\theta`$. A fixed number-conserving receiver reports the affected set or abstains. The idle qubits are the only bypass systems that can store compensating excitation during the pulse.

The best worst-label success at angle $`\theta`$ is $`\eta^*_{N,k}(r)\sin^{2k}(\theta/2)`$, where $`r`$ counts idle reference qubits. For $`1\leq k\leq\lfloor N/2\rfloor`$ the coefficient obeys

| Reference budget | Optimal normalized success |
| --- | --- |
| $`r<k`$ | $`0`$ |
| $`r=k`$ | $`1/\binom Nk`$ |
| Single affected site, any $`r`$ | $`\min\{1,r^2/N\}`$ |
| Fixed $`k`$, arbitrary $`r\geq k`$ | $`\Theta_k\!\left(\min\{1,(r^2/N)^k\}\right)`$ |

A local Bell receiver uses one idle partner per sensor. The single-site receiver matches its success with $`\lceil\sqrt N\rceil`$ idle qubits. At the minimum budget $`r=k`$, randomly partitioning the sensors into $`k`$ groups and running a one-idle-qubit receiver on each achieves the exact optimum. This construction uses a linear number of ideal pair rotations and neutral controls given configurable routing. Both comparisons count bypass qubits separately from control and routing cost.

## Read the science

- [Model](science/model.md) defines the one-round promise and proves why a conclusive label isolates the degree-$`k`$ signal.
- [Exact endpoints](science/exact-endpoints.md) derives the minimum width, the reciprocal detector and the exact single-site curve.
- [Finite frontier](science/exact-frontier.md) constructs the finite-budget optimum through a moment Gram matrix and proves the fixed-$`k`$ scaling.
- [Implementation](science/implementation.md) describes the number-conserving circuits and their cost assumptions.
- [Related work](references/RELATED_WORK.md) traces the known-angle trajectory model, the uncalibrated Bell benchmark and the reference-frame setting.

The [learning route](learning/README.md) starts from a two-sensor calculation and develops all four parts in [30 milestones](learning/MILESTONES.md). To reproduce the finite receiver checks, run `python3 verification/check_receiver.py` using the requirements in [Verification](verification/README.md).

## Purpose and contact

This repository records the science and supports self-directed learning. For discussion or potential collaboration, contact Ruge Lin at [gogoko699`gmail.com](mailto:gogoko699`gmail.com).

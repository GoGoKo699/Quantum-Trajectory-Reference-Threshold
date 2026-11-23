#!/usr/bin/env python3
"""Finite checks of the one-shot, charge-conserving trajectory receiver."""

from itertools import combinations
from math import comb, isqrt, prod, sin

import numpy as np
from numpy.testing import assert_allclose


TOL = 1e-10
ANGLES = (-0.63, 0.0, 0.31, 1.02, 2.1)


def subsets(sites, size):
    return tuple(frozenset(sites_chosen) for sites_chosen in combinations(sites, size))


def parity_table(support, moments):
    return np.array([[(-1.0) ** len(x & a) for a in moments] for x in support])


def pulse_amplitudes(probe, support, hit_sets, theta):
    phase = np.array([[theta * len(x & hit) for hit in hit_sets] for x in support])
    return probe[:, None] * np.exp(1j * phase)


def check_effects(vectors, gain):
    # The nonzero eigenvalues of the conclusive frame operator are those of
    # this smaller Gram matrix. Eigenvalues at most one leave a positive residual.
    frame_eigenvalues = np.linalg.eigvalsh(gain * (vectors.T @ vectors))
    assert_allclose(frame_eigenvalues[-1], 1, atol=TOL, rtol=0)


def check_frontier_certificate(moment_vectors, first_label, gain):
    # Project off the lower-degree columns. This Schur complement gives the
    # best symmetric gain for this fixed preparation in the moment-Gram SDP.
    lower = moment_vectors[:, :first_label]
    labels = moment_vectors[:, first_label:]
    coefficients, *_ = np.linalg.lstsq(lower, labels, rcond=None)
    residuals = labels - lower @ coefficients
    spectrum = np.linalg.eigvalsh(residuals.T @ residuals)
    assert_allclose(spectrum[0], gain, atol=TOL, rtol=0)

    gram = moment_vectors.T @ moment_vectors
    label_diagonal = np.arange(first_label, gram.shape[0])
    gram[label_diagonal, label_diagonal] -= gain
    assert_allclose(min(0.0, np.linalg.eigvalsh(gram)[0]), 0, atol=TOL, rtol=0)


def check_minimal_receiver(n, k):
    lower = sum((subsets(range(n), degree) for degree in range(k)), ())
    labels = subsets(range(n), k)
    support = lower + labels
    count = comb(n, k)
    probabilities = np.array(
        [comb(k, len(x)) / (2**k * comb(n, len(x))) for x in support]
    )
    assert_allclose(probabilities.sum(), 1, atol=TOL, rtol=0)
    probe = np.sqrt(probabilities)
    moment_vectors = probe[:, None] * parity_table(support, support)
    signed_weights = [(-2.0) ** (-k) * (-1) ** (k - len(x)) / amplitude
                      for x, amplitude in zip(support, probe)]
    reciprocals = np.array(
        [[weight if x <= hit else 0.0 for hit in labels]
         for x, weight in zip(support, signed_weights)]
    )
    expected_overlaps = np.concatenate((np.zeros((count, len(lower))), np.eye(count)), axis=1)
    assert_allclose(reciprocals.T @ moment_vectors, expected_overlaps, atol=TOL, rtol=0)
    reciprocal_gram = reciprocals.T @ reciprocals
    assert_allclose(reciprocal_gram.sum(axis=1), count, atol=TOL, rtol=0)
    gain = 1 / count
    check_effects(reciprocals, gain)
    check_frontier_certificate(moment_vectors, len(lower), gain)

    for theta in ANGLES:
        output_states = pulse_amplitudes(probe, support, labels, theta)
        probabilities = gain * abs(reciprocals.T @ output_states)**2
        expected = gain * sin(theta / 2)**(2 * k) * np.eye(count)
        assert_allclose(probabilities, expected, atol=TOL, rtol=0)


def check_subthreshold_sector(n, k, width, total_charge):
    low_charge = max(0, total_charge - width)
    high_charge = min(n, total_charge)
    support = sum((subsets(range(n), q) for q in range(low_charge, high_charge + 1)), ())
    lower = sum((subsets(range(n), d) for d in range(k)), ())
    labels = subsets(range(n), k)
    # Work on the entire allowed support. A column identity here also holds
    # after selecting fewer rows or multiplying the rows by probe amplitudes.
    moments = parity_table(support, lower + labels)
    lower_moments = moments[:, :len(lower)]
    radial_top = moments[:, len(lower):].sum(axis=1)
    coefficients, *_ = np.linalg.lstsq(lower_moments, radial_top, rcond=None)
    assert_allclose(lower_moments @ coefficients, radial_top, atol=TOL, rtol=0)


def one_hit_charge_law(n, width):
    if width == 0:
        return {n // 2: 1.0}

    root = isqrt(n)
    span = min(width, root if root * root == n else root + 1)
    left = (n - span) // 2
    right = left + span
    if span * span <= n:
        return {left: 0.5, right: 0.5}

    # Match mean N/2 and variance N/4 by mixing two fixed-mean laws.
    center = n / 2
    middle_lo = n // 2
    middle_hi = (n + 1) // 2
    outer_variance = (center - left) * (right - center)
    inner_variance = (center - middle_lo) * (middle_hi - center)
    mix = (n / 4 - inner_variance) / (outer_variance - inner_variance)
    law = {left: mix * (right - center) / span, right: mix * (center - left) / span}
    for middle in (middle_lo, middle_hi):
        law[middle] = law.get(middle, 0.0) + (1 - mix) / 2
    return law


def check_one_hit_gain(n, width):
    law = one_hit_charge_law(n, width)
    assert_allclose(sum(law.values()), 1, atol=TOL, rtol=0)
    if not all(0 <= q <= n and p >= 0 for q, p in law.items()):
        raise AssertionError(f"Invalid charge distribution for N={n}, width={width}: {law}")
    if max(law) - min(law) > width:
        raise AssertionError(f"Charge distribution exceeds width {width}: {law}")
    mean = sum(q * p for q, p in law.items())
    variance = sum((q - mean)**2 * p for q, p in law.items())
    parallel = 4 * variance / n
    transverse = 4 * sum(q * (n - q) * p for q, p in law.items()) / (n * (n - 1))
    assert_allclose(min(parallel, transverse), min(1, width * width / n), atol=TOL, rtol=0)


def check_one_hit_measurement(n, width):
    law = one_hit_charge_law(n, width)
    support = sum((subsets(range(n), q) for q, p in law.items() if p > 0), ())
    labels = subsets(range(n), 1)
    probe = np.sqrt([law[len(x)] / comb(n, len(x)) for x in support])
    moment_vectors = probe[:, None] * parity_table(support, (frozenset(),) + labels)
    gain = min(1, width * width / n)
    check_frontier_certificate(moment_vectors, 1, gain)
    if gain == 0:
        return

    columns = moment_vectors[:, 1:]
    centered_columns = columns - probe[:, None] * (probe @ columns)
    gram = centered_columns.T @ centered_columns
    reciprocals = np.linalg.solve(gram, centered_columns.T).T
    check_effects(reciprocals, gain)
    assert_allclose(reciprocals.T @ probe, 0, atol=TOL, rtol=0)
    for theta in ANGLES:
        output_states = pulse_amplitudes(probe, support, labels, theta)
        probabilities = gain * abs(reciprocals.T @ output_states)**2
        expected = gain * sin(theta / 2)**2 * np.eye(n)
        assert_allclose(probabilities, expected, atol=TOL, rtol=0)


def check_single_excitation_circuit(n):
    # Index zero denotes the bypass rail and the remaining indices the sensors.
    collapse = np.eye(n + 1)
    for rail in range(2, n + 1):
        cosine = np.sqrt((rail - 1) / rail)
        sine = 1 / np.sqrt(rail)
        exchange = np.array([[cosine, sine], [-sine, cosine]])
        collapse[[1, rail], :] = exchange @ collapse[[1, rail], :]

    pair_transform = np.eye(n + 1)
    pair_transform[:2, :2] = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    probe = collapse.T @ pair_transform.T[:, 0]
    expected_probe = np.full(n + 1, 1 / np.sqrt(2 * n))
    expected_probe[0] = 1 / np.sqrt(2)
    assert_allclose(probe, expected_probe, atol=TOL, rtol=0)

    accept_amplitudes = np.full(n + 1, -1 / np.sqrt(2))
    accept_amplitudes[:2] = (0.0, 1.0)
    accepted_map = collapse.T @ np.diag(accept_amplitudes) @ pair_transform @ collapse
    assert_allclose(accepted_map @ probe, 0.0, atol=TOL, rtol=0)
    labels = subsets(range(n), 1)
    support = (frozenset(),) + labels
    for theta in ANGLES:
        output_states = pulse_amplitudes(probe, support, labels, theta)
        accepted = accepted_map @ output_states
        assert_allclose(accepted[0, :], 0.0, atol=TOL, rtol=0)
        expected = sin(theta / 2)**2 * np.eye(n) / n
        assert_allclose(abs(accepted[1:, :])**2, expected, atol=TOL, rtol=0)


def check_balanced_split(n, k):
    # A uniformly random partition for a fixed hit set has the same occupancy
    # count distribution as a uniform hit set for a fixed balanced partition.
    groups = tuple(frozenset(range(group, n, k)) for group in range(k))
    labels = subsets(range(n), k)
    split_count = sum(all(len(group & hit) == 1 for group in groups) for hit in labels)
    conditional_gain = prod(1 / len(group) for group in groups)
    assert_allclose(split_count / len(labels) * conditional_gain, 1 / comb(n, k), atol=TOL, rtol=0)


def main():
    minimal_width_cases = 0
    for n in range(2, 11):
        for k in range(1, min(n // 2, 3) + 1):
            check_minimal_receiver(n, k)
            check_balanced_split(n, k)
            minimal_width_cases += 1

    subthreshold_sectors = 0
    for n in range(2, 9):
        for k in range(1, n // 2 + 1):
            for width in range(k):
                for total_charge in range(n + width + 1):
                    check_subthreshold_sector(n, k, width, total_charge)
                    subthreshold_sectors += 1

    one_hit_gains = 0
    one_hit_measurements = 0
    for n in range(2, 91):
        for width in range(n + 2):
            check_one_hit_gain(n, width)
            one_hit_gains += 1
            if n <= 9:
                check_one_hit_measurement(n, width)
                one_hit_measurements += 1

    for n in range(2, 33):
        check_single_excitation_circuit(n)

    print(f"Minimal-width receiver and balanced-split checks: {minimal_width_cases} (N, k) cases")
    print(f"Below-threshold radial identities: {subthreshold_sectors} charge sectors")
    print(f"One-hit covariance tradeoff: {one_hit_gains} (N, r) cases")
    print(f"One-hit full measurement and Gram certificates: {one_hit_measurements} (N, r) cases")
    print("Single-excitation circuit checks: 31 arrays, five angles each")


if __name__ == "__main__":
    main()

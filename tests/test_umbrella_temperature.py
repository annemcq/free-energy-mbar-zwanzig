"""Regression tests for temperature scaling of umbrella reduced potentials."""

import numpy as np

from src.umbrella_data import build_reduced_potentials


def test_reduced_potentials_scale_inversely_with_temperature():
    """For fixed harmonic biases, u = U/(kB*T) must scale as 1/T."""
    centers = np.array([0.0, 90.0])
    springs = np.array([10.0, 20.0])  # kJ/mol/rad^2
    samples = [np.array([10.0, 20.0]), np.array([80.0])]

    u_300, n_300 = build_reduced_potentials(centers, springs, samples, temperature=300.0)
    u_600, n_600 = build_reduced_potentials(centers, springs, samples, temperature=600.0)

    np.testing.assert_array_equal(n_300, [2, 1])
    np.testing.assert_array_equal(n_600, n_300)

    for k, count in enumerate(n_300):
        # Only compare populated samples; the remaining tensor entries are padding.
        np.testing.assert_allclose(u_600[k, :, :count], u_300[k, :, :count] / 2)
        # Verify the physical scale, not only the ratio between temperatures.
        for l in range(len(centers)):
            delta_deg = np.abs(samples[k] - centers[l])
            delta_deg = np.minimum(delta_deg, 360.0 - delta_deg)
            bias_kj_mol = 0.5 * springs[l] * np.deg2rad(delta_deg) ** 2
            expected = bias_kj_mol / (0.008314462618 * 300.0)
            np.testing.assert_allclose(u_300[k, l, :count], expected, rtol=1e-12, atol=1e-12)

    assert np.any(u_300 > 0), "Test inputs must produce nonzero energies"

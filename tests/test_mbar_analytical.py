"""
Validation test: MBAREstimator vs. the known analytical free energies
of a set of harmonic oscillators (pymbar.testsystems).

This is Part 1 of the project -- before trusting MBAR on real
umbrella-sampling data (Part 2, T4 lysozyme), confirm the
wrapper reproduces free energies we can compute in closed form.
"""

import numpy as np
import pytest
from pymbar.testsystems import HarmonicOscillatorsTestCase

from src.mbar_estimator import MBAREstimator


@pytest.fixture(scope="module")
def harmonic_case():
    return HarmonicOscillatorsTestCase(O_k=[0, 1, 2, 3, 4], K_k=[1, 2, 4, 8, 16], beta=1.0)


def test_mbar_recovers_analytical_free_energies(harmonic_case):
    N_k = np.array([200, 200, 200, 200, 200])
    x_kn, u_kln, s_n = harmonic_case.sample(N_k=N_k, mode="u_kln")

    estimator = MBAREstimator.from_u_kln(u_kln, N_k)
    f_k, df_k = estimator.free_energies_relative_to_state_0()

    analytical_f_k = harmonic_case.analytical_free_energies()
    analytical_rel = analytical_f_k - analytical_f_k[0]

    # MBAR estimate should fall within ~5 standard errors of the
    # analytical value for every state (generous tolerance to keep the
    # test robust to RNG noise across runs while still catching real bugs).
    for k in range(len(N_k)):
        assert abs(f_k[k] - analytical_rel[k]) < 5 * max(df_k[k], 1e-8), (
            f"state {k}: MBAR={f_k[k]:.4f} +/- {df_k[k]:.4f}, "
            f"analytical={analytical_rel[k]:.4f}"
        )


def test_mbar_zero_self_difference(harmonic_case):
    """Delta_f[i, i] must always be exactly zero."""
    N_k = np.array([100, 100, 100, 100, 100])
    x_kn, u_kln, s_n = harmonic_case.sample(N_k=N_k, mode="u_kln")
    estimator = MBAREstimator.from_u_kln(u_kln, N_k)
    result = estimator.free_energy_differences()
    np.testing.assert_allclose(np.diag(result["Delta_f"]), 0.0, atol=1e-10)

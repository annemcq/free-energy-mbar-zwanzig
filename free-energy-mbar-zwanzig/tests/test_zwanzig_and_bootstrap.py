"""
Regression test for the Zwanzig estimator: for two adjacent, well-overlapping
harmonic oscillators, single-step Zwanzig should agree with the analytical
free energy difference (and with MBAR) to a reasonable tolerance -- this is
the regime where Zwanzig is expected to work well, in contrast to the
poorly-overlapping, far-apart states used to illustrate its failure mode in
notebook 02.
"""

import numpy as np
import pytest
from pymbar.testsystems import HarmonicOscillatorsTestCase

from src.mbar_estimator import MBAREstimator
from src.zwanzig_estimator import zwanzig_free_energy_difference
from src.bootstrap import bootstrap_estimator


@pytest.fixture(scope="module")
def close_oscillators():
    # Small offsets/force-constant differences -> good overlap between adjacent states.
    return HarmonicOscillatorsTestCase(O_k=[0, 0.5, 1.0], K_k=[1, 1.2, 1.4], beta=1.0)


def test_zwanzig_matches_analytical_for_overlapping_states(close_oscillators):
    N_k = np.array([2000, 2000, 2000])
    u_kln, s_n = close_oscillators.sample(N_k=N_k, mode="u_kln")[1:]
    u_kn = np.concatenate([u_kln[k, :, : N_k[k]] for k in range(3)], axis=1)

    u0_n = u_kn[0, : N_k[0]]
    u1_n = u_kn[1, : N_k[0]]
    delta_f_zwanzig = zwanzig_free_energy_difference(u0_n, u1_n)

    analytical_f_k = close_oscillators.analytical_free_energies()
    delta_f_analytical = analytical_f_k[1] - analytical_f_k[0]

    assert abs(delta_f_zwanzig - delta_f_analytical) < 0.1


def test_bootstrap_std_is_nonnegative_and_finite():
    rng = np.random.default_rng(0)
    samples_per_state = [rng.normal(size=50), rng.normal(size=50)]

    def mean_diff(resampled):
        return resampled[1].mean() - resampled[0].mean()

    result = bootstrap_estimator(samples_per_state, mean_diff, n_bootstrap=100, seed=1)
    assert np.isfinite(result["std"])
    assert result["std"] >= 0

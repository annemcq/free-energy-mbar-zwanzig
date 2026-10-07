"""
Simple resampling-based bootstrap for putting error bars on estimators
that don't provide their own uncertainty (e.g. Zwanzig), and for
cross-checking MBAR's analytical standard errors.

This is not a general-purpose library, but a small, readable 
implementation for the per-state sample layout used
throughout this project (a list of per-state 1D arrays).
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np


def bootstrap_estimator(
    samples_per_state: Sequence[np.ndarray],
    estimator_fn: Callable[[list[np.ndarray]], float],
    n_bootstrap: int = 200,
    seed: int | None = 0,
) -> dict:
    """
    Resample each state's samples with replacement (keeping each
    state's sample count fixed) and re-run `estimator_fn` on each
    resampled dataset.

    Parameters
    ----------
    samples_per_state : sequence of np.ndarray
        One array per thermodynamic state; arrays may have different
        lengths and different shapes per element (e.g. each entry can
        itself be a small tuple of arrays -- estimator_fn decides what
        to do with a resampled state's data).
    estimator_fn : callable
        Takes a list of resampled per-state arrays (same structure as
        samples_per_state) and returns a scalar (or array) estimate.
    n_bootstrap : int
        Number of bootstrap replicates.
    seed : int or None
        RNG seed for reproducibility.

    Returns
    -------
    dict with 'estimates' (array of length n_bootstrap), 'mean', 'std'
    """
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(n_bootstrap):
        resampled = []
        for state_samples in samples_per_state:
            n = len(state_samples)
            idx = rng.integers(0, n, size=n)
            resampled.append(np.asarray(state_samples)[idx])
        estimates.append(estimator_fn(resampled))
    estimates = np.array(estimates)
    return {
        "estimates": estimates,
        "mean": estimates.mean(axis=0),
        "std": estimates.std(axis=0, ddof=1),
    }

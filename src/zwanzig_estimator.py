"""
Zwanzig free-energy perturbation estimators used for comparison with MBAR.

For samples drawn from state 0, the reduced free-energy difference to
state 1 is estimated as

    Delta_f = -log < exp(-(u_1 - u_0)) >_0

where u_0 and u_1 are reduced potential energies. The estimator is most
reliable when the sampled configurations have sufficient overlap between
the two states.
"""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp


def zwanzig_free_energy_difference(u0_n: np.ndarray, u1_n: np.ndarray) -> float:
    """
    Estimate Delta_f = f_1 - f_0 from samples drawn from state 0 only.

    Parameters
    ----------
    u0_n : np.ndarray, shape (N,)
        Reduced potential energy of each sample, evaluated in state 0
        (the state they were drawn from).
    u1_n : np.ndarray, shape (N,)
        Reduced potential energy of the same samples, evaluated in
        state 1 (the perturbed state).

    Returns
    -------
    float
        Delta_f in units of kT.
    """
    u0_n = np.asarray(u0_n)
    u1_n = np.asarray(u1_n)
    n = len(u0_n)
    # -log( (1/n) sum exp(-(u1-u0)) ) = -[ logsumexp(-(u1-u0)) - log(n) ]
    return -(logsumexp(-(u1_n - u0_n)) - np.log(n))


def zwanzig_pairwise_chain(u_kn: np.ndarray, N_k: np.ndarray) -> np.ndarray:
    """
    Apply Zwanzig sequentially along a chain of adjacent states
    0 -> 1 -> 2 -> ... -> K-1.

    Each pairwise estimate uses only samples drawn from the earlier state.
    The returned values are cumulative free-energy differences relative
    to state 0.

    Parameters
    ----------
    u_kn : np.ndarray, shape (K, N)
        Reduced potential energies for all samples evaluated in each state.
        Samples are stored contiguously by originating state.
    N_k : np.ndarray, shape (K,)
        Number of samples drawn from each state.

    Returns
    -------
    np.ndarray, shape (K,)
        Cumulative Delta_f values relative to state 0.
    """
    K = len(N_k)
    f_rel = np.zeros(K)
    offset = 0
    for k in range(K - 1):
        n_k = int(N_k[k])
        u_k_in_k = u_kn[k, offset:offset + n_k]
        u_kp1_in_k = u_kn[k + 1, offset:offset + n_k]
        f_rel[k + 1] = f_rel[k] + zwanzig_free_energy_difference(
            u_k_in_k, u_kp1_in_k
        )
        offset += n_k
    return f_rel

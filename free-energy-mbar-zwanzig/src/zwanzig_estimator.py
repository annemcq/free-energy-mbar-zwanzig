"""
Single-step free energy perturbation or Zwanzig is utilised here as
a less costly comparison against MBAR on the same umbrella-sampling
data (Part 2 of the project).

Zwanzig's equation for the free energy difference between two states
0 -> 1, estimated from samples drawn from state 0:

    Delta_f = -log < exp(-(u_1 - u_0)) >_0

where the average is over samples n drawn from state 0, and u_0, u_1
are reduced potential energies (already divided by kT). This is
numerically unstable when the two states barely overlap -- that
instability, and where it shows up, is exactly what the notebook
compares against MBAR's more robust multi-state estimate.
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
    0 -> 1 -> 2 -> ... -> K-1, each leg estimated only from the samples
    drawn from the earlier state of the pair. Returns cumulative
    Delta_f relative to state 0, shape (K,).

    This mirrors how a rescoring pipeline that only has 'forward'
    samples (e.g. WT-only trajectories reweighted onto mutants) would
    have to chain estimates, and is deliberately less data-efficient
    than MBAR, which pools all states at once.
    """
    K = len(N_k)
    f_rel = np.zeros(K)
    offset = 0
    for k in range(K - 1):
        n_k = int(N_k[k])
        u_k_in_k = u_kn[k, offset:offset + n_k]
        u_kp1_in_k = u_kn[k + 1, offset:offset + n_k]
        f_rel[k + 1] = f_rel[k] + zwanzig_free_energy_difference(u_k_in_k, u_kp1_in_k)
        offset += n_k
    return f_rel

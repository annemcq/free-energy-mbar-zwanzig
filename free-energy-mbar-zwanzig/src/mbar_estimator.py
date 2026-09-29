"""
Thin wrapper around pymbar.MBAR for computing free energy differences
between a set of thermodynamic states from reduced-potential-energy
matrices.

This is kept deliberately simple, as pymbar handles the core of the work.
The wrapper provides the rest of the project (the tests, the two
notebooks) with one stable and documented entry point, instead of
calling the pymbar API directly in several different places. 
"""

from __future__ import annotations

import numpy as np
from pymbar import MBAR, FES


class MBAREstimator:
    """
    Free energy estimator based on the Multistate Bennett Acceptance
    Ratio (MBAR) method.

    Parameters
    ----------
    u_kn : np.ndarray, shape (K, N)
        Reduced potential energies. u_kn[k, n] is the reduced potential
        of sample n (drawn from *some* state) evaluated in state k.
    N_k : np.ndarray, shape (K,)
        Number of samples drawn from each state k. Samples are assumed
        to be stored contiguously per originating state, matching the
        pymbar convention.
    """

    def __init__(self, u_kn: np.ndarray, N_k: np.ndarray, **mbar_kwargs):
        self.u_kn = np.asarray(u_kn)
        self.N_k = np.asarray(N_k)
        self._mbar = MBAR(self.u_kn, self.N_k, **mbar_kwargs)

    @classmethod
    def from_u_kln(cls, u_kln: np.ndarray, N_k: np.ndarray, **mbar_kwargs) -> "MBAREstimator":
        """
        Build an estimator from the (K, L, N_max) 'u_kln' layout that
        pymbar's testsystems produce, converting it to the flat (K, N)
        'u_kn' layout MBAR expects.
        """
        K, L, N_max = u_kln.shape
        assert K == L, "u_kln must have matching K and L (same number of states sampled/evaluated)"
        u_kn = np.zeros((K, N_max * K))
        col = 0
        for k in range(K):
            n_k = int(N_k[k])
            u_kn[:, col:col + n_k] = u_kln[k, :, :n_k]
            col += n_k
        u_kn = u_kn[:, :col]
        return cls(u_kn, N_k, **mbar_kwargs)

    def free_energy_differences(self) -> dict:
        """
        Returns
        -------
        dict with keys:
            'Delta_f' : (K, K) matrix, Delta_f[i, j] = f_j - f_i (in kT)
            'dDelta_f': (K, K) matrix of estimated standard errors
        """
        result = self._mbar.compute_free_energy_differences()
        return {"Delta_f": result["Delta_f"], "dDelta_f": result["dDelta_f"]}

    def free_energies_relative_to_state_0(self) -> tuple[np.ndarray, np.ndarray]:
        """Convenience accessor: f_k - f_0 and its standard error, for all k."""
        result = self.free_energy_differences()
        return result["Delta_f"][0, :], result["dDelta_f"][0, :]

    def compute_pmf(
        self,
        u_n: np.ndarray,
        x_n: np.ndarray,
        bin_edges: np.ndarray,
        reference_point: str = "from-lowest",
        n_bootstraps: int = 0,
        seed: int = -1,
    ) -> dict:
        """
        Compute a 1D potential of mean force (PMF/FES) over bins defined by
        bin_edges, using the samples' collective coordinate x_n.

        Thin wrapper around pymbar's `FES` class (the modern replacement
        for the old `MBAR.compute_pmf` / `computePMF` API removed in
        pymbar >= 4).

        Parameters
        ----------
        u_n : np.ndarray, shape (N,)
            Reduced potential energy of every sample evaluated in the
            (unbiased) state the PMF should be reported in -- typically
            u_kn[0, :] concatenated in the same sample order as x_n.
        x_n : np.ndarray, shape (N,)
            Collective variable value for every sample (same order as u_n).
        bin_edges : np.ndarray, shape (n_bins + 1,)
        n_bootstraps : int
            If > 0, also compute bootstrap uncertainties (see pymbar FES docs).

        Returns
        -------
        dict with 'f_i' (PMF in kT per bin, relative to reference_point),
        'df_i' (uncertainty per bin), and 'bin_centers'.
        """
        fes = FES(self.u_kn, self.N_k)
        histogram_parameters = {"bin_edges": [np.asarray(bin_edges)]}
        fes.generate_fes(
            u_n,
            x_n.reshape(-1, 1),
            fes_type="histogram",
            histogram_parameters=histogram_parameters,
            n_bootstraps=n_bootstraps,
        )
        bin_centers = 0.5 * (np.asarray(bin_edges)[:-1] + np.asarray(bin_edges)[1:])
        uncertainty_method = "bootstrap" if n_bootstraps > 0 else "analytical"
        result = fes.get_fes(
            bin_centers.reshape(-1, 1),
            reference_point=reference_point,
            uncertainty_method=uncertainty_method,
        )
        return {
            "f_i": result["f_i"],
            "df_i": result.get("df_i"),
            "bin_centers": bin_centers,
        }

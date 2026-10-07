"""
Utilities for loading the T4 lysozyme umbrella-sampling dataset and
constructing reduced potentials.

The dataset contains umbrella simulations of the chi torsion of a valine
sidechain in T4 lysozyme L99A with benzene bound in the cavity. It is
distributed with the pymbar umbrella-sampling examples and is based on
Mobley et al., J. Mol. Biol. 371(4):1118-1134 (2007).

Input files follow the original centers.dat / prod{k}_dihed.xvg layout.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def load_umbrella_dataset(data_dir: str | Path) -> dict:
    """
    Load umbrella centers, spring constants, and torsion-angle trajectories.

    Parameters
    ----------
    data_dir : str or Path
        Folder containing centers.dat and prod{k}_dihed.xvg files.

    Returns
    -------
    dict
        K : int
            Number of umbrella windows.
        chi0_k : np.ndarray, shape (K,)
            Umbrella centers in degrees.
        Kspring_k : np.ndarray, shape (K,)
            Spring constants in kJ/mol/rad^2.
        chi_kn : list of np.ndarray
            Torsion-angle samples for each window, in degrees.
    """
    data_dir = Path(data_dir)
    centers = np.loadtxt(data_dir / "centers.dat")
    chi0_k = centers[:, 0]
    Kspring_k = centers[:, 1]
    K = len(chi0_k)

    chi_kn = []
    for k in range(K):
        fname = data_dir / f"prod{k}_dihed.xvg"
        chi_vals = []
        with open(fname) as f:
            for line in f:
                if line.startswith(("#", "@")):
                    continue
                tokens = line.split()
                if len(tokens) < 2:
                    continue
                chi = float(tokens[1])

                # Wrap into [-180, 180)
                while chi < -180.0:
                    chi += 360.0
                while chi >= 180.0:
                    chi -= 360.0

                chi_vals.append(chi)

        chi_kn.append(np.array(chi_vals))

    return {
        "K": K,
        "chi0_k": chi0_k,
        "Kspring_k": Kspring_k,
        "chi_kn": chi_kn,
    }


def build_reduced_potentials(
    chi0_k,
    Kspring_k,
    chi_kn,
    temperature: float = 300.0,
):
    """
    Build the reduced-potential tensor for the umbrella windows.

    The spring constants are given in kJ/mol/rad^2 and converted to
    reduced units with beta = 1 / (kB * T). The underlying unbiased
    potential is common to all windows, so only the harmonic bias is
    required for the free-energy comparison.

    Parameters
    ----------
    chi0_k : array-like
        Umbrella centers in degrees.
    Kspring_k : array-like
        Harmonic spring constants in kJ/mol/rad^2.
    chi_kn : list of arrays
        Sampled torsion angles for each umbrella window, in degrees.
    temperature : float, optional
        Simulation temperature in kelvin. Default is 300 K.

    Returns
    -------
    u_kln : np.ndarray, shape (K, K, N_max)
        Reduced potential energies.
    N_k : np.ndarray, shape (K,)
        Number of samples from each umbrella window.
    """
    kB = 0.008314462618  # kJ mol^-1 K^-1
    beta = 1.0 / (kB * temperature)

    K = len(chi0_k)
    N_k = np.array([len(c) for c in chi_kn])
    N_max = int(N_k.max())
    u_kln = np.zeros((K, K, N_max))

    for k in range(K):
        n_k = N_k[k]
        chi_samples = chi_kn[k]

        dchi = chi_samples[None, :] - chi0_k[:, None]

        # Minimum-image convention for a periodic angle in degrees
        dchi = np.where(
            np.abs(dchi) > 180.0,
            360.0 - np.abs(dchi),
            dchi,
        )

        u_kln[k, :, :n_k] = (
            beta
            * (Kspring_k[:, None] / 2.0)
            * np.deg2rad(dchi) ** 2
        )

    return u_kln, N_k


def flatten_u_kln(u_kln: np.ndarray, N_k: np.ndarray) -> np.ndarray:
    """
    Convert (K, K, N_max) u_kln into the flat (K, N) u_kn layout
    expected by MBAR/FES.
    """
    K = u_kln.shape[0]
    u_kn = np.concatenate(
        [u_kln[k, :, : N_k[k]] for k in range(K)],
        axis=1,
    )
    return u_kn

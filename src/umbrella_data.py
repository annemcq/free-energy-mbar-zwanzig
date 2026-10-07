"""
Loader for the umbrella-sampling dataset bundled in `data/`.

The chosen dataset is distributed with pymbar's own examples
(choderalab/pymbar-examples, umbrella-sampling-pmf/): umbrella
sampling simulation of the chi torsion of a valine sidechain in T4
lysozyme L99A with benzene bound in the cavity (Mobley et al., J. Mol.
Biol. 2007). Note that it is not the original alanine dipeptide system.
The latter was discarded as pymbar's own repository does not actually
ship umbrella-sampling data for alanine dipeptide (its alanine dipeptide
example uses parallel tempering instead). This lysozyme dataset is the
real umbrella-sampling system pymbar ships, which is what is subsequently
used in Part 2; the methodology (MBAR, Zwanzig, bootstrap) is identical
regardless of which system it's applied to.

Format, per the original `centers.dat` / `prod{k}_dihed.xvg` layout:

    data/centers.dat
        one line per umbrella window k: "<center_deg> <spring_const>"
        (spring constant in kJ/mol/rad^2; a 3rd column, if present, is
        a per-window temperature -- unused here since all windows in
        this dataset share one temperature)

    data/prod{k}_dihed.xvg
        Gromacs xvg time series for window k; comment lines start with
        '#' or '@'. Data lines are "<time> <chi_degrees>".
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def load_umbrella_dataset(data_dir: str | Path) -> dict:
    """
    Parameters
    ----------
    data_dir : path to the folder containing centers.dat and prod{k}_dihed.xvg

    Returns
    -------
    dict with:
        'K'        : int, number of umbrella windows
        'chi0_k'   : (K,) spring center for each window, degrees
        'Kspring_k': (K,) spring constant for each window, kJ/mol/rad^2
        'chi_kn'   : list of length K, each a 1D array of chi values (degrees)
                     for that window (ragged -- windows have different N_k)
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
    Build the (K, K, N_max) reduced-potential-energy tensor u_kln for a set
    of harmonic umbrella windows on a periodic angular coordinate.

    The umbrella spring constants are given in kJ/mol/rad^2. They are
    converted to reduced units using beta = 1 / (kB * T), with a default
    temperature of 300 K matching the original pymbar umbrella-sampling
    example.

    u_kln[k, l, n] =
        beta * (Kspring_l / 2) * dchi(chi_kn[k][n], chi0_l)^2

    The underlying unbiased potential is identical for all umbrella
    windows, so only the harmonic bias is needed when constructing the
    reduced potentials used for the free-energy comparison.

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

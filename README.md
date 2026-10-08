# Free Energy Estimation: MBAR & Zwanzig

[![Tests](https://github.com/annemcq/free-energy-mbar-zwanzig/actions/workflows/tests.yml/badge.svg)](https://github.com/annemcq/free-energy-mbar-zwanzig/actions/workflows/tests.yml)

The current work carries out free energy estimation using biased/multistate simulation data, via single-step Zwanzig perturbation and Multistate Bennett Acceptance Ratio (MBAR). This methodology is first validated against an analytical solution and subsequently applied to real public umbrella-sampling data. In practice, it applies the free-energy-estimation methodology used in molecular simulation research (MBAR, Zwanzig reweighting, bootstrap error estimation) to public benchmark systems, independent of any specific research project's data.

## Structure

```text
free-energy-mbar-zwanzig/
├── README.md
├── environment.yml
├── pytest.ini
├── data/
│   └── umbrella_lysozyme_chi/     # real umbrella-sampling data (see below)
├── src/
│   ├── mbar_estimator.py          # documented wrapper around pymbar.MBAR / FES
│   ├── zwanzig_estimator.py       # single-step FEP estimator + chained variant
│   ├── bootstrap.py               # generic per-state resampling bootstrap
│   └── umbrella_data.py           # parses the umbrella-sampling dataset, builds u_kln
├── notebooks/
│   ├── 01_validate_analytical.ipynb   # MBAR vs. harmonic-oscillator analytical solution
│   └── 02_lysozyme_pmf.ipynb          # real PMF, MBAR vs. Zwanzig, bootstrap error bars
├── tests/
│   ├── test_mbar_analytical.py
│   └── test_zwanzig_and_bootstrap.py
└── results/
    └── figures/
```

Run everything:

```bash
conda env create -f environment.yml
conda activate free-energy-mbar-zwanzig
pytest tests/
jupyter nbconvert --to notebook --execute --inplace notebooks/*.ipynb
```

## Part 1 — Validation against an analytical solution

`HarmonicOscillatorsTestCase` from `pymbar.testsystems` provides a set of harmonic oscillators with known closed-form free energies. `notebooks/01_validate_analytical.ipynb` and `tests/test_mbar_analytical.py` both check that `MBAREstimator` reproduces those free energies to within a few standard errors, before applying it to real data.

![MBAR vs. analytical free energies](results/figures/mbar_vs_analytical_harmonic.png)

*MBAR's estimate (with standard error) tracks the analytical solution for every state.*

## Sampling overlap and estimator reliability

Before comparing free-energy estimates, the analysis now checks the statistical overlap between the umbrella windows using MBAR's state-overlap matrix. This provides a direct diagnostic of whether neighboring windows form a sufficiently connected sampling network.

![MBAR state-overlap matrix](results/figures/mbar_state_overlap.png)

The overlap analysis also provides context for the chained Zwanzig calculation: weak overlap between adjacent windows is a regime where single-step perturbation estimates become more sensitive to sampling noise, and uncertainty can accumulate along the chain.

## Part 2 — Real umbrella-sampling PMF: MBAR vs. Zwanzig

The real-data example uses the chi torsion of a valine sidechain in T4 lysozyme L99A with benzene bound in the cavity (Mobley et al., *J. Mol. Biol.* 371(4):1118-1134, 2007). The data come from the umbrella-sampling example distributed with `pymbar`.

`notebooks/02_lysozyme_pmf.ipynb`:

1. Loads the 26-window umbrella-sampling dataset (`src/umbrella_data.py`)
2. Subsamples each window for statistical independence (`pymbar.timeseries`)
3. Computes the 1D PMF with MBAR
4. Computes relative window free energies with a chained single-step Zwanzig estimator and compares them with MBAR
5. Estimates uncertainty in the Zwanzig chain by bootstrap resampling

![1D PMF of the lysozyme chi torsion](results/figures/lysozyme_chi_pmf_mbar.png)

*MBAR reconstruction of the one-dimensional free-energy profile along the chi torsion.*

![MBAR vs. chained Zwanzig](results/figures/mbar_vs_zwanzig_chain.png)

*MBAR and the chained Zwanzig estimates agree closely over much of the umbrella-sampling range, particularly near the reference window. The bootstrap uncertainty of the Zwanzig chain increases as pairwise steps are accumulated, and differences of a few $k_BT$ appear for some later windows. MBAR retains smaller uncertainties by combining information from all sampled states in a single multistate estimate.*

## Notes

- Built with AI assistance.
- Uses only public data and public benchmark systems — no unpublished research data.

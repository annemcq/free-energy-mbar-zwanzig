# Free Energy Estimation: MBAR & Zwanzig

[![Tests](https://github.com/annemcq/free-energy-mbar-zwanzig/actions/workflows/tests.yml/badge.svg)](https://github.com/annemcq/free-energy-mbar-zwanzig/actions/workflows/tests.yml)

The current work carries out free energy estimation using biased/multistate simulation data, via single-step Zwanzig perturbation and Multistate Bennett Acceptance Ratio (MBAR). This methodology is first validated against an analytical solution and subsequently applied to real public umbrella-sampling data. In practice, it applies the free-energy-estimation methodology used in molecular simulation research (MBAR, Zwanzig reweighting, bootstrap error estimation) to public benchmark systems, independent of any specific
research project's data.

## Structure

```
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

`HarmonicOscillatorsTestCase` from `pymbar.testsystems` provides a set of harmonic
oscillators with known closed-form free energies. `notebooks/01_validate_analytical.ipynb`
and `tests/test_mbar_analytical.py` both check that `MBAREstimator` reproduces those free
energies to within a few standard errors, before trusting it on real data.

![MBAR vs. analytical free energies](results/figures/mbar_vs_analytical_harmonic.png)

*MBAR's estimate (with standard error) tracks the analytical solution for every state.*

## Part 2 — Real umbrella-sampling PMF: MBAR vs. Zwanzig

**System of choice.** The original system was alanine dipeptide. However, upon inspection of  `pymbar`'s own example repository (`choderalab/pymbar-examples`), I discovered its *umbrella-sampling* dataset is not alanine dipeptide, but of the **chi torsion of a valine sidechain in T4 lysozyme L99A**, with benzene bound in the cavity (Mobley et al., *J. Mol. Biol.* 371(4):1118-1134, 2007). It's worth noting that alanine dipeptide is present in said repository but only in the parallel-tempering 2D-PMF example. Given that the goal is genuine 1D umbrella-sampling PMF, the system of choice hereafter is the lysozyme dataset.

`notebooks/02_lysozyme_pmf.ipynb`:
1. Loads the 26-window umbrella-sampling dataset (`src/umbrella_data.py`)
2. Subsamples each window for statistical independence (`pymbar.timeseries`)
3. Computes the 1D PMF with MBAR — recovers the expected three-rotamer-well structure
4. Computes the same relative window free energies with a **chained single-step Zwanzig**
   estimator and compares against MBAR: they agree near the reference window and diverge
   increasingly for windows further along the chain, illustrating the known failure mode of
   pairwise FEP chaining vs. a true multistate estimator
5. Adds bootstrap error bars for the Zwanzig chain (no analytical uncertainty available there)

![1D PMF of the lysozyme chi torsion](results/figures/lysozyme_chi_pmf_mbar.png)

*The recovered PMF shows the expected three-rotamer-well structure of the valine chi torsion.*

![MBAR vs. chained Zwanzig](results/figures/mbar_vs_zwanzig_chain.png)

*MBAR and the chained Zwanzig estimator agree near the reference window, but Zwanzig's error
compounds and diverges for windows further along the chain — exactly the failure mode a
multistate estimator like MBAR is designed to avoid.*

## Notes

- Built with AI assistance (as disclosed for all projects in this portfolio); the
  methodology (MBAR/Zwanzig, coarse-graining GNNs used elsewhere in the portfolio) mirrors
  standard practice in molecular-simulation free-energy estimation.
- Uses only public data and public benchmark systems — no unpublished research data.

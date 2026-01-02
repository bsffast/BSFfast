# BSFfast

BSFfast is a lightweight numerical tool providing precomputed effective
bound-state formation (BSF) cross sections for annihilation processes
in the early Universe.

This repository accompanies the paper:

**BSFfast: Rapid computation of bound-state effects on annihilation in the early Universe**
[arXiv:2512.23812](https://arxiv.org/abs/2512.23812)

The code provides tabulated effective BSF cross sections and interpolation
routines for use in Boltzmann solvers.

If you use BSFfast in a publication, please cite the paper above.

## Contact

For scientific questions, usage issues, or bug reports, please contact:

- Stefan Lederer <stefan.lederer@tum.de>
- Jan Heisig <heisig@physik.rwth-aachen.de>

## Quick start (Python)

Place the directory `BSFfast_DataCSV/` next to `BSFfast.py`, then:

```python
from BSFfast import fastXS

# choose model parameters
x = 200.0
m = 1000.0  # GeV

# rescaled models (require coupling strength alpha)
xs = fastXS("dQCD-S", x, m, alpha=0.1)

# SM QCD tables support "cutoff" (default) and "plateau"
xs_cut = fastXS("QCD-SU", x, m)
xs_plat = fastXS("QCD-SU", x, m, "plateau")

```

A more comprehensive Python example is provided in `example_use.py`.

Corresponding example files are also provided for:

- **Mathematica**: `example_use.nb`
- **C**: `example_use.c`

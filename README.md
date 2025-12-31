# BSFfast

BSFfast is a lightweight numerical tool providing precomputed effective
bound-state formation (BSF) cross sections for annihilation processes
in the early Universe.

This repository accompanies the paper:

**BSFfast: Rapid computation of bound-state effects on annihilation in the early Universe**

The code provides tabulated effective BSF cross sections and interpolation
routines for use in Boltzmann solvers.

If you use BSFfast in a publication, please cite the paper above.

## Contact

For questions, comments, or bug reports, please contact:

- stefan.lederer@tum.de  
- heisig@physik.rwth-aachen.de

## Quick start (Python)

Place the directory `BSFfast_DataCSV/` next to `BSFfast.py`, then:

```python
from BSFfast import fastXS

# choose model parameters
x = 20.0
m = 1000.0  # GeV

# rescaled models (require coupling strength alpha)
xs = fastXS("dQCD-S", x, m, alpha=0.2)

# SM QCD tables support "cutoff" (default) and "plateau"
xs_cut = fastXS("QCD-SU", x, m)
xs_plat = fastXS("QCD-SU", x, m, "plateau")

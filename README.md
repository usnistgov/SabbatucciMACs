# SabbatucciMACs.jl

This package facilitates use of the atomic photoelectric effect data described
in [Theory and calculation of the atomic photoeffect](https://doi.org/10.1016/j.radphyschem.2015.10.021)
by Lorenzo Sabbatucci and FrancescSalvat.  It is particularly designed for
use in X-ray spectroscopy calculations in which the mass absorption coefficient
(MAC) (in cm²/g) is the most common representation.

Because the MAC is very sensitive with respect to energy around the absorption
edges, it is important to maintain a consistent set.  For this reason, the
library also makes available a database of atomic edge energies (in eV).

Currently, this library is not available through the Julia package system and
must be install through the GitHub URL directly.
```
pkg> add ....  
```

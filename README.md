# PyGFN0

[![Pypi version](https://img.shields.io/pypi/v/pygfn0)](https://pypi.org/project/pygfn0/)
[![PyPI Downloads](https://static.pepy.tech/badge/pygfn0)](https://pepy.tech/projects/pygfn0)

This is a Python version for [GFN0-xTB](https://github.com/pprcht/gfn0) based on F2PY.

---

## Usage

There is an `ase.Calculator` subclass termed `GFN0` for non-PBC system. And a more low function called `gfn0` can be found in [the code](https://github.com/LiuGaoyong/PyGFNFF/blob/main/pygfn0/_pygfn0.py).

```python
from ase import Atoms
from ase.build import molecule
from ase.optimize import BFGS
from scipy.spatial.distance import pdist

from pygfn0 import GFN0

atoms = Atoms(molecule("CO"), calculator=GFN0())
opt = BFGS(atoms, logfile="-", trajectory=None)
opt.run(fmax=0.03, steps=50)
print(f"Energy: {e:.3f}eV")
print(f"C=O Length: {bl:.3f}\u212b")

# Output:
#       Step     Time          Energy          fmax
# BFGS:    0 21:48:56     -145.434618        3.546420
# BFGS:    1 21:48:56     -145.184501        9.984024
# BFGS:    2 21:48:56     -145.493902        0.833465
# BFGS:    3 21:48:56     -145.496822        0.175982
# BFGS:    4 21:48:56     -145.496955        0.004404
# Energy: -145.497eV
# C=O Length: 1.116Å
```



### Reference

1. P. Pracht, S. Grimme, et.al. A Robust Non-Self-Consistent Tight-Binding Quantum Chemistry Method for large Molecules (2019), DOI: https://doi.org/10.26434/chemrxiv.8326202.v1
2. A standalone library of the GFN0-xTB method. https://github.com/pprcht/gfn0/

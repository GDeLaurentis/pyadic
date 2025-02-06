# pyAdic

[![CI Lint](https://github.com/GDeLaurentis/pyadic/actions/workflows/ci_lint.yml/badge.svg)](https://github.com/GDeLaurentis/pyadic/actions/workflows/ci_lint.yml)
[![CI Test](https://github.com/GDeLaurentis/pyadic/actions/workflows/ci_test.yml/badge.svg)](https://github.com/GDeLaurentis/pyadic/actions/workflows/ci_test.yml)
[![Coverage](https://img.shields.io/badge/Coverage-90%25-green?labelColor=2a2f35)](https://github.com/GDeLaurentis/pyadic/actions)
[![Docs](https://github.com/GDeLaurentis/pyadic/actions/workflows/cd_docs.yml/badge.svg?label=Docs)](https://gdelaurentis.github.io/pyadic/)
[![PyPI](https://img.shields.io/pypi/v/pyadic?label=PyPI)](https://pypi.org/project/pyadic/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/pyadic.svg?label=PyPI%20downloads)](https://pypistats.org/packages/pyadic)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GDeLaurentis/pyadic/HEAD)
[![DOI](https://zenodo.org/badge/519290930.svg)](https://zenodo.org/doi/10.5281/zenodo.11114230)
[![Python](https://img.shields.io/pypi/pyversions/pyadic?label=Python)](https://pypi.org/project/pyadic/)


The `pyadic` library is Python 3 package that provides number types for finite fields $\mathbb{F}_p$ (`ModP`) and $p$-adic numbers $\mathbb{Q}_p$ (`PAdic`). The goal is to mimic the flexible behavior of built-in types, such as `int`, `float` and `complex`. Thus, one can mix-and-match the different number types, as long as the operations are consistent. In particular, `ModP` and `PAdic` are compatible with `fractions.Fraction`.

In addition to arithmetic operations, the pyadic library also provides the following functions:

- `rationalise` to perform rationalization ($\mathbb{F}_p\rightarrow \mathbb{Q}$ and $\mathbb{Q}_p \rightarrow \mathbb{Q}$);
- `finite_field_sqrt` and `padic_sqrt` to compute square roots (which may involve `FieldExtension`);
- `padic_log` to compute the $p$-adic logarithm.
- polynomial and rational function interpolation, see `interpolation.py` module.

A shout-out to [galois](https://github.com/mhostetter/galois) for a very nice tool. It is recommented for vectorized finite field operations, unless type compatibility is an issue. For scalar operation this repo is recommended. See performance comparison below.

## Installation
The package is available on the [Python Package Index](https://pypi.org/project/pyadic/)
```console
pip install pyadic
```
Alternativelty, it can be installed by cloning the repo
```console
git clone https://github.com/GDeLaurentis/pyadic.git path/to/repo
pip install -e path/to/repo
```

### Requirements
`pip` will automatically install the required packages, which are
```
numpy, sympy
```
Additionally, `pytest` is needed for testing.

### Testing
Extensive tests are implemented with [pytest](https://github.com/pytest-dev/pytest)

```console
pytest --cov pyadic/ --cov-report html tests/ --verbose
```

## Quick Start

```python
In [1]: from pyadic import PAdic, ModP
In [2]: from fractions import Fraction as Q

# 7/13 as a 12-digit 2147483647-adic number
In [3]: PAdic(Q(7, 13), 2147483647, 12)  
Out [3]: 1817101548 + 825955248*2147483647 + 1156337348*2147483647^2 + 330382099*2147483647^3 + 1321528398*2147483647^4 + 991146298*2147483647^5 + 1817101547*2147483647^6 + 825955248*2147483647^7 + 1156337348*2147483647^8 + 330382099*2147483647^9 + 1321528398*2147483647^10 + 991146298*2147483647^11 + O(2147483647^12)

# 7/13 in F_2147483647
In [4]: ModP(Q(7, 13), 2147483647)
Out [4]: 1817101548 % 2147483647

# Mapping back to rational numbers
In [5]: from pyadic.finite_field import rationalise
In [6]: rationalise(ModP(Q(7, 13), 2147483647))
Out [6]: Fraction(7, 13)
In [7]: rationalise(PAdic(Q(7, 13), 2147483647, 12))
Out [7]: Fraction(7, 13)
```

## Perfomance comparison with [galois](https://github.com/mhostetter/galois) for finite fields

Scalar instantiation and operations are faster in pyadic
```python
import numpy
from galois import GF
from pyadic import ModP
from random import randint

GFp = GF(2 ** 31 - 1)
x = randint(0, 2 ** 31 - 1)

%timeit GFp(x)
2.84 µs ± 63.5 ns

%timeit ModP(x, 2 ** 31 - 1)
297 ns ± 0.876 ns

%timeit GFp(x) ** 2
30.1 µs ± 20.6 µs 

%timeit ModP(x, 2 ** 31 - 1) ** 2
2.23 µs ± 91.8 ns
```

while galois is faster for vectorized operations (the bigger the array the bigger the gain)
```python
%timeit numpy.array([randint(0, 2 ** 31 - 1) for i in range(100)]).view(GFp) ** 2
65.6 µs ± 1.86 µs 

%timeit numpy.array([ModP(randint(0, 2 ** 31 - 1), 2 ** 31 - 1) for i in range(100)]) ** 2
351 µs ± 9.28 µs
```

However, galois requires everything to be appropriately typed, while pyadic performs type-casting on-the-fly
```python
numpy.array([randint(0, 2 ** 31 - 1) for i in range(100)]).view(GFp) / 2
TypeError

numpy.array([ModP(randint(0, 2 ** 31 - 1), 2 ** 31 - 1) for i in range(100)]) / 2
array([...], dtype=object)
```

## Citation

If you found this library useful, please consider citing it


```bibtex
@inproceedings{DeLaurentis:2023qhd,
    author = "De Laurentis, Giuseppe",
    title = "{Lips: $p$-adic and singular phase space}",
    booktitle = "{21th International Workshop on Advanced Computing and Analysis Techniques in Physics Research}: {AI meets Reality}",
    eprint = "2305.14075",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    reportNumber = "PSI-PR-23-14",
    month = "5",
    year = "2023"
}
```

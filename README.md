# pyAdic

[![Continuous Integration Status](https://github.com/GDeLaurentis/pyadic/actions/workflows/continuous_integration.yml/badge.svg)](https://github.com/GDeLaurentis/pyadic/actions)
[![Coverage](https://img.shields.io/badge/Coverage-91%25-green?labelColor=2a2f35)](https://github.com/GDeLaurentis/pyadic/actions)
[![PyPI Downloads](https://img.shields.io/pypi/dm/pyadic.svg?label=PyPI%20downloads)](https://pypi.org/project/pyadic/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GDeLaurentis/pyadic/HEAD)


## Requirements
```
numpy, sympy
```

## Installation
```
pip install -e path/to/repo
```

## Testing

```
pytest --cov pyadic/ --cov-report html tests/ --verbose
```

## Quick Start

```
from pyadic import PAdic
from fractions import Fraction as Q
# 7/13 as a 12-digit 2147483647-adic number
PAdic(Q(7, 13), 2147483647, 12)
```

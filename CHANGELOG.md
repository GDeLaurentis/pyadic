# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2024-01-02

### Added

- Univariate Newton and Thiele interpolation algorithms, `Newton_polynomial_interpolation` and `Thiele_rational_interpolation`.
- Gaussian rationals, `GaussianRational`, moved from [lips](https://github.com/GDeLaurentis/lips).
- This changelog.

### Changed

- `vec_chained_FF_rationalize` optimized for sparse tensors. New keyword `optimize_for_sparse_arrays` defaults to `True`.

### Fixed

- Precision of `PAdic` when instantiated from negative integers proportional to the prime [Issue 3](https://github.com/GDeLaurentis/pyadic/issues/3).
- Recursion issue in rationalizaton ($\mathbb{F}_p \rightarrow \mathbb{Q}$) of tensors using `numpy.vectorize`.
- Compatibility with `numpy.uint32` and `numpy.uint64`.


## [0.1.2] - 2023-04-09

### Added

- Project description in `README.md`


## [0.1.1] - 2023-04-04

### Added

- $p$-adic numbers. `PAdic`.
- Finite fields, `ModP`.


[unreleased]: https://github.com/GDeLaurentis/pyadic/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/GDeLaurentis/pyadic/releases/tag/v0.2.0
[0.1.2]: https://github.com/GDeLaurentis/pyadic/releases/tag/v0.1.2
[0.1.1]: https://github.com/GDeLaurentis/pyadic/releases/tag/v0.1.1

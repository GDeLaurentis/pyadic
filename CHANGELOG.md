# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

### Deprecated


## [0.2.4] - 2025-04-23

### Added

- Support for `FieldExtension` of an arbitrary number of square roots. Does not check for relations among square roots. With finite fields, a single square root is sufficient, assuming a branch. This approach does not assume branches, but hides relations.

### Changed

- Finite-field object `ModP` can be instantiated from a rational string, e.g. `2/3`, when the prime is also specified.

- Improved parsing of `GaussianRational` and its string representation.

### Fixed

- Fixed issue with `hash` of `ModP` and `PAdic` causing clashes in caches with `functools.lru_cache` due to hash(integer) = integer. Hashing function now hashes the string representation of the numbers.

- Fixed issue with `PAdic` instantiation from string, where if the prime `p` and the number of digits `k` were supplied it would fail to call the `__rstr__` parser expecting a rational number even if the string was already an expansion in p.

- CI doc workflow should now fail if Sphinx autodoc fails.


## [0.2.3] - 2025-02-11

### Changed

- Python 3.13 in CI.

### Fixed

- Fixed issue with arithmetic operations between `FieldExtension` object and e.g. `numpy.ndarrays`.


## [0.2.2] - 2024-06-07

### Added

- Continuous deployment of documentation via github pages.
- DOI.

### Changed

- Python 3.10, 3.11 and 3.12 are also tested in CI.


## [0.2.1] - 2024-05-04

### Added

- Multivariate Newton interpolation algorithm, `multivariate_Newton_polynomial_interpolation`.

### Changed

- Improved compatibility of `extended_euclidean_algorithm`: output is of same type as input.
- Improved `ModP` and `PAdic` constructors to handle a wider variety of inputs, e.g. `ModP('+1', 2 ** 31 - 1)` is now valid.
- Splitting CI Test and Lint, adding automatic PyPI release workflow.

### Fixed

- Rationalisation of $p$-adic non integers [Issue 4](https://github.com/GDeLaurentis/pyadic/issues/4).
- Fixed naming of `extended_euclidean_algorithm` (was `extended_euclideal_algorithm`).


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


[unreleased]: https://github.com/GDeLaurentis/pyadic/compare/v0.2.4...HEAD
[0.2.4]: https://github.com/GDeLaurentis/pyadic/compare/v0.2.3...v0.2.4
[0.2.3]: https://github.com/GDeLaurentis/pyadic/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/GDeLaurentis/pyadic/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/GDeLaurentis/pyadic/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/GDeLaurentis/pyadic/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/GDeLaurentis/pyadic/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/GDeLaurentis/pyadic/releases/tag/v0.1.1

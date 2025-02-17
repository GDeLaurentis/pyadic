import copy
import fractions
import functools
import math
import numbers

import numpy
import sympy

from . import padic
from .field_extension import FieldExtension


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def ModPfy(func):
    @functools.wraps(func)
    def wrapper_ModPfy(self, other):
        if isinstance(other, ModP):
            if self.p != other.p:
                raise ValueError(f"Can't cast numbers between different finite fields: FF{self.p} and FF{other.p}")
            return func(self, other)
        elif isinteger(other) or isinstance(other, fractions.Fraction) or (hasattr(other, "imag") and not isinstance(other, numpy.ndarray)) or isinstance(other, padic.PAdic):
            return func(self, ModP(other, self.p))
        else:
            return NotImplemented
    return wrapper_ModPfy


def isinteger(x):
    return (isinstance(x, int) or
            type(x) in [numpy.int32, numpy.int64, numpy.uint32, numpy.uint64, sympy.Integer, sympy.core.numbers.Zero] or
            (hasattr(x, "is_integer") and callable(x.is_integer) and x.is_integer()))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class ModP(object):
    """Finite field with p elements ($\\mathbb{FF}_p$), i.e. integers modulo p, with p prime."""

    __slots__ = 'n', 'p'

    def __init__(self, n, p=None):
        if p is not None and isinteger(n) and isinteger(p):
            self.n = int(n) % int(p)
            self.p = int(p)
        elif p is not None and isinstance(n, fractions.Fraction):
            self_ = ModP(n.numerator, p) / ModP(n.denominator, p)
            self.n = self_.n
            self.p = self_.p
        elif p is not None and hasattr(n, "imag"):
            res = ModP(n.real, p)
            if n.imag != 0:
                b = ModP(n.imag, p)
                i = finite_field_sqrt(ModP(-1, p))  # this fails if p is a prime power
                res += i * b
            if isinstance(res, FieldExtension):
                raise ValueError("Complex to ModP conversion requires √-1 in field or zero imaginary part.")
            self.n = res.n
            self.p = res.p
        elif isinstance(n, ModP):
            if p is not None:
                assert p == n.p
            self.n = n.n
            self.p = n.p
        elif p is None and isinstance(n, padic.PAdic):
            self.n = int(n)
            self.p = n.p ** n._k
        elif isinstance(n, str) and (n.isnumeric() or n.lstrip("+-").isnumeric()) and p is not None:
            self.n, self.p = int(n) % int(p), p
        elif isinstance(n, str):
            self.n, self.p = self.__rstr__(n)
        else:
            raise TypeError('Bad finite field constructor, (n, p) of  value:({}, {}) and type:({}, {}).'.format(n, p, type(n), type(p)))

    def __getstate__(self):
        return (int(self), self.p)

    def __setstate__(self, state):
        self.__init__(*state)

    def __int__(self):
        return self.n

    def __abs__(self):
        return 0 if self.n == 0 else 1

    def __str__(self):
        return "%d %% %d" % (self, self.p)

    @staticmethod
    def __rstr__(string):
        if "%" in string:
            return tuple(map(int, string.replace(" ", "").split("%")))
        elif "mod" in string:
            return tuple(map(int, string.replace(" ", "").split("mod")))
        else:
            raise Exception(f"String {string} not understood")

    def __repr__(self):
        return str(self)

    def __neg__(self):
        """Unary '-' operation"""
        return ModP(self.p - int(self), self.p)

    def __pos__(self):
        """Unary '+' operation"""
        return self

    @ModPfy
    def __eq__(self, other):
        return self.n == other.n and self.p == other.p

    @ModPfy
    def __add__(self, other):
        return ModP(int(self) + int(other), self.p)

    @ModPfy
    def __radd__(self, other):
        return ModP(int(other) + int(self), self.p)

    @ModPfy
    def __sub__(self, other):
        return ModP(int(self) - int(other), self.p)

    @ModPfy
    def __rsub__(self, other):
        return ModP(int(other) - int(self), self.p)

    @ModPfy
    def __mul__(self, other):
        return ModP(int(self) * int(other), self.p)

    @ModPfy
    def __rmul__(self, other):
        return ModP(int(other) * int(self), self.p)

    @ModPfy
    def __truediv__(self, other):
        # is this really needed given the decorator?
        # if not isinstance(other, ModP):
        #     other = ModP(other, self.p)
        return self * other._inv()

    @ModPfy
    def __rtruediv__(self, other):
        return other * self._inv()

    def __pow__(self, n):
        assert (isinstance(n, int) or n.is_integer())
        if n < 0:
            return 1 / self ** -n
        elif n == 0:
            return ModP(1, self.p)
        elif n % 2 == 0:
            root_2_res = self ** (n / 2)
            return root_2_res * root_2_res
        else:
            return self * (self ** (n - 1))

    def _inv(self):
        """Find multiplicative inverse of self in Z_p (Z mod p) using the extended Euclidean algorithm."""

        s, t, gcd = extended_euclidean_algorithm(int(self), self.p)

        if gcd != 1:
            raise ZeroDivisionError("Inverse of {} mod {} does not exist. Are you sure {} is prime?".format(self, self.p, self.p))

        return ModP(s, self.p)

    def __hash__(self):
        return hash(self.n) + hash(self.p)

    def sqrt(self):
        return finite_field_sqrt(self)


numbers.Number.register(ModP)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def vec_ModP(prime, optimize_for_sparse_arrays=True):
    """Vectorized version of ModP."""
    if optimize_for_sparse_arrays is False:
        return numpy.vectorize(functools.partial(ModP, p=prime), otypes="O")
    else:
        def _vec_ModP_optimized_(tensor):
            tensor_non_zero_mask = (numpy.array(tensor) != 0)
            values_to_insert = vec_ModP(prime, optimize_for_sparse_arrays=False)(tensor[tensor_non_zero_mask])
            ModPtensor = numpy.full(tensor.shape, ModP(0, prime))
            ModPtensor[tensor_non_zero_mask] = values_to_insert
            return ModPtensor
        return _vec_ModP_optimized_


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def extended_euclidean_algorithm(a, b):
    """Returns Bezout coefficients (s,t) and gcd(a,b) such that: as+bt=gcd(a,b). - Pseudocode from https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm"""

    # This ensures that the output is of the same type as the input,
    # e.g. for working with a, b in sympy.polys.rings.PolyElement
    zero = a * 0
    one = zero + 1

    (old_r, r) = (a, b)
    (old_s, s) = (one, zero)
    (old_t, t) = (zero, one)

    while r != 0:
        quotient = old_r // r
        (old_r, r) = (r, old_r - quotient * r)
        (old_s, s) = (s, old_s - quotient * s)
        (old_t, t) = (t, old_t - quotient * t)

    # output "Bézout coefficients:", (old_s, old_t)
    # output "greatest common divisor:", old_r
    # output "quotients by the gcd:", (t, s)

    return (old_s, old_t, old_r)


def gcd(a, b):
    _, _, gcd = extended_euclidean_algorithm(a, b)
    return gcd


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def MQRR(u, m, T=None):
    """Maximal Quotient Rational Reconstruction (M. B. Monagan)"""
    if T is None:
        c = 1
        T = 2 ** c * math.ceil(math.log(m, 2))
    if u == 0:
        if m > T:
            return 0
        else:
            return False
    (n, d) = (0, 0)
    (t0, r0) = (0, m)
    (t1, r1) = (1, u)
    while r1 != 0 and r0 > T:
        q = math.floor(r0 / r1)
        if q > T:
            (n, d, T) = (r1, t1, q)
        (r0, r1) = (r1, r0 - q * r1)
        (t0, t1) = (t1, t0 - q * t1)
    if d < 0:
        (n, d) = (-n, -d)
    if d == 0 or gcd(n, d) not in (1, d):
        return False
    return fractions.Fraction(n, d)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def LGreduction(u, v):
    """Lattice Gaussian Reduction - 2D version of LLL (Lenstra–Lenstra–Lovász)"""
    u, v = numpy.array(u, dtype=object), numpy.array(v, dtype=object)
    if v @ v > u @ u:
        return LGreduction(v, u)
    while v @ v < u @ u:
        (u, v) = (v, u)
        q = round(fractions.Fraction(u @ v, u @ u))
        v = v - q * u
    return (u, v)


def LGRR(a, n):
    """Lattice Gaussian Rational Reconstruction"""
    return fractions.Fraction(*LGreduction((a, 1), (n, 0))[0])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def rationalise(a, n=None, algorithm=(LGRR, MQRR)[1]):
    """Given (a, n) returns a fraction r / s such that r/s % n = a, by lattice reduction. r = sa + mn  <-> r/s % n = a"""
    if n is None:  # for FF argument
        if isinstance(a, int):
            return fractions.Fraction(a, 1)
        elif isinstance(a, ModP):
            return rationalise(int(a), a.p, algorithm)
        elif isinstance(a, padic.PAdic):
            return rationalise(int(a), a.p ** a.k, algorithm) * fractions.Fraction(a.p) ** a.n
    return algorithm(a, n)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def chinese_remainder(a1, a2):
    """Given a1 = a % n1 and a2 = a % n2 and assuming gcd(n1,n2)=1 (i.e. n1, n2 co-prime), returns a12 = a % (n1*n2)"""
    a1, n1 = int(a1), a1.p
    a2, n2 = int(a2), a2.p
    q1, q2, gcd = extended_euclidean_algorithm(n1, n2)
    assert gcd == 1
    return ModP(a1 * (q2 * n2) + a2 * (q1 * n1), n1 * n2)


def chained_chinese_remainder(*vals, primes=None):
    """Vectorized (concatenated) version of chinese remainder function."""
    if not (all([isinstance(val, ModP) for val in vals]) or len(vals) == len(primes)):
        raise Exception("Unrecognized input.")
    if primes is not None:
        vals = tuple(map(lambda x: ModP(*x), list(zip(vals, primes))))
    if len(vals) == 1:
        return vals[0]
    res = chinese_remainder(*vals[:2])
    for val in vals[2:]:
        res = chinese_remainder(res, val)
    return res


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def vec_chained_FF_rationalize(tensors, primes, factor=1, algorithm=LGRR, optimize_for_sparse_arrays=True):
    """Given list of tensors respectively mod primes, returns the rationalized tensor.
       Keyword argument factor is pre-multiplied and then devided w.r.t. reconstruction,
       can be used to aid reconstruction or verify its stability."""
    if len(tensors) != len(primes):
        raise AssertionError("As many tensors as primes must be supplied for rationalisation.")
    if not len(set([tensor.shape for tensor in tensors])) == 1:
        raise AssertionError("All tensors must have the same shape for rationalisation.")
    if optimize_for_sparse_arrays:  # optimization for sparse arrays: run the algorithm only on non-zero values.
        tensors_non_zero_mask = numpy.all(numpy.array([tensor != 0 for tensor in tensors]), axis=0)
        values_to_insert = vec_chained_FF_rationalize([tensor[tensors_non_zero_mask] for tensor in tensors], primes,
                                                      factor=factor, algorithm=algorithm, optimize_for_sparse_arrays=False)
        Qtensor = numpy.zeros_like(tensors[0], dtype='O')
        Qtensor[tensors_non_zero_mask] = values_to_insert
    else:
        vec_chained_chinese_remainder = numpy.vectorize(functools.partial(chained_chinese_remainder, primes=primes), otypes=[object])
        chained_tensors = vec_chained_chinese_remainder(*tensors)
        vec_rationalize = numpy.vectorize(functools.partial(rationalise, algorithm=algorithm), otypes="O")
        Qtensor = vec_rationalize(factor * chained_tensors) / fractions.Fraction(factor)
    return Qtensor


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@functools.lru_cache
def finite_field_sqrt(x):
    """Returns either False or the first digit of the root in the field."""
    assert isinstance(x, ModP)
    root = univariate_finite_field_solver(f"x^2-{int(x)}", dict(), x.p)
    if root is False:
        return FieldExtension(x)
    else:
        return ModP(int(root[0][sympy.symbols('x')]), x.p)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def univariate_finite_field_solver(equation, root_dict, prime):
    # !!! DO NOT MODIFY HERE !!! This function is from lips.algebraic_geometry.tools
    """Returns all possible solutions of 'equation' over a finite field of cardinality 'prime'.
       If already satisfied returns True, if no solution exists returns False."""
    equation = sympy.sympify(equation).subs(root_dict)
    if isinstance(equation, sympy.core.numbers.Integer) and equation % prime == 0:
        return True
    free_symbols = list(equation.free_symbols)
    if len(free_symbols) < 1:
        return False
    if len(free_symbols) > 1:
        raise Exception("Too many free parameters.")
    symbol = free_symbols[0]
    equation = sympy.poly(equation, modulus=prime)
    if equation == 0:
        return True
    pre_factor, factors = sympy.factor_list(sympy.factor(equation, modulus=prime))
    factors = [factor[0] for factor in factors]
    if pre_factor % prime == 0:
        return True
    linear_factors = [factor for factor in factors if sympy.diff(factor, symbol) == 1]
    if linear_factors == []:
        return False
    solutions = [ModP(int(sympy.solve(factor)[0]), prime) for factor in linear_factors]
    return update_root_dict(symbol, solutions, root_dict)


def update_root_dict(symbol, solutions, root_dict):
    # !!! DO NOT MODIFY HERE !!! this function is from lips.algebraic_geometry.tools
    """Given solutions and root_dict returns updated root_dicts."""
    root_dicts = [copy.deepcopy(root_dict)]
    root_dicts[0].update({symbol: solutions[0]})
    for solution in solutions[1:]:
        new_root_dict = copy.deepcopy(root_dict)
        new_root_dict.update({symbol: solution})
        root_dicts.append(new_root_dict)
    return root_dicts

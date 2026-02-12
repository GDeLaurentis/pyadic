import functools
import math
import numbers
import numpy
import random
import re

from fractions import Fraction as Q

from .finite_field import ModP, finite_field_sqrt, isinteger
from .field_extension import FieldExtension

fixed_relative_precision = False
all_precision_loss_warning = False


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def to_base(num, p):
    if num < p:
        return (num, )
    else:
        return (num % p, ) + (to_base(num // p, p))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def cached_function(func):
    @functools.wraps(func)
    def wrapper_cached_function(*args, **kwargs):
        call_key = (args, tuple(sorted(kwargs.items())))
        if not hasattr(func, 'cache_dict'):
            func.cache_dict = {}
        if call_key not in func.cache_dict.keys():
            call_val = func(*args, **kwargs)
            func.cache_dict[call_key] = call_val
        return func.cache_dict[call_key]
    return wrapper_cached_function


@cached_function
def full_range_random_padic_filling(p, k):
    """Returns a random number which translates in a PAdic without any zero digit."""
    return sum([random.randint(1, p - 1) * p ** i for i in range(k)])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def padicfy(func):
    @functools.wraps(func)
    def wrapper_padicfy(self, other):
        if isinstance(other, PAdic):
            if self.p != other.p:
                raise ValueError(f"Can't cast a {other.p}-adic to a {self.p}-adic.")
            return func(self, other)
        elif isinteger(other) or isinstance(other, ModP) or isinstance(other, Q) or (hasattr(other, "imag") and not isinstance(other, numpy.ndarray)):
            if func.__name__ in ["__mul__", "__rmul__"] and (isinteger(other) or isinstance(other, Q)) and other == 0:
                return 0
            return func(self, PAdic(other, self.p, max((self.n + self.k, self.k))))
        else:
            return NotImplemented
    return wrapper_padicfy


def check_orderable(func):
    @functools.wraps(func)
    def wrapper_check_orderable(self, other):
        if not isinstance(other, type(self)):
            raise TypeError("unorderable types: {} < {}".format(type(self), type(other)))
        if self.p != other.p:
            raise Exception("unorderable padics over different primes: {}, {}.".format(self.p, other.p))
        return func(self, other)
    return wrapper_check_orderable


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class PAdic(object):

    """PAdic numbers, with p prime, k digits, n valuation."""

    # This allows the intantiation step to return a field extension,
    # but doubles instantiation time, even when no field extension is required.
    # def __new__(cls, num, p=None, k=None, n=0, from_addition=False):
    #     try:
    #         PAdic_obj = super(PAdic, cls).__new__(cls)
    #         PAdic_obj.__init__(num, p, k, n, from_addition)
    #         return PAdic_obj
    #     except ValueError:
    #         if hasattr(num, "imag"):
    #             a = PAdic(num.real, p, k, n, from_addition)
    #             b = PAdic(num.imag, p, k, n, from_addition)
    #             i = padic_sqrt(PAdic(-1, p, k, n, from_addition))
    #             return a + i * b  # this is in a FieldExtension now
    #         else:
    #             raise Exception("Invalid p-adic initialisation")

    def __init__(self, num, p=None, k=None, n=0, from_addition=False):
        """0 ≤ num ≤ p ^ k - 1; p: prime; k: significant digits; n: power of prefactors of p (valuation)."""
        if isinteger(num) or isinstance(num, ModP):
            num = int(num)  # might get passed as FF number
            self.p = p
            factors_of_p = next((i for i, j in enumerate(to_base(abs(num), p)) if j != 0), None)
            if factors_of_p is None:  # leading zeros are not significant digits under any situation
                factors_of_p = k
                from_addition = True
            self.k = k - from_addition * factors_of_p  # this is the guaranteed precision
            num = int(num // p ** factors_of_p) % (p ** self.k)
            self.n = factors_of_p + n
            self.num = num
            if fixed_relative_precision is True and factors_of_p > 0:
                # this emulates the behaviour of floating point numbers,
                # where precision loss actually means random digits get added.
                self.num = self.num + p ** self.k * full_range_random_padic_filling(self.p, factors_of_p)
                self.k = self.k + factors_of_p
            if all_precision_loss_warning and self.k == 0:
                print("Lost all precision @", self)
        elif isinstance(num, Q):
            res = PAdic(num.numerator, p, k, n, from_addition) / PAdic(num.denominator, p, k, n, from_addition)
            self.num, self.p, self.k, self.n = res.num, res.p, res.k, res.n
        elif isinstance(num, PAdic):
            self.num, self.p, self.k, self.n = num.num, num.p, num.k, num.n
        elif isinstance(num, float):
            if not (Q(num).limit_denominator(10 ** 3) == Q(num).limit_denominator(10 ** 5) == Q(num).limit_denominator(10 ** 8)):
                raise ValueError(f"PAdic instantiation from float failed. It's unclear to which rational {num} corresponds to.")
            res = PAdic(Q(num).limit_denominator(10 ** 3), p, k, n, from_addition)
            self.p = res.p
            self.k = res.k
            self.n = res.n
            self.num = res.num
        elif hasattr(num, "imag"):
            res = PAdic(num.real, p, k, n, from_addition)
            if num.imag != 0:
                b = PAdic(num.imag, p, k, n, from_addition)
                i = padic_sqrt(PAdic(-1, p, k, n, from_addition))
                res += + i * b
            if isinstance(res, FieldExtension):
                raise ValueError(f"Can't create {p}-adic from {num}. A field extension is required.")
            self.p = res.p
            self.k = res.k
            self.n = res.n
            self.num = res.num
        elif isinstance(num, str) and ('O(' in num or (p is None and k is None)):
            self.num, self.p, self.k, self.n = self.__rstr__(num)
        elif isinstance(num, str):
            num = Q(num)
            res = PAdic(num.numerator, p, k, n, from_addition) / PAdic(num.denominator, p, k, n, from_addition)
            self.num, self.p, self.k, self.n = res.num, res.p, res.k, res.n
        else:
            raise Exception(f"Invalid p-adic initialisation: {num}, {p}, {k}, {n}, {from_addition}.")

    # GETTERS and SETTERS

    @property
    def num(self):
        """0 ≤ num ≤ p ^ k - 1"""
        return self._num

    @num.setter
    def num(self, value):
        if value < 0:
            raise ValueError("Padic num should be non-negative")
        self._num = value

    @property
    def p(self):
        """p: prime for the padic."""
        return self._p

    @p.setter
    def p(self, value):
        if value <= 0:
            raise ValueError("Padic p should be positive, got: {}.".format(value))
        self._p = value

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, value):
        """k: number of significant digits."""
        if value < 0:
            value = 0
        self._k = value

    @property
    def as_tuple(self):
        """Tuple reprensentation of the mantissa."""
        return (to_base(self.num, self.p) + tuple([0 for i in range(self.k)]))[:self.k]

    @property
    def as_dict(self):
        """Dict reprensentation: {valuation: mantissa digit}."""
        return dict(zip(range(self.n, self.k), self.as_tuple))

    @property
    def as_tuple_from_zero(self):
        """Tuple representation of the mantissa, shifted to start from p^0 term."""
        if self.n < 0:
            raise ValueError("PAdic tuple from zero representation is only defined for PAdic integers (non-negative valuation).")
        return (0, ) * self.n + self.as_tuple

    def __getstate__(self):
        return (self.num, self.p, self.k, self.n)

    def __setstate__(self, state):
        self.__init__(*state)

    def __str__(self):
        if self.k == 0:
            if self.n == 0:
                return "O(1)"
            elif self.n == 1:
                return f"O({self.p})"
            else:
                return f"O({self.p}^{self.n})"
        else:
            return (" + ".join(filter(lambda x: x is not None,
                                      ["{}".format(i) if (j == 0 and i != 0) else
                                       "{}*{}".format(i, self.p) if (j == 1 and i != 0) else
                                       "{}*{}^{}".format(i, self.p, j) if (i != 0) else None
                                       for i, j in zip(self.as_tuple, range(self.n, self.n + self.k))])) +
                    (" + O(1)" if self.n + self.k == 0 else f" + O({self.p})" if self.n + self.k == 1 else f" + O({self.p}^{self.n + self.k})"))

    def __repr__(self):
        return str(self)

    @staticmethod
    def __rstr__(string):
        """Constructor from string (inverse method to __str__ or __repr__)."""
        string = string.replace(" ", "")
        # get the prime
        prime = int(re.findall(r"O\((\d+)", string)[0])
        if prime == 1:  # for the case with + O(1)
            match = re.findall(r"\*(\d+)\^", string)
            if match != []:
                prime = int(match[0])
        # get the valuation
        valuation = string.split("+")[0]
        valuation = re.findall(rf"{prime}[\^\-\d+]*", valuation)
        if valuation == []:
            valuation = 0
        else:
            if "^" in valuation[0]:
                valuation = int(valuation[0].split("^")[1])
            else:
                valuation = 1
        # get the mantissa
        mantissa = [int(re.sub(r"\*[\^\-\d+]{0,}", "", entry.replace(f"{prime}", ""))) for entry in string.split("+")[:-1]]
        significant_digits = len(mantissa)
        mantissa = sum([entry * prime ** i for i, entry in enumerate(mantissa)])
        return (mantissa, prime, significant_digits, valuation)

    # COMPARISON

    @padicfy
    def __eq__(self, other):
        return all([self.num == other.num, self.p == other.p, self.k == other.k, self.n == other.n])

    @padicfy
    def __le__(self, other):
        return self.n >= other.n

    @padicfy
    def __lt__(self, other):
        return self.n > other.n

    @padicfy
    def __ge__(self, other):
        return self.n <= other.n

    @padicfy
    def __gt__(self, other):
        return self.n < other.n

    # ALGEBRA

    def __int__(self):
        if self.n < 0:
            raise ValueError(f"Only p-adic integers can be converted to int, not p-adic numbers with negative valuation. Received valuation {self.n}.")
        return self.num * self.p ** self.n

    def __abs__(self):
        return PAdic(0, self.p, 0, self.n)

    @padicfy
    def __add__(self, other):
        if self.n > other.n:
            return other + self
        else:
            return PAdic((self.num + other.num * self.p ** (other.n - self.n)), self.p,
                         self.k if self.k < (other.n - self.n) + other.k else (other.n - self.n) + other.k, self.n, from_addition=True)

    @padicfy
    def __radd__(self, other):
        return other + self

    @padicfy
    def __sub__(self, other):
        return self + (- other)

    @padicfy
    def __rsub__(self, other):
        return - (self - other)

    @padicfy
    def __mul__(self, other):
        return PAdic((self.num * other.num) % self.p ** self.k, self.p, min([self.k, other.k]), self.n + other.n)

    @padicfy
    def __rmul__(self, other):
        return self * other

    @padicfy
    def __truediv__(self, other):
        return PAdic(int(self.num * ModP(other.num, other.p ** other.k)._inv()) % self.p ** self.k, self.p, min([self.k, other.k]), self.n - other.n)

    @padicfy
    def __div__(self, other):
        return self.__truediv__(other)

    @padicfy
    def __rtruediv__(self, other):
        return other / self

    @padicfy
    def __rdiv__(self, other):
        return self.__rtruediv__(other)

    def __neg__(self):
        """Unary '-' operation"""
        return PAdic((-1 * self.num) % self.p ** self.k, self.p, self.k, self.n)

    def __pos__(self):
        """Unary '+' operation"""
        return self

    def __pow__(self, n):
        assert (isinstance(n, int) or n.is_integer())
        if n < 0:
            return 1 / self ** -n
        elif n == 0:
            return PAdic(1, self.p, self.k)
        elif n % 2 == 0:
            root_2_res = self ** (n / 2)
            return root_2_res * root_2_res
        else:
            return self * (self ** (n - 1))

    def __hash__(self):
        return hash(str(self))


numbers.Number.register(PAdic)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def padic_log(w, base=None):
    """
    If the valuation is 0, then bring in radius of convergence using Fermat’s little theorem (w ** (p - 1) = 1 mod p),
    i.e. log_p(w) = log_p(w^(p - 1)) / (p - 1), with w^(p - 1) = 1 + x, x ~ O(p).
    Otherwise, factor out p ^ valuation, use that log(a * b) = log(a) + log(b) and that log_p(p ^ valuation) = valuation
    What about branch? See sage docs.
    """
    if w == 1:    # for compatibility with integers / exact (infinite-precision) p-adic
        return 0
    if base is None:
        return padic_log(w, w.p)
    if base is w.p:
        if w.n == 0:
            x = w ** (w.p - 1) - 1
            return sum([(-1) ** (n + 1) * x ** n / n for n in range(1, w.k)]) / (w.p - 1)
        elif w.n > 0:
            return w.n + padic_log(w / w.p ** w.n, w.p)
        else:
            return w.n + padic_log(w * w.p ** -w.n, w.p)
    else:
        return padic_log(w, base=w.p) / padic_log(PAdic(base, w.p, w.k), base=w.p)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def refine_sqrt_precision(x, s):
    """Given (s | s^2 - x << 1), makes s^2 closer to x."""
    return s + (x - s ** 2) / (2 * s)


@functools.lru_cache
def padic_sqrt(x):
    """Working precision padic sqrt."""
    assert isinstance(x, PAdic)
    ffx = ModP(x.as_tuple[0], x.p)
    root = finite_field_sqrt(ffx)
    if isinstance(root, FieldExtension):
        return FieldExtension(x)
    if x.n % 2 != 0:  # sqrt(x) with x ~ O(p^(odd power))
        raise NotImplementedError("Unramified field extension")
    root = PAdic(int(root), x.p, x.k, x.n // 2)
    for i in range(math.ceil(math.log(x.k, 2))):
        root = refine_sqrt_precision(x, root)
    return root

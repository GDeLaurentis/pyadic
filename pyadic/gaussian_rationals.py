import random
import mpmath
import functools

from fractions import Fraction


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def rand_rat_frac():
    return Fraction(random.randrange(-100, 101), random.randrange(1, 201))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def to_gaussian_rational(other):
    if type(other) is GaussianRational:
        return other
    elif type(other) is int:
        return GaussianRational(Fraction(other, 1), Fraction(0, 1))
    elif type(other) is Fraction:
        return GaussianRational(other, Fraction(0, 1))
    elif type(other) is complex:
        assert other.real.is_integer() and other.imag.is_integer()
        return GaussianRational(Fraction(int(other.real), 1), Fraction(int(other.imag), 1))
    elif type(other) is mpmath.mpc:
        assert mpmath.isint(other, gaussian=True)
        return GaussianRational(Fraction(int(other.real), 1), Fraction(int(other.imag), 1))
    else:
        raise Exception("Can't convert {} of type {} to GaussianRational.".format(other, type(other)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def gaussian_rational(func):
    @functools.wraps(func)
    def wrapper_gaussian_rational(arg1, other):
        other = to_gaussian_rational(other)
        return func(arg1, other)
    return wrapper_gaussian_rational


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class GaussianRational(object):

    __slots__ = ('real', 'imag')

    def __init__(self, fraction_real, fraction_imag=0):
        if isinstance(fraction_real, (int, str)):
            fraction_real = Fraction(fraction_real)
        if isinstance(fraction_imag, (int, str)):
            fraction_imag = Fraction(fraction_imag)
        if not (isinstance(fraction_real, Fraction) and isinstance(fraction_imag, Fraction)):
            raise ValueError(f"Gaussian Rational invalid initialisation {fraction_real} {fraction_imag}")
        self.real = fraction_real
        self.imag = fraction_imag

    def __str__(self):
        if self.imag != 0:
            return f"({self.real}+{self.imag}j)"
        else:
            return f"{self.real}"

    def __repr__(self):
        return f"GaussianRational({self.real}, {self.imag})"

    @gaussian_rational
    def __eq__(self, other):
        return self.real == other.real and self.imag == other.imag

    @gaussian_rational
    def __add__(self, other):
        return GaussianRational(self.real + other.real, self.imag + other.imag)

    @gaussian_rational
    def __radd__(self, other):
        return other + self

    @gaussian_rational
    def __sub__(self, other):
        return GaussianRational(self.real - other.real, self.imag - other.imag)

    @gaussian_rational
    def __rsub__(self, other):
        return other - self

    @gaussian_rational
    def __mul__(self, other):
        return GaussianRational(self.real * other.real - self.imag * other.imag, self.real * other.imag + self.imag * other.real)

    @gaussian_rational
    def __rmul__(self, other):
        return other * self

    @gaussian_rational
    def __truediv__(self, other):
        mod_other = other * other.conjugate()
        return self * other.conjugate() * (1 / mod_other.real)

    @gaussian_rational
    def __rtruediv__(self, other):
        return other / self

    def __neg__(self):
        return -1 * self

    def __pow__(self, n):
        assert isinstance(n, int) or n.is_integer()
        if n == 0:
            return 1
        elif n % 2 == 0:
            root_2_res = self ** (n / 2)
            return root_2_res * root_2_res
        else:
            return self * (self ** (n - 1))

    def __complex__(self):
        return complex(self.real, self.imag)

    def __abs__(self):
        abs_self = self * self.conjugate()
        assert abs_self.imag == 0
        return abs_self.real

    def conjugate(self):
        return GaussianRational(self.real, - self.imag)

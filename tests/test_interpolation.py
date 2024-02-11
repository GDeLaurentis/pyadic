import sympy

from pyadic.interpolation import Newton_polynomial_interpolation, Thiele_rational_interpolation, \
    multivariate_Newton_polynomial_interpolation

t = sympy.symbols('t')
t1, t2, t3 = sympy.symbols('t1:4')


def Ptest0(tval):
    return (1)


def Ptest1(tval):
    return (tval ** 20 + tval - 1)


def Ptest2(t1, t2, t3):
    return (t1 + 2 * t2) + 3 * (t1 * t2) + t3 ** 5


def Rtest0(tval):
    return (1)


def Rtest1(tval):
    return (tval ** 20 + tval - 1) / (5 - tval)


def test_Newton_polynomial_interpolation_trivial():
    assert Newton_polynomial_interpolation(Ptest0, 2 ** 31 - 1, verbose=True) == 1


def test_Newton_polynomial_interpolation():
    assert Newton_polynomial_interpolation(Ptest1, 2 ** 31 - 1, verbose=True) == t ** 20 + t - 1


def test_Newton_polynomial_interpolation_multivariate():
    assert multivariate_Newton_polynomial_interpolation(Ptest2, 2 ** 31 - 1, verbose=True) == (t1 + 2 * t2) + 3 * (t1 * t2) + t3 ** 5


def test_Thiele_rational_interpolation_trivial():
    assert Thiele_rational_interpolation(Rtest0, 2 ** 31 - 61, verbose=True) == 1


def test_Thiele_rational_interpolation():
    assert Thiele_rational_interpolation(Rtest1, 2 ** 31 - 61, verbose=True) == (t ** 20 + t - 1) / (5 - t)


def test_Thiele_rational_interpolation_as_continued_fractionl():
    Thiele_rational_interpolation(Rtest1, 2 ** 31 - 1, as_continued_fraction=True, verbose=True)

import pickle
import random

from fractions import Fraction as Q

import pyadic

from pyadic import PAdic
from pyadic.padic import padic_sqrt, padic_log
from pyadic.finite_field import rationalise, LGRR, MQRR


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_picklable():
    obj = PAdic(3 + 2 * 7, 7, 3)
    mydump = pickle.dumps(obj, protocol=2)
    loaded = pickle.loads(mydump)
    assert obj == loaded


def test_instantiation_from_complex_when_in_fied():
    assert PAdic(1j, 2 ** 31 - 19, 5)


def test_instantiation_from_string():
    assert PAdic("3 + O(7)") == PAdic(3, 7, 1)


def test_str_non_zero():
    p, k = 10007, 3
    assert str(PAdic(1 + 2 * p + 3 * p ** 2, p, k)) == "1 + 2*{p} + 3*{p}^2 + O({p}^{k})".format(p=p, k=k)


def test_addition():
    p, k = 10007, 3
    a, b = random.randrange(0, 10000), random.randrange(0, 10000)
    assert PAdic(a + b, p, k) == PAdic(a, p, k) + PAdic(b, p, k)
    assert PAdic(a + b, p, k) == a + PAdic(b, p, k)
    assert PAdic(a + b, p, k) == PAdic(a, p, k) + b


def test_addition_with_zero():
    # positive valuation
    p, k, n = 10007, 3, 2
    a = PAdic(random.randrange(0, 10000), p, k, n)
    assert a + 0 == a
    # zero valuation
    p, k, n = 10007, 3, 0
    a = PAdic(random.randrange(0, 10000), p, k, n)
    assert a + 0 == a
    # negative valuation
    p, k, n = 10007, 3, -2
    a = PAdic(random.randrange(0, 10000), p, k, n)
    assert a + 0 == a


def test_multiplication_by_one():
    # positive valuation
    p, k, n = 10007, 3, 2
    a = PAdic(random.randrange(0, 10000), p, k, n)
    assert a * 1 == a
    # zero valuation
    p, k, n = 10007, 3, 0
    a = PAdic(random.randrange(0, 10000), p, k, n)
    assert a * 1 == a
    # negative valuation
    p, k, n = 10007, 3, -2
    a = PAdic(random.randrange(0, 10000), p, k, n)
    assert a * 1 == a


def test_unary_neg():
    p, k, n = 10007, 3, 2
    a = PAdic(random.randrange(0, 10000), p, k, n)
    assert - a == -1 * a


def test_unary_pos():
    p, k, n = 10007, 3, 2
    a = PAdic(random.randrange(0, 10000), p, k, n)
    assert + a == 1 * a


def test_sqrt_in_field_and_negative_pow():
    p, k, n = 10007, 10, 0
    a = PAdic(random.randrange(0, 10000), p, k, n)
    b = a ** -2
    c = 1 / padic_sqrt(b)
    assert a == c or a == -c


def test_absolute_value():
    p = 2 ** 31 - 1
    assert PAdic(1, p, 3, -1) > PAdic(1, p, 3, 0) > PAdic(1, p, 3, 1)
    assert PAdic(1, p, 3, 1) < PAdic(1, p, 3, 0) < PAdic(1, p, 3, -1)
    assert abs(PAdic(1, p, 3, 0)) == abs(PAdic(2, p, 3, 0))


def test_rationalisation_padic():
    assert rationalise(PAdic(Q(7, 13), 2147483647, 12), algorithm=LGRR) == Q(7, 13)
    assert rationalise(PAdic(Q(7, 13), 2147483647, 12), algorithm=MQRR) == Q(7, 13)


def test_fixed_relative_precision_is_more_stable_but_not_O_guaranteed():
    p = 2 ** 31 - 1
    x = PAdic(1, p, 3)
    y = PAdic(1 - p, p, 3)
    pyadic.padic.fixed_relative_precision = False
    assert (1 / (x - y) - 1 / (x - y)).n == 1                   # fewer digits in result
    assert str(1 / (x - y) - 1 / (x - y)) == "O(2147483647)"    # but precision of result is correctly tracked to O(p)
    pyadic.padic.fixed_relative_precision = True
    assert (1 / (x - y) - 1 / (x - y)).n == 2                   # this has one more correct digit O(p^2)
    assert str(1 / (x - y) - 1 / (x - y)) != "O(2147483647^2)"  # but from O(p^2) included it is filled with random numbers
    pyadic.padic.fixed_relative_precision = False               # reset to default value (would be better with context manager in case of assert failue)


def test_padic_log():
    p = 2 ** 31 - 1
    w = PAdic(2, p, 5, 0)
    # compare to value obtained from sage (this also tests __rstr__)
    assert padic_log(w) == PAdic("2078209981*2147483647 + 2112846813*2147483647^2 + 2089755591*2147483647^3 + 1033332184*2147483647^4 + O(2147483647^5)")
    w = PAdic(2, p, 5, -3)
    padic_log(w)  # branch info?
    w = PAdic(2, p, 5, +3)
    padic_log(w)  # branch info?

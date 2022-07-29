import pickle
import random

from pyadic import PAdic
from pyadic.padic import padic_sqrt


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_picklable():
    obj = PAdic(3 + 2 * 7, 7, 3)
    mydump = pickle.dumps(obj, protocol=2)
    loaded = pickle.loads(mydump)
    assert obj == loaded


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


def test_sqrt_in_field():
    p, k, n = 10007, 10, 0
    a = PAdic(random.randrange(0, 10000), p, k, n)
    b = a ** 2
    c = padic_sqrt(b)
    assert a == c or a == -c

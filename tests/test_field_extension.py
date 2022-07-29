import pickle

from fractions import Fraction as Q

from pyadic.padic import PAdic, padic_sqrt
from pyadic.field_extension import FieldExtension


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_picklable():
    square = PAdic(8926681701808096183360073083313707792198548824, 2 ** 31 - 1, 5)
    sqrt = padic_sqrt(square)
    assert isinstance(sqrt, FieldExtension)
    mydump = pickle.dumps(sqrt, protocol=2)
    loaded = pickle.loads(mydump)
    assert sqrt == loaded


def test_inverse1():
    square = PAdic(8926681701808096183360073083313707792198548824, 2 ** 31 - 1, 5)
    inverse = 1 / FieldExtension(square, (square * 0 + 1, square * 0 - 1))
    assert inverse.tuple[0] == inverse.tuple[1]


def test_inverse2():
    square = PAdic(8926681701808096183360073083313707792198548824, 2 ** 31 - 1, 5)
    sqrt = padic_sqrt(square)
    assert 1 / sqrt * sqrt == square * 0 + 1


def test_multiplication_and_addition():
    square = PAdic(8926681701808096183360073083313707792198548824, 2 ** 31 - 1, 5)
    sqrt = padic_sqrt(square)
    assert 4 * square == (sqrt + sqrt) ** 2


def test_inverse_power():
    square = PAdic(8926681701808096183360073083313707792198548824, 2 ** 31 - 1, 5)
    sqrt = padic_sqrt(square)
    assert 1 / sqrt == sqrt ** - 1


def test_power():
    square = PAdic(8926681701808096183360073083313707792198548824, 2 ** 31 - 1, 5)
    sqrt = padic_sqrt(square)
    assert sqrt ** 3 == sqrt * sqrt * sqrt


def test_arithmetics():
    j = padic_sqrt(PAdic(-1, 2 ** 31 - 1, 5))
    a = PAdic(Q(3, 5), 2 ** 31 - 1, 5)
    b = PAdic(Q(7, 13), 2 ** 31 - 1, 5)
    c = PAdic(Q(2, 3), 2 ** 31 - 1, 5)
    assert (a * j - j * b) / c / (+ j) == (a - b) / c

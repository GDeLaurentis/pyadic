import pickle
import random
import pytest

from fractions import Fraction

from pyadic import ModP, PAdic
from pyadic.finite_field import extended_euclideal_algorithm, rationalise, MQRR, LGRR, finite_field_sqrt
from pyadic.field_extension import FieldExtension


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_picklable():
    obj = ModP('2 % 10007')
    mydump = pickle.dumps(obj, protocol=2)
    loaded = pickle.loads(mydump)
    assert obj == loaded

def test_same_class_instantiation_and_unary_plus():
    assert ModP(ModP(123, 2 ** 31 - 19)) == + ModP(123, 2 ** 31 - 19)


def test_instantiation_from_padic():
    assert ModP(PAdic(123, 2 ** 31 - 19, 3)) == ModP(123, (2 ** 31 - 19) ** 3)


def test_instantiation_with_complex():
    p = 2 ** 31 - 19
    assert ModP(1, p) * (2 + 3j) == ModP(2, p) + finite_field_sqrt(ModP(-1, p)) * ModP(3, p) == ModP(2 + 3j, p)


def test_instantiation_with_complex_requires_field_extension():
    p = 2 ** 31 - 1
    with pytest.raises(ValueError):
        ModP(2 + 3j, p)

def test_nonsense_instantiation():
    with pytest.raises(TypeError):
        ModP((), [])


def test_addition():
    p = 10007
    a, b = random.randrange(0, 10000), random.randrange(0, 10000)
    assert (a + b) % p == ModP(a, p) + ModP(b, p)
    assert (a + b) % p == a + ModP(b, p)
    assert (a + b) % p == ModP(a, p) + b


def test_failed_operation_different_FFs():
    a = ModP('2 % 3')
    b = ModP('2 % 5')
    with pytest.raises(ValueError):
        a + b


def test_addition_with_fraction_is_symmetric():
    a = ModP('2 % 5')
    b = Fraction(1, 3)
    assert a + b == b + a
    assert isinstance(a + b, ModP) and isinstance(b + a, ModP)


def test_subtraction():
    p = 10007
    a, b = random.randrange(0, 10000), random.randrange(0, 10000)
    assert (a - b) % p == ModP(a, p) - ModP(b, p)
    assert (a - b) % p == a - ModP(b, p)
    assert (a - b) % p == ModP(a, p) - b


def test_multiplication():
    p = 10007
    a, b = random.randrange(0, 10000), random.randrange(0, 10000)
    assert (a * b) % p == ModP(a, p) * ModP(b, p)
    assert (a * b) % p == a * ModP(b, p)
    assert (a * b) % p == ModP(a, p) * b

def test_division():
    p = 10007
    a, b = random.randrange(1, p - 1), random.randrange(1, p - 1)
    assert ModP(a, p) / b == ModP(Fraction(a, b), p)


def test_negation():
    p = 10007
    a = random.randrange(0, 10000)
    assert (-a) % p == -ModP(a, p)


def test_inverse():
    p = 10007
    a = ModP(random.randrange(1, p), p)
    assert a * a._inv() == 1


def test_rstr():
    assert ModP('2 mod 4') == ModP('2 % 4')


def test_failed_inverse():
    a = ModP('2 % 4')
    with pytest.raises(ZeroDivisionError):
        a * a._inv()

def test_pow():
    p = 10007
    x = random.randrange(1, p)
    assert ModP(x, p) ** -13 == 1 / ModP(x ** 13, p)


def test_hash():
    hash(ModP(12345, 10007))


def test_extended_euclideal_algorithm():
    for i in range(10):
        a, b = random.randint(1, 1000), random.randint(1, 1000)
        s, t, gcd = extended_euclideal_algorithm(a, b)
        assert(a * s + b * t == gcd)


def test_reconstruction_MQRR_2147483647():
    assert rationalise(298260199, 2147483647, algorithm=MQRR) == Fraction(-51071, 36)


def test_reconstruction_LGRR_2147483647():
    assert rationalise(298260199, 2147483647, algorithm=LGRR) == Fraction(11326, 42041)


def test_reconstruction_LGRR_2147483647_pow12():
    assert rationalise(4479489461410435237106627746985045825552416200368757674235908629372482248377113245126581341636163272702,
                       9619630365287747226839050681966839463919428531629782475127367001763589500187642982976435876178980851951693987841, algorithm=MQRR) == -1
    assert rationalise(4479489461410435237106627746985045825552416200368757674235908629372482248377113245126581341636163272702,
                       9619630365287747226839050681966839463919428531629782475127367001763589500187642982976435876178980851951693987841, algorithm=LGRR) == -1


def test_sqrt_in_field_extension():
    assert isinstance(finite_field_sqrt(ModP(-1, 2 ** 31 - 1)), FieldExtension)

def test_sqrt_in_fintie_field():
    assert not isinstance(finite_field_sqrt(ModP(-1, 2 ** 31 - 19)), FieldExtension)

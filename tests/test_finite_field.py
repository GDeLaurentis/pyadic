import pickle
import random
import pytest
import sympy
import numpy

from fractions import Fraction as Q

from pyadic import ModP, PAdic
from pyadic.finite_field import vec_ModP, extended_euclideal_algorithm, rationalise, MQRR, LGRR, \
    finite_field_sqrt, chained_chinese_remainder, vec_chained_FF_rationalize
from pyadic.field_extension import FieldExtension
from pyadic.primes import primes


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


def test_isscalar():
    assert numpy.isscalar(ModP(3, 7))


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
    b = Q(1, 3)
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


def test_multiplication_with_array():
    p = 10007
    a = random.randrange(0, 10000)
    array = [0, 1, 2, 3]
    assert numpy.all(ModP(a, p) * numpy.array(array) == numpy.array([ModP(a, p) * entry for entry in array]))


def test_division():
    p = 10007
    a, b = random.randrange(1, p - 1), random.randrange(1, p - 1)
    assert ModP(a, p) / b == ModP(Q(a, b), p)


def test_negation():
    p = 10007
    a = random.randrange(0, 10000)
    assert (-a) % p == -ModP(a, p)


def test_inverse():
    p = 10007
    a = ModP(random.randrange(1, p), p)
    assert a * a._inv() == 1


def test_str_eq_repr():
    p = 2 ** 31 - 19
    a = ModP(random.randrange(1, p), p)
    assert str(a) == repr(a)


def test_rstr():
    assert ModP('2 mod 4') == ModP('2 % 4')


def test_rstr_invalid():
    with pytest.raises(Exception):
        ModP('2 : 4')


def test_failed_inverse():
    a = ModP('2 % 4')
    with pytest.raises(ZeroDivisionError):
        a * a._inv()


def test_pow():
    p = 10007
    x = random.randrange(1, p)
    assert ModP(x, p) ** -13 == 1 / ModP(x ** 13, p)


def test_operation_not_understood_by_ModP():
    x = sympy.symbols("x")
    assert ModP("2 % 3") * x == x * ModP("2 % 3")


def test_trivial_abs():
    assert abs(ModP(1, 10007)) == abs(ModP(2, 10007))
    assert abs(ModP(-1, 10007)) > abs(ModP(0, 10007))


def test_hash():
    hash(ModP(12345, 10007))


def test_extended_euclideal_algorithm():
    for i in range(10):
        a, b = random.randint(1, 1000), random.randint(1, 1000)
        s, t, gcd = extended_euclideal_algorithm(a, b)
        assert a * s + b * t == gcd


def test_reconstruction_MQRR_2147483647():
    assert rationalise(298260199, 2147483647, algorithm=MQRR) == Q(-51071, 36)


def test_reconstruction_LGRR_2147483647():
    assert rationalise(298260199, 2147483647, algorithm=LGRR) == Q(11326, 42041)


def test_trivial_rationalisation():
    assert rationalise(1234) == 1234


def test_rationalisation_large_power():
    assert rationalise(ModP(Q(7, 13), 2147483647 ** 12), algorithm=LGRR) == Q(7, 13)
    assert rationalise(ModP(Q(7, 13), 2147483647 ** 12), algorithm=MQRR) == Q(7, 13)


def test_chained_chinese_remainder():
    assert rationalise(chained_chinese_remainder(ModP(Q(-7, 3), 2 ** 31 - 1), ModP(Q(-7, 3), 2 ** 31 - 19))) == Q(-7, 3)


def test_sqrt_in_field_extension():
    assert isinstance(finite_field_sqrt(ModP(-1, 2 ** 31 - 1)), FieldExtension)


def test_sqrt_in_fintie_field():
    assert not isinstance(finite_field_sqrt(ModP(-1, 2 ** 31 - 19)), FieldExtension)


def test_vec_chained_FF_rationalize():
    Qmatrix = numpy.array([[Q(16, 9973), Q(10007)], [Q(99991, 4), Q(100003)]])
    used_primes = primes[:4]
    FF_matrices = [vec_ModP(prime)(Qmatrix).astype(int) for prime in used_primes]
    assert numpy.all(vec_chained_FF_rationalize(FF_matrices, used_primes) == Qmatrix)

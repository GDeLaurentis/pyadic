import numpy
import pytest
import pickle
import hashlib

from fractions import Fraction as Q

from syngular import Field
from pyadic.gaussian_rationals import GaussianRational, rand_rat_frac


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_arithmetic_operations():

    a = GaussianRational(rand_rat_frac(), rand_rat_frac())
    b = GaussianRational(rand_rat_frac(), rand_rat_frac())

    assert numpy.isclose(complex(a + b), complex(a) + complex(b))
    assert numpy.isclose(complex(a - b), complex(a) - complex(b))
    assert numpy.isclose(complex(a * b), complex(a) * complex(b))
    assert numpy.isclose(complex(a / b), complex(a) / complex(b))
    assert numpy.isclose(complex(a ** 17), complex(a) ** 17)


def test_instantiation_from_strings():
    assert GaussianRational("2/3+1j") == GaussianRational(Q("2/3"), Q("1"))
    assert GaussianRational("-3+3/128j") == GaussianRational(Q("-3"), Q("3/128"))
    assert GaussianRational("51/6-7/8j") == GaussianRational(Q("51/6"), Q("-7/8"))
    assert GaussianRational("-4") == GaussianRational(Q("-4"), Q("0"))
    assert GaussianRational("3/4") == GaussianRational(Q("3/4"), Q("0"))
    assert GaussianRational("+2j") == GaussianRational(Q("0"), Q("2"))
    assert GaussianRational("-33/4j") == GaussianRational(Q("0"), Q("-33/4"))
    assert GaussianRational("0") == GaussianRational(Q("0"), Q("0"))


def test_instantiation_from_complex():
    assert GaussianRational(3.3 + 1.4j) == GaussianRational('33/10', '7/5')


@pytest.mark.parametrize(
    'original', [
        GaussianRational("0"),
        GaussianRational("3/4"),
        GaussianRational("-33/4j"),
        GaussianRational("-3+3/128j")
    ]
)
def test_serializable_and_hash_stable(original):
    dumped = pickle.dumps(original)
    loaded = pickle.loads(dumped)

    assert original == loaded

    hash1 = hashlib.sha256(pickle.dumps(original)).hexdigest()
    hash2 = hashlib.sha256(pickle.dumps(loaded)).hexdigest()

    assert hash1 == hash2


def test_arithmetics_with_mpc():
    a, b = Field("mpc", 0, 300).random(), Field("gaussian rational", 0, 0).random()
    assert a + b == b + a
    assert a * b == b * a

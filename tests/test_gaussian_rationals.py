import numpy

from fractions import Fraction as Q

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

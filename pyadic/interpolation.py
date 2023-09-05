import functools
import numpy
import random
import sympy

from collections.abc import Iterator

from pyadic import ModP


def Newton_polynomial_interpolation(f, prime, seed=0, verbose=False):
    """Univariate polynomial interpolation of f(t), samples taken modulo prime. See arXiv:1608.01902 section 3.2"""
    t_sequence_generator = FFSequenceGenerator(prime, seed)
    t = sympy.symbols('t')
    avals, tvals, subtracted = [], [], [f]

    def fsubtracted(i, t):
        return (subtracted[i - 1](t) - avals[i - 1]) / (t - tvals[i - 1])

    while avals[-2:] != [0, 0]:
        if verbose:
            print(f"\r@ {len(avals)}, {avals}", end="")
        tvals += [next(t_sequence_generator)]
        avals += [subtracted[-1](tvals[-1])]
        fsubtracted_partial = functools.partial(fsubtracted, len(avals))
        subtracted += [fsubtracted_partial]
    if verbose:
        print(f"\rFinished after {len(avals)} samples: {avals}.", end=" ")
    FFGF = sympy.GF(prime).frac_field(t)
    tpoly = FFGF(0)
    for aval, tval in zip(avals[:-2][::-1], tvals[:-2][::-1]):
        tpoly = int(aval) + (FFGF(t) - int(tval)) * tpoly
    return tpoly.as_expr()


def Thiele_rational_interpolation(f, prime, as_continued_fraction=False, seed=0, verbose=False):
    """Univariate rational interpolation of f(t), samples taken modulo prime. See arXiv:1608.01902 section 3.3"""
    t_sequence_generator = FFSequenceGenerator(prime, seed)
    t = sympy.symbols('t')
    avals, tvals, subtracted = [], [], [f]

    def fsubtracted(i, t):
        denom = (subtracted[i - 1](t) - avals[i - 1])
        if (numpy.isscalar(denom) and denom == 0) or (numpy.all(denom == 0)):
            raise ZeroDivisionError
        return (t - tvals[i - 1]) / denom

    while True:
        if verbose:
            print(f"\r@ {len(avals)}, {avals}", end="")
        tvals += [next(t_sequence_generator)]
        try:
            avals += [subtracted[-1](tvals[-1])]
        except ZeroDivisionError:
            tvals = tvals[:-1]
            break
        fsubtracted_partial = functools.partial(fsubtracted, len(avals))
        subtracted += [fsubtracted_partial]
    if verbose:
        print(f"\rFinished after {len(avals)} samples, {avals}.", end=" ")
    if as_continued_fraction:
        tpoly = "oo"
        for aval, tval in zip(avals[:][::-1], tvals[:][::-1]):
            tpoly = f"({int(aval)}+(t-{int(tval)})/(" + tpoly + "))"
        oo = sympy.oo  # noqa, used in eval
        return eval(tpoly)
    # get the single fraction (i.e. simplify the continued fraction)
    FFGF = sympy.GF(prime).frac_field(t)
    tpoly = FFGF(int(avals[-1]))
    for aval, tval in zip(avals[-2::-1], tvals[-2::-1]):
        tpoly = int(aval) + (FFGF(t) - int(tval)) / tpoly
    return tpoly.as_expr()


class FFSequenceGenerator(Iterator):
    """Random number generator decoupled from global state."""

    def __init__(self, prime, seed=None):
        self.prime = prime
        self.seed = seed
        random.seed(seed)
        self.local_state = random.getstate()

    def __next__(self):
        global_state = random.getstate()
        random.setstate(self.local_state)
        next_value = ModP(random.randint(1, self.prime - 1), self.prime)
        self.local_state = random.getstate()
        random.setstate(global_state)
        return next_value


def interpolation_t_sequence(length, prime, seed=0):
    """Returns the pseudo-random sequence of t-values (mod prime) used in Newton polynomial intrpolation and in Thiele rational interpolation."""
    FF_generator = FFSequenceGenerator(prime, seed)
    return [next(FF_generator) for i in range(length)]

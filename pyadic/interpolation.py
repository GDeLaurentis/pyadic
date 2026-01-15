import functools
import numpy
import random
import sympy
import inspect

from collections.abc import Iterator

from pyadic import ModP


def Newton_polynomial_interpolation(f, prime, seed=0, verbose=False, as_nested_sum=False, as_expr=True):
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
    if as_nested_sum:
        tpoly = "0"
        for aval, tval in zip(avals[:-2][::-1], tvals[:-2][::-1]):
            tpoly = f"({int(aval)}+(t-{int(tval)})*({tpoly}))"
        if verbose:
            print(f"\nNested sum representation: {tpoly}")
        return eval(tpoly)
    FFGF = sympy.GF(prime).frac_field(t)
    tpoly = FFGF(0)
    for aval, tval in zip(avals[:-2][::-1], tvals[:-2][::-1]):
        tpoly = int(aval) + (FFGF(t) - int(tval)) * tpoly
    if as_expr:
        return tpoly.as_expr()
    return tpoly


def Thiele_rational_interpolation(f, prime, as_continued_fraction=False, seed=0, verbose=False, as_expr=True):
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
    if len(avals) == 0:
        return sympy.oo
    FFGF = sympy.GF(prime).frac_field(t)
    tpoly = FFGF(int(avals[-1]))
    for aval, tval in zip(avals[-2::-1], tvals[-2::-1]):
        tpoly = int(aval) + (FFGF(t) - int(tval)) / tpoly
    if as_expr:
        return tpoly.as_expr()
    return tpoly


def multivariate_Newton_polynomial_interpolation(f, prime, seed=0, depth=0, verbose=False):
    """Recursive multivariate polynomial interpolation of f(ts), samples taken modulo prime."""

    # End of recursion condition: univariate function
    function_signature = inspect.signature(f)
    num_args = len(function_signature.parameters)
    if num_args == 1:
        tpoly = Newton_polynomial_interpolation(f, prime, seed=seed, verbose=verbose)
        tpoly = tpoly.subs({'t': 't1', })
        return tpoly

    t_sequence_generator = FFSequenceGenerator(prime, seed=seed + depth)  # do not use same sequence for all variables!!!
    avals, tvals, subtracted = [], [], [f]

    dynamic_fsubtracted_string = f"""def fsubtracted(i, {''.join([f't{i}, ' for i in range(1, num_args + 1)])}):
    return (subtracted[i - 1]({''.join([f't{i}, ' for i in range(1, num_args + 1)])}) - \
    ModP(int(avals[i - 1]({''.join([f'int(t{i}), ' for i in range(2, num_args + 1)])})), {prime})) / (t1 - tvals[i - 1])
    """
    local_dict = {'ModP': ModP, 'subtracted': subtracted, 'avals': avals, 'tvals': tvals}  # Include any necessary objects in the local dictionary
    exec(dynamic_fsubtracted_string, local_dict)
    fsubtracted = local_dict['fsubtracted']

    MAX_SAMPLES = 1000  # sets iteration limit
    for counter in range(MAX_SAMPLES):
        if avals[-2:] == [0, 0]:
            break
        if verbose:
            if depth != 0:
                print()
            print(f"@ depth: {depth} - samples: {len(avals)}, {(avals, tvals)}", end="\n")
        tvals += [next(t_sequence_generator)]
        local_dict = {'ModP': ModP, 'subtracted': subtracted}  # Include any necessary objects in the local dictionary
        function_string = f"""def rest_function({''.join([f't{i}, ' for i in range(2, num_args + 1)])}):
        return subtracted[-1](ModP('{tvals[-1]}'), {''.join([f't{i}, ' for i in range(2, num_args + 1)])})"""
        exec(function_string, local_dict)
        rest_function = local_dict['rest_function']
        avals += [multivariate_Newton_polynomial_interpolation(rest_function, prime, seed=seed, depth=depth + 1, verbose=verbose)]
        avals[-1] = avals[-1].subs({f't{i}': f't{i+1}' for i in range(1, num_args)}, simultaneous=True)
        avals[-1] = sympy.poly(avals[-1], sympy.symbols([f't{i+1}' for i in range(1, num_args)]), modulus=prime)
        fsubtracted_partial = functools.partial(fsubtracted, len(avals))
        subtracted += [fsubtracted_partial]

    if verbose:
        print(f"\nFinished after {len(avals)} samples: {avals}.", end="\n")

    FFGF = sympy.GF(prime).frac_field(*sympy.symbols([f't{i}' for i in range(1, num_args + 1)]))
    tpoly = FFGF(0)
    for aval, tval in zip(avals[:-2][::-1], tvals[:-2][::-1]):
        tpoly = FFGF(aval.as_expr()) + (FFGF(sympy.symbols('t1')) - int(tval)) * tpoly
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

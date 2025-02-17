# -*- coding: utf-8 -*-

import fractions
import numbers
import numpy
import operator
import itertools

from functools import reduce


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class FieldExtension(object):

    """Field extension by multiple sqrts."""

    def get_sorting_key(self):
        from .finite_field import ModP
        from .padic import PAdic
        return int if self.base_type in (ModP, PAdic) else abs

    @staticmethod
    def field_extension_basis(squares, key=abs):
        squares = [square for square in squares if square != 1]
        squares = sorted(set(squares), key=key)
        basis = {()}
        for r in range(1, len(squares) + 1):
            for subset in itertools.combinations(squares, r):
                basis.add(subset)
        return sorted(basis, key=lambda x: (len(x), tuple(key(elem) for elem in x)))

    def __new__(cls, squares, dict_=None):
        fe_obj = super(FieldExtension, cls).__new__(cls)
        fe_obj.__init__(squares, dict_)
        if all(val == 0 or key == () for key, val in fe_obj.dict_.items()):
            return fe_obj.dict_[()]
        return fe_obj

    def __init__(self, squares, dict_=None):
        try:
            self.base_type = type(squares[-1])
        except (TypeError, KeyError, IndexError):
            self.base_type = type(squares)
            squares = (squares, )
        sorting_key = self.get_sorting_key()
        self.squares = tuple(sorted(set(squares), key=sorting_key))
        full_basis = self.field_extension_basis(self.squares, key=sorting_key)

        if dict_ is None:
            zero = self.squares[0] * 0
            one = zero + 1
            self.dict_ = {basis_elem: zero for basis_elem in full_basis}
            # Set the coefficients of single-entry square roots to 1
            for basis_elem in full_basis:
                if len(basis_elem) == 1 and basis_elem[0] != 1:
                    self.dict_[basis_elem] = one
        else:
            self.dict_ = {key: value for key, value in dict_.items() if value != 0}

        # Determine which squares are actually used
        used_squares = set()
        for basis_elem in self.dict_.keys():
            used_squares.update(basis_elem)

        # Update squares and recompute basis
        self.squares = tuple(sorted(used_squares, key=sorting_key))
        self.basis = self.field_extension_basis(self.squares, key=sorting_key)

        # Project the dict onto the new basis
        self.dict_ = {key: self.dict_.get(key, 0) for key in self.basis}

        # print(self.squares, self.basis, self.dict_)

    def __repr__(self):
        terms = []
        for basis_entry in self.basis:
            if basis_entry == ():
                if self.dict_.get(basis_entry, 0) != 0:
                    terms.append(f"({self.dict_.get(basis_entry, 0)})")
            else:
                if self.dict_.get(basis_entry, 0) != 0:
                    sqrt_terms = "·".join(f"√({s})" for s in basis_entry if s != 1)
                    terms.append(f"({self.dict_.get(basis_entry, 0)})·{sqrt_terms}")
        return " + ".join(terms) if terms else "0"

    def __getnewargs__(self):
        return self.__getstate__()

    def __getstate__(self):
        return (self.squares, self.dict_)

    def __setstate__(self, state):
        self.__init__(*state)

    def __eq__(self, other):
        if not isinstance(other, FieldExtension):
            return NotImplemented
        return self.squares == other.squares and self.dict_ == other.dict_

    def __add__(self, other):
        if isinstance(other, (self.base_type, int, complex, numpy.integer, fractions.Fraction)):
            new_dict = self.dict_.copy()
            new_dict[()] = new_dict.get((), 0) + other
            return FieldExtension(self.squares, new_dict)
        elif isinstance(other, FieldExtension):
            new_squares = sorted(set(self.squares).union(other.squares), key=self.get_sorting_key())
            # Perform the addition by combining the coefficients for each basis element
            new_dict = self.dict_.copy()
            for key, value in other.dict_.items():
                new_dict[key] = new_dict.get(key, 0) + value
            # Return a new FieldExtension with the merged squares
            return FieldExtension(new_squares, new_dict)
        else:
            # for debugging only - otherwise allow to use __radd__ from other
            # raise TypeError(f"Unsupported {other} of type {type(other)} for addition with {self} of type {type(self)}.")
            return NotImplemented

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        if isinstance(other, (self.base_type, int, complex, numpy.integer, fractions.Fraction)):
            new_dict = self.dict_.copy()
            for key in new_dict:
                new_dict[key] *= other
            return FieldExtension(self.squares, new_dict)
        elif isinstance(other, FieldExtension):
            new_squares = sorted(set(self.squares).union(other.squares), key=self.get_sorting_key())
            new_dict = {}
            for self_key, self_value in self.dict_.items():
                for other_key, other_value in other.dict_.items():
                    common_keys = set(self_key).intersection(set(other_key))
                    combined_value = self_value * other_value * reduce(operator.mul, common_keys, 1)
                    combined_key = tuple(sorted(set(self_key).union(other_key) - common_keys, key=self.get_sorting_key()))
                    if combined_key in new_dict:
                        new_dict[combined_key] += combined_value
                    else:
                        new_dict[combined_key] = combined_value
            return FieldExtension(new_squares, new_dict)
        else:
            # for debugging only - otherwise allow to use __rmull__ from other
            # raise TypeError(f"Unsupported {other} of type {type(other)} for multiplication with {self} of type {type(self)}.")
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, (self.base_type, int, complex, numpy.integer, fractions.Fraction)):
            return FieldExtension(self.squares, {key: value / other for key, value in self.dict_.items()})
        elif isinstance(other, FieldExtension):
            return self * other._inverse()
        else:
            # for debugging only - otherwise allow to use __rtruediv__ from other
            # raise TypeError(f"Unsupported {other} of type {type(other)} for division with {self} of type {type(self)}.")
            return NotImplemented

    def __div__(self, other):
        return self.__truediv__(other)

    def __rtruediv__(self, other):
        return other * self._inverse()

    def __rdiv__(self, other):
        return self.__rtruediv__(other)

    def __sub__(self, other):
        return self + (- other)

    def __rsub__(self, other):
        return - (self - other)

    def __neg__(self):
        """Unary '-' operation"""
        return -1 * self

    def __pos__(self):
        """Unary '+' operation"""
        return self

    def naive_conjugate(self, sqrt):
        """Conjugates on the i-th root."""
        return FieldExtension(self.squares, {key: (-1) ** (-1 if sqrt in key else 0) * value for key, value in self.dict_.items()})

    @property
    def conjugate(self):
        conjugate = 1
        result = self
        for i, sqrt in enumerate(self.squares):
            new_conjugate = result.naive_conjugate(sqrt)
            result *= new_conjugate
            conjugate *= new_conjugate
            if not isinstance(result, FieldExtension):
                break
        else:
            raise Exception(f"Failed to compute conjugate of {self}.")
        assert isinstance(self * conjugate, self.base_type)
        return conjugate

    def _inverse(self):
        conjugate = self.conjugate
        return conjugate / (self * conjugate)

    def __pow__(self, n):
        assert (isinstance(n, int) or n.is_integer())
        if n < 0:
            return 1 / self ** -n
        elif n == 0:
            return self * 0 + 1
        elif n % 2 == 0:
            root_2_res = self ** (n / 2)
            return root_2_res * root_2_res
        else:
            return self * (self ** (n - 1))


numbers.Number.register(FieldExtension)

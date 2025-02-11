# -*- coding: utf-8 -*-

import fractions
import numbers
import numpy


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class FieldExtension(object):

    """Field extension by sqrt."""

    def __new__(cls, square, tuple_=None):
        fe_obj = super(FieldExtension, cls).__new__(cls)
        fe_obj.__init__(square, tuple_)
        if tuple_ is not None and tuple_[1] == square * 0:
            return fe_obj.tuple[0]
        else:
            return fe_obj

    def __init__(self, square, tuple_=None):
        if tuple_ is None:
            zero = square * 0
            one = zero + 1
            self.tuple = (zero, one)
        else:
            self.tuple = tuple_
        self.square = square
        self.base_type = type(square)

    def __repr__(self):
        return "(" + str(self.tuple[0]) + ") + (" + str(self.tuple[1]) + ")·√(" + str(self.square) + ")"

    def __getnewargs__(self):
        return self.__getstate__()

    def __getstate__(self):
        return (self.square, self.tuple)

    def __setstate__(self, state):
        self.__init__(*state)

    def __eq__(self, other):
        return self.square == other.square and self.tuple[0] == other.tuple[0] and self.tuple[1] == other.tuple[1]

    def __add__(self, other):
        if isinstance(other, (self.base_type, int, complex, numpy.integer, fractions.Fraction)):
            return FieldExtension(self.square, (self.tuple[0] + other, self.tuple[1]))
        elif isinstance(other, FieldExtension):
            assert self.square == other.square
            return FieldExtension(self.square, (self.tuple[0] + other.tuple[0], self.tuple[1] + other.tuple[1]))
        else:
            # for debugging only - otherwise allow to use __radd__ from other
            # raise TypeError(f"Unsupported {other} of type {type(other)} for addition with {self} of type {type(self)}.")
            return NotImplemented

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        if isinstance(other, (self.base_type, int, complex, numpy.integer, fractions.Fraction)):
            return FieldExtension(self.square, (self.tuple[0] * other, self.tuple[1] * other))
        elif isinstance(other, FieldExtension):
            assert self.square == other.square
            return FieldExtension(self.square, (self.tuple[0] * other.tuple[0] + self.square * self.tuple[1] * other.tuple[1],
                                                self.tuple[0] * other.tuple[1] + self.tuple[1] * other.tuple[0]))
        else:
            # for debugging only - otherwise allow to use __rmull__ from other
            # raise TypeError(f"Unsupported {other} of type {type(other)} for multiplication with {self} of type {type(self)}.")
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, (self.base_type, int, complex, numpy.integer, fractions.Fraction)):
            return FieldExtension(self.square, (self.tuple[0] / other, self.tuple[1] / other))
        elif isinstance(other, FieldExtension):
            assert self.square == other.square
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

    def _inverse(self):
        return FieldExtension(self.square, (self.tuple[0], -self.tuple[1])) / (self.tuple[0] ** 2 - self.square * self.tuple[1] ** 2)

    def __pow__(self, n):
        assert (isinstance(n, int) or n.is_integer())
        if n < 0:
            return 1 / self ** -n
        elif n == 0:
            return FieldExtension(self.square, (self.square * 0 + 1, self.square * 0))
        elif n % 2 == 0:
            root_2_res = self ** (n / 2)
            return root_2_res * root_2_res
        else:
            return self * (self ** (n - 1))


numbers.Number.register(FieldExtension)

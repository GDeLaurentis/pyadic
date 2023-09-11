.. seampy documentation master file, created by
   sphinx-quickstart on Tue Sep 24 12:11:24 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyAdic
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

The pyadic library is Python 3 package that provides number types for finite fields (ModP) and $p$-adic numbers (PAdic). 
The goal is mimic the flexible behavior of built-in types, such as int, float and complex. 
Thus, one can mix-and-match the different number types, as long as the operations are consistent. 
In particular, ModP and PAdic are compatible with fractions.Fraction.


Installation
------------

Installation is easy with pip::

  pip install pyadic

alternatively the package can be cloned from github at https://github.com/GDeLaurentis/pyadic.


Quick start
-----------

In [1]: from pyadic import PAdic, ModP
In [2]: from fractions import Fraction as Q

# 7/13 as a 12-digit 2147483647-adic number
In [3]: PAdic(Q(7, 13), 2147483647, 12)  
Out [3]: 1817101548 + 825955248*2147483647 + 1156337348*2147483647^2 + 330382099*2147483647^3 + 1321528398*2147483647^4 + 991146298*2147483647^5 + 1817101547*2147483647^6 + 825955248*2147483647^7 + 1156337348*2147483647^8 + 330382099*2147483647^9 + 1321528398*2147483647^10 + 991146298*2147483647^11 + O(2147483647^12)

# 7/13 in F_2147483647
In [4]: ModP(Q(7, 13), 2147483647)
Out [4]: 1817101548 % 2147483647

# Mapping back to rational numbers
In [5]: from pyadic.finite_field import rationalise
In [6]: rationalise(ModP(Q(7, 13), 2147483647))
Out [6]: Fraction(7, 13)
In [7]: rationalise(PAdic(Q(7, 13), 2147483647, 12))
Out [7]: Fraction(7, 13)

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. Hidden TOCs

.. toctree::
   :caption: Modules Documentation
   :maxdepth: 2

   modules

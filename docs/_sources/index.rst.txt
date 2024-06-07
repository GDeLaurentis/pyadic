pyAdic
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

The pyadic library is Python 3 package that provides number types for finite fields (ModP) and *p*-adic numbers (PAdic). 
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

.. code-block:: python
   :caption: Instantiating number pyadic types
   :linenos:
   
   from pyadic import PAdic, ModP
   from fractions import Fraction as Q

   # 7/13 as a 12-digit 2147483647-adic number
   PAdic(Q(7, 13), 2147483647, 12)  
   > 1817101548 + 825955248*2147483647 + 1156337348*2147483647^2 + 330382099*2147483647^3 + 1321528398*2147483647^4 + 991146298*2147483647^5 + 1817101547*2147483647^6 + 825955248*2147483647^7 + 1156337348*2147483647^8 + 330382099*2147483647^9 + 1321528398*2147483647^10 + 991146298*2147483647^11 + O(2147483647^12)

   # 7/13 in F_2147483647
   ModP(Q(7, 13), 2147483647)
   > 1817101548 % 2147483647


.. code-block:: python
   :caption: Mapping back to rational numbers
   :linenos:
   :lineno-start: 12

   from pyadic.finite_field import rationalise

   rationalise(ModP(Q(7, 13), 2147483647))
   > Fraction(7, 13)
   rationalise(PAdic(Q(7, 13), 2147483647, 12))
   > Fraction(7, 13)


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

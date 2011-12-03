.. Description of pymc21 module
.. _pymc21-page:

=========================================
Bringing Nonzero Elements to the Diagonal
=========================================

The :mod:`pymc21` Module
========================

.. automodule:: pyorder.pymc21.pymc21

Functions available
-------------------

.. autofunction:: nonzerodiag


Examples
========

Basic Usage
-----------

This first example calls the Fortran subroutine directly with the matrix data
in compressed sparse column format.

.. warning::

   Keep in mind that when calling the Fortran subroutines directly, all indices
   in the matrix data must be 1-based, i.e., row indices range from 1 through
   nrow and column indices range from 1 through ncol.

.. literalinclude:: ../../examples/demo_mc21.py
   :linenos:

Python Interface
----------------

In this second example, Python usage is intended. Matrix indices must therefore
be 0-based. The permutation vector returned by :func:`nonzerodiag` is now also
0-based.

.. literalinclude:: ../../examples/demo_hb_mc21.py
   :linenos:

Reorganizing a square matrix so it has a nonzero diagonal is useful for
instance when factorizing a submatrix of a basis in the Simplex method for
linear programming.

The effect of :func:`nonzerodiag` is illustrated on the following example,
`shl_0 <http://www.cise.ufl.edu/research/sparse/matrices/HB/shl_0.html>`_ from
the `University of Florida Sparse Matrix Collection
<http://www.cise.ufl.edu/research/sparse/matrices>`_. I used the matrix in
Rutherford-Boeing format.

.. image:: shl0_mc21.png

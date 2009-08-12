.. Introduction to PyOrder

============
Introduction
============

PyOrder is a Python package that provides access to several means to reorder
sparse matrices. Typically, sparse matrices are reordered prior to
factorization so as to preserve sparsity of the factors, ensure that there are
as many nonzero elements as possible on the diagonal, or decrease the envelope,
wavefront, profile or semi-bandwidth.

.. note::

   Orderings designed to promote sparsity of the factors of a matrix are
   different in nature from the other types of orderings mentioned in the
   previous paragraph and are not considered in PyOrder. If you would like to
   see such orderings (e.g., AMD, CAMD, COLAMD, etc.) included in future
   releases of PyOrder, please let me know.

A Few Concepts about Symmetric Sparse Matrices
==============================================

Although not all ordering functions in this package assume symmetry, most
concepts are defined for symmetric matrices only.

If :math:`A` is a square :math:`n \times n` symmetric sparse matrix,
let :math:`a_{ij}` be the element at the intersection of the :math:`i`-th row
and the :math:`j`-th column. The following concepts are important for
reordering purposes.

The :math:`i`-th *wavefront* of :math:`A` is the number of nonzero rows in the
submatrix :math:`A(i:n,1:i)`, :math:`i = 1, \ldots, n`. It is
denoted :math:`f_i`.

The maximum and mean-square wavefronts are

.. math::

   \max_{i=1,\ldots,n} f_i
   \quad \text{and} \quad
   \frac{1}{n} \sum_{i=1}^{n} f_i^2.

If :math:`m_i` is the column index of the first nonzero element on
the :math:`i`-th row of :math:`A`, the *envelope* of :math:`A` is the set of
nonzero elements lying between :math:`(i,m_i)` and the diagonal, exclusive of
the diagonal, i.e.

.. math::

   \text{Env}(A) = \left\{
   (i,j) \mid 1 \leq i \leq n, \ m_i \leq j < i
   \right\}.

The *profile* :math:`P(A)` of the matrix is the number of elements
in :math:`\text{Env}(A)` plus the number of nonzero elements on the diagonal,
i.e.,

.. math::

   P = \sum_{i=1}^{n} (i - m_i + 1).

It is possible to show that those quantities are linked via the relationship

.. math::

   P = \sum_{i=1}^{n} f_i.


Availability
============

PyOrder is essentially a set of interfaces to subroutines from the HSL
`<http://www.cse.scitech.ac.uk/nag/hsl>`_ (formerly the Harwell Subroutine
Library). Subroutines from the HSL may be obtained free of charge under certain
conditions. See `<http://hsl.rl.ac.uk/hsl2007/hsl20074researchers.html>`_ to
decide whether the terms apply to you. The HSL subroutines relevant to PyOrder
are not packaged together with the Python interfaces and should be obtained
separately.

.. note::

   At this time, only double precision real data is supported. Please let me
   know if you would like support for single precision and/or complex data.

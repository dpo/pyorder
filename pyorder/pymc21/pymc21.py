"""
A Python interface to the HSL subroutine MC21AD.

**References**

.. [Duff81a] I. S. Duff, *On Algorithms for Obtaining a Maximum Transversal*,
             ACM Trans. Math. Software, **7**, pp. 315-330, 1981.

.. [Duff81b] I. S. Duff, *Algorithm 575: Permutations for a zero-free
             diagonal*, ACM Trans. Math. Software, **7**, pp. 387-390, 1981.

.. moduleauthor: dominique.orban@gerad.ca
"""

import numpy as np
from mc21module import mc21ad

__docformat__ = 'restructuredtext'

def nonzerodiag(nrow, colind, rowptr):
    """
    Given the sparsity pattern of a square sparse matrix in compressed row
    (csr) format, attempt to find a *row* permutation so the row-permuted
    matrix has a nonzero diagonal, if this is possible. This function assumes
    that the matrix indexing is 0-based. Note that in general, this function
    does not preserve symmetry. The method used is a depth-first search with
    lookahead described in [Duff81a]_ and [Duff81b]_.

    :parameters:

        :nrow: The number of rows of the input matrix.

        :colind: An integer array (or list) of length nnz giving the column
                 indices of the nonzero elements in each row.

        :rowptr: An integer array (or list) of length nrow+1 giving the indices
                 of the first element of each row in colind.

    :returns:

        :perm: An integer array of length nrow giving the variable permutation.
               If irow and jcol are two integer arrays describing the pattern of
               the input matrix in triple format, perm[irow] and jcol
               describe the permuted matrix.

        :nzdiag: The number of nonzeros on the diagonal of the permuted matrix.
    """
    lenrows = rowptr[1:] - rowptr[:-1]
    perm, nzdiag = mc21ad(colind+1, rowptr[:-1]+1, lenrows)
    return (perm-1, nzdiag)

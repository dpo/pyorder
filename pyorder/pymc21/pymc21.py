"""
A Python interface to the HSL subroutine MC21AD.
dominique.orban@gerad.ca
"""

__docformat__ = 'restructuredtext'

import numpy as np
from mc21module import mc21ad

def nonzerodiag(nrow, colind, rowptr):
    """
    Given the sparsity pattern of a square sparse matrix in compressed row (csr)
    format, attempt to find a *row* permutation so the row-permuted matrix has a nonzero
    diagonal, if this is possible. This function assumes that the matrix
    indexing is 0-based.

    :parameters:

        :nrow: The number of rows of the input matrix.

        :colind: An integer array (or list) of length nnz giving the column
                 indices of the nonzero elements in each row.

        :rowptr: An integer array (or list) of length nrow+1 giving the indices
                 of the first element of each row in colind.

    :return values:

        :perm: An integer array of length nrow giving the variable permutation.
               If irow and jcol are two integer arrays describing the pattern of
               the input matrix in triple format, perm[irow] and jcol
               describe the permuted matrix.

        :nzdiag: The number of nonzeros on the diagonal of the permuted matrix.
    """
    lenrows = rowptr[1:] - rowptr[:-1]
    perm, nzdiag = mc21ad(colind+1, rowptr[:-1]+1, lenrows)
    return (perm-1, nzdiag)

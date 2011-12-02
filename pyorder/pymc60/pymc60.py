"""
A Python interface to the HSL subroutine MC60AD.

The functions in this module compute a symmetric permutation of a sparse
symmetric matrix so as to reduce its profile, wavefront, or bandwidth via
Sloan's method [SLO]_ or the reverse Cuthill-McKee method [CM]_.

References

.. [CM] E. Cuthill and J. McKee, *Reducing the bandwidth of sparse symmetric
        matrices* In Proc. 24th Nat. Conf. ACM, pages 157-172, 1969.

.. [RS] J. K. Reid and J. A. Scott, *Ordering symmetric sparse matrices for
        small profile and wavefront*, International Journal for Numerical
        Methods in Engineering, **45** (12), pp. 1737--1755, 1999.

.. [SLO] S. W. Sloan, *An algorithm for profile and wavefront reduction of
         sparse matrices*, International Journal of Numerical Methods in
         Engineering, **23**, pp. 239--251, 1986.

.. moduleauthor:: dominique.orban@gerad.ca
"""

__docformat__ = 'restructuredtext'

import numpy as np
import mc60module as mc60

def sloan(n, rowind, colptr, icntl=[0,6], weight=[2,1]):
    """
    Apply Sloan's algorithm to reduce the profile and wavefront of a sparse
    symmetric matrix. Either the lower or the upper triangle of the input
    matrix should be given in compressed sparse column (csc) or compressed
    sparse row (csr) format. This includes the diagonal of the matrix. A set of
    weights can be supplied to define the priority function in Sloan's method.

    :parameters:

        :n: The order of the input matrix.

        :rowind: An integer array (or list) of length nnz giving the row
                 indices of the nonzero elements in each column.

        :colptr: An integer array (or list) of length n+1 giving the indices of
                 the first element of each column in rowind.

    Note that since either triangle can be given in either csc or csr format,
    the words 'row' and 'column' may be swapped in the description above.
    The indexing in rowind and colptr should be zero-based.

    :keywords:

        :icntl: An integer array (or list) of length two of control parameters
                used during the first phase, where the input data is checked.
                The method terminates if duplicates of out-of-range indices are
                discovered (icntl[0]=0) or ignores them (icntl[0]=1). No
                diagnostic messages will be output if icntl[1]=0. If icntl[1]
                is > 0, it gives the unit number (in the Fortran sense) where
                diagonostic messages are output.

        :weight: An integer array (or list) of length two giving the weights
                 in Sloan's priority function. Reid and Scott (1999) recommend
                 to apply the method twice, with either [2,1] and [16,1], or
                 with [1,2] and [16,1], and to retain the best result.

    :returns:

        :perm: An integer array of length n giving the variable permutation. If
               irow and jcol are two integer arrays describing the pattern of
               the input matrix in triple format, perm[irow] and perm[jrow]
               describe the permuted matrix.

        :rinfo: A real array of length 4 giving statistics on the permuted
                matrix.
                rinfo[0] = profile
                rinfo[1] = maximum wavefront
                rinfo[2] = semi-bandwidth
                rinfo[3] = root-mean-square wavefront.
    """
    return reorder_matrix(n, rowind, colptr, icntl, jcntl=[0,0], weight=[2,1])

def rcmk(n, rowind, colptr, icntl=[0,6]):
    """
    Apply the reverse Cuthill-McKee algorithm to reduce the bandwidth of
    a sparse symmetric matrix. Either the lower or the upper triangle of the
    input matrix should be given in compressed sparse column (csc) or
    compressed sparse row (csr) format. This includes the diagonal of the
    matrix.

    :parameters:

        :n: The order of the input matrix.

        :rowind: An integer array (or list) of length nnz giving the row
                 indices of the nonzero elements in each column.

        :colptr: An integer array (or list) of length n+1 giving the indices of
                 the first element of each column in rowind.

    Note that since either triangle can be given in either csc or csr format,
    the words 'row' and 'column' may be swapped in the description above.
    The indexing in rowind and colptr should be zero-based.

    :keywords:

        :icntl: An integer array (or list) of length two of control parameters
                used during the first phase, where the input data is checked.
                The method terminates if duplicates of out-of-range indices are
                discovered (icntl[0]=0) or ignores them (icntl[0]=1). No
                diagnostic messages will be output if icntl[1]=0. If icntl[1]
                is > 0, it gives the unit number (in the Fortran sense) where
                diagonostic messages are output.

    :returns:

        :perm: An integer array of length n giving the variable permutation. If
               irow and jcol are two integer arrays describing the pattern of
               the input matrix in triple format, perm[irow] and perm[jcol]
               describe the permuted matrix.

        :rinfo: A real array of length 4 giving statistics on the permuted
                matrix.
                rinfo[0] = profile
                rinfo[1] = maximum wavefront
                rinfo[2] = semi-bandwidth
                rinfo[3] = root-mean-square wavefront.
    """
    return reorder_matrix(n, rowind, colptr, icntl, jcntl=[1,0])

def reorder_matrix(n, rowind, colptr, icntl=[0,6], jcntl=[0,0], weight=[2,1]):
    """
    Helper function called by `sloan` and `rcm` performing the bulk of the
    work when applying Sloan's method or the reverse Cuthill-McKee algorithm to
    a symmetric sparse matrix.
    """
    # Make room for pattern of the whole matrix and adjust to be 1-based
    icptr = colptr.copy() + 1
    irn = np.empty(2*(icptr[-1] - 1), dtype=np.int32)
    irn[:icptr[-1]-1] = rowind.copy() + 1

    # Check data
    info = mc60.mc60ad(irn, icptr, icntl)

    # Compute supervariables
    nsup, svar, vars = mc60.mc60bd(irn, icptr)

    # Permute reduced matrix
    permsv = np.empty(nsup, dtype=np.int32)
    pair = np.empty((2,nsup/2), dtype=np.int32)
    info = mc60.mc60cd(n, irn, icptr[:nsup+1], vars[:nsup],
                       jcntl, permsv, weight, pair)

    # Compute profile, maximum wavefront, semi-bandwidth and root-mean-square
    # wavefront of the permuted matrix
    rinfo = mc60.mc60fd(n, irn, icptr[:nsup+1], vars[:nsup], permsv)

    # Obtain variable permutation from supervariable permutation
    perm, possv = mc60.mc60dd(svar, vars[:nsup], permsv)

    # Adjust permutation to make it 0-based
    perm -= 1

    # There are other things to return!
    return (perm, rinfo)

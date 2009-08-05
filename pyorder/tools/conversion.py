import numpy as np

def coord2csc(ncol, irow, jcol, val=None):
    """
    Convert a sparse matrix in triple format into compressed sparse column
    (csc) format. To convert to compressed sparse row (csr) format, substitute
    the number of rows for ncol and swap irow and jcol.

    :arguments:

        :ncol: Number of columns in the input matrix.
        :irow: Integer Numpy array of length nnz containing the row indices of
               the nonzeros.
        :jcol: Integer Numpy array of length nnz containing the column indices
               of the nonzeros.

    :keywords:

        :val: Numpy array of length nnz containing the numerical values.

    :returns:

        :rowind: Integer Numpy array of length nnz containing the row indices
                 of the nonzeros in each column.
        :colptr: Integer Numpy array of length ncol+1 containing the positions
                 in rowind of the start of each column.
        :values: Numpy array of length nnz containing the numerical values.
    """
    nnz = len(irow)
    if len(jcol) != nnz:
        raise ValueError, 'irow and jcol must have the same length'

    rowind = np.zeros(nnz, dtype=np.int)
    colptr = np.zeros(ncol+1, dtype=np.int)
    if val is not None: values = np.zeros(nnz)

    # Store number of nonzeros in each column
    for k in range(nnz):
        j = jcol[k]
        colptr[j] += 1
    colptr[ncol] = nnz

    # Go backwards to find the row index of the first nonzero in each column
    for j in range(ncol-1, -1, -1):
        colptr[j] = colptr[j+1] - colptr[j]

    # Copy entries
    for k in range(nnz):
        j = jcol[k]
        elem = colptr[j]
        rowind[elem] = irow[k]
        if val is not None: values[elem] = val[k]
        colptr[j] = elem + 1

    # Restore colptr
    for j in range(ncol-1, 0, -1):
        colptr[j] = colptr[j-1]
    colptr[0] = 0

    if val is None:
        return (rowind, colptr)
    else:
        return (rowind, colptr, values)


def csc2coord(rowind, colptr):
    """
    Convert a sparse matrix in compressed sparse column (csc) format into
    triple format. To convert to compressed sparse row (csr) format, swap irow
    and jcol.

    :arguments:

        :rowind: Integer Numpy array of length nnz containing the row indices
                 of the nonzeros in each column.
        :colptr: Integer Numpy array of length ncol+1 containing the positions
                 in rowind of the start of each column.

    :returns:

        :irow: Integer Numpy array of length nnz containing the row indices of
               the nonzeros.
        :jcol: Integer Numpy array of length nnz containing the column indices
               of the nonzeros.
    """
    nnz = len(rowind)
    ncol = len(colptr)-1
    irow = np.empty(nnz, np.int)
    jcol = np.empty(nnz, np.int)
    for i in range(ncol):
        for j in range(colptr[i], colptr[i+1]):
            jcol[j] = i
            irow[j] = rowind[j]
    return (irow,jcol)

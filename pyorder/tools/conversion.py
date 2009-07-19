import numpy as np

def coord2csc(ncol, irow, jcol, val=None):
    """
    Convert a sparse matrix in triple format into compressed sparse column
    (csc) format. To convert to compressed sparse row (csr) format, substitute
    the number of rows for ncol and swap irow and jcol.
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

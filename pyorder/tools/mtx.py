"""
Read a sparse matrix in Matrix Market format

.. moduleauthor: D. Orban <dominique.orban@gerad.ca>
"""

import numpy as np
from string import atoi, atof

class MatrixMarketMatrix:
    """
    A MatrixMarketMatrix object represents a sparse matrix read from a file.
    This file must describe a sparse matrix in the MatrixMarket file format.

    See http://math.nist.gov/MatrixMarket for more information.

    Example: mat = MatrixMarketMatrix('1138bus.mtx')
    """

    def __init__(self, fname, **kwargs):
        self.comments= ''
        self.dtype = None
        self.irow = None
        self.jcol = None
        self.values = None
        self.nrow = self.ncol = self.nnz = 0
        self.shape = (0,0)
        self.symmetric = self.Hermitian = self.skewsym = False

        fp = open(fname)
        pos = self.readHeader(fp)
        nrecs = self.readData(fp,pos)
        fp.close()

        if nrecs != self.nnz:
            raise ValueError, 'Read %d records. Expected %d.' % (nrecs,self.nnz)

    def readHeader(self, fp):
        fp.seek(0)
        hdr = fp.readline().split()
        if hdr[1] != 'matrix' or hdr[2] != 'coordinate':
            raise TypeError, 'Type not supported: %s' % hdr[1:3]

        # Determine entries type
        dtypes = {'real': np.float,
                  'complex': np.complex,
                  'integer': np.int,
                  'pattern': None}

        self.dtype = dtypes[hdr[3]]

        # Determine symmetry
        if hdr[4] == 'symmetric':
            self.symmetric = True
        elif hdr[4] == 'Hermitian':
            self.Hermitian = True
        elif hdr[4] == 'skew-symmetric':
            self.skewsym = True

        # Read comments
        line = fp.readline()
        while line[0] == '%':
            self.comments += line
            line = fp.readline()

        # Return current position
        return fp.tell() - len(line)

    def readData(self, fp, pos):
        fp.seek(pos)
        size = fp.readline().split()
        self.nrow = atoi(size[0])
        self.ncol = atoi(size[1])
        self.shape = (self.nrow, self.ncol)
        self.nnz  = atoi(size[2])
        self.irow = np.empty(self.nnz, dtype=np.int)
        self.jcol = np.empty(self.nnz, dtype=np.int)

        if self.dtype is not None:
            self.values = np.empty(self.nnz, dtype=self.dtype)

        # Read in data
        k = 0
        for line in fp.readlines():
            line = line.split()
            self.irow[k] = atoi(line[0])-1
            self.jcol[k] = atoi(line[1])-1
            if self.dtype == np.int:
                self.values[k] = atoi(line[2])
            elif self.dtype == np.float:
                self.values[k] = atof(line[2])
            elif self.dtype == np.complex:
                self.values[k] = complex(atof(line[2]),atof(line[3]))
            k += 1

        return k

    def find(self):
        """
        Return the sparse matrix in triple format (val,irow,jcol). If the
        matrix data type is `None`, i.e., only the matrix sparsity pattern
        is available, this method returns (irow,jcol).
        """
        if self.dtype is not None:
            return (self.values,self.irow,self.jcol)
        return (self.irow,self.jcol)

if __name__ == '__main__':

    import sys
    fname = sys.argv[1]
    A = MatrixMarketMatrix(fname)
    print 'type: ', A.dtype
    print 'symmetric, Hermitian, skew: ', A.symmetric, A.Hermitian, A.skewsym
    print 'comments:' ; print A.comments
    
    from pyorder.tools.spy import FastSpy
    import matplotlib.pyplot as plt
    fig = plt.figure()
    FastSpy(A.nrow,A.ncol,A.irow,A.jcol,
            sym=(A.symmetric or A.Hermitian or A.skewsym),
            #val=A.values,
            ax=fig.gca(),
            )
    plt.show()

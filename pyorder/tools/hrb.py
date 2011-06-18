# Copyright (c) 2008 dominique.orban@gmail.com
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""
Provides access to sparse linear systems described in Harwell-Boeing or
Rutherford-Boeing format. This module exposes the two classes
HarwellBoeingMatrix and RutherfordBoeingData. For more information, see
the references below.

**References**

.. [DGL] I.S. Duff, R.G. Grimes and J.G. Lewis,
         Sparse Matrix Test Problems, ACM Transactions on
         Mathematical Software, 15(1), p.1-14, 1989

.. [HBUG] `<ftp://ftp.cerfacs.fr/pub/algo/matrices/harwell_boeing/userguide.ps.Z>`_

.. [HB] `<http://math.nist.gov/MatrixMarket/data/Harwell-Boeing>`_

.. [RBC] The Rutherford-Boeing Sparse Matrix Collection,
        I.S. Duff, R.G. Grimes and J.G. Lewis, Technical Report RAL-TR-97-031,
        Rutherford Appleton Laboratory, Chilton, OX, UK, 1997.
        (`<ftp://ftp.numerical.rl.ac.uk/pub/reports/duglRAL97031.pdf>`_)

.. [RB] `<http://www.cerfacs.fr/algor/Softs/RB>`_

.. moduleauthor:: dominique.orban@gerad.ca
"""
import numpy
import FortranFormat as F

class HarwellBoeingMatrix:
    """
    Imports a sparse matrix from a file in Harwell-Boeing format. The matrix is
    stored in compressed sparse row format in (self.ind, self.ip, self.val).
    Right-hand sides, if any, are stored in self.rhs. Right-hand sides can be
    stored as dense vectors, in which case self.rhs has shape (nrow, nrhs), as
    sparse vectors, in which case they are stored in compressed sparse column
    format in (self.rhsptr, self.rhsind, self.rhs), or in elemental format
    (typically when the matrix itself is stored in finite-element format), in
    which case self.rhs has shape (nnzero, nrhs).

    Note that the matrix indices are zero-based, i.e., row indices range
    from 0 through nrow-1 and column indices range from 0 through ncol-1.

    The matrix can be subsequently converted to triple format with
        (row, col) = self.find()

    :keywords:

        :patternOnly:  do not read matrix element values (False)
        :readRhs:      read right-hand sides, if any (False)
        :readGuess:    read starting guess, if any (False)
        :realSol:      read solution vector, if any (False)
    """
    def __init__(self, fname, **kwargs):

        self.title = ''
        self.key = ''
        self.ip = self.ind = self.val = None
        self.nrhs   = 0
        self.rhsptr = self.rhsind = self.rhs = None
        self.guess = self.sol = None

        fp = open(fname, 'r')
        self.readMatrix(fp, **kwargs)
        fp.close()

    def find(self):
        if self.ip is None or self.ind is None: return None
        row = numpy.empty(self.nnzero, numpy.int)
        col = numpy.empty(self.nnzero, numpy.int)
        # Adjust indices to 0-based scheme
        for i in range(self.ncol):
            for j in range(self.ip[i], self.ip[i+1]):
                col[j] = i
                row[j] = self.ind[j]
        return(row,col)

    def readArray(self, fp, which, nelm, format):
        fmt = F.FortranFormat(str(format))
        ind = 0
        while ind < nelm-1:
            fdata = F.FortranLine(fp.readline(), fmt)
            ind2 = min(ind + len(fdata.data), nelm)
            which[ind:ind2] = fdata.data[0:ind2-ind]
            ind = ind2
        return

    def fortranRead(self, stream, format):
        fmt = F.FortranFormat(format)
        fdata = F.FortranLine(stream, fmt)
        return fdata.data

    def readMatrix(self, fp, **kwargs):
        self.patternOnly = kwargs.get('patternOnly', False)
        self.readRhs = kwargs.get('readRhs', False)
        self.readGuess = kwargs.get('readGuess', False)
        self.readSol = kwargs.get('readSol', False)
        fRead = self.fortranRead

        # Read matrix header
        fp.seek(0)
        (self.title, self.key) = fRead(fp.readline(), 'A72,A8')
        (totcrd, ptrcdr, indcrd, valcrd, rhscrd) = fRead(fp.readline(), '5I14')
        (self.mxtype, self.nrow, self.ncol, self.nnzero, neltvl) = \
            fRead(fp.readline(), 'A3,11X,4I14')
        (ptrfmt, indfmt, valfmt, rhsfmt) = fRead(fp.readline(), '2A16,2A20')

        # Decide whether matrix is symmetric and real or complex
        dataType = numpy.float
        if self.mxtype[0]=='C': dataType = numpy.complex
        self.issym = (self.mxtype[1]=='S')

        # Read right-hand side info if present
        if rhscrd > 0:
            (rhstyp, self.nrhs, nrhsix) = fRead(fp.readline(), 'A3,11X,2I14')

        # Set up arrays to hold matrix pattern
        self.ip = numpy.empty(self.ncol+1, dtype=numpy.int)
        self.ind = numpy.empty(self.nnzero, dtype=numpy.int)

        # Read matrix pattern
        self.readArray(fp, self.ip, self.ncol+1, ptrfmt)
        self.readArray(fp, self.ind, self.nnzero, indfmt)

        # Adjust indices to be 0-based
        self.ip -= 1
        self.ind -= 1

        if self.patternOnly or self.mxtype[0] == 'P': return

        # Read matrix values
        if self.mxtype[2] == 'A': # Matrix is assembled
            vallen = self.nnzero
        else:                     # Finite-element format, not assembled
            vallen = neltvl

        self.val = numpy.empty(vallen, dtype=dataType)
        self.readArray(fp, self.val, vallen, valfmt)

        # Read right-hand sides, if any
        if not self.readRhs or self.nrhs == 0: return
        if rhstyp[0] == 'F':
            # Read dense right-hand sides
            self.rhs = numpy.empty((self.nrow, self.nrhs), dataType)
            for i in range(self.nrhs):
                self.readArray(fp, self.rhs[:,i], self.nrow, rhsfmt)
        elif self.mxtype[2] == 'A':
            # Read sparse right-hand sides
            self.rhsptr = numpy.empty(self.nrhs+1, numpy.int)
            self.rhsind = numpy.empty(self.nrhsix, numpy.int)
            self.rhs    = numpy.empty(self.nrhsix, dataType)
            self.readArray(fp, self.rhsptr, self.nrhs+1, ptrfmt)
            self.readArray(fp, self.rhsind, self.nrhsix, indfmt)
            self.readArray(fp, self.rhs,    self.nrhsix, rhsfmt)
            self.rhsind -= 1
            self.rhsptr -= 1
        else:
            # Read elemental right-hand sides
            self.rhs = numpy.empty((self.nnzero, self.nrhs), dataType)
            for i in range(self.nrhs):
                self.readArray(fp, self.rhs[:,i], self.nnzero, rhsfmt)

        # Read initial guesses if requested (always dense)
        if self.readGuess and rhstyp[1] == 'G':
            self.guess = numpy.empty((self.nrow, self.nrhs), dataType)
            for i in range(self.nrhs):
                self.readArray(fp, self.guess[:,i], self.nrow, rhsfmt)

        # Read solution vectors if requested (always dense)
        if self.readSol and rhstyp[2] == 'X':
            self.sol = numpy.empty((self.nrow, self.nrhs), dataType)
            for i in range(self.nrhs):
                self.readArray(fp, self.sol[:,i], self.nrow, rhsfmt)

        return


class RutherfordBoeingData(HarwellBoeingMatrix):
    """
    Imports data from a file in Rutherford-Boeing format. The data is held in
    (self.ind, self.ip, self.val). If the data represents a sparse matrix, the
    three arrays represent the matrix stored in compressed sparse row format.
    Otherwise, the three arrays represent the supplementary data. Refer to the
    Rutherford-Boeing documentation for more information (reference [4] in the
    docstring for the present module.)

    Note that the matrix indices are zero-based, i.e., row indices range from
    0 through nrow-1 and column indices range from 0 through ncol-1.

    The data can be subsequently converted to triple format with
        (row, col) = self.find()

    :keywords:

        :patternOnly:  do not read data values (False)
    """

    def readMatrix(self, fp, **kwargs):
        self.patternOnly = kwargs.get('patternOnly', False)
        fRead = self.fortranRead

        # Read matrix header
        fp.seek(0)
        (self.title, self.key) = fRead(fp.readline(), 'A72,A8')
        (buffer1,) = fRead(fp.readline(), 'A80')
        (buffer2,) = fRead(fp.readline(), 'A80')

        self.issym = False

        if buffer2[2] in ['a', 'e']:
            (self.mxtype, self.m, self.nvec, self.ne, self.neltvl) = \
                fRead(buffer2, 'A3,11X,4(1X,I13)')

            self.nrow = self.m
            self.ncol = self.nvec
            self.nnzero = self.ne
            self.issym = (self.mxtype[1]=='s')
            dataType = numpy.float
            if self.mxtype[0]=='c':
                dataType = numpy.complex
            elif self.mxtype[0]=='i':
                dataType = numpy.int

            (ptrfmt, indfmt, valfmt) = fRead(fp.readline(), '2A16,A20')

            np1 = self.nvec + 1
            if self.mxtype[1:2] == 're': np1 = 2*self.nvec + 1
            self.ip = numpy.empty(np1, numpy.int)
            self.readArray(fp, self.ip, np1, ptrfmt)

            self.ind = numpy.empty(self.ne, numpy.int)
            self.readArray(fp, self.ind, self.ne, indfmt)

            # Adjust indices to be 0-based
            self.ip -= 1
            self.ind -= 1

            # Stop here if pattern only is requested/available
            if self.patternOnly or self.mxtype[0]=='p' or self.mxtype[0]=='x':
                return

            # Read values
            nreal = self.ne
            if self.neltvl > 0: nreal = self.neltvl
            self.val = numpy.empty(nreal, dtype=dataType)
            self.readArray(fp, self.val, nreal, valfmt)

        else:

            # Read supplementary data
            self.dattyp, positn, orgniz, caseid, numerf, m, nvec, nauxd = \
                fRead(buffer1, 'A3, 2A1, 1X, A8, 2X, A1, 3(2X, I13)')
            auxfm1, auxfm2, auxfm3 = fRead(buffer2, '3A20')

            # Read integer data
            if (self.dattyp=='rhs' and orgniz=='s') or \
                    self.dattyp in ['ipt','icv','ord']:
                if self.dattyp=='ord':
                    nauxd = m * nvec
                    auxfm = auxfm1
                else:
                    self.ip = numpy.empty(nvec+1, numpy.int)
                    self.readArray(fp, self.ip, nvec+1, auxfm1)
                    auxfm = auxfm2
                    self.ip -= 1 # Adjust to 0-based indexing
                self.ind = numpy.empty(nauxd, numpy.int)
                self.readArray(fp, self.ind, nauxd, auxfm)
                self.ind -= 1 # Adjust to 0-based indexing
                if self.dattyp != 'rhs': return

            if self.patternOnly: return

            # Read values
            dataType = numpy.float
            if numerf=='c':
                dataType = numpy.complex
            elif numerf=='i':
                dataType = numpy.int
            if self.dattyp != 'rhs': nauxd = m*nvec
            if self.dattyp == 'rhs' and orgniz == 's':
                auxfm = auxfm3
            else:
                auxfm = auxfm1
            self.val = numpy.empty(nauxd, dataType)
            self.readArray(fp, self.val, nauxd, auxfm)

        return


if __name__ == '__main__':

    # Demo the HarwellBoeingMatrix and RutherfordBoeingData classes
    # In case of Rutherford-Boeing supplementary data, some of the following
    # may not make sense.

    import sys
    fname = sys.argv[1]
    #M = HarwellBoeingMatrix(fname, patternOnly=False, readRhs=True)
    M = RutherfordBoeingData(fname, patternOnly=False)
    print 'Data of order (%-d,%-d) with %-d nonzeros' % (M.nrow,M.ncol,M.nnzero)

    numpy.set_printoptions(precision=8,
                           threshold=10,
                           edgeitems=3,
                           linewidth=80,
                           suppress=False)

    # Comment this out for Rutherford-Boeing data
    #if M.readRhs:
    #    for i in range(M.nrhs):
    #        print M.rhs[:,i]

    # Plot sparsity pattern
    try:
        import pylab
    except:
        sys.stderr.write('Pylab is required for the demo\n')
        sys.exit(1)

    (row,col) = M.find()
    fig = pylab.figure()
    ax = fig.gca()
    ax.plot(col, row, 'ks', markersize=1, linestyle='None')
    if M.issym: ax.plot(row, col, 'ks', markersize=1, linestyle='None')
    ax.set_xlim(xmin=-1, xmax=M.ncol)
    ax.set_ylim(ymin=M.nrow, ymax=-1)
    if M.nrow == M.ncol:
        ax.set_aspect('equal')
    pylab.title(M.title, fontsize='small')
    pylab.show()

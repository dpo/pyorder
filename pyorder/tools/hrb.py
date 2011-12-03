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
from conversion import csc2coord, coord2csc
from fortranformat import FortranRecordReader, FortranRecordWriter
from fortranformat import config as ff_config

# Possibility to set ff_config.RET_WRITTEN_VARS_ONLY = True
# See:
# https://bitbucket.org/brendanarnold/py-fortranformat/issue/1/reading-less-records-than-specified-by#comment-531697

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
        (val, row, col) = self.find()

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
        self._readMatrix(fp, **kwargs)
        fp.close()

    def data(self):
        """
        Return matrix data in compressed sparse row format.

        :returns:
            :val:   array of values of nonzero elements
            :ip:    array of pointer to rows
            :ind:   array of indices of nonzero elements in each row
        """
        #if self.ip is None or self.ind is None: return (None,None,None)
        return (self.val, self.ip, self.ind)

    def find(self):
        """
        Return matrix data in coordinate format.

        :returns:
            :val:   array of values of nonzero elements
            :irow:   array of indices of nonzero elements in each row
            :jcol:   array of indices of nonzero elements in each column
        """
        if self.ip is None or self.ind is None: return (None,None,None)
        irow, jcol = csc2coord(self.ind, self.ip)
        return (self.val, irow, jcol)

    def _readArray(self, fp, which, nelm, format):
        #print 'Reading %d values with format %s' % (nelm, format)
        fmt = FortranRecordReader(format)
        ind = 0
        while ind < nelm-1:
            fdata = fmt.read(fp.readline())
            ind2 = min(ind + len(fdata), nelm)
            which[ind:ind2] = fdata[:ind2-ind]
            ind = ind2
        # Read last line, if any.
        if ind < nelm:
            fdata = fmt.read(fp.readline())
            which[ind:] = fdata[:nelm-ind]
        return

    def _fortranRead(self, stream, format):
        fmt = FortranRecordReader(format)
        fdata = fmt.read(stream)
        return fdata

    def _readMatrix(self, fp, **kwargs):
        self.patternOnly = kwargs.get('patternOnly', False)
        self.readRhs = kwargs.get('readRhs', False)
        self.readGuess = kwargs.get('readGuess', False)
        self.readSol = kwargs.get('readSol', False)
        fRead = self._fortranRead

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
        self._readArray(fp, self.ip, self.ncol+1, ptrfmt)
        self._readArray(fp, self.ind, self.nnzero, indfmt)

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
        self._readArray(fp, self.val, vallen, valfmt)

        # Read right-hand sides, if any
        if not self.readRhs or self.nrhs == 0: return
        if rhstyp[0] == 'F':
            # Read dense right-hand sides
            self.rhs = numpy.empty((self.nrow, self.nrhs), dataType)
            for i in range(self.nrhs):
                self._readArray(fp, self.rhs[:,i], self.nrow, rhsfmt)
        elif self.mxtype[2] == 'A':
            # Read sparse right-hand sides
            self.rhsptr = numpy.empty(self.nrhs+1, numpy.int)
            self.rhsind = numpy.empty(self.nrhsix, numpy.int)
            self.rhs    = numpy.empty(self.nrhsix, dataType)
            self._readArray(fp, self.rhsptr, self.nrhs+1, ptrfmt)
            self._readArray(fp, self.rhsind, self.nrhsix, indfmt)
            self._readArray(fp, self.rhs,    self.nrhsix, rhsfmt)
            self.rhsind -= 1
            self.rhsptr -= 1
        else:
            # Read elemental right-hand sides
            self.rhs = numpy.empty((self.nnzero, self.nrhs), dataType)
            for i in range(self.nrhs):
                self._readArray(fp, self.rhs[:,i], self.nnzero, rhsfmt)

        # Read initial guesses if requested (always dense)
        if self.readGuess and rhstyp[1] == 'G':
            self.guess = numpy.empty((self.nrow, self.nrhs), dataType)
            for i in range(self.nrhs):
                self._readArray(fp, self.guess[:,i], self.nrow, rhsfmt)

        # Read solution vectors if requested (always dense)
        if self.readSol and rhstyp[2] == 'X':
            self.sol = numpy.empty((self.nrow, self.nrhs), dataType)
            for i in range(self.nrhs):
                self._readArray(fp, self.sol[:,i], self.nrow, rhsfmt)

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
        (val, row, col) = self.find()

    :keywords:

        :patternOnly:  do not read data values (False)
    """

    def __init__(self, fname, **kwargs):

        HarwellBoeingMatrix.__init__(self, fname, **kwargs)
        self.transposed = kwargs.get('transposed', False)

    def _readMatrix(self, fp, **kwargs):
        self.patternOnly = kwargs.get('patternOnly', False)
        fRead = self._fortranRead

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
            self._readArray(fp, self.ip, np1, ptrfmt)

            self.ind = numpy.empty(self.ne, numpy.int)
            self._readArray(fp, self.ind, self.ne, indfmt)

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
            self._readArray(fp, self.val, nreal, valfmt)

        else:

            # Read supplementary data
            (self.dattyp, self.positn, self.orgniz, self.caseid, numerf, m,
            nvec, self.nauxd) = fRead(buffer1,
                                      'A3,2A1,1X,A8,2X,A1,3(2X,I13)')
            auxfm1, auxfm2, auxfm3 = fRead(buffer2, '3A20')
            self.nrow = m
            self.nvec = self.ncol = nvec

            # Read integer data
            if (self.dattyp=='rhs' and self.orgniz=='s') or \
                    self.dattyp in ['ipt','icv','ord']:
                if self.dattyp=='ord':
                    self.nauxd = m * nvec
                    auxfm = auxfm1
                else:
                    self.ip = numpy.empty(nvec+1, numpy.int)
                    self._readArray(fp, self.ip, nvec+1, auxfm1)
                    auxfm = auxfm2
                    self.ip -= 1 # Adjust to 0-based indexing
                self.ind = numpy.empty(self.nauxd, numpy.int)
                self._readArray(fp, self.ind, self.nauxd, auxfm)
                self.ind -= 1 # Adjust to 0-based indexing
                if self.dattyp != 'rhs': return

            if self.patternOnly: return

            # Read values
            dataType = numpy.float
            if numerf=='c':
                dataType = numpy.complex
            elif numerf=='i':
                dataType = numpy.int
            if self.dattyp != 'rhs': self.nauxd = m*nvec
            if self.dattyp == 'rhs' and self.orgniz == 's':
                auxfm = auxfm3
            else:
                auxfm = auxfm1
            self.val = numpy.empty(self.nauxd, dataType)
            self._readArray(fp, self.val, self.nauxd, auxfm)
            self.nnzero = self.nauxd

        return

    def _mul(self, other):
        # No type or dimension checking for now...
        y = numpy.zeros(self.nrow)
        for col in range(self.ncol):
            for k in range(self.ip[col],self.ip[col+1]):
                row = self.ind[k]
                val = self.val[k]
                y[row] += val * other[col]
                if self.issym and row != col:
                    y[col] += val * other[row]
        return y

    def _rmul(self, other):
        # No type or dimension checking for now...
        y = numpy.zeros(self.ncol)
        for col in range(self.ncol):
            for k in range(self.ip[col],self.ip[col+1]):
                row = self.ind[k]
                val = self.val[k]
                y[col] += val * other[row]
                if self.issym and row != col:
                    y[row] += val * other[col]
        return y

    def __mul__(self, other):
        if self.transposed:
            return self._rmul(other)
        else:
            return self._mul(other)


# Helper functions.

def get_int_fmt(n):
    "Return appropriate format for integer arrays."
    fmts = ['(40I2)', '(26I3)', '(20I4)', '(16I5)', '(13I6)', '(11I7)',
            '(10I8)', '(8I9)', '(8I10)', '(7I11)', '(4I20)']
    nlines = [40,26,20,16,13,11,10,8,8,7,4]
    nn = n
    for k in range(1,n+1):
        if nn < 10: break
        nn /= 10
    if k <= 10: return (fmts[k+1], nlines[k+1])
    return (fmts[10], nlines[10])

def get_real_fmt(p):
    "Return appropriate format for real array (1 <= p <= 17)."
    fmt = ['(8E10.1E3)', '(7E11.2E3)', '(6E12.3E3)', '(6E13.4E3)',
            '(5E14.5E3)', '(5E15.6E3)', '(5E16.7E3)', '(4E17.8E3)',
            '(4E18.9E3)', '(4E19.10E3)', '(4E20.11E3)', '(3E21.12E3)',
            '(3E22.13E3)', '(3E23.14E3)', '(3E24.15E3)', '(3E25.16E3)']
    fmt1 = ['(1P,8E10.1E3)', '(1P,7E11.2E3)', '(1P,6E12.3E3)',
            '(1P,6E13.4E3)', '(1P,5E14.5E3)', '(1P,5E15.6E3)',
            '(1P,5E16.7E3)', '(1P,4E17.8E3)', '(1P,4E18.9E3)',
            '(1P,4E19.10E3)', '(1P,4E20.11E3)', '(1P,3E21.12E3)',
            '(1P,3E22.13E3)', '(1P,3E23.14E3)', '(1P,3E24.15E3)',
            '(1P,3E25.16E3)']
    lens = [8,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3]
    return (fmt[p-2], fmt1[p-2], lens[p-2])

def fortranWriteLine(data, stream, fformat):
    "Write `data` to `stream` according to Fortran format `fformat`."
    fmt = FortranRecordWriter(fformat)
    stream.write(fmt.write(data))
    stream.write('\n')
    return

def fortranWriteArray(data, chunk_size, stream, fformat):
    """
    Write array `data` to `stream`, possibly using multiple lines,
    according to Fortran format `fformat`.
    """
    nelts = len(data)
    nelts_per_line = nelts/chunk_size
    #print 'Writing %d elements %d per line...' % (nelts, chunk_size)
    for k in range(nelts_per_line):
        chunk = data[k*chunk_size:(k+1)*chunk_size]
        fortranWriteLine(chunk, stream, fformat)
    if nelts_per_line*chunk_size < nelts:
        fortranWriteLine(data[nelts_per_line*chunk_size:], stream, fformat)
    return

# End of helper functions.


def write_rb(fname, nrow, ncol, ip, ind,
             val=None, precision=17, symmetric=False, skew=False,
             title='Generic', key='Generic'):
    """
    Write a sparse matrix to file in Rutherford-Boeing format. The input matrix
    must be described in compressed column format by the arrays `ip` and `ind`.
    Rows must be ordered in each column.
    If numerical values are to be written to file, they should be specified in
    the array `val`.

    Currently only supports assembled matrices in compressed column format.
    """

    patternOnly = (val == None)
    rectangular = (nrow != ncol)
    ne = len(ind)

    # Check that columns are in order.
#    for j in range(ncol):
#        if ip[j+1] <= ip[j]:
#            raise ValueError, 'Columns must be ordered.'

    # Check that rows are in order in each column.
#    for j in range(ncol):
#        for k in range(ip[j], ip[j+1]-1):
#            if ind[k] >= ind[k+1]:
#                raise ValueError, 'Rows must be ordered in each column.'

    # Set mxtype.
    mxtype0 = 'r'
    if patternOnly: mxtype0 = 'p'
    if symmetric: mxtype1 = 's'
    if skew: mxtype1 = 'z'
    if rectangular: mxtype1 = 'r'

    mxtype = mxtype0 + mxtype1 + 'a'

    # Set format and number card images for pointer array.
    (ptrfmt, ptrn) = get_int_fmt(ne+1)
    ptrcrd = ncol/ptrn + 1

    # Set format and number card images for index array.
    (indfmt, indn) = get_int_fmt(nrow)
    indcrd = (ne-1)/indn + 1

    # Set number of card images for numerical entries.
    if patternOnly:
        valcrd = 0 ; valfmi = ' '
    else:
        (valfmi, valfmo, valn) = get_real_fmt(precision)
        valcrd = (ne-1)/valn + 1

    totcrd = ptrcrd + indcrd + valcrd
    neltvl = 0

    fp = open(fname, 'w')

    lt = len(title)
    if lt < 72: title = title + (72-lt)*' '
    lk = len(key)
    if lk < 8: key = key + (8-lk)*' '

    # Write header.
    fortranWriteLine([title, key], fp, 'A72,A8')
    fortranWriteLine([totcrd, ptrcrd, indcrd, valcrd], fp, 'I14,3(1X,I13)')
    fortranWriteLine([mxtype, nrow, ncol, ne, neltvl], fp, 'A3,11X,4(1X,I13)')
    fortranWriteLine([ptrfmt, indfmt, valfmi], fp, '2A16,A20')

    # Write pointer and index arrays. Ajust for 1-based indexing.
    #print 'Writing pointer array...'
    fortranWriteArray(ip+1, ptrn, fp, ptrfmt)
    #print 'Writing index array...'
    fortranWriteArray(ind+1, indn, fp, indfmt)

    # Write matrix entries.
    neltvl = ne
    if not patternOnly:
        #print 'Writing matrix entries...'
        fortranWriteArray(val, valn, fp, valfmo)

    fp.close()
    return


def write_rb_from_rb(fname, matrix):
    """
    Convenience function to write a matrix in RB format to file.
    """
    (ip, ind, val) = matrix.data()
    write_rb(fname, matrix.nrow, matrix.ncol, ip, ind, val,
             symmetric=matrix.issym)
    return


def write_rb_from_coord(fname, nrow, ncol, irow, jcol, val=None, **kwargs):
    """
    Convenience function to write a matrix in coordinate format to file
    in RB format.
    """
    (ind, ip, values) = coord2csc(ncol, irow, jcol, val)
    write_rb(fname, nrow, ncol, ip, ind, val, **kwargs)
    return


def write_aux(fname, nrow, nvec, precision=17, title='Generic', key='Generic',
        caseid='Generic', dattyp='rhs', positn='r', orgniz='d', nauxd=None,
        ip=None, ind=None, val=None):
    """
    Write supplementary data to file in Rutherford-Boeing format.

    Only real data is supported for now.
    """

    data_types = ['ord','rhs','sln','est','evl','svl','evc','svc','sbv','sbm',
                  'sbp','ipt','icv','lvl','geo','avl']
    organizations = ['s','d','e']
    positions = ['r','l','s']

    if dattyp not in data_types:
        raise ValueError, 'Unknown data type: %s' % dattyp
    if positn not in positions:
        raise ValueError, 'Unknown position: %s' % positn
    if orgniz not in organizations:
        raise ValueError, 'Unknown organization: %s' % orgniz

    if dattyp in ['evl','svl','lvl','sbp']: nvec = 1
    if dattyp in ['evl','svl','lvl','sbv','sbm','sbp','avl']: positn = ' '
    if dattyp != 'rhs': orgniz = ' '

    numerf = 'r'
    if dattyp == 'ord': numerf = 'i'
    if dattyp in ['ipt','icv']: numerf = 'p'

    if orgniz != 'e':
        nauxd = 0
        if orgniz == 's' or dattyp in ['icv','ipt']:
            if ip is None:
                raise ValueError, 'Need pointer array for data type %s' % dattyp
            nauxd = ip(nvec)-1
        if orgniz == 'd': nauxd = nrow*nvec

    # Set data formats.
    auxfm1 = auxfm2 = auxfm3 = ' '
    fm1 = fm3 = ' '

    if dattyp in ['ipt','icv']:
        (auxfm1, n1) = get_int_fmt(nauxd+1)
        (auxfm2, n2) = get_int_fmt(nrow)
    elif dattyp == 'ord':
        (auxfm1, n1) = get_int_fmt(nrow)
    else:
        if precision < 2 or precision > 17: precision = 17
        if dattyp == 'rhs' and orgniz == 's':
            (auxfm1, n1) = get_int_fmt(nauxd+1)
            (auxfm2, n2) = get_int_fmt(nrow)
            (auxfm3, fm3, n3) = get_real_fmt(precision)
        else:
            (auxfm1, fm1, n1) = get_real_fmt(precision)

    fp = open(fname, 'w')

    lt = len(title)
    if lt < 72: title = title + (72-lt)*' '
    lk = len(key)
    if lk < 8: key = key + (8-lk)*' '

    # Write header.
    fortranWriteLine([title, key], fp, 'A72,A8')
    fortranWriteLine([dattyp, positn, orgniz, caseid, numerf, nrow, nvec,
        nauxd], fp, 'A3,2A1,1X,A8,1X,A1,3(1X,I13)')
    fortranWriteLine([auxfm1, auxfm2, auxfm3], fp, '3A20')

    # Write pointer and index arrays. Ajust for 1-based indexing.
    if (dattyp == 'rhs' and orgniz == 's') or dattyp in ['ipt','icv']:
        #print 'Writing pointer array...'
        fortranWriteArray(ip+1, n1, fp, auxfm1)
        #print 'Writing index array...'
        fortranWriteArray(ind+1, n2, fp, auxfm2)

    # Write entries.
    #print 'Writing entries...'
    if dattyp == 'rhs' and orgniz == 's':
        fortranWriteArray(val, n3, fp, fm3)
    else:
        fortranWriteArray(val, n1, fp, fm1)

    fp.close()
    return


def write_aux_from_rb(fname, rbdata):
    """
    Convenience function to write supplementary data in RB format to file.
    """
    (ip, ind, val) = rbdata.data()
    write_aux(fname, rbdata.nrow, ip.shape[0], dattyp=rbdata.dattyp,
              title=rbdata.title, key=rbdata.key, caseid=rbdata.caseid,
              positn=rbdata.positn, orgniz=rbdata.orgniz, nauxd=rbdata.nauxd,
              ip=ip, ind=ind, val=val)

def write_aux_from_numpy(fname, array, **kwargs):
    """
    Convenience function write a numpy array to file in RB format.
    """
    if len(array.shape) == 1:
        nvec = 1
    else:
        nvec = array.shape[1]
    write_aux(fname, array.shape[0], nvec, val=array, **kwargs)


if __name__ == '__main__':

    # Demo the HarwellBoeingMatrix and RutherfordBoeingData classes
    # In case of Rutherford-Boeing supplementary data, some of the following
    # may not make sense.

    import sys
    fname = sys.argv[1]
    plot = False

    numpy.set_printoptions(precision=8,
                           threshold=10,
                           edgeitems=3,
                           linewidth=80,
                           suppress=False)

    #M = HarwellBoeingMatrix(fname, patternOnly=False, readRhs=True)
    # Comment this out for Rutherford-Boeing data
    #if M.readRhs:
    #    for i in range(M.nrhs):
    #        print M.rhs[:,i]

    M = RutherfordBoeingData(fname, patternOnly=False)
    print 'Data of order (%-d,%-d) with %-d nonzeros' % (M.nrow,M.ncol,M.nnzero)

    # Plot sparsity pattern
    if plot:
        try:
            import pylab
        except:
            sys.stderr.write('Pylab is required for the demo\n')
            sys.exit(1)

        (val, row,col) = M.find()
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

    # Write data back to file
    #(ip, ind, val) = M.data()
    #print val
    #write_rb('newmat.rb', M.nrow, M.ncol, ip, ind, val, symmetric=M.issym)
    x = numpy.ones(M.ncol)
    #y = M*x
    #write_aux_from_numpy('rhs.rb', y)
    write_aux_from_numpy('sol.rb', x)

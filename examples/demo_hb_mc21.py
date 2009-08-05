"""
Illustrate usage of the pymc21 module, using an input matrix in Harwell-Boeing
or Rutherford-Boeing format. Supply a file name as input argument on the
command line and uncomment below as appropriate.
"""

import sys
import numpy as np
from pyorder.pymc21.pymc21 import nonzerodiag
from pyorder.tools.hrb import HarwellBoeingMatrix, RutherfordBoeingData
from pyorder.tools.spy import FastSpy
import pylab

if len(sys.argv) < 2:
    sys.stderr.write('Supply input matrix as argument\n')
    sys.exit(1)

fname = sys.argv[1]
#M = HarwellBoeingMatrix(fname, patternOnly=True, readRhs=False)
M = RutherfordBoeingData(fname, patternOnly=True, readRhs=False)

if M.nrow != M.ncol:
    sys.stderr.write('Input matrix must be square\n')
    sys.exit(1)

perm, nzdiag = nonzerodiag(M.nrow, M.ind, M.ip)

(irow, jcol) = M.find()
left = pylab.subplot(121)
FastSpy(M.nrow, M.ncol, irow, jcol, sym=M.issym, ax=left.get_axes())

right = pylab.subplot(122)
FastSpy(M.nrow, M.ncol, perm[irow], jcol, sym=M.issym, ax=right.get_axes())
pylab.show()

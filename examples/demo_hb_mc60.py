"MC60 demo with input matrix in HB or RB format"

import numpy as np
from pyorder.pymc60 import sloan, rcmk
from pyorder.tools.hrb import HarwellBoeingMatrix, RutherfordBoeingData
from pyorder.tools.spy import FastSpy
import pylab
import sys

if len(sys.argv) < 2:
    sys.stderr.write('Data file name must be supplied\n')
    sys.exit(1)

fname = sys.argv[1]
#M = HarwellBoeingMatrix(fname, patternOnly=True, readRhs=False)
M = RutherfordBoeingData(fname, patternOnly=True, readRhs=False)

if M.nrow != M.ncol or not M.issym:
    sys.stderr.write('Input matrix must be square and symmetric\n')
    sys.exit(1)

# Compute reverse Cuthill-McKee ordering
perm, rinfo = rcmk(M.nrow, M.ind, M.ip)

# Or: Compute Sloan's ordering
#perm, rinfo = sloan(M.nrow, M.ind, M.ip)


# Plot original matrix
(_, irow, jcol) = M.find()
left = pylab.subplot(121)
FastSpy(M.nrow, M.ncol, irow, jcol, sym=M.issym,
        ax=left.get_axes(), title='Original')

# Apply permutation and plot reordered matrix
right = pylab.subplot(122)
FastSpy(M.nrow, M.ncol, perm[irow], perm[jcol], sym=M.issym,
        ax=right.get_axes(), title='Reordered')
pylab.show()

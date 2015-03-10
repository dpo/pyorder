"""
MC60 demo with input matrix in HB or RB format. We use the mc60module directly.
This allows the user to specify a supervariable ordering and pseudoperipheral
pairs directly. See the MC60 documentation.
"""

import sys
from pyorder.pymc60 import mc60module
import numpy as np
from pyorder.tools.hrb import HarwellBoeingMatrix, RutherfordBoeingData

if len(sys.argv) < 2:
    sys.stderr.write('Data file name must be supplied\n')
    sys.exit(1)

fname = sys.argv[1]
#M = HarwellBoeingMatrix(fname, patternOnly=True, readRhs=False)
M = RutherfordBoeingData(fname, patternOnly=True, readRhs=False)

if M.nrow != M.ncol:
    sys.stderr.write('Input matrix must be square\n')
    sys.exit(1)

if not M.issym:
    sys.stderr.write('Input matrix must be symmetric\n')
    sys.exit(1)

icntl = np.array([0,6], dtype=np.int32) # Abort on error
#jcntl = np.array([0,0], dtype=np.int32) # Sloan's alg with auto choice
jcntl = np.array([1,0], dtype=np.int32) # RCM alg with auto choice
weight = np.array([2,1])                # Weights in Sloan's alg

# Store lower triangle of symmetric matrix in csr format
# We make copies of M.ip and M.ind because MC60 will modify them
n = M.nrow
icptr = M.ip.copy() + 1  # Convert to 1-based indexing (for Fortran)

# Make room for pattern of the whole matrix
irn = np.empty(2*(icptr[-1] - 1), dtype=np.int32)
irn[:icptr[-1]-1] = M.ind.copy() + 1

# Check data
info = mc60module.mc60ad(irn, icptr, icntl)

# Compute supervariables
nsup, svar, vars = mc60module.mc60bd(irn, icptr)
print 'The number of supervariables is ', nsup

# Permute reduced matrix
permsv = np.empty(nsup, dtype=np.int32)
pair = np.empty((2,nsup/2), dtype=np.int32)
info = mc60module.mc60cd(n, irn, icptr[:nsup+1], vars[:nsup],
                         jcntl, permsv, weight, pair)

# Compute profile and wavefront
rinfo = mc60module.mc60fd(n, irn, icptr[:nsup+1], vars[:nsup], permsv)

# Obtain variable permutation from supervariable permutation
perm, possv = mc60module.mc60dd(svar, vars[:nsup], permsv)

np.set_printoptions(precision=3, linewidth=80, threshold=10, edgeitems=3)
print 'The variable permutation is ', perm
print 'The profile is ', rinfo[0]
print 'The maximum wavefront is ', rinfo[1]
print 'The semibandwidth is ', rinfo[2]
print 'The root-mean-square wavefront is ', rinfo[3]

try:
    import pylab
    from pyorder.tools.spy import FastSpy
    # Plot original matrix
    (_, irow, jcol) = M.find()
    left = pylab.subplot(121)
    right = pylab.subplot(122)
    FastSpy(M.nrow, M.ncol, irow, jcol, sym=M.issym,
            ax=left.get_axes(), title='Original')

    # Apply permutation and plot reordered matrix
    perm -= 1   # Convert to 0-based indexing
    FastSpy(M.nrow, M.ncol, perm[irow], perm[jcol], sym=M.issym,
            ax=right.get_axes(), title='Reordered')
    pylab.show()
except:
    pass

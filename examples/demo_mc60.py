"MC60 demo from the HSL MC60 spec sheet"
import numpy as np
from pyorder.pymc60 import mc60module

icntl = np.array([0,6], dtype=np.int32) # Abort on error
jcntl = np.array([0,0], dtype=np.int32) # Sloan's alg with auto choice
weight = np.array([2,1])                # Weights in Sloan's alg

# Store lower triangle of symmetric matrix in csr format (1-based)
n = 5
icptr = np.array([1,6,8,9,10,11], dtype=np.int32)  # nnz = 10
irn = np.empty(2*(icptr[-1]-1), dtype=np.int32)
irn[:icptr[-1]-1] = np.array([1,2,3,4,5,2,3,3,4,5], dtype=np.int32)

# Check data
info = mc60module.mc60ad(irn, icptr, icntl)

# Compute supervariables
nsup, svar, vars = mc60module.mc60bd(irn, icptr)
print 'The number of supervariables is ', nsup

# Permute reduced matrix
permsv = np.empty(nsup, dtype=np.int32)
pair = np.empty((2, nsup/2), dtype=np.int32)
info = mc60module.mc60cd(n,irn,icptr[:nsup+1],vars[:nsup],jcntl,permsv,weight,pair)

# Compute profile and wavefront
rinfo = mc60module.mc60fd(n, irn, icptr[:nsup+1], vars[:nsup], permsv)

# Obtain variable permutation from supervariable permutation
perm, possv = mc60module.mc60dd(svar, vars[:nsup], permsv)
print 'The variable permutation is ', perm    # 1-based
print 'The profile is ', rinfo[0]
print 'The maximum wavefront is ', rinfo[1]
print 'The semibandwidth is ', rinfo[2]
print 'The root-mean-square wavefront is ', rinfo[3]

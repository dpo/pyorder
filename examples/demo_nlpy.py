
try:
    from nlpy.model import AmplModel
except:
    msg='NLPy is required to run this demo. See http://nlpy.sf.net'
    raise RuntimeError, msg

from pyorder.tools import coord2csc
from pyorder.pymc60 import sloan, rcmk
from pyorder.tools.spy import FastSpy
import numpy as np
import matplotlib.pyplot as plt

nlp = AmplModel('truss18bars.nl')
x = np.random.random(nlp.n)
y = np.random.random(nlp.m)
H = nlp.hess(x,y)
(val,irow,jcol) = H.find()
(rowind, colptr, values) = coord2csc(nlp.n, irow, jcol, val)  # Convert to CSC

perm1, rinfo1 = rcmk(nlp.n, rowind, colptr)   # Reverse Cuthill-McKee
perm2, rinfo2 = sloan(nlp.n, rowind, colptr)  # Sloan's method

left = plt.subplot(131)
FastSpy(nlp.n, nlp.n, irow, jcol, sym=True,
        ax=left.get_axes(), title='Original')

# Apply permutation 1 and plot reordered matrix
middle = plt.subplot(132)
FastSpy(nlp.n, nlp.n, perm1[irow], perm1[jcol], sym=True,
        ax=middle.get_axes(), title='Rev. Cuthill-McKee (semibandwidth=%d)' % rinfo1[2])

# Apply permutation 2 and plot reordered matrix
right = plt.subplot(133)
FastSpy(nlp.n, nlp.n, perm2[irow], perm2[jcol], sym=True,
        ax=right.get_axes(), title='Sloan (semibandwidth=%d)' % rinfo2[2])

#plt.savefig('mpvc.pdf',bbox_inches='tight')
plt.show()

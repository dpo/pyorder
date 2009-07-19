import numpy as np
from pyorder.pymc21 import mc21module
n = 4
icn = np.array([1,4,3,4,1,4,2,4], dtype=np.int32)
ip = np.array([1,3,5,7], dtype=np.int32)
lenr = np.array([2,2,2,2], dtype=np.int32)
iperm,numnz = mc21module.mc21ad(icn,ip,lenr)
print 'iperm = ', iperm, ' (1-based)'
print 'numnz = ', numnz

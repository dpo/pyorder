"""
Compute a row permutation so as to bring nonzero elements on the diagonal of a
square sparse matrix in compressed sparse row (csr) format.
This module provides an interface to the HSL subroutine MC21.
"""

from pymc21 import *

__all__ = filter(lambda s:not s.startswith('_'), dir())

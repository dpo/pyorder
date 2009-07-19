"""
Profile/wavefront/semi-bandwidth reduction by way of Sloan's and the reverse
Cuthill-McKee algorithms. This module provides an interface to the HSL
subroutine MC60.
"""

from pymc60 import *

__all__ = filter(lambda s:not s.startswith('_'), dir())

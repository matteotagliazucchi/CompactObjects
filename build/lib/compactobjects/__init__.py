import sys

from .compact_stars import CompactStars
from .eos import ImplicitEos, PolytropicEos
from .constants import *

if sys.version_info < (3,):
    raise ImportError(
"""You are running 'compactobjects' on Python 2.

'compactobjects' is designed to work with Python 3
Please upgrade to Python 3.
""")
import sys

from .compact_star import CompactStar
from .eos import ImplicitEos, PressureEdenPolytropic, PressureDensityPolytropic, PressureDensityPiecewise
from .utils import conversion_dict, eos_lib
from .read_eos import glue_crust_core_eos

if sys.version_info < (3,):
    raise ImportError(
"""You are running 'compactobjects' on Python 2.

'compactobjects' is designed to work with Python 3
Please upgrade to Python 3.
""")
"""

ViennaLS
========

ViennaLS is a level set library developed for high performance
topography simulations. The main design goals are simplicity and efficiency,
tailored towards scientific simulations. ViennaLS can also be used for
visualisation applications, although this is not the main design target.
"""

from __future__ import annotations
import sys as sys
from viennals._core import setNumThreads
from . import _core
from . import common
from . import d2
from . import d3

__all__: list = ["d2", "d3", "common", "setNumThreads", "__version__"]

def _windows_dll_path(): ...

__version__: str = "VIENNALS_VERSION"

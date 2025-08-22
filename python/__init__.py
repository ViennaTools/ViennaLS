"""
ViennaLS
========

ViennaLS is a level set library developed for high performance
topography simulations. The main design goals are simplicity and efficiency,
tailored towards scientific simulations. ViennaLS can also be used for
visualisation applications, although this is not the main design target.
"""

import sys


def _windows_dll_path():

    import os

    additional_paths = [
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "viennals.libs")
    ]

    for path in additional_paths:
        os.add_dll_directory(path)
        os.environ["PATH"] = path + os.pathsep + os.environ["PATH"]


if sys.platform == "win32":
    _windows_dll_path()

from . import _core as _C  # the binary inside the package
import sys as _sys

d2 = _C.d2
d3 = _C.d3
common = _C.common
_sys.modules[__name__ + ".d2"] = d2
_sys.modules[__name__ + ".d3"] = d3
_sys.modules[__name__ + ".common"] = common

# top-level API
setNumThreads = _C.setNumThreads
__version__ = getattr(_C, "__version__", None)
__all__ = ["d2", "d3", "common", "setNumThreads", "__version__"]
del _C, _sys

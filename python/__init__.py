"""
ViennaLS
========

ViennaLS is a level set library developed for high performance
topography simulations. The main design goals are simplicity and efficiency,
tailored towards scientific simulations. ViennaLS can also be used for
visualisation applications, although this is not the main design target.
"""


def _windows_dll_path():

    import os

    additional_paths = [
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "viennals.libs")
    ]

    for path in additional_paths:
        os.add_dll_directory(path)
        os.environ["PATH"] = path + os.pathsep + os.environ["PATH"]


import sys as _sys
from . import _core as _C  # the binary inside the package


if _sys.platform == "win32":
    _windows_dll_path()

# bring d2 and d3 into the top-level namespace
d2 = _C.d2
d3 = _C.d3
_sys.modules[__name__ + ".d2"] = d2
_sys.modules[__name__ + ".d3"] = d3


# forward any other (common) names to _core (PEP 562)
def __getattr__(name):
    return getattr(_C, name)


def __dir__():
    return sorted(set(globals()) | set(dir(_C)))

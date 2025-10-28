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

if _sys.platform == "win32":
    _windows_dll_path()

from . import _core as _C  # the binary inside the package

# bring d2 and d3 into the top-level namespace
d2 = _C.d2
d3 = _C.d3
_sys.modules[__name__ + ".d2"] = d2
_sys.modules[__name__ + ".d3"] = d3
PROXY_DIM = 2  # default dimension for proxy classes


def setDimension(d: int):
    """Set the dimension of the simulation (2 or 3).

    Parameters
    ----------
    d: int
        Dimension of the simulation (2 or 3).
    """
    global PROXY_DIM
    if d == 2 or d == 3:
        PROXY_DIM = d
    else:
        raise ValueError("Dimension must be 2 or 3.")


def __getattr__(name):
    # 1) common/top-level from _core
    try:
        return getattr(_C, name)
    except AttributeError as e_core:
        pass
    # 2) fallback to current default dimension
    m = d2 if PROXY_DIM == 2 else d3
    try:
        return getattr(m, name)
    except AttributeError:
        raise AttributeError(
            f"module {__name__!r} has no attribute {name!r}"
        ) from e_core


def __dir__():
    return sorted(set(globals()) | set(dir(_C)) | set(dir(d2)) | set(dir(d3)))

"""
ViennaLS is a header-only C++ level set library developed for high performance topography simulations. The main design goals are simplicity and efficiency, tailored towards scientific simulations. ViennaLS can also be used for visualization applications, although this is not the main design target.
"""

from __future__ import annotations
import typing
from viennals import common
from viennals import d2
from viennals import d3

__all__: list[str] = ["common", "d2", "d3", "setNumThreads", "version"]

def setNumThreads(arg0: typing.SupportsInt) -> None: ...

__version__: str = "VIENNALS_VERSION"
version: str = "VIENNALS_VERSION"

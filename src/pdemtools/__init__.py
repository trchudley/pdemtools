from importlib.metadata import version

# import submodules
from . import load
from . import data
from ._index_search import search

# import xarray accessor
from ._accessor import DemAccessor

# import hidden modules for test purposes
from . import _coreg
from . import _geomorphometry
from . import _utils

__version__ = "1.1.0"

__all__ = ["search", "DemAccessor"]

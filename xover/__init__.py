"""
xover: crossover analysis tools for GLODAP
==========================================

The current version of xover does not contain code for running an initial
crossover analysis, but rather for doing the inversion to find a set of
suggested adjustments for a given network of crossovers, using the furthest-
first approach (new for GLODAP v3; earlier versions used a different method).

The main function needed is `xover.inversion.furthest_first`; suggested import
as follows:

```python
import xover.inversion as xinv

ff = xinv.furthest_first(xovers, **kwargs)
```

Consult the function docstring for more information about the format of
`xovers` and other possible `kwargs`, plus the contents of the output `ff`.
"""

from .meta import __version__


__all__ = ["__version__"]

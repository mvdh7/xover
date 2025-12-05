# xover: crossover analysis tools for GLODAP

[![Run tests](https://github.com/mvdh7/xover/actions/workflows/run-tests.yml/badge.svg)](https://github.com/mvdh7/xover/actions/workflows/run-tests.yml)
[![PyPI version](https://img.shields.io/pypi/v/xover.svg?style=popout)](https://pypi.org/project/xover/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/xover.svg)](https://anaconda.org/conda-forge/xover)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2600228.svg)](https://doi.org/10.5281/zenodo.2600228)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Installation

    pip install xover
    conda install conda-forge::xover

## Use

The current version of xover does not contain code for running an initial crossover analysis, but rather for doing the inversion to find a set of suggested adjustments for a given network of crossovers, using the furthest-first approach (new for GLODAP v3; earlier versions used a different method).

The main function needed is `xover.inversion.furthest_first`; suggested import as follows:

```python
import xover.inversion as xinv

ff = xinv.furthest_first(xovers, **kwargs)
```

Consult the function docstring for more information about the format of `xovers` and other possible `kwargs`, plus the contents of the output `ff`.

## Version 1 (MATLAB)

Version 1 of xover contained some tools for performing the initial crossover analysis in MATLAB.  The code can be found in the [matlab](https://github.com/mvdh7/xover/tree/matlab) branch, on the [v1.0 release](https://github.com/mvdh7/xover/releases/tag/v1.0) and/or on Zenodo with the DOI [10.5281/zenodo.2600229](https://doi.org/10.5281/zenodo.2600229).

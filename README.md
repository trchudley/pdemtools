# pDEMtools

__Conveniently search, download, and process ArcticDEM and REMA products__


[![conda-forge version](https://anaconda.org/conda-forge/pdemtools/badges/version.svg)](https://anaconda.org/conda-forge/pdemtools) [![PyPI version](https://badge.fury.io/py/pdemtools.svg)](https://pypi.org/project/pdemtools/) [![Documentation Status](https://readthedocs.org/projects/pdemtools/badge/?version=latest)](https://pdemtools.readthedocs.io/en/latest/?badge=latest) [![Unit Tests](https://github.com/trchudley/pdemtools/actions/workflows/unit_test.yml/badge.svg)](https://github.com/trchudley/pdemtools/actions/workflows/unit_test.yml) [![JOSS paper](https://joss.theoj.org/papers/10.21105/joss.07149/status.svg)](https://doi.org/10.21105/joss.07149)

pDEMtools provides a convenient set of functions to explore, download, and preprocess high-resolution DEMs of the polar regions from the ArcticDEM (Porter _et al._  2022; 2023) and Reference Elevation Model of Antarctica (REMA; Howat _et al._ 2022a, b) products, courtesy of the Polar Geospatial Center (PGC).

The first aim of pDEMtools is to enable access to ArcticDEM and REMA mosaics and multitemporal strips using the `search()` function and `load` module:

 - **`search()`**: This function aims to allow users to easily query the PGC STAC API to find relevant ArcticDEM and REMA strips for their areas of interest. 
 - **`load`**: This module provides simple one-line functions to preview and download strips and mosaics from the relevant AWS bucket to an `xarray` Dataset.

The second aim is to provide (pre)processing functions _specific_ to the sort of uses that ArcticDEM and REMA users might want (e.g. a focus on ice sheet and cryosphere work), as well as the particular _strengths_ of ArcticDEM and REMA datasets (high-resolution and multitemporal). Tools include:

 - Terrain attribute derivation (hillshade, slope, aspect, various curvatures) using a 5x5 polynomial fit suited for high-resolution data.
 - Quick geoid correction using BedMachine source data.
 - Simple coregistration for quick elevation change analysis.
 - Identifying/masking sea level and icebergs.

Rather than introducing custom classes, pDEMtools will always try and return DEM data as an [`xarray`](https://docs.xarray.dev/en/stable/) DataArray with geospatial metadata via the [`rioxarray`](https://corteva.github.io/rioxarray/stable/) extension. The aim is to allow the user to quickly move beyond pDEMtools into their own analysis in whatever format they desire, be that `xarray`, `numpy` or `dask` datasets, DEM-specific Python packages such as [`xdem`](https://github.com/GlacioHack/xdem) for advanced coregistration or [`richdem`](https://github.com/r-barnes/richdem) for flow analysis, or exporting to geospatial file formats for analysis beyond Python.

Contact: thomas.r.chudley@durham.ac.uk

## Quick Install

The latest release of pdemtools can installed using `conda`:

```
$ conda install pdemtools -c conda-forge
```

Please visit the [pDEMtools readthedocs](https://pdemtools.readthedocs.io/) for more information on installing, using, and contributing to pDEMtools.

## Cite

A software paper for `pdemtools` is published in the [Journal of Open Source Software](https://joss.theoj.org/), and can be cited as follows:

> Chudley, T. R., and Howat, I. M. (2024). pDEMtools: conveniently search, download, and process ArcticDEM and REMA products. _Journal of Open Source Software_, _9_(102), 7149, doi.org/10.21105/joss.07149

or by using `bibtex`:

```
@article{Chudley2024,
  title = {pDEMtools: conveniently search,  download,  and process ArcticDEM and REMA products},
  volume = {9},
  ISSN = {2475-9066},
  url = {http://dx.doi.org/10.21105/joss.07149},
  DOI = {10.21105/joss.07149},
  number = {102},
  journal = {Journal of Open Source Software},
  publisher = {The Open Journal},
  author = {Chudley,  Thomas R. and Howat,  Ian M.},
  year = {2024},
  pages = {7149}
}
```

When using ArcticDEM and REMA products, please [cite](#refererences) the datasets appropriately and [acknowledge](#acknowledgements) the PGC.

Several algorithms implemented in the library were developed by others. These will be highlighted in the documentation, and the original authors should be properly cited when used. For example:

> We masked sea ice and melange following the method of Shiggins _et al._ (2023) as implemented in pDEMtools (Chudley and Howat, 2024).

<!-- # To do

The tool is presented _as-is_, but requests/contributions to functionality are welcome (thomas.r.chudley@durham.ac.uk). Avenues for future work include the following:

 - Quicker preview downloads of hillshades and DEMs through use of the GeoTIFF overviews and the `rxr.open_rasterio()` `overview_level` function. This can result in uneven x/y resolutions though, so perhaps another option for upsampling may be useful as an accessor utility.
 - Implement Ian Howat's blunder filter algorithm. -->



# Refererences

Howat, I., _et al._ (2022a). The Reference Elevation Model of Antarctica – Strips, Version 4.1. _Harvard Dataverse_ https://doi.org/10.7910/DVN/X7NDNY

Howat, I., _et al._ (2022b). The Reference Elevation Model of Antarctica – Mosaics, Version 2, _Harvard Dataverse_ https://doi.org/10.7910/DVN/EBW8UC

Porter, C., _et al._ (2022). ArcticDEM - Strips, Version 4.1. _Harvard Dataverse_. https://doi.org/10.7910/DVN/OHHUKH

Porter, C., _et al._ (2023), ArcticDEM, Version 4.1, _Harvard Dataverse_. https://doi.org/10.7910/DVN/3VDC4W


# Acknowledgements

**ArcticDEM:** DEMs are provided by the Polar Geospatial Center under NSF-OPP awards 1043681, 1559691, and 1542736.

**REMA:** DEMs are provided by the Byrd Polar and Climate Research Center and the Polar Geospatial Center under NSF-OPP awards 1543501, 1810976, 1542736, 1559691, 1043681, 1541332, 0753663, 1548562, 1238993 and NASA award NNX10AN61G. Computer time provided through a Blue Waters Innovation Initiative. DEMs produced using data from Maxar.


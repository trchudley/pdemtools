.. pdemtools documentation master file, created by
   sphinx-quickstart on Fri Feb 16 10:00:00 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pDEMtools
=========

**Conveniently search, download, and preprocess ArcticDEM and REMA products.**

.. image:: _static/arcticdem_header.jpg
  :width: 800
  :alt: A hillshaded DEM of Helheim Glacier

pDEMtool aims to provide a convenient set of functions to explore, download, and preprocess high-resolution DEMs of the polar regions from the ArcticDEM (Porter _et al._  2022; 2023) and Reference Elevation Model of Antarctica (REMA; Howat _et al._ 2022a, b) products, courtesy of the Polar Geospatial Center (PGC).

The primary aim of pDEMtools is to enable access to ArcticDEM and REMA mosaics and multitemporal strips using the ``search()`` function and ``load`` module:

 - ``search()``: This function aims to replicate the kind of convenient catalogue searching available when querying a dynamic STAC catalogue (e.g. ``pystac_client``), allowing users to easily find relevant ArcticDEM and REMA strips for their areas of interest. 
 - ``load``: This module provides simple one-line functions to preview and download strips and mosaics from the relevant AWS bucket to an ``xarray`` Dataset.

pDEMtools does not aim to provide a comprehensive DEM analysis environment (there are other excellent tools for that, such as `xdem <https://github.com/GlacioHack/xdem>`_ for advanced coregistration or `richdem <https://github.com/r-barnes/richdem>`_ for flow analysis), but some simple (pre)processing functions are provided. It is hoped that they are *specific* to the sort of uses that ArcticDEM and REMA users might want (e.g. a focus on ice sheet and cryosphere work), as well as the particular *strengths* of ArcticDEM and REMA datasets (high-resolution and multitemporal). Tools include:

 - Identifying/masking sea level and icebergs.
 - Terrain attribute derivation (hillshade, slope, aspect, various curvatures) using a 5x5 polynomial fit suited for high-resolution data.
 - Quick geoid correction using BedMachine source data.
 - Simple coregistration for quick elevation change analysis.

Rather than introducing custom classes, pDEMtools will always try and return DEM data as an `xarray <https://docs.xarray.dev/en/stable/>`_ DataArray with geospatial metadata via the `rioxarray <https://corteva.github.io/rioxarray/stable/>`_ extension. The aim is to allow the user to quickly move beyond pDEMtools into their own analysis in whatever format they desire, be that `xarray` datasets, `numpy` or `dask` arrays, or exporting to geospatial file formats for analysis beyond Python.

**Contact me:** Tom Chudley, thomas.r.chudley@durham.ac.uk

.. toctree::
   :maxdepth: 2
   :caption: Getting Started:

   index
   getting_started/install.md
   getting_started/supplementary_datasets.md
   getting_started/cite.md

.. toctree::
   :maxdepth: 2
   :caption: Examples:

   index
   notebooks/mosaic_and_terrain
   notebooks/strip_search_and_dem_difference
   notebooks/get_icebergs

.. toctree::
   :maxdepth: 3
   :glob:
   :caption: API reference:

   index
   api/search.rst
   api/load.rst
   api/data.rst
   api/pdt_accessor.rst

.. toctree::
   :maxdepth: 2
   :glob:
   :caption: Appendix:

   index
   appendix/version_updates.md
   appendix/references.md
   appendix/acknowledgements.md


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

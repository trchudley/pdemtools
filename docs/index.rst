.. pdemtools documentation master file, created by
   sphinx-quickstart on Fri Feb 16 10:00:00 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |version badge| image:: https://anaconda.org/conda-forge/pdemtools/badges/version.svg
    :alt: conda-forge package link and version
    :target: https://anaconda.org/conda-forge/pdemtools

.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.07149/status.svg
    :alt: JOSS publication link and DOI
    :target: https://doi.org/10.21105/joss.07149

pDEMtools     |version badge|  |JOSS|
=====================================

.. .. image:: _static/pdemtools_logo_1600px.png
..    :align: center
..    :alt: pDEMtools logo
..    :width: 1600
..    :height: 349

**Conveniently search, download, and process ArcticDEM and REMA products.**

.. figure:: _static/arcticdem_header.jpg
  :width: 800
  :align: center
  :alt: A hillshaded DEM of Helheim Glacier

``pdemtools`` provides a convenient set of functions to explore, download, and preprocess high-resolution DEMs of the polar regions from the ArcticDEM (Porter *et al.*  2022; 2023) and Reference Elevation Model of Antarctica (REMA; Howat *et al.* 2022a, b) products, courtesy of the Polar Geospatial Center (PGC).

The first aim of ``pdemtools`` is to enable access to ArcticDEM and REMA mosaics and multitemporal strips using the ``search()`` function and ``load`` module:

 - ``search()``: This function aims to allow users to easily query the PGC STAC API to find relevant ArcticDEM and REMA strips for their areas of interest. 
 - ``load``: This module provides simple one-line functions to preview and download strips and mosaics from the relevant AWS bucket to an ``xarray`` Dataset.

The second aim is to provide (pre)processing functions *specific* to the sort of uses that ArcticDEM and REMA users might want (e.g. a focus on ice sheet and cryosphere work), as well as the particular *strengths* of ArcticDEM and REMA datasets (high-resolution and multitemporal). Tools include:

 - Terrain attribute derivation (hillshade, slope, aspect, various curvatures) using a 5x5 polynomial fit suited for high-resolution data.
 - Quick geoid correction using BedMachine source data.
 - Simple coregistration for quick elevation change analysis.
 - Identifying/masking sea level and icebergs.

Rather than introducing custom classes, pDEMtools will always try and return DEM data as an `xarray <https://docs.xarray.dev/en/stable/>`_ DataArray with geospatial metadata via the `rioxarray <https://corteva.github.io/rioxarray/stable/>`_ extension. The aim is to allow the user to quickly move beyond pDEMtools into their own analysis in whatever format they desire, be that `xarray`, `numpy` or `dask` datasets, DEM-specific Python packages such as `xdem <https://github.com/GlacioHack/xdem>`_ for advanced coregistration or `richdem <https://github.com/r-barnes/richdem>`_ for flow analysis, or exporting to geospatial file formats for analysis beyond Python.

**Cite pDEMtools:** 

   Chudley, T. R., and Howat, I. M. (2024). pDEMtools: conveniently search, download, and process ArcticDEM and REMA products. *Journal of Open Source Software*, 9(102), 7149, doi.org/10.21105/joss.07149

.. toctree::
   :maxdepth: 2
   :hidden: 
   :caption: Getting Started:

   getting_started/why_pdemtools.md
   getting_started/install.md
   getting_started/supplementary_datasets.md
   getting_started/cite.md

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Examples:

   examples/mosaic
   examples/strip_search
   examples/terrain_parameters
   examples/coregistration
   examples/get_icebergs
   examples/batch_processing.md

.. toctree::
   :maxdepth: 3
   :hidden:
   :caption: API reference:

   api/search.rst
   api/load.rst
   api/data.rst
   api/pdt_accessor.rst

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Appendix:

   appendix/faq.md
   appendix/community_guidelines.md
   appendix/version_updates.md
   appendix/references.md
   appendix/acknowledgements.md

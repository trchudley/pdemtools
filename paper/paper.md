---
title: 'pDEMtools: conveniently search, download, and process ArcticDEM and REMA products'
tags:
  - Python
  - digital-elevation-models
  - remote-sensing
  - geomorphometry
  - glaciology
  - Greenland
  - Antarctica
  - arctic
authors:
  - name: Thomas R. Chudley
    orcid: 0000-0001-8547-1132
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Ian M. Howat
    orcid: 0000-0002-8072-6260
    affiliation: "2, 3"
affiliations:
 - name: Department of Geography, Durham University, Durham, UK
   index: 1
 - name:  Byrd Polar and Climate Research Center, Ohio State University, Columbus, OH, USA
   index: 2
 - name: School of Earth Sciences, Ohio State University, Columbus, OH, USA
   index: 3

date: 10 May 2024
bibliography: paper.bib

---


# Summary

`pdemtools` is a Python package designed for accessing, processing, and handling high-resolution Digital Elevation Models (DEMs) of the polar regions from the ArcticDEM and Reference Elevation Model of Antarctica (REMA) projects. Tools are provided to search, filter, and download ArcticDEM and REMA data, as well as to fulfill common preprocessing requirements such as geoid correction, coregistration, and calculating terrain attributes. The aim is to allow users to quickly move beyond basic DEM data management into their own analyses.


# Statement of need

[ArcticDEM](https://www.pgc.umn.edu/data/arcticdem/) and [REMA](https://www.pgc.umn.edu/data/rema/) are high-resolution, time-stamped 2-metre-resolution DEMs of the polar regions provided by the Polar Geospatial Center (PGC). They are extracted by applying stereo auto-correlation techniques [@noh_surface_2017] to pairs of submetre Maxar satellite imagery. The data includes Worldview-1, Worldview-2, Worldview-3, and GeoEye-1, beginning in 2007 (ArcticDEM) or 2009 (REMA) and ongoing to the present day. Products are available as tens of thousands of time-stamped 'strips' [@porter_arcticdem_2022; @howat_remastrips_2022] constructed from individual scene pairs, or as a single mosaic [@porter_arcticdem_2023; @howat_remamosaic_2022] compiled from the combined stack of strips. Strips allow users to perform change detection by comparing data from different seasons or years, whilst mosaics provide a consistent and comprehensive product over the entire polar regions.

As Earth Science has moved into the 'big data' era, increasing amounts of Arctic- and Antarctic-focused resources are available as public, cloud-optimised datasets. New approaches are providing Python tools to act as combined API and processing tools, such as `icepyx` [@scheick_icepyx_2023] or `pypromice` [@how_pypromice_2023]. From 2022 (ArcticDEM v4.1 and REMA v2), the PGC DEM products are [hosted](https://polargeospatialcenter.github.io/stac-browser/#/external/pgc-opendata-dems.s3.us-west-2.amazonaws.com/pgc-data-stac.json) as Cloud Optimised GeoTIFFs (CoGs) in a SpatioTemporal Asset Catalog (STAC), a standardised structure for cataloguing spatiotemporal data. However, the PGC STAC is currently a 'static' rather than 'dynamic' STAC, which means there is no convenient Application Programming Interface (API) for searching the datasets in response to user queries. This limits the ability of users to programmatically interact with ArcticDEM and REMA data in a quick and efficient manner. The `pdemtools` package has two aims: the first is to provide a Python-focussed alternative for searching and downloading ArcticDEM and REMA data, emulating dynamic STAC query tools such as `pystac` [@radiant_pystac_2024]; whilst the second is to provide commonly used processing functions specific to the needs of ArcticDEM and REMA users (a focus on ice sheet and cryosphere work), as well as the particular strengths of the ArcticDEM and REMA datasets (high-resolution and multitemporal).

The `pdemtools` `search()` tool and `load` module allow for convenient access to the ArcticDEM and REMA datasets. Mosaics can be downloaded from a one-line `load.mosaic()` function, whilst the `search()` function allows for convenient filtering of a locally downloading ArcticDEM/REMA strip index according to variables such as date, region of interest, spatial coverage, temporal baseline, source sensors, accuracy, and cross-track data. The results of searches are returned as a `geopandas` dataframe [@jordahl_geopandas_2024], and can be downloaded using the `load.from_search()` function. Elevation models are returned as `xarray` DataArrays [@hoyer_xarray_2017] with geospatial metadata via the `rioxarray` extension [@corteva_rioxarray_2024] - a standard format for storing and processing n-dimensional geospatial data within the geospatial Python community. By utilising standardised formats, the aim is to allow the user to quickly move beyond `pdemtools` into their own analysis in whatever format they desire, be that `xarray`, `numpy` or `dask` datasets, DEM analysis Python packages such as `xdem` [@glaciohack_xdem_2020] for advanced coregistration or `richdem` [@barnes_richdem_2016] for flow analysis, or exporting to geospatial file formats for analysis beyond Python.

After download, there exist a number of (pre-)processing steps that are near universally common in topographic analyses. These include geoid-correction, co-registration of time-series data, and/or the construction of terrain parameters such as hillshade, slope, aspect, and curvature. `pdemtools` contains pre-built functions to perform these processing steps, as well as further functionality specific to ArcticDEM and REMA use cases. For instance, we include functions to quickly extract the EIGEN-6C4 geoid [@foerste_geoid_2014] and Greenland/Antarctic bedrock masks directly from local versions of the Greenland and Antarctic BedMachine datasets [@morlighem_icebridge_2020; @morlighem_measures_2020], reprojecting and resampling the data to match the target DEM. Options for ingesting user-provided mask and geoid data are also provided. Additionally, partial derivatives of the surface used to calculate terrain parameters ($\frac{\partial z}{\partial x}$, $\frac{\partial z}{\partial y}$, $\frac{\partial^2 z}{\partial x^2}$, $\frac{\partial^2 z}{\partial y^2}$, $\frac{\partial^2 z}{\partial x \partial y}$) are calculated following @florinsky_computation_2009, as opposed to more common methods such as @zevenbergen_quantitative_1987. The newer approach computes partial derivatives of elevation based on fitting a third-order polynomial, by the least-squares approach, to a 5 $\times$ 5 window as opposed to the more common 3 $\times$ 3 window. This is more appropriate for high-resolution DEMs: curvature over a 10 m window for the 2 m resolution ArcticDEM/REMA strips will lead to a local denoising effect that limits the impact of noise common in high-resolution photogrammetric products. These methods are also adapted into a co-registration routine, which otherwise follows the commonly used approach of @nuth_coregistration_2011.

We aim to grow `pdemtools` by implementing new methods developed by the ArcticDEM and REMA research community. For instance, we currently include sea-level-filtering and iceberg detection routines outlined by @shiggins_automated_2023, and invite community contributions or requests of other routines that will be of use to users of `pdemtools`. Ongoing research projects making use of `pdemtools` are applying ArcticDEM and REMA data to the mapping of crevasses, ice cliff heights, and subglacial lakes, as well as the initiation of ice sheet models. It has also been used within training exercises at the 2024 Polar Geospatial Center Data Workshop, contributing to a growing international network of `pdemtools` users.


# Documentation

Package documentation is available at the `pdemtools` [readthedocs](https://pdemtools.readthedocs.io). 


# Acknowledgements

Chudley was supported by a Leverhulme Early Career Fellowship (ECF-2022-589). ArcticDEM data were provided by the Polar Geospatial Center under NSF-OPP awards 1043681, 1559691, and 1542736. REMA data were provided by the Byrd Polar and Climate Research Center and the Polar Geospatial Center under NSF-OPP awards 1543501, 1810976, 1542736, 1559691, 1043681, 1541332, 0753663, 1548562, 1238993 and NASA award NNX10AN61G. Computer time was provided through a Blue Waters Innovation Initiative. DEMs were produced using data from Maxar. We are grateful to the PGC for making index strips in the parquet format available from their website; to Connor Shiggins for his input on the iceberg detection function; to Andrew Sole and Adrien Wehrle for identifying and helping to resolve bugs in early versions of `pdemtools`; to the reviewers for their helpful comments; and to the colleagues who helped to test this package in the early phases of development.


# References

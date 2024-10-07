# Supplementary datasets

To make the most of pDEMtools, two supplementary datasets must be available locally.

## ArcticDEM/REMA strip index

The first dataset is the ArcticDEM or REMA strip index made available by the PGC  (in GeoParquet format), used by the `search` function. Strip index GeoParquet files can be downloaded from the PGC ([Greenland](https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/), [Antarctica](https://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/)). 

The appriate files to download are the `ArcticDEM_Strip_Index_s2s041_gpqt.zip` and `REMA_Strip_Index_s2s041_gpkg.zip` - note the __GeoParquet file format__ (ending `_gpqt.zip`), *not* the similarly named GeoPackage (ending `_gpkg.zip`). No mosaic files are necessary - these are included with `pdemtools`. 

Unzip the files before use. They can be placed anywhere in your file system - the filepath, as a string, is provided to the relevant `pdemtools` functions.

> __NOTE:__ These index files are maintained by the PGC, so we cannot guarantee that newly released updates to ArcticDEM/REMA index files will not temporarily break `pdemtools` functions until we can fix and update the software. From version 0.8.3 onwards, these files are known to work up to and including the latest PGC data update (ArcticDEM Strip Release 2023, released Aug 2024).

## BedMachine

The second dataset is the Greenland BedMachine (v5; Morlighem _et al._ 2022a) or Antarctica BedMachine (v3; Morlighem _et al._ 2022b) in netcdf format. This dataset provides the default EIGEN-6C4 geoid and ice/bedrock mask used by the geoid correction and coregistration functions (NB: for applications outside of the ice sheets, it is possible to use your own geoid/bedrock mask). 

The appropriate versions of BedMachine can be downloaded from the NSIDC ([Greenland](https://nsidc.org/data/idbmg4/versions/5), [Antarctica](https://nsidc.org/data/nsidc-0756/versions/3)). Follow these links to the NSIDC website. Once you're there, the most convenient way of downloading is via the 'Data Access Tool'. __Ensure you download the netcdf version__ (ending `*.nc`), not the tif version.

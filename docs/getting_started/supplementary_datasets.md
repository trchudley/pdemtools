# Supplementary datasets

Two supplementary datasets can be relied upon locally to enhance the capability of `pdemtools`.

## BedMachine

Greenland BedMachine (v5; Morlighem _et al._ 2022a) or Antarctica BedMachine (v3; Morlighem _et al._ 2022b) provide the default EIGEN-6C4 geoid and ice/bedrock mask used by the geoid correction and coregistration functions (NB: for applications outside of the ice sheets, it is possible to use your own geoid/bedrock mask). 

The appropriate versions of BedMachine, in netcdf format, can be downloaded from the NSIDC ([Greenland](https://nsidc.org/data/idbmg4/versions/5), [Antarctica](https://nsidc.org/data/nsidc-0756/versions/3)). Follow these links to the NSIDC website. Once you're there, the most convenient way of downloading is via the 'Data Access Tool'. __Ensure you download the netcdf version__ (ending `*.nc`), not the tif version.


## ArcticDEM/REMA strip index

As of `pdemtools` `v1.0.0`, the `pdt.search()` functions search the ArcticDEM and REMA strips by seamlessly quering to the online PGC dynamic STAC API. Prior to this, a local ArcticDEM or REMA strip index was required to be downloaded. The functionality to search a local copy of the index is still available within the `pdt.search()` function, and may be useful for those that do not wish to rely on an online connection to search the PGC catalogue.

ArcticDEM or REMA strip index GeoParquet files can be downloaded from the PGC ([Greenland](https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/), [Antarctica](https://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/)). The appriate files to download are the `ArcticDEM_Strip_Index_s2s041_gpqt.zip` and `REMA_Strip_Index_s2s041_gpkg.zip` - note the __GeoParquet file format__ (ending `_gpqt.zip`), *not* the similarly named GeoPackage (ending `_gpkg.zip`). No mosaic files are necessary - these are included with `pdemtools`. 

Unzip the files before use. They can be placed anywhere in your file system - the filepath, as a string, is provided to the relevant `pdemtools` functions.

> __NOTE:__ These index files are maintained by the PGC, so we cannot guarantee that newly released updates to ArcticDEM/REMA index files will not temporarily break `pdemtools` functions until we can fix and update the software. From version 0.8.3 onwards, these files are known to work up to and including the latest PGC data update (ArcticDEM Strip Release 2023, released Aug 2024).

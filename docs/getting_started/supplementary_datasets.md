# Supplementary datasets

To make the most of pDEMtools, two supplementary datasets must available locally:

## ArcticDEM/REMA strip index

The first is the ArcticDEM or REMA strip index made available by the PGC  (in GeoParquet format), used by the `search` function. Strip index GeoParquet files can be downloaded from the PGC ([Greenland](https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/), [Antarctica](https://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/)). To enable rapid searching, *please download the GeoParquet file format*: these files end in `_gpqt.zip`. Unzip them before use.

## BedMachine

The second is the Greenland BedMachine (v5; Morlighem _et al._ 2022a) or Antarctica BedMachine (v3; Morlighem _et al._ 2022b), which is the default geoid/bedrock mask used by the geoid correction and coregistration functions (NB: for applications outside of the ice sheets, it is possible to use your own geoid/bedrock mask). BedMachine can be downloaded from the NSIDC ([Greenland](https://nsidc.org/data/idbmg4/versions/5), [Antarctica](https://nsidc.org/data/nsidc-0756/versions/3)).

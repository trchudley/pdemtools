# Supplementary datasets

To make the most of pDEMtools, two supplementary datasets must available locally. The first is the ArcticDEM or REMA strip index made available by the PGC, used by the `search` function. The second is the Greenland BedMachine (v5; Morlighem _et al._ 2022a) or Antarctica BedMachine (v3; Morlighem _et al._ 2022b), which is the default geoid/bedrock mask used by the geoid correction and coregistration functions (NB: for applications outside of the ice sheets, functions exist for using your own geoid/bedrock mask).

Users have two options:

**Manual download:** Strip indexes can be downloaded from the PGC ([Greenland](https://www.pgc.umn.edu/data/arcticdem/), [Antarctica](https://www.pgc.umn.edu/data/rema/)) and BedMachine from the NSIDC ([Greenland](https://nsidc.org/data/idbmg4/versions/5), [Antarctica](https://nsidc.org/data/nsidc-0756/versions/3)). It is highly recommended PGC strip directory can converted to `.parquet` using, for example, the GeoPandas `to_parquet()` function (see ['Storing strip index files'](#storing-strip-index-files) section below).

**Download scripts:** For convenience the `supp_data` directory contains scripts to download the relevant files directly. The files will be downloaded into the directory the scripts are run. 

```bash
cd supp_data
python download_bedmachine_antarctica_v3.py
python download_bedmachine_greenland_v5.py
python download_index_REMA.py
python download_index_ArcticDEM.py
```

The BedMachine download scripts are provided by the NSIDC, and require an Earthdata user account and password to be provided. The index download scripts are sometimes blocked by the PGC to prevent scraping, so you may have to revert to downloading them manually.

## Storing strip index files

It is **highly recommended** that you store the strip index files as a `*.feather` or `*.parquet` format. You can export a geopandas GeoDataFrame as these formats using the `to_feather()` or `to_parquet()` options - e.g:

```python
import geopandas as gpd
gdf.read_file('.../ArcticDEM_Strip_Index_s2s041_gpkg.gpkg')
gdf.to_parquet("ArcticDEM_Strip_Index_s2s041.parquet")
```

Both of these formats are column-based storage formats designed for 'big data': they are both more compressed and more efficient to read/filter than shapefiles, geopackages, or the like. There are advantages and disadvantages to both (I use parquet) but it's worth at least using one. My own testing for the strip index files showed that the parquet and feather formats are a quarter of the size of the default PGC shapefile and, on my personal laptop, using either increased the speed of loading using geopandas from nearly two minutes with a shapefile to only a few seconds:

```python
%%time
gdf_shp = geopandas.read_file(fpath_shp, intersects=geometry)

# CPU times: user 1min 41s, sys: 4.49 s, total: 1min 45s
# Wall time: 1min 48s

%%time
gdf_feather = geopandas.read_feather(fpath_feather)
gdf_feather = gdf_feather[gdf_feather.intersects(geometry)]

# CPU times: user 1.98 s, sys: 1.62 s, total: 3.6 s
# Wall time: 3.96 s

%%time
gdf_parquet = geopandas.read_parquet(fpath_parquet)
gdf_parquet = gdf_parquet[gdf_parquet.intersects(geometry)]

# CPU times: user 2.1 s, sys: 1.54 s, total: 3.63 s
# Wall time: 3.89 s
```


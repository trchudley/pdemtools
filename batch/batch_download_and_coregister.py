"""
This script automatically downloads and coregisters a complete history of ArcticDEM or
REMA strips within a given AOI. The strips will be created within a named directory in 
the current working directory.

The script will download the 2 m ArcticDEM or REMA mosaic for coregistration, and then
loop through the strip record, downloading and coregistering against stable ground
identified within the BedMachine mask. If coregistration is succesful, the file will 
be saved ending *_coreg.tif. If coregistration fails (most usually because no stable
bedrock is available in the strip), the filename will not have `_coreg` appended.

Note: this will only work in regions of Greenland and Antarctica where bare rock
is exposed. If you require coregistration in regions outside of Greenland or Antarctica, please 
drop me an email - I would be happy to help.

Tom Chudley | thomas.r.chudley@durham.ac.uk
Durham University

v1 | 2023-11 | Initial script
v2 | 2024-07 | Updated to make *_coreg.tif fname dependent on sucessful coregistration.

"""

import os

import rioxarray as rxr
import numpy as np
import pdemtools as pdt

from glob import glob
from math import isnan

# ----------------------------------------------------------------------------------- #
# EDIT THIS SECTION TO SPECIFY DEPENDENT FILES, THE REGION, AND BOUNDS OF INTEREST
# ----------------------------------------------------------------------------------- #

# Filepath to ArcticDEM parquet file -- download from https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ArcticDEM_Strip_Index_latest_gpqt.zip
index_fpath = ".../ArcticDEM_Strip_Index_s2s041_gpqt.parquet"

# Filepath to BedMachine v5 netcdf file -- download from  https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ArcticDEM_Strip_Index_latest_gpqt.zip
bm_fpath = ".../data/bedmachine_5/BedMachineGreenland-v5.nc"

# The name of your study area - this will be used to name output files and directories
region = "isunguatasermia"

# "arcticdem" or "rema"
dataset = "arcticdem"

# the bounds of your study area, in EPSG:3413 for ArcticDEM or EPSG:3031 for REMA
xmin, ymin, xmax, ymax = -247000, -2500000, -212000, -2487000

# Parameters with which to filter ArcticDEM dataset. Note that you may wish to add
# more or less - feel free to modify the search function at line 99.
dates = "20010101/20231231"
baseline_max_hours = 24
min_aoi_frac = 0.1

# ----------------------------------------------------------------------------------- #
# END OF EDITABLE PARAMETER SECTION
# ----------------------------------------------------------------------------------- #

# define AOI
bounds = (xmin, ymin, xmax, ymax)

print(f"Downloading data for {region}:")

# Create output directory
outdir = f"{region}_data"
if not os.path.exists(outdir):
    os.mkdir(outdir)

# Get reference DEM (ArcticDEM mosaic) if it doesn't already exist
reference_dem_fpath = os.path.join(outdir, f"{region}_arcticdem_mosaic_2m.tif")

if not os.path.exists(reference_dem_fpath):
    print("\nDownloading reference DEM...")

    reference_dem = pdt.load.mosaic(
        dataset=dataset,  # must be `arcticdem` or `rema`
        resolution=2,     # must be 2, 10, or 32
        bounds=bounds,    # (xmin, ymin, xmax, ymax) or shapely geometry
        version="v4.1",   # optional: desired version (defaults to most recent)
    )
    reference_dem.rio.to_raster(
        reference_dem_fpath, compress="ZSTD", predictor=3, zlevel=1
    )

else:
    print("\nLoading reference DEM...")
    reference_dem = pdt.load.from_fpath(
        os.path.join(outdir, f"{region}_arcticdem_mosaic_2m.tif"), bounds=bounds
    )

reference_dem = reference_dem.squeeze()
bedrock_mask = pdt.data.bedrock_mask_from_bedmachine(bm_fpath, reference_dem)

# Search for DEM strips
print("\nSearching for DEM strips...")
gdf = pdt.search(
    index_fpath,
    bounds,
    dates=dates,
    # months=[6, 7, 8, 9],
    # years=[2019],
    baseline_max_hours=baseline_max_hours,
    # sensors=["WV03", "WV02", "WV01"],
    # accuracy=2,
    min_aoi_frac=min_aoi_frac,
)
gdf = gdf.sort_values("acqdate1")

n_strips = len(gdf)

print(f"{n_strips} strips found")

i = 1

print("\nDownloading DEM strips...")
for _, row in gdf.iterrows():
    date = row.acqdate1.date()
    date_str = date.strftime("%Y%m%d")
    dem_id = row.dem_id

    out_fname = os.path.join(outdir, f"{date_str}_{dem_id}")

    # If the file doesn't yet exist, download.
    if len(glob(f'{out_fname}*')) == 0:

        print(f"\nDownloading {i}/{n_strips} {os.path.basename(out_fpath)}...")

        # download DEM
        dem = pdt.load.from_search(row, bounds=bounds, bitmask=True)
        dem.compute()  # rioxarray uses lazy evaluation, so we can force the download using the `.compute()` function.
        
        # pad to full size of AOI
        dem = dem.rio.pad_box(*bounds, constant_values=np.nan)

        # coregister DEM, with return_stats=True.
        dem = dem.pdt.coregister(reference_dem, bedrock_mask, max_horiz_offset=50, return_stats=True)
        
        # return_stats=True, so extract DEM as well as RMSE (which will be NaN if the coregisteration failed)
        rmse = dem[-1]
        dem = dem[0]

        # check whether coreg worked, and construct filename appropriately
        if isnan(rmse) == True:
            out_fpath = out_fname + '.tif'
        else:
            out_fpath = out_fname + '_coreg.tif'

        # Export to geotiff
        dem.rio.to_raster(out_fpath, compress="ZSTD", predictor=3, zlevel=1)
        del dem

    i += 1

print("Finished")

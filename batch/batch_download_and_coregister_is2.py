"""
This script automatically downloads and coregisters a complete history of ArcticDEM or
REMA strips within a given AOI. Coregistration will take place against the contemporaneous
ICESat-2 data, which will be downloaded automatically. The strips will be created within 
a new named directory in working directory within which the script will run.

The script will loop through the strip record, downloading the strip DEM as well as the 
contemporaneous ICESat-2 data, and then coregistering the data. Options to geoid-correct
the data or coregister only against stable ground identified within the Greenland/Antarctic 
BedMachine mask (or a custom bedrock mask if provided as a filepath to suitable shapefile). 
If coregistration is succesful, the file will be saved ending *_coreg.tif, or *_coreg_dz.tif 
if only vertical coregistration could be applied. If coregistration fails (most usually 
there is not overlapping reference data available within the strip), the filename will 
not have `_coreg` appended.

If any further assistance or custom additions are required, please drop me an email. I 
would be happy to help.

Tom Chudley | thomas.r.chudley@durham.ac.uk
Durham University

v1 | 2025-01 | Initial script, based upon the `batch_download_and_coregister_dem.py`
               version 3 script.

"""

import os, json

import pdemtools as pdt
import rioxarray as rxr
import numpy as np

from glob import glob
from math import isnan
from warnings import warn

from pandas import to_datetime

# ----------------------------------------------------------------------------------- #
# EDIT THIS SECTION TO SPECIFY DEPENDENT FILES, THE REGION, AND BOUNDS OF INTEREST
# ----------------------------------------------------------------------------------- #

# The name of your study area - this will be used to name output files and directories
region = "isunguatasermia"

# "arcticdem" or "rema"
dataset = "arcticdem"

# the bounds of your study area, in EPSG:3413 for ArcticDEM or EPSG:3031 for REMA
xmin, ymin, xmax, ymax = -247000, -2500000, -212000, -2487000

# Parameters with which to filter ArcticDEM/REMA dataset. Note that you may wish to
# further refine parameters - feel free to modify the `pdt.search()`` function at line
# 153, in consultation with the pdemtools documentation.
dates = "20010101/20231231"
baseline_max_hours = 24
min_aoi_frac = 0.1

# Filepath to BedMachine v5 netcdf file, if operating in Greenland/Antarctica. This
# will allow for automated extraction of a stable ground mask for coregistration,
# as well as a geoid for geoid correction. Leave blank if operating outside of
# Greenland/Antarctica and/or providing your own mask and geoid data.
#  - Download for Greenland: https://nsidc.org/data/idbmg4/versions/5
#  - Download for Antarctica: https://nsidc.org/data/nsidc-0756/versions/3
bedmachine_fpath = ".../data/bedmachine_5/BedMachineGreenland-v5.nc"

# option to geoid-correct data using the bedmachine geoid, defaults to False
geoid_correct = False

# filepath to custom geoid geotiff. Default is `None`, which will extract from BedMachine
# Greenland/Antarctica instead.
custom_geoid_fpath = None

# option to coregister only against stable ground, defaults to False
mask_coreg = False

# filepath of custom mask of stable ground, if operating outside of Antarctica/Greenland
# This should be provided as a shapefile or equivalent geopandas-readable vector.
# Default is `None`, which will extract from BedMachine Greenland/Antarctica instead.
custom_mask_fpath = None

# this script will check
override_datecheck = False

# ----------------------------------------------------------------------------------- #
# END OF EDITABLE PARAMETER SECTION
# ----------------------------------------------------------------------------------- #

# sanity check date string:
first_date_str = dates.split("/")[0]
first_date = to_datetime(first_date_str, format="%Y%m%d")
cutoff_date = to_datetime("2018-10-04")
if first_date < cutoff_date:
    if override_datecheck:
        warn(
            f"You have searched for data beginning {first_date_str}. No ICESat-2 data "
            "exists prior to 2018-10-04. As `override_datecheck` is set to True, "
            "any data prior to this date will be downloaded (uncorrected) anyway."
        )
    else:
        raise ValueError(
            f"You have searched for data beginning {first_date_str}. No ICESat-2 data "
            "exists prior to 2018-10-04. To override this check and download the (uncorrected) "
            "data anyway, set `override_datecheck = True`"
        )

# define AOI
bounds = (xmin, ymin, xmax, ymax)

print(f"\nDownloading data for {region}:")

# Create output directory
outdir = f"{region}_data"
if not os.path.exists(outdir):
    os.mkdir(outdir)

# Search for DEM strips
print("\nSearching for DEM strips...")
gdf = pdt.search(
    dataset=dataset,
    bounds=bounds,
    dates=dates,
    # months=[6, 7, 8, 9],
    # years=[2019],
    baseline_max_hours=baseline_max_hours,
    # sensors=["WV03", "WV02", "WV01"],
    # accuracy=2,
    min_aoi_frac=min_aoi_frac,
)

n_strips = len(gdf)

print(f"{n_strips} strips found")

i = 1

print("\nDownloading DEM strips...")

first_loop = True

for _, row in gdf.iterrows():
    date = row.pdt_datetime1.date()
    date_str = date.strftime("%Y%m%d")
    dem_id = row.pdt_id

    out_fname = os.path.join(outdir, f"{date_str}_{dem_id}")

    # If the file doesn't yet exist, download.
    if len(glob(f"{out_fname}*")) == 0:

        print(f"\nDownloading {i}/{n_strips} {os.path.basename(out_fname)}...")

        # download DEM
        dem = pdt.load.from_search(row, bounds=bounds, bitmask=True)
        dem.compute()  # rioxarray uses lazy evaluation, so we can force the download using the `.compute()` function.

        # pad to full size of AOI
        dem = dem.rio.pad_box(*bounds, constant_values=np.nan)

        # If in the first loop, download the optional bedrock mask and geoid if requested
        if first_loop:

            # get bedrock mask
            if mask_coreg:
                if custom_mask_fpath is None:
                    bedrock_mask = pdt.data.bedrock_mask_from_bedmachine(
                        bedmachine_fpath, dem
                    )
                else:
                    bedrock_mask = pdt.data.bedrock_mask_from_vector(
                        custom_mask_fpath, dem
                    )
            else:
                bedrock_mask = None

            if geoid_correct:
                # Download and save geoid
                if custom_geoid_fpath is None:
                    geoid = pdt.data.geoid_from_bedmachine(bedmachine_fpath, dem)
                else:
                    geoid = pdt.data.geoid_from_raster(custom_geoid_fpath, dem)
            else:
                geoid = None

            first_loop = False

        # get the IS2 coregistration data
        print("Getting IS2 coregistration data...")
        is2_gdf = pdt.data.icesat2_atl06(dem, date, stable_mask=bedrock_mask)

        # coregister DEM, with return_stats=True.
        dem, metadata = dem.pdt.coregister_is2(
            is2_gdf,
            stable_mask=bedrock_mask,
            max_horiz_offset=50,
            return_stats=True,
        )

        # geoid correct if necessary
        if geoid_correct:
            print("Geoid-correcting...")
            dem = dem.pdt.geoid_correct(geoid)

        # check whether coreg worked, and construct filename appropriately
        if metadata["coreg_status"] == "failed":
            out_fpath = out_fname + ".tif"
        if metadata["coreg_status"] == "coregistered":
            out_fpath = out_fname + "_coreg.tif"
        elif metadata["coreg_status"] == "dz_only":
            out_fpath = out_fname + "_coreg_dz.tif"
        else:
            warn(
                f"Unknown coregistration status {metadata['coreg_status']}. Saving wihtout extension.",
                UserWarning,
                stacklevel=2,
            )
            out_fpath = out_fname + ".tif"

        # Export to geotiff
        dem.rio.to_raster(out_fpath, compress="ZSTD", predictor=3, zlevel=1)
        del dem

        # Export metadata as json
        json_fpath = out_fpath.rsplit(".", 1)[0] + ".json"
        with open(json_fpath, "w") as json_file:
            json.dump(metadata, json_file, indent=4)

    i += 1

print("Finished")

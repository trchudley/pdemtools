"""This module contains functions necessary extracting relevant datasets for DEM
processing (geoids, masks, etc), resampled to match the DEM xarray object.
"""

import datetime
from typing import Optional, Literal
from warnings import warn

import rioxarray as rxr
import geopandas as gpd
import sliderule

from shapely.geometry import box

from rasterio.enums import Resampling
from xarray import DataArray, Dataset
from geopandas.geodataframe import GeoDataFrame
from sliderule import sliderule, icesat2
from pandas import to_datetime

from ._utils import clip, get_resolution

# import os
# from shapely.geometry.polygon import Polygon


def geoid_from_bedmachine(bm_fpath: str, target_rxd: DataArray) -> DataArray:
    """Extracts the BedMachine geoid (EIGEN-6C4), bilinearly resampled to match the
    target dataset.

    :param bm_fpath: Filepath to BedMachine dataset, defaults to None
    :type bm_fpath: str
    :param target_rxd: (rio)xarray dataset that BedMachine will be resampled to match
    :type target_rxd: DataArray

    :returns: geoid for the target_rxd region as an xarray DataArray
    :rtype: DataArray"""

    geoid = rxr.open_rasterio(f"{bm_fpath}")["geoid"]
    geoid_crs = geoid.rio.crs
    geoid = geoid.squeeze().astype("float32").rio.write_crs(geoid_crs)
    geoid = geoid.rio.reproject_match(target_rxd, Resampling.bilinear)

    return geoid.squeeze()


def geoid_from_raster(fpath: str, target_rxd: DataArray = None) -> DataArray:
    """Extracts an arbritary geoid stored as a raster dataset, bilinearly resampled to
    match the target dataset.

    :param fpath: Filepath to geoid raster
    :type fpath: str
    :param target_rxd: (rio)xarray dataset/array that the raster will be resampled to
        match. Optional, defaults to None
    :type target_rxd: DataArray

    :returns: geoid for the target_rxd region as an xarray DataArray
    :rtype: DataArray
    """

    geoid = rxr.open_rasterio(f"{fpath}")
    geoid_crs = geoid.rio.crs
    geoid = geoid.squeeze().astype("float32").rio.write_crs(geoid_crs)

    if target_rxd != None:
        geoid = geoid.rio.reproject_match(target_rxd, Resampling.bilinear)

    return geoid.squeeze()


def bedrock_mask_from_vector(
    vector: str | GeoDataFrame, target_rxd: DataArray
) -> DataArray:
    """Construct boolean bedrock mask from a Geopandas vector file of bedrock areas and
    a given target rioxarray dataset. Returns mask where bedrock values are 1 and
    outside are 0.

    :param vector: either a GeoPandas GeoDataFrame of bedrock areas, or a filepath to
        a Geopandas-readable vector file (e.g. shapefile, geopackage, etc.).
    :type vector_fpath: str | GeoDataFrame
    :param target_rxd: (rio)xarray dataset that BedMachine will be resampled to match
    :type target_rxd: DataArray

    :returns: bedrock mask for the target_rxd region as a (rio)xarray DataArray
    :rtype: DataArray
    """

    if type(vector) == str:
        gdf_clip = gpd.read_file(vector)
    elif type(vector) == GeoDataFrame:
        gdf_clip = vector
    else:
        raise ValueError(
            "Input `vector` must be either a filepath string or GeoPandas GeoDataFrame"
        )
    target_rxd = target_rxd.rio.write_nodata(-9999)  # Enforce -9999 as nodata value
    target_clip = target_rxd.rio.clip(gdf_clip.geometry.values, drop=False)
    return (target_clip.where(target_clip != -9999) * 0 + 1).fillna(0).squeeze()


def bedrock_mask_from_bedmachine(bm_fpath: str, target_rxd: DataArray) -> DataArray:
    """Construct boolean bedrock mask from bedmachine and a given target rioxarray
    dataset. Returns mask where bedrock values are 1 and outside are 0.

    :param bm_fpath: Filepath to BedMachine dataset, defaults to None
    :type bm_fpath: str
    :param target_rxd: (rio)xarray dataset that BedMachine will be resampled to match
    :type target_rxd: DataArray

    :returns: bedrock mask for the target_rxd region as a (rio)xarray DataArray
    :rtype: DataArray
    """

    # Open geoid
    mask = rxr.open_rasterio(f"{bm_fpath}")["mask"]
    mask_crs = mask.rio.crs

    # Get geoid-projected geometry of the extent of the target dataset,
    # with a bit of a buffer for safety
    clip_bounds = (
        gpd.GeoDataFrame(
            geometry=[box(*target_rxd.rio.bounds())], crs=target_rxd.rio.crs
        )
        .to_crs(mask_crs)
        .buffer(1000)
        .total_bounds
    )

    # Clip geoid to smaller extent
    mask = clip(mask, clip_bounds)

    # bedmachine mask sets ice-free land to 1. Set everything else to zero.
    mask = mask.where(mask == 1, other=0)

    # Reample using bilinear resampling rather than nearest neighbour, as mask
    # resolution is much larger than DEM resolution
    mask = mask.squeeze().astype("float32").rio.write_crs(mask_crs)
    mask = mask.rio.reproject_match(
        match_data_array=target_rxd, resampling=Resampling.bilinear
    )

    # Set interpolated values to boolean (1/0) values again
    mask = mask.round()

    return mask.squeeze()


def icesat2_atl06(
    target_rxd: DataArray | tuple,
    date: str | datetime.date,
    stable_mask: Optional[DataArray] = None,
    days_r: Optional[int] = [7, 14, 28, 45],
    min_is2_points: Optional[int] = 100,
    filter_points: Optional[bool] = True,
    include_dt_column: Optional[bool] = True,
    epsg: Optional[int] = None,
    sr_srt: Optional[int] = 3,
    sr_cnf: Optional[int] = 1,
    sr_ats: Optional[float] = 20.0,
    sr_cnt: Optional[int] = 10,
    sr_len: Optional[float] = 20.0,
    sr_res: Optional[float] = 10.0,
    sr_track: Optional[int] = 0,
    sr_sigma_r_max: Optional[int] = 5,
):
    """Get ICESat-2 ATL06 data for a given target rioxarray dataset. Returns
    geopandas dataframe of ATL06-SR, processed and downloaded using the
    `sliderule` package.

    Parameters beginning with `sr_` are the sliderule parameters for
    `sliderule.icesat2.atl06p()`. All parameters are left to default, apart from
    the length `len` (20 m rather than 40 m) and step distance `res` (10 m rather
    than 20 m), in order to increase the resolution of the returned data to
    approximate the 10 m resolution of the default coregsitration DEM resolution.
    More information on parameters  are provided here:
    https://slideruleearth.io/web/rtd/user_guide/ICESat-2.html

    :param target_rxd: (rio)xarray dataset that ICESat-2 will be downloaded for.
        Also accepts a tuple of bounds in format [xmin, ymin, xmax, ymax], if the
        `epsg` parameter is provided.
    :type target_rxd: DataArray
    :param date: Date of dataset, as an interpretable string or datetime.date object.
    :type date: str | datetime.date
    :param stable_mask: (rio)xarray mask dataset where 1=valid, 0=masked. Will be
        used to spatially filter the ICESat-2 data. Defaults to None.
    :type stable_mask: DataArray
    :param days_r: Day range within which to search for ICESat-2 data. Should be
        ascending list. Defaults to [7, 14, 28, 45], meaning that all ICESat-2 data
        within 7 days will be returned, and then, if not returning more points
        than `min_is2_points` (default 100), 28 days, 14 days, etc.
    :type days_r: List[int]
    :param min_is2_points: Minimum number of ICESat-2 points. If there are fewer,
        then the param `days_r` will be increased. If, after the maximum `days_r`
        has been searched and there are still insufficient points, then the
        function will return the available results with a warning. Defaults to 50.
    :type min_is2_points: int
    :param filter_points: If True, then only ICESat-2 points that cover valid
        pixels in the `target_rxd` and `stable_mask` will be returned.
        Defaults to True.
    :type filter_points: bool
    :param include_dt_column: If True, then the datetime column `dt` will be included
        in the returned datafram that shows the temporal offset between the ICESat-2
        point and the provided `date`. Defaults to True.
    :type include_dt_column: bool
    :param epsg: EPSG code for the target_rxd dataset. If `target_rxd` is a tuple,
        this must be provided.
    :type epsg: int
    :param sr_srt: sliderule surface type parameter. 0-land, 1-ocean, 2-seaice,
        3-landice, 4-inlandwater. Defaults to 3.
    :type sr_srt: int
    :param sr_cnf: sliderule confidence threshold.  Default 1 (within 10 m).
    :type sr_cnf: int
    :param sr_ats: sliderule minimum along track spread. Defaults to 20.0 m.
    :type sr_ats: float
    :param sr_cnt: sliderule minimum photon count in segment. Defaults to 10.
    :type sr_cnt: int
    :param sr_len: sliderule length of each extent in metres. Defaults to 20.0.
    :type sr_len: float
    :param sr_res: sliderule step distance for extents. Defaults to 10.0.
    :type sr_res: float
    :param sr_track: sliderule track. 0: all tracks, 1: gt1, 2: gt2, 3: gt3
        Defaults to 0.
    :type sr_track: int
    :param sr_sigma_r_max: sliderule maximum robust dispersion in metres.
        Defaults to 5.
    :type sr_sigma_r_max: int

    :returns: geopandas dataframe of ICESat-2 data
    :rtype: geopandas.DataFrame
    """

    # Sanity check date
    dem_date = to_datetime(date)
    cutoff_date = to_datetime("2018-10-04")
    if dem_date < cutoff_date:
        warn(
            f"You have searched for data beginning {dem_date}. No ICESat-2 data "
            "exists prior to 2018-10-04. This function will return no usable data."
        )

    # Get bounds in right format
    if isinstance(target_rxd, (list, tuple)) and len(target_rxd) == 4:
        if epsg is None:
            raise ValueError(
                "If querying with a tuple of bounds, `epsg` must be defined."
            )
        bounds = target_rxd

    elif isinstance(target_rxd, DataArray):
        bounds = target_rxd.rio.bounds()
        epsg = target_rxd.rio.crs
        # resolution = get_resolution(target_rxd)
    else:
        raise ValueError("`target_rxd` must be xarray object or tuple")

    if isinstance(days_r, int):
        days_r = [days_r]

    # Get search region as shapely geometry in epsg:4326
    gdf_4326 = gpd.GeoDataFrame(geometry=[box(*bounds)], crs=epsg).to_crs(
        4326
    )  # .geometry.values[0]

    # connect to sliderule
    icesat2.init("slideruleearth.io")
    sr_region = sliderule.toregion(gdf_4326)

    # set search parameters (apart from dates)
    params = {
        "poly": sr_region["poly"],
        "srt": sr_srt,  # Surface. 0-land, 1-ocean, 2-seaice, 3-landice (default), 4-inlandwater
        "cnf": sr_cnf,  # Confidence. Default 1 (within 10 m). 2: Low. 3: Medium. 4: High.
        "ats": sr_ats,  # Mininum along track spread. SR Default: 20. pDEMtoold default: 10
        "cnt": sr_cnt,  # Minimum photon count in segment. Default 10.
        "len": sr_len,  # Extent length. ATL06 default is 40 metres.
        "res": sr_res,  # Step distance. ATL06 default is 20 metres.
        "track": sr_track,  # Integer: 0: all tracks, 1: gt1, 2: gt2, 3: gt3. Default 0.
        "sigma_r_max": sr_sigma_r_max,  # Max robust dispersion [m]. Default 5.
    }

    first_iteration = True

    for r in days_r:

        if first_iteration:
            first_iteration = False
        else:
            print(f"Expanding temporal search.")

        date_min = to_datetime(date) - datetime.timedelta(days=r)
        date_max = to_datetime(date) + datetime.timedelta(days=r)
        date_min = date_min.strftime("%Y-%m-%d")
        date_max = date_max.strftime("%Y-%m-%d")

        params["t0"] = date_min
        params["t1"] = date_max

        print(f"Querying points within {r} days... ", end="")

        gdf = icesat2.atl06p(params).to_crs(epsg)

        n_points = len(gdf)

        # Determine number of points that intersect with valid data (according
        # to NaN values and mask)
        if n_points > 0:
            # point coordinates dataset
            coords_ds = Dataset(
                {
                    "x": (["points"], gdf.geometry.x.values),
                    "y": (["points"], gdf.geometry.y.values),
                }
            )

            gdf_flt = gdf.copy()

            # Determine raster value and mask at each point
            gdf_flt["values"] = target_rxd.interp(coords_ds, method="nearest").values
            if stable_mask is not None:
                gdf_flt["mask"] = stable_mask.interp(coords_ds, method="nearest").values

            # Filter values at each point
            gdf_flt = gdf_flt[~gdf_flt["values"].isnull()]
            if stable_mask is not None:
                gdf_flt = gdf_flt[gdf_flt["mask"] > 0]

            # Determine remaining number of points
            n_points_flt = len(gdf_flt)
        else:
            gdf_flt = gdf
            n_points_flt = n_points

        # Print statement
        if n_points == 0:
            print(f"No points found intersecting.", end=" ")
        elif n_points_flt < min_is2_points:
            print(
                f"{n_points} points found with {n_points_flt} intersecting dataset. Below minimum threshold ({min_is2_points}).",
                end=" ",
            )
        else:
            print(f"{n_points} points found with {n_points_flt} intersecting dataset.")
            break

    # Return filtered points if requested
    if n_points > 0:

        if filter_points:
            gdf = gdf_flt

            # Drop filter column values
            gdf.drop(["values"], axis=1, inplace=True)
            if stable_mask is not None:
                gdf.drop(["mask"], axis=1, inplace=True)

        if include_dt_column is True:
            # include 'dt' column
            gdf["request_date_dt"] = gdf.index - to_datetime(date)

    if n_points < min_is2_points:
        warn("Returning GeoDataFrame with number of points below minimum threshold.")
        return gdf
    else:
        return gdf

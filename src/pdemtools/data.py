"""Functions for extracting relevant datasets (geoids, masks, etc) resampled to match
a rioxarray object."""

import os

import rioxarray as rxr
import geopandas as gpd

from shapely.geometry import box
from shapely.geometry.polygon import Polygon
from rasterio.enums import Resampling
from xarray import DataArray

from ._utils import clip


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


def bedrock_mask_from_vector(vector_fpath: str, target_rxd: DataArray) -> DataArray:
    """Construct boolean bedrock mask from a Geopandas-readable vector file (e.g.
    shapefile, geopackage, etc.) of bedrock areas and a given target rioxarray dataset.
    Returns mask where bedrock values are 1 and outside are 0.

    Args:
        vector_fpath (str): file path to a Geopandas-readable vector file (e.g.
            shapefile, geopackage, etc.) of bedrock areas.
        target_rxd (DataArray): _description_

    Returns:
        DataArray: (rio)xarray dataset that BedMachine will be resampled to match
    """

    gdf_clip = gpd.read_file(vector_fpath)
    target_clip = target_rxd.rio.clip(gdf_clip.geometry.values, drop=False)
    return (target_clip.where(target_clip >= 0) * 0 + 1).fillna(0).squeeze()


def bedrock_mask_from_bedmachine(bm_fpath: str, target_rxd: DataArray) -> DataArray:
    """Construct boolean bedrock mask from bedmachine and a given target rioxarray
    dataset. Returns mask where bedrock values are 1 and outside are 0.

    :param bm_fpath: Filepath to BedMachine dataset, defaults to None
    :type bm_fpath: str
    :param target_rxd: (rio)xarray dataset that BedMachine will be resampled to match
    :type target_rxd: DataArray

    :returns: bedrock mask for the target_rxd region as an xarray DataArray
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
    mask = mask.rio.reproject_match(target_rxd, Resampling.bilinear)

    # Set interpolated values to boolean (1/0) values again
    mask = mask.round()

    return mask.squeeze()


def mask_from_geometry(
    target_xds: DataArray,
    geometry: str | Polygon | gpd.GeoDataFrame,
) -> DataArray:
    """Construct boolean mask from a given geometry and target rioxarray dataset.
    Returns mask where values within geometry are 1 and outside are 0.

    :param dataarray: dataarray to construct a mask for.
    :param geometry: geometry (or filepath to geometry).

    :returns mask_dataarray: xarray DataArray of mask.
    """

    # if geometry is a filepath string, (os.path.exists to check), then load a geometry

    # fillna() da object, multiply by one, and then clip to geometry. fillna() again
    # to zero

    # return object

    raise NotImplementedError()

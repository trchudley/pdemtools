"""This module contains functions necessary to open an ArcticDEM/REMA strip as an
xarray DataArray, from either local or AWS sources.
"""

import os
from typing import Optional, Literal, Union
from importlib import resources

import rioxarray as rxr
import geopandas as gpd
import numpy as np

from rioxarray.merge import merge_arrays
from xarray import DataArray
from shapely.geometry.polygon import Polygon
from geopandas.geodataframe import GeoDataFrame

from ._utils import clip

# arctic and rema valid versions for STAC retrival
VERSIONS = {"arcticdem": ["v3.0", "v4.1"], "rema": ["v2.0"]}

# filenames of mosaic indexes in ./src/pdemtools/mosaic_index directory
ARCTICDEM_V3_INDEX_FNAME = "ArcticDEM_Mosaic_Index_v3_gpkg.gpkg"
ARCTICDEM_V3_INDEX_2M_LAYER_NAME = "ArcticDEM_Mosaic_Index_v3_2m"
ARCTICDEM_V4_INDEX_FNAME = "ArcticDEM_Mosaic_Index_v4_1_gpkg.gpkg"
ARCTICDEM_V4_INDEX_2M_LAYER_NAME = "ArcticDEM_Mosaic_Index_v4_1_2m"
REMA_V2_INDEX_FNAME = "REMA_Mosaic_Index_v2_gpkg.gpkg"
REMA_V2_INDEX_2M_LAYER_NAME = "REMA_Mosaic_Index_v2_2m"

# aws location
PREFIX = "https://pgc-opendata-dems.s3.us-west-2.amazonaws.com"

# valid mosaic resolutions
VALID_MOSAIC_RES = ["2m", "10m", "32m"]


def from_fpath(
    dem_fpath: str,
    bounds: Optional[Union[tuple, Polygon]] = None,
    bitmask_fpath: Optional[str] = None,
    pad: Optional[bool] = False,
    chunks: Optional[Union[int, tuple, dict]] = None,
) -> DataArray:
    """Given a filepath (local or an AWS link), loads the desired ArcticDEM/REMA DEM
    strip as an ``xarray`` ``DataArray``. Option to filter to bounds and bitmask, if
    provided.

    If AWS link is provided, strip will be downloaded from the relevant AWS bucket. 2 m
    DEM strips are large in size and loading remotely from AWS may take some time.

    :param dem_fpath: Filepath of DEM strip
    :type dem_fpath: str
    :param bounds: Clip to bounds [xmin, ymin, xmax, ymax], in EPSG:3413 (ArcticDEM) or
        EPSG:3031 (REMA). Will accept a shapely geometry to extract bounds from.
        Defaults to None
    :type bounds: tuple | Polygon, optional
    :param mask_fpath: Path to _bitmask.tif file used to mask the DEM, defaults to None
    :type mask_fpath: str, optional
    :param pad: If the DEM strip is not the full extent of the given bounds,
        pad with NaNs to match the full bounds. Defaults to False.
    :type pad: bool
    :param chunks: Chunk size for `rioxarray.open_rasterio`, triggering `dask`
        parallelisation and lazy evaluation/loading. Defaults to `None`.
    :type chunks: int | tuple | dict

    :returns: xarray DataArray of DEM strip
    :rtype: DataArray
    """

    # Sanitise input
    if not isinstance(pad, bool):
        raise ValueError("pad must be True or False")

    # Open dataarray using rioxarray
    dem = rxr.open_rasterio(dem_fpath, chunks=chunks)

    # Convert shapely geometry to bounds
    if isinstance(bounds, Polygon):
        bounds = bounds.bounds

    # Clip if requested, or get whole bounds if not
    if bounds is not None:
        dem = clip(dem, bounds)
    else:
        bounds = dem.rio.bounds()
        if pad is not False:
            raise ValueError("If `bounds` is not provided, `pad` must be False")

    # Filter -9999.0 values
    dem = dem.where(dem > -9999.0)

    # Mask using bitmask if requested
    if bitmask_fpath is not None:
        mask = rxr.open_rasterio(bitmask_fpath)
        if bounds is not None:
            mask = clip(mask, bounds)
        dem = dem.where(mask == 0)
        del mask

    # Remove `band` dim
    dem = dem.squeeze(drop=True)

    # Enforce CF-compliant names
    dem["x"].attrs["axis"] = "X"
    dem["x"].attrs["long_name"] = "x coordinate of projection"
    dem["x"].attrs["standard_name"] = "projection_x_coordinate"
    dem["x"].attrs["units"] = "metre"

    dem["y"].attrs["axis"] = "Y"
    dem["y"].attrs["long_name"] = "y coordinate of projection"
    dem["y"].attrs["standard_name"] = "projection_y_coordinate"
    dem["y"].attrs["units"] = "metre"

    if pad is True:
        dem = dem.rio.pad_box(*bounds, constant_values=np.nan)

    return dem


def preview(
    row: GeoDataFrame,
    bounds: Optional[Union[tuple, Polygon]] = None,
    bitmask: Optional[bool] = True,
    pad: Optional[bool] = False,
):
    """Loads a 10 m hillshade preview of the desired ArcticDEM/REMA DEM strip as an
    ``xarray`` ``DataArray``, for preliminary plotting and assessment. Option to filter
    to bounds.

    :param row: A selected row from the GeoDataFrame output of pdemtools.search. Can
        either select a row manually using gdf.isel[[i]] where `i` is the desired row,
        or provide the entire GeoDataFrame,
    :type row: GeoDataFrame
    :param bounds: Clip to bounds [xmin, ymin, xmax, ymax], in EPSG:3413 (ArcticDEM) or
        EPSG:3031 (REMA). Will accept a shapely geometry to extract bounds from.
        Defaults to None
    :type bounds: tuple | Polygon, optional
    :param bitmask: If True, loads the masked hillshade according to the PGC bitmask.
        If False, loads the unmasked hillshade. Defaults to True
    :type bitmask: bool
    :param pad: If the DEM strip is not the full extent of the given bounds,
        pad with NaNs to match the full bounds. Defaults to True.
    :type pad: bool

    :returns: xarray DataArray of DEM strip 10 m hillshade
    :rtype: DataArray
    """

    if bitmask is True:
        try:
            hillshade_10m_url = row.href_hillshade_masked.values[0]
        except:
            hillshade_10m_url = row.href_hillshade_masked
    elif bitmask is False:
        try:
            hillshade_10m_url = row.href_hillshade.values[0]
        except:
            hillshade_10m_url = row.href_hillshade
    else:
        raise ValueError("masked must be True or False")

    preview = from_fpath(hillshade_10m_url, bounds)

    if bounds is not None and pad is True:
        preview = preview.rio.pad_box(*bounds, constant_values=np.nan)

    return preview.where(preview > 0)


def from_search(
    row: GeoDataFrame,
    bounds: Optional[Union[tuple, Polygon]] = None,
    bitmask: Optional[bool] = True,
    pad: Optional[bool] = False,
    chunks: Optional[Union[int, tuple, dict]] = None,
):
    """Given a row from the GeoDataFrame output of ``pdemtools.search()``, loads the 2
    m DEM strip of the desired ArcticDEM/REMA DEM strip as an xarray DataArray.

    Downloads from the relevant AWS bucket, as an xarray ``DataArray``. 2 m DEM strips
    are large in size and loading remotely from AWS may take some time. Option to
    filter to bounds and bitmask.

    :param row: A selected row from the GeoDataFrame output of pdemtools.search. Select
        the row manually using ``row = gdf.isel[[i]]``, or loop through the rows using
        ``for _, row in gdf.iterrows():``.
    :type row: GeoDataFrame
    :param bounds: Clip to bounds [xmin, ymin, xmax, ymax], in EPSG:3413 (ArcticDEM) or
        EPSG:3031 (REMA). Will accept a ``shapely.Polygon`` geometry to extract bounds
        from. Defaults to None
    :type bounds: tuple | Polygon, optional
    :param bitmask: Choose whether apply the associated bitmask. Defaults to True
    :type bitmask: bool, optional
    :param pad: If the DEM strip is not the full extent of the given bounds,
        pad with NaNs to match the full bounds. Defaults to False.
    :type pad: bool
    :param chunks: Chunk size for `rioxarray.open_rasterio`, triggering `dask`
        parallelisation and lazy evaluation/loading. Defaults to `None`.
    :type chunks: int | tuple | dict

    :returns: xarray DataArray of DEM strip
    :rtype: DataArray
    """

    try:
        dem_url = row.href_dem.values[0]
    except:
        dem_url = row.href_dem

    # try:
    #     s3url = row.s3url.values[0]
    # except:
    #     s3url = row.s3url

    # json_url = "https://" + s3url.split("/external/")[1]
    # dem_url = json_url.replace(".json", "_dem.tif")

    # Construct bitmask fpath, if required
    if bitmask is True:
        try:
            bitmask_url = row.href_mask.values[0]
        except:
            bitmask_url = row.href_mask
    else:
        bitmask_url = None

    # Pass AWS URL locations to load_local command
    return from_fpath(
        dem_url,
        bounds,
        bitmask_url,
        pad,
        chunks,
    )


def from_id(
    dataset: Literal["arcticdem", "rema"],
    geocell: str,
    dem_id: str,
    bounds: Optional[Union[tuple, Polygon]] = None,
    bitmask: Optional[bool] = True,
    bucket: Optional[str] = "https://pgc-opendata-dems.s3.us-west-2.amazonaws.com",
    version: Optional[str] = "s2s041",
    preview: Optional[bool] = False,
    pad: Optional[bool] = False,
    chunks: Optional[Union[int, tuple, dict]] = None,
) -> DataArray:
    """An alternative method of loading the selected ArcticDEM/REMA strip, which
    requires only the geocell and the dem_id (e.g. ``geocell = 'n70w051'``, ``dem_id =
    'SETSM_s2s041_WV01_20200709_102001009A689B00_102001009B63B200_2m_lsf_seg2'``).

    Downloads from the relevant AWS bucket, as an xarray ``DataArray``. Option to
    filter to bounds and bitmask. 2 m DEM strips are large in size and loading
    remotely from AWS may take some time.

    NOTE: This function is kept for vestigial purposes, and may be removed or
    significantly altered in a future release.

    :param dataset: Either 'arcticdem' or 'rema'. Case-insensitive.
    :type dataset: str
    :param geocell: Geographic grouping of ArcticDEM / REMA strip. e.g. 'n70w051'.
    :type geocell: str
    :param dem_id: ArcticDEM/REMA strip ID. e.g.
        'SETSM_s2s041_WV01_20200709_102001009A689B00_102001009B63B200_2m_lsf_seg2'
    :type dem_id: str
    :param bounds: Clip to bounds [xmin, ymin, xmax, ymax], in EPSG:3413 (ArcticDEM) or
        EPSG:3031 (REMA), defaults to None
    :type bounds: tuple, optional
    :param bitmask: Choose whether apply the associated bitmask, defaults to True
    :type bitmask: bool, optional
    :param bucket: AWS buck link, defaults to
        'https://pgc-opendata-dems.s3.us-west-2.amazonaws.com'
    :type bucket: str
    :param version: Version string, defaults to 's2s041'
    :type version: str
    :param preview: Return just a link to the STAC preview page, defaults to False
    :type preview: bool, optional
    :param pad: If the DEM strip is not the full extent of the given bounds,
        pad with NaNs to match the full bounds. Defaults to False.
    :type pad: bool
    :param chunks: Chunk size for `rioxarray.open_rasterio`, triggering `dask`
        parallelisation and lazy evaluation/loading. Defaults to `None`.
    :type chunks: int | tuple | dict

    :return: xarray DataArray of DEM strip
    :rtype: DataArray
    """

    # Sanitise data
    dataset = dataset.lower()
    geocell = geocell.lower()

    if preview is True:
        browser_prefix = "https://polargeospatialcenter.github.io/stac-browser/#/external/pgc-opendata-dems.s3.us-west-2.amazonaws.com"
        preview_fpath = (
            f"{browser_prefix}/{dataset}/strips/{version}/2m/{geocell}/{dem_id}.json"
        )
        return preview_fpath

    # Construct DEM fpath
    dem_fpath = f'{bucket}/{dataset}/"strips"/{version}/2m/{geocell}/{dem_id}_dem.tif'
    # dem_fpath = os.path.join(
    #     bucket, dataset, "strips", version, "2m", geocell, f"{dem_id}_dem.tif"
    # )

    # Construct bitmask fpath, if required
    if bitmask is True:
        bitmask_fpath = (
            f'{bucket}/{dataset}/"strips"/{version}/2m/{geocell}/{dem_id}_bitmask.tif'
        )
        # bitmask_fpath = os.path.join(
        #     bucket, dataset, "strips", version, "2m", geocell, f"{dem_id}_bitmask.tif"
        # )
    else:
        bitmask_fpath = None

    # Pass AWS URL locations to load_local command
    return from_fpath(
        dem_fpath,
        bounds,
        bitmask_fpath,
        pad,
        chunks,
    )


def mosaic(
    dataset: Literal["arcticdem", "rema"],
    resolution: Literal["2m", "10m", "32m"],
    bounds: Union[tuple, Polygon] = None,
    version: Optional[Literal["v2.0", "v3.0", "v4.1"]] = None,
    chunks: Optional[Union[int, tuple, dict]] = None,
):
    """Given a dataset, resolution, and bounding box, downloads the ArcticDEM or REMA
    mosiac from AWS.

    :param dataset: The desired dataset, either 'arcticdem' or 'rema'.
        Case-instensitive.
    :type datasat: str
    :param resolution: The desired mosaic resolution to download - must be either ``2m``,
        ``10m``, or ``32m`` (will also accept ``2``, ``10``, and ``32`` as ``int`` types)
    :type resolutions: str | int
    :param version: Desired ArcticDEM or REMA version. Must be a valid version available
        from the PGC STAC API (e.g. ``v3.0`` or ``v4.1`` for ArcticDEM, or ``v2.0`` for REMA).
    :type version: str
    :param bounds: Clip to bounds [xmin, ymin, xmax, ymax], in EPSG:3413 (ArcticDEM) or
        EPSG:3031 (REMA). Will accept a shapely geometry to extract bounds from.
    :type bounds: tuple | Polygon, optional
    :param chunks: Chunk size for `rioxarray.open_rasterio`, triggering `dask`
        parallelisation and lazy evaluation/loading. Defaults to `None`.
    :type chunks: int | tuple | dict

    :return: xarray DataArray of DEM mosaic
    :rtype: DataArray
    """

    # sanity check that datset and versioning is correct versioning is valid for selected dataset
    dataset = dataset.lower()
    if dataset not in VERSIONS.keys():
        raise ValueError(
            f"Dataset must be one of {VERSIONS.keys}. Currently `{dataset}`."
        )
    if version is None:  # pick the most recent dataset
        version = VERSIONS[dataset][-1]
    else:
        if version not in VERSIONS[dataset]:
            raise ValueError(
                f"Version of {dataset} must be one of {VERSIONS[dataset]}. Currently `{version}`"
            )

    # get resolution as str
    if isinstance(resolution, int):
        resolution = f"{resolution}m"

    # check if valid resolution
    if resolution not in VALID_MOSAIC_RES:
        raise ValueError(
            f"Resolution must be one of {VALID_MOSAIC_RES}. Currently `{resolution}`"
        )

    # Sanitise shapely geometry to bounds tuple
    if isinstance(bounds, Polygon):
        bounds = bounds.bounds

    # get dataset version
    if dataset == "arcticdem" and version == "v3.0":
        layer = ARCTICDEM_V3_INDEX_2M_LAYER_NAME
    elif dataset == "arcticdem" and version == "v4.1":
        layer = ARCTICDEM_V4_INDEX_2M_LAYER_NAME
    elif dataset == "rema" and version == "v2.0":
        layer = REMA_V2_INDEX_2M_LAYER_NAME
    else:
        raise ValueError(
            "Cannot retrive internal index filepath for specified dataset and version."
        )

    # Load tiles that intersect with AOI
    index_fpath = _get_index_fpath(dataset, version=version)
    with open(index_fpath, "rb") as file:
        tiles = gpd.read_file(index_fpath, layer=layer, bbox=bounds, engine="fiona")

        if len(tiles) < 1:
            raise ValueError(
                f"No {dataset} mosaic tiles found to intersect with bounds {aoi}"
            )

        # get aws filepaths from the tiles dataframe
        fpaths = []
        for _, row in tiles.iterrows():
            fpath = _aws_link(
                row, dataset=dataset, version=version, resolution=resolution
            )
            fpaths.append(fpath)

        tiles = None  # release tiles object

    # remove duplicates in 10m and 32m (which load supertiles, not tiles)
    fpaths = list(set(fpaths))

    # load dem(s)
    dems = []
    for fpath in fpaths:
        dem = rxr.open_rasterio(fpath, chunks=chunks).rio.clip_box(*bounds)
        dems.append(dem)

    if len(fpaths) == 1:
        dem = rxr.open_rasterio(fpaths[0], chunks=chunks).rio.clip_box(*bounds)

    # If multiple dems, merge them - NB I don't know whether this breaks lazy
    # evaluation for chunked data
    if len(dems) > 1:
        dem = merge_arrays(dems)
    else:
        dem = dems[0]

    dems = None  # release dem objects for memory management

    # Filter -9999.0 values to np.nan
    dem = dem.where(dem > -9999.0)

    # Remove `band` dim
    dem = dem.squeeze(drop=True)

    # Enforce CF-compliant names
    dem["x"].attrs["axis"] = "X"
    dem["x"].attrs["long_name"] = "x coordinate of projection"
    dem["x"].attrs["standard_name"] = "projection_x_coordinate"
    dem["x"].attrs["units"] = "metre"

    dem["y"].attrs["axis"] = "Y"
    dem["y"].attrs["long_name"] = "y coordinate of projection"
    dem["y"].attrs["standard_name"] = "projection_y_coordinate"
    dem["y"].attrs["units"] = "metre"

    return dem


def _get_index_fpath(
    dataset: Literal["arcticdem", "rema"],
    version: Literal["v2.0", "v3.0", "v4.0"],
):
    """Given ``arcticdem`` or ``rema``, gets the filepath of the package dataset using the
    ``importlib`` library. ARCTICDEM and REMA global variables necessary.
    """

    # get dataset version
    if dataset == "arcticdem" and version == "v3.0":
        fname = ARCTICDEM_V3_INDEX_FNAME
    elif dataset == "arcticdem" and version == "v4.1":
        fname = ARCTICDEM_V4_INDEX_FNAME
    elif dataset == "rema" and version == "v2.0":
        fname = REMA_V2_INDEX_FNAME
    else:
        raise ValueError(
            "Cannot retrive internal index filepath for specified dataset and version."
        )

    return resources.files("pdemtools.mosaic_index").joinpath(fname)


def _aws_link(
    row,
    dataset: Literal["arcticdem", "rema"],
    version: str,
    resolution: Literal["2m", "10m", "32m"],
    prefix: Optional[str] = PREFIX,
):
    """Using inputs from mosaic() function and the AWS location from the global variable
    `PREFIX`, construct the filepath of the relevant ArcticDEM or REMA mosaic tile.
    """
    # Construct appropriate suffix, considering ArcticDEM v3.0's alternate naming scheme
    if dataset == "arcticdem" and version == "v3.0":
        suffix = f"_{resolution}_{version}_reg_dem.tif"
    else:
        suffix = f"_{resolution}_{version}_dem.tif"

    # Construct appropriate filename given resolution
    if resolution == "2m":
        fname = f"{row.tile}{suffix}"
    elif resolution in ["10m", "32m"]:
        fname = f"{row.supertile}{suffix}"
    else:
        raise ValueError(f"Input `resolution` must be one of ['2m', '10m', '32m']")

    # Return appropriate filepath.
    return f"{prefix}/{dataset}/mosaics/{version}/{resolution}/{row.supertile}/{fname}"

"""
This module contains simple geospatial utility functions used across the other modules.
"""

import rioxarray as rxr

from xarray import DataArray

def get_resolution(xds: DataArray) -> float:
    """Retrive resolution from rioxarray Dataset. Must be the same in `x` and `y`
    direction.

    :param dem_xds: xarray Dataset of DEM strip
    :type dem_xds: Dataset

    :returns: Resolution of DEM strip
    :rtype: float
    """

    resolutions = xds.rio.resolution()
    if abs(resolutions[0]) != abs(resolutions[1]):
        raise ValueError("Dataset resolution not the same in `x` and `y` directions")
    else:
        return abs(resolutions[0])


def geospatial_match(
    rxd_1: DataArray, rxd_2: DataArray, return_info: bool = False
) -> bool:
    """Check whether two (rio)xarray DataArrays or Datasets have the same geospatial
    information.

    :param rxd_1: First (rio)xarray Dataset to compare
    :type rxd_1: DataArray
    :param rxd_2: Second (rio)xarray Dataset to compare
    :type rxd_2: DataArray

    :returns: True if Datasets match, False if not
    :rtype: bool
    """

    failed = []

    if not rxd_1.rio.shape == rxd_2.rio.shape:
        failed.append("shape")
    if not rxd_1.rio.resolution() == rxd_2.rio.resolution():
        failed.append("resolution")
    if not rxd_1.rio.bounds() == rxd_2.rio.bounds():
        failed.append("bounds")
    if not rxd_1.rio.crs == rxd_2.rio.crs:
        failed.append("crs")

    if len(failed) == 0:
        return True
    else:
        if return_info == True:
            return failed
        else:
            return False


def clip(
    rxd: DataArray,
    bounds: tuple,
) -> DataArray:
    """Clips (rio)xarray DataArray to bounds in format [xmin, ymin, xmax, ymax] using
    the rioxarray clip_box() function. Provide bounds in matching CRS (e.g. EPSG:3413
    for ArcticDEM strips and EPSG:3031 for REMA strips).

    :param rxd: (rio)xarray object to clip
    :type rxd: DataArray
    :param bounds: List of AOI bounds in format [xmin, ymin, xmax, ymax]
    :type bounds: tuple

    :returns: Clipped xarray Dataset
    :rtype: Dataset
    """

    # TODO: More informative error handling of intersection issues than rioxarray
    # currently provides

    return rxd.rio.clip_box(*bounds)

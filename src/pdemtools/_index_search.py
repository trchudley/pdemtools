"""Function to search ArcticDEM and REMA index files.
"""

import os

from typing import Optional

import pandas as pd
import geopandas as gpd

from shapely.geometry import box, polygon

SENSORS = ["WV03", "WV02", "WV01", "GE01"]
MONTHS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]


def search(
    index_fpath: str,
    bounds: Optional[tuple | polygon.Polygon] = None,
    min_aoi_frac: Optional[float] = None,
    dates: Optional[str | tuple] = None,
    months: Optional[int | tuple] = None,
    years: Optional[int | tuple] = None,
    baseline_max_hours: Optional[int] = None,
    sensors: Optional[str | tuple] = None,
    is_xtrack: Optional[bool] = None,
    accuracy: Optional[float | tuple] = None,
):
    """Efficiently search the ArcticDEM and REMA strip index files provided by the
    Polar Geospatial Center.

    :param index_fpath: Filepath to a local copy of the ArcticDEM or REMA strip index
        file, available from the Polar Geospatial Center. For speed, it is tecommended
        that the index is `.parquet` format.
    :type index_fpath: str
    :param bounds: Filter to strips that intersect with bounds [xmin, ymin, xmax, ymax],
        in EPSG:3413 (ArcticDEM) or EPSG:3031 (REMA). Accepts a tuple or a
        ``shapely.Polygon`` geometry to extract bounds from.
    :type bounds: tuple | Polygon, optional
    :param min_aoi_frac: Filter to strips than cover more the defined fraction of the
        area of interest (defined by the `bounds` variable). Must be between 0 and 1,
        and ``bounds`` must be defined.
    :type min_aoi_frac: float, optional
    :param dates: Filter strips to a date range. Dates can be provided as a tuple of two
        strings, or a single string with a `/` seperator. Date strings must be
        interpetable by the pandas.to_datetime() tool.
    :type dates: str | tuple, optional
    :param months: Filter strips to only certain months. Provide as a tuple of integers
        (e.g. for June, July, August strips only, set ``months = [6,7,8]``).
    :type months: tuple, optional
    :param years: Filter strips to only certain yeara. Provide as a tuple of integers
        (e.g. for 2011 and 2021  only, set `years = [2011,2021]`).
    :type years: tuple, optional
    :param baseline_max_hours: Filter strips to only those constructed from stereopairs
        acquired less than the provided number of hours apart.
    :type baseline_max_hours: int, optional
    :param sensors: Filter scenes to only those consrtructed from the provided
        satellite sensors. Full list is ["WV03", "WV02", "WV01", "GE01"]
    :type sensors: tuple, optional
    :param is_xtrack: Filter based on whether stereopairs are cross-track imagery.
        True = return only cross-track. False = return only non-cross-track.
    :type is_xtrack: bool, optional
    :param accuracy: Filter to strip accuracies based on the provided average height
        accuracy in metres (`rmse` in the strip index). NB that this column
        included NaN values (-1.0) so the option is provided to include only a single
        value as an upper range (e.g. 2), or a tuple of two values in order to include
        a lower bound and filter NaN values (e.g [0, 2]).
    :type accuracy: float | tuple, optional

    :returns: Strip index filtered to desired variables.
    :rtype: GeoDataFrame
    """

    # Determine whether index is for ArcticDEM or REMA (strangely, quite hard to do -
    # current approach is simply to check whether 'arcticdem' or 'rema' strings occur
    # in the filename of the index file)
    if "arcticdem" in os.path.basename(index_fpath).lower():
        epsg = 3413
    elif "rema" in os.path.basename(index_fpath).lower():
        epsg = 3031
    else:
        raise ValueError(
            "Cannot determine whether index file is ArcticDEM or REMA (please ensure 'arcticdem' or 'rema' is in filename)"
        )

    # Sanitise inputs first (prior to loading, to catch mistakes before time is wasted
    # loading index files)

    # Sanitise input: bounds
    if bounds != None:
        if type(bounds) != polygon.Polygon:
            geom = box(*bounds)
        else:
            geom = bounds
        # convert geom to EPSG:4326
        geom_4326 = gpd.GeoDataFrame(geometry=[geom], crs=epsg).to_crs(4326).geometry[0]

    # Sanitise input: min_aoi_frac
    if min_aoi_frac != None:
        if bounds == None:
            raise ValueError("`bounds` variable must be provided to use `min_aoi_frac`")
        if (min_aoi_frac < 0) or (min_aoi_frac > 1):
            raise ValueError(
                f"`min_aoi_frac` must be a value between 0 and 1. Currently {min_aoi_frac}"
            )

    # Sanitise input: dates
    if dates != None:
        if type(dates) == str:
            dates = dates.split("/")

        if len(dates) != 2:
            raise ValueError(
                "Date range must be tuple of two strings/None values, or a string of two dates seperated by `/` seperators"
            )

    # Sanitise input: months
    if months != None:
        if type(months) is int:
            months = [months]

        if any(m not in MONTHS for m in months):
            raise ValueError(f"`months` variables must be in range {MONTHS}")

    # Sanitise input: sensors
    if sensors != None:
        if type(sensors) is str:
            sensors = [sensors]

        if any(s not in SENSORS for s in sensors):
            raise ValueError(f"`sensors` variables must be in list {SENSORS}")

    # Sanitise input: accuracy
    if accuracy != None:
        if (type(accuracy) == int) or (type(accuracy) == float):
            accuracy = [accuracy]
        if len(accuracy) > 2:
            raise ValueError("`accuracy` must be single value or tuple of length 2")

    # Open the index geometry file, according to file type
    _, extension = os.path.splitext(index_fpath)

    if extension == ".parquet":
        # print("reading")
        gdf = gpd.read_parquet(index_fpath)
        if bounds != None:
            gdf = gdf[gdf.intersects(geom_4326)]

    else:
        print(
            "For quicker searching, it is highly recommended to store strip index files as a `.parquet` or `.feather` format."
        )
        if bounds == None:
            gdf = gpd.read_file(index_fpath)
        else:
            gdf = gpd.read_file(index_fpath, intersects=geom_4326)

    # Construct necessary datetime columns if necessary
    if any([dates, months, baseline_max_hours]):
        gdf["time1"] = pd.to_datetime(gdf.acqdate1)
        gdf["time2"] = pd.to_datetime(gdf.acqdate2)
        if baseline_max_hours != None:
            gdf["dem_baseline_hours"] = (
                (gdf.time2 - gdf.time1)
                .values.astype("timedelta64[h]")
                .astype("float32")
            )
            gdf["dem_baseline_hours"] = abs(
                gdf["dem_baseline_hours"].values
            )  # Absolute d_t, in hours
        if (months != None) or (years != None):
            gdf["time_mean"] = gdf.time1 + (gdf.time2 - gdf.time1) / 2
            gdf["year"] = gdf.time_mean.dt.year
            gdf["month"] = gdf.time_mean.dt.month

   # Filter according to date (acqdate)
    if dates != None:
        if dates[0] != None:
            datetime1 = pd.to_datetime(dates[0])
            try:
                gdf = gdf[gdf["time1"] > datetime1]
            except:  # account for UTC-aware datetime storage in newer PGC index files.
                datetime1 = pd.to_datetime(dates[0], utc=True)
                gdf = gdf[gdf["time1"] > datetime1]

        if dates[1] != None:
            datetime2 = pd.to_datetime(dates[1])
            try:
                gdf = gdf[gdf["time2"] < datetime2]
            except:  # account for UTC-aware datetime storage in newer PGC index files.
                datetime2 = pd.to_datetime(dates[1], utc=True)
                gdf = gdf[gdf["time2"] < datetime2]
            # can probably make this more efficient than just try-except - e.g. check if
            # date column is utc-aware and if so, make datetime utc-aware too.

    # Filter to months
    if months != None:
        gdf = gdf[(gdf["month"].isin(months))]

    # Filter to years
    if years != None:
        gdf = gdf[(gdf["year"].isin(years))]

    # Filter to time seperations
    if baseline_max_hours != None:
        gdf = gdf[(gdf["dem_baseline_hours"] <= baseline_max_hours)]

    # Filter to only selected sensors
    if sensors != None:
        gdf = gdf[(gdf["sensor1"].isin(sensors)) & (gdf["sensor2"].isin(sensors))]

    # Filter to crosstrack
    if is_xtrack != None:
        gdf = gdf[(gdf["is_xtrack"] == int(is_xtrack))]

    # Filter to accuracy range
    if accuracy != None:
        # print(accuracy, len(accuracy))
        if len(accuracy) == 1:
            # print(f"accuracy == 1: {len(accuracy)==1}")
            gdf = gdf[(gdf["rmse"] <= accuracy[0])]
        else:
            gdf = gdf[(gdf["rmse"] >= accuracy[0]) & (gdf["rmse"] <= accuracy[1])]

    # Convert from EPSG:4326 to appropriate ArcticDEM/REMA crs:
    gdf = gdf.to_crs(epsg)

    # Clip geometries to bounds of AOI
    if bounds != None:
        gdf = gpd.clip(gdf, geom)

    # Filter to geometry
    if min_aoi_frac != None:
        aoi_area = geom.area
        gdf["aoi_frac"] = gdf.area / aoi_area
        gdf = gdf[gdf["aoi_frac"] > min_aoi_frac]

    # Drop created columns
    gdf = gdf.drop(
        [
            "time1",
            "time2",
            "dem_baseline_hrs",
            "time_mean",
            "month",
            "year",
            "aoi_frac",
        ],
        axis=1,
        errors="ignore",
    )

    return gdf

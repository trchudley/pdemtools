"""Custom 'accessor' for xarray objects, to allow for filtering, geomorphological
functions, and simple coregistration.
"""

from warnings import warn
from typing import Optional, Literal

import datetime

import xarray as xr
import numpy as np

from xarray import DataArray, Dataset
from geopandas import GeoDataFrame
from pandas import to_datetime
from rasterio.enums import Resampling
from cv2 import connectedComponentsWithStats, CV_32S

from ._geomorphometry import (
    p,
    q,
    r,
    s,
    t,
    slope,
    aspect,
    hillshade,
    plan_curvature,
    horizontal_curvature,
    vertical_curvature,
    mean_curvature,
    gaussian_curvature,
    unsphericity_curvature,
    minimal_curvature,
    maximal_curvature,
)
from ._coreg import coregister, shift_dem  # coregisterdems, coregisterpoints
from ._utils import get_resolution, geospatial_match


LIST_METHODS = ["Florinsky", "ZevenbergThorne"]
LIST_ATTRIBUTES = [
    "slope",
    "aspect",
    "hillshade",
    "plan_curvature",
    "horizontal_curvature",
    "vertical_curvature",
    "horizontal_excess",
    "vertical_excess",
    "mean_curvature",
    "gaussian_curvature",
    "unsphericity_curvature",
    "minimal_curvature",
    "maximal_curvature",
]
LIST_REQUIRING_SECOND_ORDER = [
    "plan_curvature",
    "horizontal_curvature",
    "vertical_curvature",
    "horizontal_excess",
    "vertical_excess",
    "mean_curvature",
    "gaussian_curvature",
    "unsphericity_curvature",
    "minimal_curvature",
    "maximal_curvature",
]


@xr.register_dataarray_accessor("pdt")
class DemAccessor:
    """This class provides ``rioaxarray``-compatible ``DataArray`` objects with access
    to custom functions and methods.

    Functions are accessible via the ``pdt`` accessor: for instance, if you have a DEM
    loaded as the variable ``dem``, you are able to access the ``terrain`` function via
    ``dem.pdt.terrain()``.
    """

    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._center = None

    def geoid_correct(self, geoid: DataArray) -> DataArray:
        """Geoid correct the DEM, given a geoid as an xarray DataArray. The function
        will attempt to reproject the geoid to match if it does not already to do.

        :param geoid: A suitable geoid, as a (rio)xarray dataarray

        :returns: geoid-corrected (rio)xarray dataarray
        :rtype: DataArray"""

        if not geospatial_match(self._obj, geoid):
            geoid = geoid.rio.reproject_match(self._obj, Resampling.bilinear)

        return self._obj - geoid

    def coregister_is2(
        self,
        points_df: GeoDataFrame,
        stable_mask: Optional[DataArray] = None,
        return_stats: Optional[bool] = False,
        max_horiz_offset: float = 50,
        rmse_step_thresh: float = -0.001,
        max_iterations: int = 5,
        verbose: Optional[bool] = True,
        resample: Optional[bool] = True,
        resample_res: Optional[float] = 10,
    ) -> DataArray:
        """Coregisters the scene against reference ICESat-2 data based on the Nuth and Kääb (2011)
        method, as implemented by the PGC.

        Original code available here:
        https://github.com/PolarGeospatialCenter/setsm_postprocessing_python/blob/master/lib/scenes2strips.py

        :points_df: GeoDataFrame of ICESat-2 points to use for coregistration, as provided
            by the `data.load.icesat2_atl06()` function.
        :type points_df: GeoDataFrame
        :param stable_mask: Mask dataarray where 1 = stable region to be used for
            coregistration. If none, coregisters based on the entire scene.
        :type stable_mask: xarray.DataArray
        :param return_stats: If true, returns a metadata dictionary alongside
            the coregistered DEM dataarray. The dictionary contains the coregistration
            status (coregistered, dz_only, failed), the transformation parameter,
            the translation error, and the RMS error. Defaults to False
        :type return_stats: bool
        :param max_horiz_offset: maximum horizontal offset, beyond which XY
            coregistration is ignored and only Z values are corrected. Defaults
            to 15.
        :type max_horiz_offset: float
        :param rmse_step_thresh: threshold for the change in RMSE between two consecutive
            iterations. Defaults to -0.001.
        :type rmse_step_thresh: float
        :param max_iterations: maximum number of iterations. Defaults to 5.
        :type max_iterations: int
        :param verbose: If True, print status updates. Defaults to True.
        :type verbose: bool
        :param resample: If True, resample the DEM to a larger resolution. Defaults to True.
        :type resample: bool
        :param resample_res: Resolution to resample the DEM to. Defaults to 10.
        :type resample_res: float
        """

        # 0. SANITY CHECK
        if stable_mask is None:
            print(
                "No `stable_mask` provided: assuming entire scene is suitable reference surface."
            )
        else:
            check_match = geospatial_match(self._obj, stable_mask, return_info=True)
            if check_match is not True:
                raise ValueError(
                    "Input DEM and stable_mask do not share geospatial information: "
                    f"{check_match}. Consider padding (`.rio.pad_box`) or reprojecting "
                    "(`.rio.repoject_match`) stable_mask."
                )

        # 1. PRE-PROCESSING: INPUT DEM

        resolution = get_resolution(self._obj)

        if (resample is True) and (resolution < resample_res):

            print(
                f"Upsampling {resolution} m DEM to {resample_res} m resolution"
                " for coregistration process..."
            )

            dem_in = self._obj.squeeze().rio.reproject(
                self._obj.rio.crs,  # Maintain the original coordinate reference system
                resolution=resample_res,  # Set the target resolution to 10 meters
                resampling=Resampling.average,  # Upscale with average resampling method
            )

            dem_in_vals = dem_in.values

            if stable_mask is not None:

                mask_in = stable_mask.squeeze().rio.reproject(
                    self._obj.rio.crs,
                    resolution=resample_res,
                    resampling=Resampling.nearest,
                )

                mask_in_vals = mask_in.values

            else:
                mask_in_vals = None

            x_in_vals = dem_in.x.values
            y_in_vals = dem_in.y.values
            resolution = resample_res
            resampled = True

        else:

            dem_in = self._obj.squeeze()
            dem_in_vals = dem_in.values

            if stable_mask is not None:
                mask_in_vals = stable_mask.squeeze().values
            else:
                mask_in_vals = None

            x_in_vals = dem_in.x.values
            y_in_vals = dem_in.y.values
            resampled = False

        # 2. PRE-PROCESSING: REFERENCE POINTS
        # Construct points array [x, y, z]. Can include all the points -
        # coregister function will filter out NaN and masked values.
        points_vals = [
            points_df.geometry.x.values,
            points_df.geometry.y.values,
            points_df.h_mean.values,
        ]

        # 3. COREGISTRATION
        new_dem_array, metadata_dict = coregister(
            dem_in_vals,
            points_vals,
            x_in_vals,
            y_in_vals,
            resolution,
            mask_in_vals,
            max_horiz_offset=max_horiz_offset,
            rmse_step_thresh=rmse_step_thresh,
            max_iterations=max_iterations,
            verbose=verbose,
        )

        # ADD ICESAT-2 SPECIFIC TO METADATA
        if "request_date_dt" in points_df.columns:

            if metadata_dict["coreg_status"] == "failed":
                metadata_dict["points_dt_days_max"] = None
                metadata_dict["points_dt_days_count"] = None

            else:
                dt_days_max = int(
                    points_df["request_date_dt"].dt.round("d").dt.days.abs().max()
                )
                metadata_dict["points_dt_days_max"] = dt_days_max

                dt_days = points_df["request_date_dt"].dt.round("d").dt.days
                dt_days_counts_dict = dt_days.value_counts().to_dict()
                metadata_dict["points_dt_days_count"] = dt_days_counts_dict

        metadata_dict["coregistration_type"] = "reference_icesat2"

        # APPLY COREGISTRATION TO NON-RESAMPLED DATASET
        if resampled:

            # If the coregistration estimation has been run on a DEM that has been
            # upsampled to 10 m, we have to apply the calculated offsets to original
            # high-resolution image

            if metadata_dict["coreg_status"] == "failed":
                new_dem_array = self._obj.squeeze().values.astype("float32")

            else:
                pn = [
                    metadata_dict["z_offset"],
                    metadata_dict["x_offset"],
                    metadata_dict["y_offset"],
                ]

                # Apply offsets
                if pn[1] != 0 and pn[2] != 0:
                    print("Applying transform to original image...")
                    new_dem_array = shift_dem(
                        self._obj.squeeze().values,
                        pn,
                        self._obj.squeeze().x.values,
                        self._obj.squeeze().y.values,
                    ).astype("float32")

                else:
                    print("Applying dz to original image...")
                    new_dem_array = (self._obj.squeeze() - pn[0]).astype("float32")

        # Return coregistered DEM
        if return_stats is True:
            return (self._obj.fillna(0) * 0 + new_dem_array, metadata_dict)

        else:
            return self._obj.fillna(0) * 0 + new_dem_array

    def coregister(
        self,
        reference: DataArray,
        stable_mask: Optional[DataArray] = None,
        return_stats: Optional[bool] = False,
        max_horiz_offset: float = 50,
        rmse_step_thresh: float = -0.001,
        max_iterations: int = 5,
        verbose: Optional[bool] = True,
    ) -> DataArray:
        """This function is deprecated, and redirects to the `coregister_dems`
        function instead. Please consult the `coregister_dems` function for
        documentation.
        """

        warn(
            "The `coregister` function is deprecated and will be removed in a "
            "future version. Please use coregister_dems instead.",
            DeprecationWarning,
            stacklevel=2,
        )

        return self.coregister_dems(
            reference,
            stable_mask,
            return_stats,
            max_horiz_offset,
            rmse_step_thresh,
            max_iterations,
            verbose,
        )

    def coregister_dems(
        self,
        reference: DataArray,
        stable_mask: Optional[DataArray] = None,
        return_stats: Optional[bool] = False,
        max_horiz_offset: float = 50,
        rmse_step_thresh: float = -0.001,
        max_iterations: int = 5,
        verbose: Optional[bool] = True,
    ) -> DataArray:
        """
        Coregisters the scene against a reference DEM based on the Nuth and Kääb (2011)
        method, as implemented by the PGC.

        Original code available here:
        https://github.com/PolarGeospatialCenter/setsm_postprocessing_python/blob/master/lib/scenes2strips.py

        :param reference: reference DEM DataArray of same extent and resolution.
        :type reference: xarray.DataArray
        :param stable_mask: Mask dataarray where 1 = stable region to be used for
            coregistration. If none, coregisters based on the entire scene.
        :type stable_mask: xarray.DataArray
        :param return_stats: If true, returns a metadata dictionary alongside
            the coregistered DEM dataarray. The dictionary contains the coregistration
            status (coregistered, dz_only, failed), the transformation parameter,
            the translation error, and the RMS error. Defaults to False
        :type return_stats: bool
        :param max_horiz_offset: maximum horizontal offset, beyond which XY
            coregistration is ignored and only Z values are corrected. Defaults
            to 15.
        :type max_horiz_offset: float
        :param rmse_step_thresh: break coregistration loop if rmse step is above
            threshold, defaults to -0.001.
        :type rmse_step_thresh: float
        :param max_iterations: max iterations to attempt, defaults to 5.
        :type max_iterations: int
        :param verbose: Print detailed progress statements, defaults to True.
        :type verbose: bool

        :returns: coregistered (rio)xarray dataarray
        :rtype: DataArray
        :returns: metadata dictionary (if return_stats=True)
        :rtype: dict
        """

        check_match = geospatial_match(self._obj, reference, return_info=True)
        if not check_match is True:
            raise ValueError(
                "Input DEM and reference DEM do not share geospatial information: "
                f"{check_match}. Consider padding (`.rio.pad_box`) or reprojecting "
                "(`.rio.repoject_match`) DEMs."
            )

        if stable_mask is None:
            print(
                "No `stable_mask` provided: assuming entire scene is suitable "
                "reference surface."
            )
            stable_mask = self._obj * 0 + 1
        else:
            check_match = geospatial_match(self._obj, stable_mask, return_info=True)
            if not check_match is True:
                raise ValueError(
                    "Input DEM and stable_mask do not share geospatial information: "
                    "{check_match}. Consider padding (`.rio.pad_box`) or reprojecting "
                    "(`.rio.repoject_match`) stable_mask."
                )

        resolution = get_resolution(self._obj)

        new_dem_array, metadata_dict = coregister(
            self._obj.values,
            reference.values,
            reference.x.values,
            reference.y.values,
            resolution,
            stable_mask.values,
            max_horiz_offset=max_horiz_offset,
            rmse_step_thresh=rmse_step_thresh,
            max_iterations=max_iterations,
            verbose=verbose,
        )

        metadata_dict["coregistration_type"] = "reference_dem"

        # print(f"Translating: {trans[0]:.2f} X, {trans[1]:.2f} Y, {trans[2]:.2f} Z")
        # new_dem_array = shift_dem(self._obj.values, trans)

        if return_stats is True:
            return (self._obj.fillna(0) * 0 + new_dem_array, metadata_dict)

        else:
            return self._obj.fillna(0) * 0 + new_dem_array

    def terrain(
        self,
        attribute: Optional[str | list[str]] = "all",
        method: Literal["Florinsky", "ZevenbergThorne"] = "Florinsky",
        degrees: Optional[bool] = True,
        resolution: Optional[float] = None,
        hillshade_multidirectional: Optional[bool] = False,
        hillshade_altitude: Optional[float] = 45.0,
        hillshade_azimuth: Optional[float] = 315.0,
        hillshade_z_factor: Optional[float] = 2.0,
    ) -> DataArray | Dataset:
        """
        Returns terrain attributes as an xarray DataSet. Available attributes are as
        follows:

         - Slope (`slope`)
         - Aspect (`aspect`)
         - Hillshade (`hillshade`)
         - Plan curvature (`plan_curvature`)
         - Horizontal curvature (`horizontal_curvature`)
         - Vertical curvature (`vertical_curvature`)
         - Horizontal excess curvature (`horizontal_excess`)
         - Vertical excess curvature (`vertical_excess`)
         - Mean curvature (`mean_curvature`)
         - Gaussian curvature (`gaussian_curvature`)
         - Unsphericity curvature (`unsphericity_curvature`)
         - Minimal curvature (`minimal_curvature`)
         - Maximal curvature (`maximal_curvature`)

        Explanations of the above attributes, and their potential applications, can
        be found in Florinsky (2017; doi.org/10.1177/0309133317733667)

        By default, computation of partial derivatives for geomorphometric variable
        calculation is performed following Florinsky (2009), who calculates the third-
        order partial derivative from a 5 x 5 window, which may be more appropriate
        for high-resolution DEMs with some amount of noise. The `method` parameter
        allows for the selection of a more traditional method using a second-order
        derivative from a 3 x 3 window (`ZevenbergThorne`), as outlined by Zevenbergen
        and Throrne (1987).

        :param attribute: The attribute(s) to calculate, as a string or list.
        :type attribute: str | list
        :param method: Method to calculate geomorphometric parameters: "Florisnky" or
            "ZevenbergThorne", defaults to "Florinsky".
        :type method: Literal["Florinsky", "ZevenbergThorne"]
        :param resolution: The DEM resolution. Calculated if empty, defaults to None
        :type resolution: float
        :param degrees: Outputs degrees if True, radians if False, defaults to True
        :type degrees: bool
        :param hillshade_multidirectional: Calculates multidirectional hillshade
            following USGS/GDAL method (Mark, 1992). Defaults to False.
        :type hillshade_multidirectional: bool
        :param hillshade_altitude: Hillshade altitude in degrees (0-90°) from vertical,
            defaults to 45˚
        :type hillshade_altitude: float
        :param hillshade_azimuth: Hillshade azimuth in degrees (0-360°) clockwise from
            North, defaults to 315
        :type hillshade_azimuth: float
        :param hillshade_z_factor: Vertical exaggeration factor, defaults to 1
        :type hillshade_z_factor: float

        :returns: selected terrain attributes as a (rio)xarray DataSet
        :rtype: DataSet
        """

        # Validate attribute(s), resolution, and method
        if isinstance(attribute, str):
            if attribute.lower() == "all":
                attribute = LIST_ATTRIBUTES
            else:
                attribute = [attribute]

        for a in attribute:
            if a.lower() not in LIST_ATTRIBUTES:
                raise ValueError(
                    f"`{a}` is not a valid attribute. Atributes be from list {LIST_ATTRIBUTES}"
                )

        if resolution is None:
            resolution = get_resolution(self._obj)

        if method.lower() not in [m.lower() for m in LIST_METHODS]:
            raise ValueError(
                f"Geomorphometric method must be one of: {LIST_METHODS}. Currently {method}."
            )

        if not isinstance(hillshade_multidirectional, bool):
            raise ValueError("`hillshade_multidirectional` must be boolean")

        # calculate p and q
        p_arr = p(self._obj.values, resolution, method).astype(np.float32)
        q_arr = q(self._obj.values, resolution, method).astype(np.float32)

        # if requested variables require second order variables, calculate them
        attrs_requiring_second_order = [
            attr for attr in attribute if attr in LIST_REQUIRING_SECOND_ORDER
        ]
        if len(attrs_requiring_second_order) >= 1:
            r_arr = r(self._obj.values, resolution, method).astype(np.float32)
            s_arr = s(self._obj.values, resolution, method).astype(np.float32)
            t_arr = t(self._obj.values, resolution, method).astype(np.float32)

        # Calculate attributes
        if ("slope" in attribute) or (
            ("hillshade" in attribute) & (hillshade_z_factor == 1)
        ):
            slope_arr = slope(p_arr, q_arr)

        if (
            ("aspect" in attribute)
            or ("hillshade" in attribute)
            or (hillshade_multidirectional is True)
        ):
            aspect_arr = aspect(p_arr, q_arr)

        if "hillshade" in attribute:
            if hillshade_multidirectional is False:
                # calculate unidirectional hillshade
                if hillshade_z_factor == 1:
                    hillshade_arr = hillshade(
                        np.rad2deg(slope_arr),
                        np.rad2deg(aspect_arr),
                        hillshade_altitude,
                        hillshade_azimuth,
                        norm=True,
                    )
                if hillshade_z_factor != 1:
                    hillshade_arr = hillshade(
                        np.rad2deg(
                            slope(
                                p(
                                    self._obj.values * hillshade_z_factor,
                                    resolution,
                                    method,
                                ),
                                q(
                                    self._obj.values * hillshade_z_factor,
                                    resolution,
                                    method,
                                ),
                            )
                        ),
                        np.rad2deg(aspect_arr),
                        hillshade_altitude,
                        hillshade_azimuth,
                        norm=True,
                    )
            else:
                # calculate multidirectional hillshade following
                # https://pubs.usgs.gov/of/1992/of92-422/of92-422.pdf

                # calculate weights
                w_225 = np.sin(aspect_arr - np.deg2rad(225)) ** 2
                w_270 = np.sin(aspect_arr - np.deg2rad(270)) ** 2
                w_315 = np.sin(aspect_arr - np.deg2rad(315)) ** 2
                w_360 = np.sin(aspect_arr - np.deg2rad(360)) ** 2

                if hillshade_z_factor == 1:
                    slope_input = slope_arr

                if hillshade_z_factor != 1:
                    slope_input = slope(
                        p(
                            self._obj.values * hillshade_z_factor,
                            resolution,
                            method,
                        ),
                        q(
                            self._obj.values * hillshade_z_factor,
                            resolution,
                            method,
                        ),
                    )

                hs = (
                    w_225
                    * hillshade(
                        np.rad2deg(slope_input),
                        np.rad2deg(aspect_arr),
                        60,
                        225,
                        norm=False,
                    )
                    + w_270
                    * hillshade(
                        np.rad2deg(slope_input),
                        np.rad2deg(aspect_arr),
                        60,
                        270,
                        norm=False,
                    )
                    + w_315
                    * hillshade(
                        np.rad2deg(slope_input),
                        np.rad2deg(aspect_arr),
                        60,
                        315,
                        norm=False,
                    )
                    + w_360
                    * hillshade(
                        np.rad2deg(slope_input),
                        np.rad2deg(aspect_arr),
                        60,
                        360,
                        norm=False,
                    )
                ) / 2

                del w_225, w_270, w_315, w_360

                # get normalised hillshade with darkest points as 0
                # hillshade_arr = 1 - (hs - np.nanmin(hs)) / (
                hillshade_arr = (hs - np.nanmin(hs)) / (np.nanmax(hs) - np.nanmin(hs))
                del hs

        if "plan_curvature" in attribute:
            plan_curvature_arr = plan_curvature(p_arr, q_arr, t_arr, r_arr, s_arr)

        if ("horizontal_curvature" in attribute) or ("horizontal_excess" in attribute):
            horizontal_curvature_arr = horizontal_curvature(
                p_arr, q_arr, t_arr, r_arr, s_arr
            )

        if ("vertical_curvature" in attribute) or ("vertical_excess" in attribute):
            vertical_curvature_arr = vertical_curvature(
                p_arr, q_arr, t_arr, r_arr, s_arr
            )

        if (
            ("mean_curvature" in attribute)
            or ("unsphericity_curvature" in attribute)
            or ("maximal_curvature" in attribute)
            or ("minimal_curvature" in attribute)
            or ("horizontal_excess" in attribute)
            or ("vertical_excess" in attribute)
        ):
            mean_curvature_arr = mean_curvature(p_arr, q_arr, t_arr, r_arr, s_arr)

        if (
            ("gaussian_curvature" in attribute)
            or ("unsphericity_curvature" in attribute)
            or ("maximal_curvature" in attribute)
            or ("minimal_curvature" in attribute)
            or ("horizontal_excess" in attribute)
            or ("vertical_excess" in attribute)
        ):
            gaussian_curvature_arr = gaussian_curvature(
                p_arr, q_arr, t_arr, r_arr, s_arr
            )

        if (
            ("unsphericity_curvature" in attribute)
            or ("maximal_curvature" in attribute)
            or ("minimal_curvature" in attribute)
            or ("horizontal_excess" in attribute)
            or ("vertical_excess" in attribute)
        ):
            unsphericity_curvature_arr = unsphericity_curvature(
                mean_curvature_arr, gaussian_curvature_arr
            )

        if (
            ("minimal_curvature" in attribute)
            or ("horizontal_excess" in attribute)
            or ("vertical_excess" in attribute)
        ):
            minimal_curvature_arr = minimal_curvature(
                mean_curvature_arr, unsphericity_curvature_arr
            )

        if "maximal_curvature" in attribute:
            maximal_curvature_arr = maximal_curvature(
                mean_curvature_arr, unsphericity_curvature_arr
            )

        # Then convert necessary values into degrees from radians...
        if degrees is True:
            if "slope" in attribute:
                slope_arr = np.rad2deg(slope_arr)
            if "aspect" in attribute:
                aspect_arr = np.rad2deg(aspect_arr)

        # Then convert into xarray dataset:
        xds = self._obj.to_dataset(name="dem", promote_attrs=True)

        if "slope" in attribute:
            xds["slope"] = self._obj * 0 + slope_arr
        if "aspect" in attribute:
            xds["aspect"] = self._obj * 0 + aspect_arr
        if "hillshade" in attribute:
            xds["hillshade"] = self._obj * 0 + hillshade_arr

        if "plan_curvature" in attribute:
            xds["plan_curvature"] = self._obj * 0 + plan_curvature_arr
        if "horizontal_curvature" in attribute:
            xds["horizontal_curvature"] = self._obj * 0 + horizontal_curvature_arr
        if "vertical_curvature" in attribute:
            xds["vertical_curvature"] = self._obj * 0 + vertical_curvature_arr
        if "horizontal_excess" in attribute:
            xds["horizontal_excess"] = self._obj * 0 + (
                horizontal_curvature_arr - minimal_curvature_arr
            )
        if "vertical_excess" in attribute:
            xds["vertical_excess"] = self._obj * 0 + (
                vertical_curvature_arr - minimal_curvature_arr
            )
        if "mean_curvature" in attribute:
            xds["mean_curvature"] = self._obj * 0 + mean_curvature_arr
        if "gaussian_curvature" in attribute:
            xds["gaussian_curvature"] = self._obj * 0 + gaussian_curvature_arr
        if "unsphericity_curvature" in attribute:
            xds["unsphericity_curvature"] = self._obj * 0 + unsphericity_curvature_arr
        if "minimal_curvature" in attribute:
            xds["minimal_curvature"] = self._obj * 0 + minimal_curvature_arr
        if "maximal_curvature" in attribute:
            xds["maximal_curvature"] = self._obj * 0 + maximal_curvature_arr

        # If you have only asked for one attribute, return a dataarray rather than a
        # dataset
        if len(attribute) == 1:
            xds = xds[attribute[0]]

        return xds

    def mask_ocean(
        self,
        candidate_height_thresh_m: Optional[float] = 10,
        candidate_area_thresh_km2: Optional[float] = 1,
        near_sealevel_thresh_m: Optional[float] = 5,
        return_sealevel_as_zero: Optional[bool] = False,
        return_mask: Optional[bool] = False,
    ) -> DataArray:
        """Returns a DEM with mélange/ocean regions filtered out, adapting the method
        of Shiggins _et al._ (2023). If no likely sea level is identified, returns the
        original DEM.

        WARNING: The DEM must be geoid-corrected for this function to work correctly.

        :param candidate_height_thresh_m: Maximum value relative to geoid to be considered
            as SL, in m, defaults to 10
        :type candidate_height_thresh_m: float
        :param candidate_area_thresh_km2: Minimum area beneath `candidate_height_thresh_m`
            to be considered for sea level assessment, in km^2. Defaults to 1
        :type candidate_area_thresh_km2: float
        :param near_sealevel_thresh_m: Filter out regions below this value, in metres above
            sea level. Defaults to 5
        :type near_sealevel_thresh_m: float
        :param return_sealevel_as_zero: If True, subtracts that estimated sea level from
            the DEM, returning new height values with the estimated sea level set to 0
            m. Necessary for calculating accurate iceberg heights. Defaults to False.
        :type return_sealevel_as_zero: bool
        :param return_mask: Return the sea level mask rather than the masked DEM.
            Defaults to False.
        :type return_mask: bool

        :returns: Filtered DEM as xarray DataArray
        :rtype: DataArray
        """

        # resolution = get_resolution(self._obj)

        # mask = self._obj.pdt.get_ocean_mask(
        #     candidate_height_thresh_m,
        #     candidate_area_thresh_km2,
        #     near_sealevel_thresh_m,
        # )

        est_sea_level = self._obj.pdt.get_sea_level(
            candidate_height_thresh_m,
            candidate_area_thresh_km2,
        )

        if est_sea_level is None:
            mask = None
        else:
            mask = self._obj > (est_sea_level + near_sealevel_thresh_m)

        if return_mask is True:
            if mask is None:
                warn("No sea level detected in image. Returning `None` as mask object.")
            return mask
        else:
            if mask is None:
                warn("No sea level detected in image. Returning original DEM")
                return self._obj
            else:
                if return_sealevel_as_zero is True:
                    return self._obj.where(mask) - est_sea_level
                elif return_sealevel_as_zero is False:
                    return self._obj.where(mask)
                else:
                    raise ValueError("return_sealevel_as_zero must be True or False")

    def get_sea_level(
        self,
        candidate_height_thresh_m: Optional[float] = 10,
        candidate_area_thresh_km2: Optional[float] = 1,
    ) -> float:
        """Returns estimated sea level following method of Shiggins _et al._ (2023). If no
        candidate sea level is identified, None is returned. DEM must be geoid-corrected.

        :param candidate_height_thresh_m: Maximum value relative to geoid to be considered
            as SL, in m, defaults to 10
        :type candidate_height_thresh_m: float
        :param candidate_area_thresh_km2: Minimum area beneath `candidate_height_thresh_m`
            to be considered for sea level assessment, in km^2, defaults to 1
        :type candidate_area_thresh_km2: float

        :returns: Estimated sea level
        :rtype: float
        """

        resolution = get_resolution(self._obj)

        # Get values close to sea level as 1D numpy array
        near_sealevel_values = self._obj.values.ravel()
        near_sealevel_values = near_sealevel_values[
            near_sealevel_values < candidate_height_thresh_m
        ]

        # Skip melange filtering if the candidate region is less than the threshold area
        thresh_px_n = int(candidate_area_thresh_km2 * 1e6 / (resolution * 2))
        if len(near_sealevel_values) < thresh_px_n:
            return None

        # Else, construct a 0.25m resolution histogram between -15 and +15:
        else:
            bin_edges = np.arange(-15.125, 15.375, 0.25)
            bin_centres = bin_edges[:-1] + 0.125
            hist, _ = np.histogram(near_sealevel_values, bins=bin_edges)

            # The estimated sea level is the centre of the modal bin
            est_sea_level = bin_centres[np.argmax(hist)]

            return est_sea_level

    def mask_icebergs(
        self,
        area_thresh_m2: Optional[float | int] = 1e6,
        retain_icebergs: Optional[bool] = False,
        return_mask: Optional[bool] = False,
        connectivity: Literal[4, 8] = 4,
    ) -> DataArray:
        """After masking the ocean using the `mask_ocean()` function, icebergs will
        remain. This function gives you the opportunity to mask them as well by
        identifying the size of connected groups of pixels and masking above/below a
        threshold.

        By default, this function masks icebergs and retains terrestrial ice/land.
        However, by setting `retain_icebergs = True`, the function will instead mask
        ice/land and retain icebergs, allowing iceberg area and volume to be assessed.
        For accurate volume assessment, the DEM 0 m value must be set to the estimated
        sea level height from the `get_sea_level()` function.  This can be automated by
        setting `return_sealevel_as_zero = True` in the mask_ocean() function.

        :param area_thresh_m2: Size threshold between icebergs and terrestrial ice/land, in
            m2. Defaults to 1e6 m2 (1 km2)
        :type area_thresh_m2: float | int
        :param retain_icerbergs: By default, this function masks icebergs and retains
            terrestrial ice/land (by masking pixel groups below the threshold area). By
            switching the parameter to True, the function will instead mask out larger
            pixel groups and retain the smaller icebergs.
        :type retain_icebergs: bool
        :param return_mask: If True, returns mask of icebergs (where icebergs = 1). If
            False, returns input DEM with icebergs masked. Defaults to False.
        :type return_mask: bool
        :type area_thresh_m2: float | int
        :param connectivity: Connectivity with which to calculate iceberg regions. Must
            be 4 or 8. Defaults to 4.
        :type connectivity: int

        :returns: Either masked DEM or mask as xarray DataArray, depending on
            the input value of `return_mask`.
        :rtype: DataArray

        """

        resolution = get_resolution(self._obj)

        # Generate boolean mask in int8 format for cv2.connectedComponentsWithStats()
        connected_input = (~np.isnan(self._obj.values)).astype(np.int8).squeeze()

        # Generate labels and statistics
        num_labels, labels, stats, _ = connectedComponentsWithStats(
            connected_input, connectivity, CV_32S
        )

        # calculate areas of connected components from number of pixels and resolution
        areas = stats[:, 4] * resolution * resolution

        # Get labels that are greater than area_thesh
        if (retain_icebergs is True) and (return_mask is False):
            labels_valid = np.where(
                areas < area_thresh_m2, np.arange(num_labels), np.nan
            )
        else:
            labels_valid = np.where(
                areas > area_thresh_m2, np.arange(num_labels), np.nan
            )

        labels_valid = labels_valid[~np.isnan(labels_valid)]

        if return_mask is True:
            return self._obj * 0 + ~np.isin(labels, labels_valid)

        else:
            return self._obj.where(np.isin(labels, labels_valid))

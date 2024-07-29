"""
This module contains the functions necessary for coregistration according to Nuth and
Kääb (2011), based on the coregisterdems() function by Erik Husby at line 720 of the 
PGC setsm_processing scenes2strips.py module.

Original available at:
https://github.com/PolarGeospatialCenter/setsm_postprocessing_python/blob/fd36fd54933ec43f587902a4fdcd1acbd90951c2/lib/scenes2strips.py
"""

import operator

import numpy as np

from scipy import stats
from osgeo import gdal, gdal_array

from ._geomorphometry import p_f, q_f

# ignore divide-by-zero error
np.seterr(divide="ignore", invalid="ignore")


def coregisterdems(
    dem1,  # Reference DEM
    dem2,  # DEM to be coregistered
    x,
    y,
    mask,
    res,
    max_horiz_offset=50,
    rmse_step_thresh=-0.001,
    max_iterations=5,
):
    """
    Simplified version of Erik Husby's coregisterdems() Python function.

    INPUTS:
    dem_1, dem_2: 2D arrays (of same shape) of dems. dem2 is the dem to be coregistered
    mask: mask of regions to be used in coregistration process (1=VALID FOR COREGISTRATION)

    OUTPUTS:
    trans: the [dz,dx,dy] transformation parameter
    trans_err: 1-sigma errors of trans
    rms: root mean square of the transformation in the vertical from the residuals

    If the registration fails due to lack of overlap, NaNs are returned in p and perr.
    If the registration fails to converge or exceeds the maximum shift, the median
    vertical offset is applied.

    """

    # initial trans and RMSE settings
    p = np.zeros((3, 1))  # p  is prior iteration trans var
    pn = p.copy()  # pn is current iteration trans var
    perr = np.zeros((3, 1))  # perr is prior iteration regression errors
    pnerr = perr.copy()  # pnerr is current iteration regression errors
    d0 = np.inf  # initial RMSE

    # Edge case markers
    meddz = None
    return_meddz = False
    critical_failure = False

    it = 0
    while True:
        it += 1

        print(f"Planimetric Correction Iteration {it}")

        print(f"Offset (z,x,y): {pn[0, 0]:.3f}, {pn[1, 0]:.3f}, {pn[2, 0]:.3f}")
        # print(f"pn: {pn}")

        # Break loop if conditions reached
        if np.any(np.abs(pn[1:]) > max_horiz_offset):
            print(
                f"Maximum horizontal offset ({max_horiz_offset}) exceeded."
                "Consider raising the threshold if offsets are large."
            )
            return_meddz = True
            break

        # Apply offsets
        if pn[1] != 0 and pn[2] != 0:
            dem2n = shift_dem(dem2, pn.T[0], x, y).astype('float32')
        else:
            dem2n = dem2 - pn[0].astype('float32')

        # # Calculate slopes - original method from PGC
        # sy, sx = np.gradient(dem2n, res)
        # sx = -sx

        # Calculate slope - using Florinsky slope method (p = sx, q = sy)
        sy = q_f(dem2n, res)
        sx = p_f(dem2n, res)
        sy = -sy
        sx = -sx

        # Difference grids.
        dz = dem2n - dem1

        # Mask (in full script, both m1 and m2 are applied)
        dz[mask == 0] = np.nan

        # If no overlap between scenes, break the loop
        if np.all(np.isnan(dz)):
            print("No overlap")
            critical_failure = True
            break

        # Filter NaNs and outliers.
        n = (
            ~np.isnan(sx)
            & ~np.isnan(sy)
            & (np.abs(dz - np.nanmedian(dz)) <= 3 * np.nanstd(dz))
        )
        n_count = np.count_nonzero(n)

        if n_count < 10:
            print("Too few ({}) registration points: 10 required".format(n_count))
            critical_failure = True
            break

        # Get RMSE
        d1 = np.sqrt(np.mean(np.power(dz[n], 2)))
        print("RMSE = {}".format(d1))

        # Keep median dz if first iteration.
        if it == 1:
            meddz = np.median(dz[n])
            meddz_err = np.std(dz[n] / np.sqrt(n_count))
            d00 = np.sqrt(np.mean(np.power(dz[n] - meddz, 2)))

        # Get improvement in RMSE
        rmse_step = d1 - d0  # initial d0 == inf

        # break if rmse above threshold
        if rmse_step > rmse_step_thresh or np.isnan(d0):
            print(
                "RMSE step in this iteration ({:.5f}) is above threshold ({}), "
                "stopping and returning values of prior iteration. ".format(
                    rmse_step, rmse_step_thresh
                )
            )
            # If fails after first registration attempt,
            # set dx and dy to zero and subtract the median offset.
            if it == 2:
                print("Second iteration regression failure")
                return_meddz = True
            break
        elif it == max_iterations:
            print("Maximum number of iterations ({}) reached".format(max_iterations))
            break

        # Keep this adjustment.
        dem2out = dem2n.copy()
        p = pn.copy()
        perr = pnerr.copy()
        d0 = d1

        # Build design matrix.
        X = np.column_stack((np.ones(n_count, dtype=np.float32), sx[n], sy[n]))
        sx, sy = None, None  # release for data amangement

        # Solve for new adjustment.
        p1 = np.reshape(np.linalg.lstsq(X, dz[n], rcond=None)[0], (-1, 1))

        # Calculate p errors.
        _, R = np.linalg.qr(X)
        RI = np.linalg.lstsq(R, np.identity(3, dtype=np.float32), rcond=None)[0]
        nu = X.shape[0] - X.shape[1]  # residual degrees of freedom
        yhat = np.matmul(X, p1)  # predicted responses at each data point
        r = dz[n] - yhat.T[0]  # residuals
        normr = np.linalg.norm(r)

        dz = None  # release for memory managment

        rmse = normr / np.sqrt(nu)
        tval = stats.t.ppf((1 - 0.32 / 2), nu)

        se = rmse * np.sqrt(np.sum(np.square(np.abs(RI)), axis=1, keepdims=True))
        p1err = tval * se

        # Update shifts.
        pn = p + p1
        pnerr = np.sqrt(np.square(perr) + np.square(p1err))

        # END OF LOOP

    if return_meddz:
        print("Returning median vertical offset: {:.3f}".format(meddz))
        dem2out = dem2 - meddz
        p = np.array([[meddz, 0, 0]]).T
        perr = np.array([[meddz_err, 0, 0]]).T
        d0 = d00

    if critical_failure:
        print("Regression critical failure, returning NaN trans and RMSE")
        dem2out = dem2
        p = np.full((3, 1), np.nan)
        perr = np.full((3, 1), np.nan)
        d0 = np.nan

    print(
        "Final offset (z,x,y): {:.3f}, {:.3f}, {:.3f}".format(p[0, 0], p[1, 0], p[2, 0])
    )
    print("Final RMSE = {}".format(d0))

    return dem2out, p.T[0], perr.T[0], d0


def shift_dem(dem, trans, x, y):
    """
    Shifts DEM according to translation factors ascertained in coregisterdems function

    INPUT
    dem: DEM to correct (dem2 in coregisterdems parlance)
    trans: (dz, dx, dy) translation parameters

    """

    # print("Translating in Z direction")
    dem_coreg = dem - trans[0]

    if (trans[1] == 0) & (trans[2] == 0):
        pass
        # print("XY translation variables == 0. Skipping horizontal translation.")
    else:
        # print("Translating in XY direction.")
        print(f"Translating: {trans[0]:.2f} Z, {trans[1]:.2f} X, {trans[2]:.2f} Y")

        # Interpolation grid
        #     x, y = np.meshgrid(np.arange(dem.shape[1]), np.arange(dem.shape[0]))
        # x = np.arange(dem.shape[1])
        # y = np.arange(dem.shape[0])
        xi = x - trans[1]
        yi = y - trans[2]
        # xi = x + trans[2]
        # yi = y + trans[1]

        # Check that uniform spacing is maintained (sometimes rounding errors).
        if len(np.unique(np.diff(xi))) > 1:
            xi = np.round(xi, 4)
        if len(np.unique(np.diff(yi))) > 1:
            yi = np.round(yi, 4)

        # Correct in xy dimensions - linearly interpolating floating data to reference grid
        dem_coreg = interp2_gdal(
            xi, yi, dem_coreg, x, y, "linear", extrapolate=False, oob_val=np.nan
        )

    return dem_coreg


def interp2_gdal(X, Y, Z, Xi, Yi, interp_str, extrapolate=False, oob_val=np.nan):
    """
    Resample array data from one set of x-y grid coordinates to another.
    From PGC's raster_array_tools.py

    INPUT:
    X, Y: Original grid coordinates (meshgrid)
    Z: Array to be redampled.
    Zi, Yi: New grid coordinated (meshgrid)
    Extrapolate: Interpolate values for pixels that fall outside range of old grid coords. If False, set values to `oob_val`
    oob_val: Fill regions of output array where new grid coords fall outside range of old grid coords.

    OUTPUT:
    Zi: The resampled array.
    """

    dtype_gdal, promote_dtype = dtype_np2gdal(Z.dtype)
    if promote_dtype is not None:
        Z = Z.astype(promote_dtype)

    interp_gdal = interp_str2gdal(interp_str)

    mem_drv = gdal.GetDriverByName("MEM")

    ds_in = mem_drv.Create("", X.size, Y.size, 1, dtype_gdal)
    ds_in.SetGeoTransform((X[0], X[1] - X[0], 0, Y[0], 0, Y[1] - Y[0]))
    ds_in.GetRasterBand(1).WriteArray(Z)

    ds_out = mem_drv.Create("", Xi.size, Yi.size, 1, dtype_gdal)
    ds_out.SetGeoTransform((Xi[0], Xi[1] - Xi[0], 0, Yi[0], 0, Yi[1] - Yi[0]))

    gdal.ReprojectImage(
        ds_in,
        ds_out,
        "",
        "",
        interp_gdal,
    )
    # gdal.ReprojectImage(
    #     ds_in,
    #     ds_out,
    #     "",
    #     "",
    #     interp_gdal,
    # )

    Zi = ds_out.GetRasterBand(1).ReadAsArray()

    if not extrapolate:
        interp2_fill_oob(X, Y, Zi, Xi, Yi, oob_val)
    
    ds_in, ds_out = None, None

    return Zi


def dtype_np2gdal(dtype_np):
    # TODO: Write docstring.

    if dtype_np == bool:  # np.bool:
        promote_dtype = np.uint8
    elif dtype_np == np.int8:
        promote_dtype = np.int16
    elif dtype_np == np.float16:
        promote_dtype = np.float32
    else:
        promote_dtype = None

    if promote_dtype is not None:
        warn(
            "NumPy array data type ({}) does not have equivalent GDAL data type and is not "
            "supported, but can be safely promoted to {}".format(
                dtype_np, promote_dtype(1).dtype
            )
        )
        dtype_np = promote_dtype

    dtype_gdal = gdal_array.NumericTypeCodeToGDALTypeCode(dtype_np)
    if dtype_gdal is None:
        raise InvalidArgumentError(
            "NumPy array data type ({}) does not have equivalent "
            "GDAL data type and is not supported".format(dtype_np)
        )

    return dtype_gdal, promote_dtype


def interp_str2gdal(interp_str):
    """
    Turrns `interp_str` string into a GDAL object for use in gdal.ReprojectImage.
    Accepts "nearest", "linear", "cubic", "spline", "lanczos", "average", "mode".
    """

    interp_choices = (
        "nearest",
        "linear",
        "cubic",
        "spline",
        "lanczos",
        "average",
        "mode",
    )

    interp_dict = {
        "nearest": gdal.GRA_NearestNeighbour,
        "linear": gdal.GRA_Bilinear,
        "bilinear": gdal.GRA_Bilinear,
        "cubic": gdal.GRA_Cubic,
        "bicubic": gdal.GRA_Cubic,
        "spline": gdal.GRA_CubicSpline,
        "lanczos": gdal.GRA_Lanczos,
        "average": gdal.GRA_Average,
        "mode": gdal.GRA_Mode,
    }

    if interp_str not in interp_dict:
        raise ValueError(
            "`interp` must be one of {}, but was '{}'".format(
                interp_choices, interp_str
            )
        )

    return interp_dict[interp_str]


def interp2_fill_oob(X, Y, Zi, Xi, Yi, fillval=np.nan, coord_grace=True):
    # Rows and columns of Zi outside the domain of Z are made NaN.
    # Assume X and Y coordinates are monotonically increasing/decreasing
    # so hopefully we only need to work a short way inwards from the edges.

    Xi_size = Xi.size
    Yi_size = Yi.size

    Xmin = min(X[0], X[-1])
    Ymax = max(Y[0], Y[-1])

    x_lfttest_val = X[0]
    x_rgttest_val = X[-1]
    y_toptest_val = Y[0]
    y_bottest_val = Y[-1]

    if x_lfttest_val == Xmin:
        # X-coords increase from left to right.
        x_lfttest_op = operator.lt
        x_rgttest_op = operator.gt
    else:
        # X-coords decrease from left to right.
        x_lfttest_op = operator.gt
        x_rgttest_op = operator.lt

    if y_toptest_val == Ymax:
        # Y-coords decrease from top to bottom.
        y_toptest_op = operator.gt
        y_bottest_op = operator.lt
    else:
        # Y-coords increase from top to bottom.
        y_toptest_op = operator.lt
        y_bottest_op = operator.gt

    if coord_grace:
        x_grace = (X[1] - X[0]) / 64
        y_grace = (Y[1] - Y[0]) / 16
        x_lfttest_val -= x_grace
        x_rgttest_val += x_grace
        y_toptest_val -= y_grace
        y_bottest_val += y_grace
        x_lfttest_op = operator.le if x_lfttest_op(0, 1) else operator.ge
        x_rgttest_op = operator.le if x_rgttest_op(0, 1) else operator.ge
        y_toptest_op = operator.le if y_toptest_op(0, 1) else operator.ge
        y_bottest_op = operator.le if y_bottest_op(0, 1) else operator.ge

    i = 0
    while x_lfttest_op(Xi[i], x_lfttest_val) and i < Xi_size:
        Zi[:, i] = fillval
        i += 1
    i = -1
    while x_rgttest_op(Xi[i], x_rgttest_val) and i >= -Xi_size:
        Zi[:, i] = fillval
        i -= 1
    j = 0
    while y_toptest_op(Yi[j], y_toptest_val) and j < Yi_size:
        Zi[j, :] = fillval
        j += 1
    j = -1
    while y_bottest_op(Yi[j], y_bottest_val) and j >= -Yi_size:
        Zi[j, :] = fillval
        j -= 1

    return Zi

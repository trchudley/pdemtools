"""
Functions for calculating geomorphometry variables.
"""

import numba
import numpy as np


"""
Numba functions for quadratic 3x3 polynomial fit following Zevenbergen and Thorne (1987)
"""


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _p_zt(z, w):
    return (-z[0, -1] + z[0, 1]) / (2 * w)


@numba.njit
def p_zt(z, w=2):
    return _p_zt(z, w).astype('float32')


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _q_zt(z, w):
    return (z[-1, 0] - z[1, 0]) / (2 * w)


@numba.njit
def q_zt(z, w=2):
    return _q_zt(z, w).astype('float32')


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _r_zt(z, w):
    return ((z[0, -1] + z[0, 1]) / 2 - z[0, 0]) / (w**2)


@numba.njit
def r_zt(z, w=2):
    return _r_zt(z, w).astype('float32')


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _t_zt(z, w):
    return ((z[-1, 0] + z[1, 0]) / 2 - z[0, 0]) / (w**2)


@numba.njit
def t_zt(z, w=2):
    return _t_zt(z, w).astype('float32')


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _s_zt(z, w):
    return (1) * (-z[-1, -1] + z[-1, 1] + z[1, -1] - z[1, 1]) / (4 * w**2)


@numba.njit
def s_zt(z, w=2):
    return _s_zt(z, w).astype('float32')


"""
 Numba functions for third order polynomial 5x5 fit following Florinsky (2009)
"""


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _p_f(z, w):
    return (1 / (420 * w)) * (
        44 * (z[-2, 1] + z[2, 1] - z[-2, -1] - z[2, -1])
        + 31
        * (
            z[-2, -2]
            + z[2, -2]
            - z[-2, 2]
            - z[2, 2]
            + 2 * (z[-1, 1] + z[1, 1] - z[-1, -1] - z[1, -1])
        )
        + 17 * (z[0, 2] - z[0, -2] + 2 * (z[0, 1] - z[0, -1]))
        + 5 * (z[-1, 2] + z[1, 2] - z[-1, -2] - z[1, -2])
    )


@numba.njit
def p_f(z, w=2):
    return _p_f(z, w).astype('float32')


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _q_f(z, w):
    return (1 / (420 * w)) * (
        44 * (z[-1, -2] + z[-1, 2] - z[1, -2] - z[1, 2])
        + 31
        * (
            z[2, -2]
            + z[2, 2]
            - z[-2, -2]
            - z[-2, 2]
            + 2 * (z[-1, -1] + z[-1, 1] - z[1, -1] - z[1, 1])
        )
        + 17 * (z[-2, 0] - z[2, 0] + 2 * (z[-1, 0] - z[1, 0]))
        + 5 * (z[-2, -1] + z[-2, 1] - z[2, -1] - z[2, 1])
    )


@numba.njit
def q_f(z, w=2):
    return _q_f(z, w).astype('float32')


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _r_f(z, w):
    return (1 / (35 * w**2)) * (
        2
        * (
            z[-2, -2]
            + z[-2, 2]
            + z[-1, -2]
            + z[-1, 2]
            + z[0, -2]
            + z[0, 2]
            + z[1, -2]
            + z[1, 2]
            + z[2, -2]
            + z[2, 2]
        )
        - 2 * (z[-2, 0] + z[-1, 0] + z[0, 0] + z[1, 0] + z[2, 0])
        - z[-2, -1]
        - z[-2, 1]
        - z[-1, -1]
        - z[-1, 1]
        - z[0, -1]
        - z[0, 1]
        - z[1, -1]
        - z[1, 1]
        - z[2, -1]
        - z[2, 1]
    )


@numba.njit
def r_f(z, w=2):
    return _r_f(z, w).astype('float32')


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _t_f(z, w):
    return (1 / (35 * w**2)) * (
        2
        * (
            z[-2, -2]
            + z[-2, -1]
            + z[-2, 0]
            + z[-2, 1]
            + z[-2, 2]
            + z[2, -2]
            + z[2, -1]
            + z[2, 0]
            + z[2, 1]
            + z[2, 2]
        )
        - 2 * (z[0, -2] + z[0, -1] + z[0, 0] + z[0, 1] + z[0, 2])
        - z[-1, -2]
        - z[-1, -1]
        - z[-1, 0]
        - z[-1, 1]
        - z[-1, 2]
        - z[1, -2]
        - z[1, -1]
        - z[1, 0]
        - z[1, 1]
        - z[1, 2]
    )


@numba.njit
def t_f(z, w=2):
    return _t_f(z, w).astype('float32')


@numba.stencil(func_or_mode="constant", cval=np.nan)
def _s_f(z, w):
    return (1 / (100 * w**2)) * (
        z[-1, 1]
        + z[1, -1]
        - z[-1, -1]
        - z[1, 1]
        + 4 * (z[-2, 2] + z[2, -2] - z[-2, -2] - z[2, 2])
        + 2
        * (
            z[-2, 1]
            + z[-1, 2]
            + z[1, -2]
            + z[2, -1]
            - z[-2, -1]
            - z[-1, -2]
            - z[1, 2]
            - z[2, 1]
        )
    )


@numba.njit
def s_f(z, w=2):
    return _s_f(z, w).astype('float32')


"""
Wrapper functions to choose between Florinsky and Zevenberg methods
"""


def p(z, w, method):
    if method.lower() == "florinsky":
        return p_f(z.squeeze(), w)
    elif method.lower() == "zevenbergthorne":
        return p_zt(z.squeeze(), w)
    else:
        raise ValueError(
            "`method` variable must be one of ['Florinsky', 'ZevenbergThorne']"
        )


def q(z, w, method):
    if method.lower() == "florinsky":
        return q_f(z.squeeze(), w)
    elif method.lower() == "zevenbergthorne":
        return q_zt(z.squeeze(), w)
    else:
        raise ValueError(
            "`method` variable must be one of ['Florinsky', 'ZevenbergThorne']"
        )


def r(z, w, method):
    if method.lower() == "florinsky":
        return r_f(z.squeeze(), w)
    elif method.lower() == "zevenbergthorne":
        return r_zt(z.squeeze(), w)
    else:
        raise ValueError(
            "`method` variable must be one of ['Florinsky', 'ZevenbergThorne']"
        )


def s(z, w, method):
    if method.lower() == "florinsky":
        return s_f(z.squeeze(), w)
    elif method.lower() == "zevenbergthorne":
        return s_zt(z.squeeze(), w)
    else:
        raise ValueError(
            "`method` variable must be one of ['Florinsky', 'ZevenbergThorne']"
        )


def t(z, w, method):
    if method.lower() == "florinsky":
        return t_f(z.squeeze(), w)
    elif method.lower() == "zevenbergthorne":
        return t_zt(z.squeeze(), w)
    else:
        raise ValueError(
            "`method` variable must be one of ['Florinsky', 'ZevenbergThorne']"
        )


"""
Functions for processing curvature from derivatives
"""


def slope(p_arr, q_arr):
    """outputs in radians"""
    return np.arctan(np.sqrt(p_arr**2 + q_arr**2))


def aspect(p_arr, q_arr):
    """outputs in radians"""

    return np.arctan2(p_arr, q_arr) + np.pi

    # return np.deg2rad(
    #     -90 * (1 - np.sign(q_arr)) * (1 - np.abs(np.sign(p_arr)))
    #     + 180 * (1 + np.sign(p_arr))
    #     - (180 / np.pi)
    #     * np.sign(p_arr)
    #     * np.arccos(-q_arr / np.sqrt(p_arr**2 + q_arr**2))
    # )


def hillshade(slope, aspect, altitude=45, azimuth=315, norm=True):
    """accepts degrees"""

    if slope.max() < 2 * np.pi:
        raise Warning(
            "Maximum slope value is < 2π (slope.max()), ensure input is degrees"
        )
    if aspect.max() < 2 * np.pi:
        raise Warning(
            "Maximum slope value is < 2π (aspect.max()), ensure input is degrees"
        )

    hs = 255.0 * (
        (np.cos(np.deg2rad(altitude)) * np.cos(np.deg2rad(slope)))
        + (
            np.sin(np.deg2rad(altitude))
            * np.sin(np.deg2rad(slope))
            * np.cos(np.deg2rad(azimuth) - np.deg2rad(aspect))
        )
    ).astype("float32")

    # hs = np.cos(np.pi * 0.5 - np.deg2rad(aspect) - np.deg2rad(azimuth)) * np.sin(
    #     np.deg2rad(slope)
    # ) * np.sin(np.pi * 0.5 - np.deg2rad(altitude)) + np.cos(np.deg2rad(slope)) * np.cos(
    #     np.pi * 0.5 - np.deg2rad(altitude)
    # )

    # return normalised values
    if norm == True:
        return (hs - np.nanmin(hs)) / (np.nanmax(hs) - np.nanmin(hs))
    else:
        return hs


def plan_curvature(p_arr, q_arr, t_arr, r_arr, s_arr):
    return -(
        ((q_arr**2 * r_arr) - (2 * p_arr * q_arr * s_arr) + (p_arr**2 * t_arr))
        / np.sqrt((p_arr**2 + q_arr**2) ** 3)
    )


def horizontal_curvature(p_arr, q_arr, t_arr, r_arr, s_arr):
    return -(
        ((q_arr**2 * r_arr) - (2 * p_arr * q_arr * s_arr) + (p_arr**2 * t_arr))
        / ((p_arr**2 + q_arr**2) * np.sqrt(1 + p_arr**2 + q_arr**2))
    )


def vertical_curvature(p_arr, q_arr, t_arr, r_arr, s_arr):
    return -(
        ((p_arr**2 * r_arr) + (2 * p_arr * q_arr * s_arr) + (q_arr**2 * t_arr))
        / ((p_arr**2 + q_arr**2) * np.sqrt((1 + p_arr**2 + q_arr**2) ** 3))
    )


def mean_curvature(p_arr, q_arr, t_arr, r_arr, s_arr):
    return -(
        (1 + q_arr**2) * r_arr - 2 * p_arr * q_arr * s_arr + (1 + p_arr**2) * t_arr
    ) / (2 * np.sqrt(1 + p_arr**2 + q_arr**2) ** 3)


def gaussian_curvature(p_arr, q_arr, t_arr, r_arr, s_arr):
    return (r_arr * t_arr - s_arr**2) / ((1 + p_arr**2 + q_arr**2) ** 2)


def unsphericity_curvature(mean, gaussian):
    return np.sqrt(mean**2 - gaussian)


def maximal_curvature(mean, unsphericity):
    return mean + unsphericity


def minimal_curvature(mean, unsphericity):
    return mean - unsphericity

import pytest

import rioxarray as rxr
import xarray as xr
import geopandas as gpd
import numpy as np

import pdemtools as pdt

from numpy.testing import assert_allclose
from shapely.geometry import box


def generate_test_data_1():
    """Generates a 2D sin pattern"""

    x = np.arange(-500, 500, 2)
    y = np.arange(-500, 500, 2)

    X, Y = np.meshgrid(x, y)
    Z = np.sin((X**2 + Y**2) / np.max(X**2) * 10) * 10 + 30

    dem = xr.DataArray(
        data=Z,
        coords=[("x", np.flip(x)), ("y", y)],
        dims=["y", "x"],  # Adjust order if "x" and "y" represent column coordinates
    )
    dem.rio.write_crs("EPSG:3413", inplace=True)
    dem.rio.write_transform(inplace=True)  # Set the georeferencing

    return dem


def generate_test_data_2(x_offset=-13.4, y_offset=+14.6, z_offset=15.3):
    """As for `generate_test_data_1` but offsets by fixed amount"""

    x = np.arange(-500, 500, 2)
    y = np.arange(-500, 500, 2)

    x2 = np.arange(-500 + x_offset, 500 + x_offset, 2)
    y2 = np.arange(-500 + y_offset, 500 + y_offset, 2)

    X, Y = np.meshgrid(x2, y2)
    Z = np.sin((X**2 + Y**2) / np.max(X**2) * 10) * 10 + 30 + z_offset

    dem = xr.DataArray(
        data=Z,
        coords=[("x", np.flip(x)), ("y", y)],
        dims=["y", "x"],  # Adjust order if "x" and "y" represent column coordinates
    )
    dem.rio.write_crs("EPSG:3413", inplace=True)
    dem.rio.write_transform(inplace=True)  # Set the georeferencing

    return dem


def generate_test_data_3(z_offset=-30):
    """Generates a 2D sin pattern"""

    x = np.arange(-10000, 10000, 32)
    y = np.arange(-10000, 10000, 32)

    X, Y = np.meshgrid(x, y)
    Z = np.sin((X**2 + Y**2) / np.max(X**2) * 10) * 10 + 30 + z_offset

    dem = xr.DataArray(
        data=Z,
        coords=[("x", np.flip(x)), ("y", y)],
        dims=["y", "x"],  # Adjust order if "x" and "y" represent column coordinates
    )
    dem.rio.write_crs("EPSG:3413", inplace=True)
    dem.rio.write_transform(inplace=True)  # Set the georeferencing

    return dem


def test_geomorphometry():

    variable_min_max = {
        "slope": (0.0, 24.721811294555664),
        "aspect": (0.0, 360.0),  # (1.3660377589985728e-05, 359.9999694824219),
        "hillshade": (0.0, 1.0),
        "horizontal_curvature": (-0.0007999989320524037, 0.0007856095326133072),
        "vertical_curvature": (-0.027519449591636658, 0.022581202909350395),
        "mean_curvature": (-0.013776100240647793, 0.011318518780171871),
        "gaussian_curvature": (-1.0310398465662729e-05, 9.49732475419296e-06),
        "unsphericity_curvature": (0.0, 0.013750489801168442),
        "minimal_curvature": (-0.02751944772899151, 0.0007855556323193014),
        "maximal_curvature": (-0.0007999998051673174, 0.022581204771995544),
    }

    dem = generate_test_data_1()
    dem = dem.pdt.terrain(variable_min_max.keys())

    for variables in variable_min_max.keys():
        min_val = dem[variables].min().item()
        max_val = dem[variables].max().item()
        print(f"{variables} | min: {min_val} | max: {max_val}")
        assert_allclose(min_val, variable_min_max[variables][0])
        assert_allclose(max_val, variable_min_max[variables][1])


def test_coregistration():

    dem_1 = generate_test_data_1()
    dem_2 = generate_test_data_2()

    mask_gdf = gpd.GeoDataFrame(geometry=[box(-200, -100, 0, 100)], crs=3413)
    mask = pdt.data.bedrock_mask_from_vector(mask_gdf, dem_1)

    dem_2_coreg = dem_2.pdt.coregister_dems(dem_1, mask)

    # diff = (dem_2 - dem_1).mean().item()
    diff_coreg = (dem_2_coreg - dem_1).mean().item()

    # expected_diff = 15.476494902013904
    expected_diff_coreg = 0.19603341170854155

    # print(diff, corrected_diff)

    assert_allclose(expected_diff_coreg, diff_coreg)


def test_geoid_correction():

    dem = generate_test_data_1()
    geoid = dem * 0 + 30

    dem_geoidcor = dem.pdt.geoid_correct(geoid)

    # expected_mean = 30.835231084697448
    expected_geoidcor_mean = 0.8352310846974458

    # mean = dem.mean().item()
    geoidcor_mean = dem_geoidcor.mean().item()

    assert_allclose(expected_geoidcor_mean, geoidcor_mean)


def test_mask_ocean():

    dem = generate_test_data_3(z_offset=-25)
    dem_masked = dem.pdt.mask_ocean(near_sealevel_thresh_m=5)

    orig_min = dem.min().item()
    new_min = dem_masked.min().item()
    expected_min = 0.0014336554957807834

    # print(dem_masked)
    # print(dem_masked.rio.resolution())
    # print(dem.pdt.get_sea_level())
    # print(orig_min)
    # print(new_min)

    assert new_min > orig_min
    assert_allclose(new_min, expected_min)


if __name__ == "__main__":
    test_geomorphometry()
    test_coregistration()
    test_geoid_correction()
    test_mask_ocean()

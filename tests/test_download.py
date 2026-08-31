import pytest

import dask.array
import pdemtools as pdt

from numpy.testing import assert_allclose


def test_download():

    # Test download a small section of Store Glacier
    bounds = (-205650, -2132203, -204722, -2131404)
    expected_mean_elev = 263.5019226074219

    dem = pdt.load.mosaic(
        dataset="arcticdem", bounds=bounds, resolution=32, version="v4.1"
    )

    mean_elev = dem.mean().item()

    assert_allclose(mean_elev, expected_mean_elev, rtol=1e-5)


def test_lazy_mosaic_merge():
    """Regression test: bounds spanning the boundary between ArcticDEM 32m
    supertiles `19_38` and `19_39`, so `mosaic()` must merge >1 tile. With
    `chunks` set, the merge must stay lazy (dask-backed) rather than forcing
    a compute, and the merged result must retain valid CRS/transform info.
    """

    bounds = (-201000, -2132203, -199000, -2131404)

    dem = pdt.load.mosaic(
        dataset="arcticdem", bounds=bounds, resolution=32, version="v4.1", chunks=True
    )

    # merge must not have triggered computation
    assert isinstance(dem.data, dask.array.Array)

    # geospatial metadata must survive the merge
    assert dem.rio.crs is not None
    assert dem.rio.transform() is not None

    # values must still be correct once computed
    computed = dem.compute()
    assert computed.notnull().any()


if __name__ == "__main__":
    test_download()
    test_lazy_mosaic_merge()

import pytest

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


if __name__ == "__main__":
    test_download()

import pytest

import pdemtools as pdt

from importlib import resources

SAMPLE_INDEX_FNAME = "test_arcticdem_index_kiv_steenstrup.parquet"
SAMPLE_INDEX_FPATH = resources.files("pdemtools.test_data").joinpath(SAMPLE_INDEX_FNAME)


def test_search():

    bounds = (459000, -2539000, 473000, -2528000)

    gdf = pdt.search(
        SAMPLE_INDEX_FPATH,
        bounds,
        dates="20170101/20221231",
        months=[6, 7, 8, 9],
        years=[2017, 2021],
        baseline_max_hours=24,
        sensors=["WV03", "WV02", "WV01"],
        min_aoi_frac=0.7,
    )

    assert len(gdf) == 6


if __name__ == "__main__":
    test_search()

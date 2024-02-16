# Install

## Installation Instructions

 `pdemtools` can be downloaded from [Github](https://github.com/) and installed using [`pip`](https://pip.pypa.io/en/stable/), with dependencies installable using [`conda`](https://docs.conda.io/en/latest/):

```bash
# Clone the pdemtools github repository
git clone git@github.com:trchudley/pdemtools.git

# Move to pdemtools directory
cd pdemtools

# Conda option 1: create a fresh conda environment
conda env create -f environment.yml
conda activate pdemtools_env

# -- OR --

# Conda option 2 : update existing conda environment
conda activate existing_env
conda env update -f environment.yml

# Install pdemtools
pip install -e .
```

It is not normally advisable to mix `pip` and `conda`, but in this use case `conda` will install all the dependencies (listed below) and `pip` sets up the appropriate paths to `pdemtools`.

 If you use [`mamba`](https://mamba.readthedocs.io/en/latest/) instead of `conda`, simple use `mamba` instead of `conda` as a drag-and-drop replacement.

## Dependencies

`pdemtools` dependencies are as follows:

 - rioxarray (and xarray)xw
 - geopandas (and pandas)
 - shapely
 - openCV
 - numba
 - numpy
 - scipy
 - GDAL
 - pyarrow
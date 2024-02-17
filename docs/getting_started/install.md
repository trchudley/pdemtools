# Install

`pdemtools` can be installed locally by cloning the GitHub repository. Releases on the [PyPi](https://pypi.org/) package index and [conda-forge](https://conda-forge.org/) repository are planned in the near future. 

## Cloning from GitHub

`pdemtools` can be cloned from the [Github repoository](https://github.com/), and installed in your Python environment by navigating to the top level of the repository and using [`pip`](https://pip.pypa.io/en/stable/).

```bash
# Clone the pdemtools github repository
git clone git@github.com:trchudley/pdemtools.git

# Move to pdemtools directory
cd pdemtools

# Initiate a new conda environment with dependencies
conda env create -f environment.yml
conda activate pdemtools_env

# Install pdemtools
pip install -e .
```

<!-- It is not normally advisable to mix `pip` and `conda`, but in this use case `conda` will install all the dependencies (listed below) and `pip` sets up the appropriate paths to `pdemtools`. -->

## Dependencies

`pdemtools` has the following dependencies:

 - rioxarray (and xarray)
 - geopandas (and pandas)
 - shapely
 - openCV
 - numba
 - numpy
 - scipy
 - GDAL
 - pyarrow

## A note on installing dependencies

There are [known errors when installing the GDAL package with pip](https://github.com/OSGeo/gdal/issues/2827), meaning that the GDAL dependency may create an error when installing through pip. As a result, please ensure GDAL is installed on your system prior to installing through pip. 

The easiest route to achieve this is to install GDAL alongside the other `pdemtools` dependencies using [`conda`](https://docs.conda.io/en/latest/) or its drop-in replacement [`mamba`](https://mamba.readthedocs.io/en/latest/) prior to installing pdemtools:

```bash
mamba create -n my_pdemtools_env "python>3.9"
mamba install rioxarray rasterio geopandas pandas shapely numpy gdal opencv-python scipy numba pyarrow
```

For ease of use, there is also [an `environment.yml` file](https://github.com/trchudley/pdemtools/blob/main/environment.yml) included within the GitHub repository to create a new environment with the required dependencies:

```bash
mamba env create -f environment.yml
mamba activate pdemtools_env
```

Or, you can use it to update a pre-existing environment:

```bash
mamba activate existing_env
mamba env update -f environment.yml
```






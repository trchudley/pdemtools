# Install

The recommended method to install `pdemtools` is with `conda`:

```bash
conda install pdemtools -c conda-forge
```

or with `mamba`:

```bash
mamba install pdemtools -c conda-forge
```

If you don't already have `conda` or `mamba` installed, the recommended method is via `micromamba`. Instructions can be found [here](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

It is also possible to install `pdemtools` with `pip`:

```bash
pip install pdemtools
```

However, there are [known errors when installing the GDAL package with pip](https://github.com/OSGeo/gdal/issues/2827), meaning that the GDAL dependency may create an error when installing through `pip`. If this occurs, please ensure GDAL is installed on your system by other means prior to installing through `pip`, or revert to using `conda`. 


## Developer install

For development purposes, `pdemtools` can be cloned from the [Github repoository](https://github.com/), and installed in your Python environment installed locally using `pip`. A conda environment configuration file is included for rapid set-up.

```bash
# Clone the pdemtools github repository
git clone git@github.com:trchudley/pdemtools.git

# Move to pdemtools directory
cd pdemtools

# Initiate a new conda environment with dependencies
conda env create -f environment.yml -n pdemtools_env
conda activate pdemtools_env

# Install pdemtools
pip install -e .
```

`pdemtools` has built-in unit tests, which can be run by installing `pytest` into the same environment and entering the commend `pytest` whilst in the top directory.

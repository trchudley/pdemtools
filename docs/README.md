
# Website/Docs

This directory contains the files necessary for sphinx to construct the readthedocs website (https://pdemtools.readthedocs.io/).

Unless you want to edit the website, you shouldn't need to touch anything in this directory. To edit the website, please refer to relevant [sphinx documentation](https://sphinx-tutorial.readthedocs.io/).

Website updates will be automatically ingested by readthedocs and updated on the https://pdemtools.readthedocs.io/en/latest/ site. At every version update, the 'latest' site will be saved as the 'stable' version.

If you wish to ensure the website compiles before uploading to readthedocs, you can do so by ensuring you have installed the necessary tools in your `pdemtools` python environment:

```bash
conda install sphinx==7.2.6 myst-parser nbsphinx nbsphinx-link sphinx-book-theme
```

The, to build the website using sphinx:

```bash
cd docs
make html
```

The website will then be available to view in `_build/html`. If it's struggling to compile, `make clean` might fix things in a pinch.

When you're happy, push to git (including `git add .` if you've got new files!).

Note that there will always be some quirks between the local sphinx compilation and the readthedocs version that might mean things aren't going to be perfect when they're online. It's worth checking any new/edited pages in the 'latest' version of the website before a new `pdemtools` version is released and commits a 'stable' version of the website for a while!
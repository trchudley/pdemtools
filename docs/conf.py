# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pdemtools"
# copyright = "2024, Tom Chudley"
author = "Tom Chudley"
# release = "0.6"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax"
    "myst_parser",  # mardown parser
    "nbsphinx",  # notebook viewer
    "nbsphinx_link",  # link to outside root directory
    # "sphinx_rtd_theme",  # readthedocs theme
    "sphinx_book_theme",  # readthedocs theme
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}


# add_module_names = False

# Add path to pdemtools
import os, sys

sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../src"))
sys.path.insert(0, os.path.abspath("../src/pdemtools"))

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "alabaster"
# html_theme = "sphinx_rtd_theme"
html_theme = "sphinx_book_theme"
html_static_path = ["_static"]

# Book theme options
html_theme_options = {
    "repository_url": "https://github.com/trchudley/pdemtools",
    "use_repository_button": True,
    "home_page_in_toc": True,
    #    "logo": {
    #       "image_light": "_static/logo-light.png",
    #       "image_dark": "_static/logo-dark.png",
    #    }
}
html_title = "pDEMtools"

# Site-wide logo!
html_logo = "_static/pdemtools_logo_1600px.png"

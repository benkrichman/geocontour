import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'geocontour/geocontour')))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'geocontour'
copyright = '2023, Benjamin Krichman'
author = 'Benjamin Krichman'
release = '1.2.2'
version = '1.2.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon','sphinx.ext.autodoc','myst_parser']
templates_path = ['_templates']
exclude_patterns = []
autodoc_member_order = 'bysource'
myst_enable_extensions=['html_image']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['source/_static']
html_logo = 'icon_geocontour.png'

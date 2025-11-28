# Configuration file for the Sphinx documentation builder.
#
# This file redirects users from quantms.readthedocs.io to docs.quantms.org

# -- Project information -----------------------------------------------------
project = "quantms"
copyright = "2024, bigbio"
author = "bigbio"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx_reredirects",
]

# -- Options for HTML output -------------------------------------------------
html_theme = "alabaster"

# -- sphinx-reredirects configuration ----------------------------------------
# Redirect the root page to the new documentation site
redirects = {
    "index": "https://docs.quantms.org/en/latest/",
}

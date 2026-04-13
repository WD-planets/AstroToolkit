import os
import sys

from bokeh.embed import file_html
from bokeh.layouts import Column, GridBox, Row
from bokeh.plotting import figure
from bokeh.resources import CDN

sys.path.insert(0, os.path.abspath("./tutorials"))
sys.path.insert(0, os.path.abspath("../../src"))

sphinx_gallery_conf = {
    "examples_dirs": "tutorials",
    "gallery_dirs": "auto_tutorials",
    "write_computation_times": False,
    "filename_pattern": r"\.py$",
    "ignore_pattern": r"^_.*\.py$",
    "reference_url": {"ATK": None},
    "run_stale_examples": False,
}

# -------------------
# Project Information
# -------------------

project = "AstroToolkit"
copyright = "2026, Ethan Moorfield"
author = "Ethan Moorfield"
release = "1.8.0"

html_favicon = "_static/logo/icon.png"

# ---------------------
# General Configuration
# ---------------------

# avoid duplicate label warnings
autosectionlabel_prefix_document = True


def setup(app):
    app.add_css_file("stylesheet.css")


extensions = [
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx_gallery.gen_gallery",
    "sphinxcontrib.video",
    "bokeh.sphinxext.bokeh_plot",
    "sphinx.ext.mathjax",
]
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# ------------
# HTML Options
# ------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {"collapse_navigation": False}

# -------------------
# InterSphinx Options
# -------------------

intersphinx_mapping = {
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "bokeh": ("https://docs.bokeh.org/en/latest/", None),
}

# --------------------------------
# Sphinx-Gallery Bokeh Integration
# --------------------------------


def _repr_html_(self):
    return file_html(self, CDN, "Bokeh Figure")


figure._repr_html_ = _repr_html_
Row._repr_html_ = _repr_html_
Column._repr_html_ = _repr_html_
GridBox._repr_html_ = _repr_html_

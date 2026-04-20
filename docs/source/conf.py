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

from docutils import nodes
from sphinx.ext.autodoc import ClassDocumenter
from sphinx.util.inspect import isfunction, ismethod


def autodoc_skip_member(app, what, name, obj, skip, options):
    """
    Skips private or default members
    """

    if name.startswith("__") and name.endswith("__"):
        return True
    if name.startswith("_"):
        return True

    return skip


def process_docstring(app, what, name, obj, options, lines):
    """
    Adds required parameters from _required attr for object initialisation via from_table() and from_dataframe()
    """

    # Only act on methods
    if what != "method":
        return

    # Get method name
    method_name = name.split(".")[-1]

    if method_name not in {"from_dataframe", "from_table"}:
        return

    # Get the owning class
    owner = getattr(obj, "__self__", None)

    if owner and hasattr(owner, "_required"):
        lines.append("")
        lines.append("The following parameters are required as keyword arguments:")
        for item in owner._required:
            lines.append("")
            lines.append(f"- {item}")


def setup(app):
    app.connect("autodoc-skip-member", autodoc_skip_member)
    app.connect("autodoc-process-docstring", process_docstring)
    app.add_css_file("stylesheet.css")


extensions = [
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx_gallery.gen_gallery",
    "sphinxcontrib.video",
    "bokeh.sphinxext.bokeh_plot",
    "sphinx.ext.mathjax",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
numpydoc_xref_param_type = False
numpydoc_attributes_as_param_list = False

autosummary_generate = True
autosummary_generate_overwrite = True
autosummary_imported_members = True

autodoc_class_signature = "separated"
autodoc_typehints = "signature"
autodoc_typehints_format = "short"
autodoc_preserve_defaults = True
autodoc_member_order = "bysource"
autodoc_inherit_docstrings = True
autodoc_default_options = {
    "show-inheritance": True,
    "members": False,
    "inherited-members": False,
    "undoc-members": False,
    "special-members": False,
}

# ------------
# HTML Options
# ------------

html_static_path = ["_static"]
html_theme = "pydata_sphinx_theme"
html_js_files = ["force_light.js"]
html_context = {"default_mode": "light"}
html_theme_options = {"navbar_end": ["navbar-icon-links"]}

# -------------------
# InterSphinx Options
# -------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
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

import os
import sys

import astropy.units as u
from bokeh.embed import file_html
from bokeh.layouts import Column, GridBox, Row
from bokeh.plotting import figure
from bokeh.resources import CDN

from ATK.utilities.docstrings import ATTR_DOCSTRINGS

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


def attach_units(app, what, name, obj, options, lines):
    if what != "class":
        return

    if not hasattr(obj, "_units"):
        return

    units_override = {u.deg: "deg", u.arcsec: "arcsec", u.arcmin: "arcmin"}

    lines.append("")
    lines.append(".. rubric:: Units")
    # not in use currently but doesn't hurt to leave
    lines.append(f".. _{obj.__name__}_Units:")
    lines.append("")
    lines.append("The following attributes are automatically converted to :class:`~astropy.units.Quantity` with a default unit unless one is explictly provided:")

    for attr, unit in obj._units.items():
        lines.append("")
        u_unit = u.Unit(unit)

        unit_str = None
        for k, v in units_override.items():
            if u_unit == k:
                unit_str = v
                break

        if unit_str is None:
            if u_unit == u.one:
                unit_str = "dimensionless"
            else:
                unit_str = u_unit.to_string("unicode")

        lines.append(f"- ``{attr}`` - {unit_str}")

    lines.append("")
    lines.append("|")


def attach_plotting_params(app, what, name, obj, options, lines):
    if what != "class":
        return

    if not hasattr(obj, "_plot_params"):
        return

    lines.append("")
    lines.append(".. rubric:: Plotting Arguments")
    # not in use currently but doesn't hurt to leave
    lines.append(f".. _{obj.__name__}_Plotting_Arguments:")
    lines.append("")
    lines.append("The following keyword arguments are accepted when plotting via :meth:`~ATK.Models.DataSet.plot()` or :meth:`~ATK.Models.DataSet.open`.")

    for param, info in obj._plot_params.items():
        lines.append("")
        lines.append(f"{param} : {info[0]}")
        info[1] = info[1].lstrip()
        for line in info[1].split("\n"):
            lines.append(f"   {line.lstrip().rstrip()}")

    lines.append("")
    lines.append("|")


def attach_data_methods(app, what, name, obj, options, lines):
    if what != "class":
        return

    if not hasattr(obj, "_data_methods") and not hasattr(obj, "_group_data_methods_doc"):
        return

    lines.append("")
    lines.append(".. rubric:: Data Methods")
    # not in use currently but doesn't hurt to leave
    lines.append(f".. _{obj.__name__}_Data_Methods:")
    lines.append("")
    lines.append(f"The following **Data Methods** are supported by :class:`~ATK.Models.{obj.__name__}` - either individually or through :meth:`DataSet.apply() <ATK.Models.DataSet.apply>`:")

    data_methods = getattr(obj, "_data_methods", ())
    group_data_methods = getattr(obj, "_group_data_methods_doc", {})
    all_data_methods = list(data_methods) + list(group_data_methods.keys())

    for method in all_data_methods:
        lines.append("")
        lines.append(f"- :meth:`~ATK.Models.{obj.__name__}.{method}`")
    lines.append("")
    lines.append("|")


def attach_plot_methods(app, what, name, obj, options, lines):
    if what != "class":
        return

    if not hasattr(obj, "_plot_methods_doc") and not hasattr(obj, "_group_plot_methods_doc"):
        return

    lines.append("")
    lines.append(".. rubric:: Plot Methods")
    # not in use currently but doesn't hurt to leave
    lines.append(f".. _{obj.__name__}_Plot_Methods:")
    lines.append("")
    lines.append(f"The following **Plot Methods** are supported by :class:`~ATK.Models.{obj.__name__}` through :meth:`DataSet.apply() <ATK.Models.DataSet.apply>`:")

    plot_methods = getattr(obj, "_plot_methods_doc", {})
    group_plot_methods = getattr(obj, "_group_plot_methods_doc", {})
    all_plot_methods = plot_methods | group_plot_methods

    for method, link in all_plot_methods.items():
        lines.append("")
        lines.append(f"- ``{method}`` (see :doc:`here <{link}>`)")
    lines.append("")
    lines.append("|")


def remove_self_from_signature(app, what, name, obj, options, signature, return_annotation):
    if what == "method" and signature:
        if signature.startswith("(self, "):
            signature = "(" + signature[len("(self, ") :]
        elif signature == "(self)":
            signature = "()"
    return signature, return_annotation


def common_attr_docstrings(app, what, name, obj, options, lines):
    """
    Replaces attribute docstring 'DOC_OVERRIDE' with a pre-existing docstring
    """

    if what != "attribute":
        return

    attr_name = name.split(".")[-1]

    if attr_name not in ATTR_DOCSTRINGS:
        return

    content = [l.strip() for l in lines if l.strip()]
    if "DOC_OVERRIDE" not in content:
        return

    doc = ATTR_DOCSTRINGS[attr_name].strip("\n")
    lines[:] = doc.splitlines()


def setup(app):
    app.connect("autodoc-process-docstring", attach_units)
    app.connect("autodoc-process-docstring", attach_data_methods)
    app.connect("autodoc-process-docstring", attach_plot_methods)
    app.connect("autodoc-process-docstring", attach_plotting_params)
    app.connect("autodoc-skip-member", autodoc_skip_member)
    app.connect("autodoc-process-signature", remove_self_from_signature)
    app.connect("autodoc-process-docstring", common_attr_docstrings)
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

nitpicky = False

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

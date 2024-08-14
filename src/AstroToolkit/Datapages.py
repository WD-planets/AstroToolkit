"""
This module assists in creating datapages, and provides additional panels for this purpose. It is greatly recommended to see the datapage examples alongside this section.
"""

from .Configuration.baseconfig import ConfigStruct

config = ConfigStruct()
config.read_config()


def buttons(radius="config", grid_size="config", pos=None, source=None):
    """buttons(source/pos, *)

    Generates Vizier and SIMBAD search buttons for a given target.

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees
    :param int, optional radius: search radius in arcseconds, default = config
    :param int, optional grid_size: datapage grid size used for scaling, default = config

    :return: Vizier and SIMBAD search buttons
    :rtype: bokeh figure

    |

    """

    from .DatapageElements.datapage_buttons import get_search_buttons
    from .Input.input_validation import check_inputs

    inputs = check_inputs(
        {"radius": radius, "pos": pos, "source": source, "grid_size": grid_size},
        "buttons",
    )
    radius, pos, source, grid_size = (
        inputs["radius"],
        inputs["pos"],
        inputs["source"],
        inputs["grid_size"],
    )

    config.read_config()
    if radius == "config":
        radius = float(config.datapage_search_button_radius)
    if grid_size == "config":
        grid_size = int(config.datapage_grid_size)

    return get_search_buttons(
        radius=radius, source=source, pos=pos, grid_size=grid_size
    )


def datatable(selection, source=None, pos=None, radius="config"):
    """datatable(selection,source/pos, *)
    Generates a datatable for a given target.

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees
    :param int, optional radius: search radius to use when fetching data from supported surveys in arcseconds, default = config
    :param dict<str> or dict<dict> selection: dict of datatable entries, with each entry taking a different format depending on the use case:

    To generate entires for supported surveys quickly and easily:

    .. code-block:: python

        survey: "default"

    where:

    :param str survey: supported survey, from:

        - gaia
        - panstarrs
        - skymapper
        - galex
        - sdss
        - wise
        - twomass (2MASS)
        - rosat

    |

    For customised entries within supported surveys:

    .. code-block:: python

        survey: {
                "parameters": parameters
                "errors": errors
                "notes": notes
                }

    where:

    :param str survey: supported survey, from:

        - gaia
        - panstarrs
        - skymapper
        - galex
        - rosat
        - sdss
        - wise
        - twomass (2MASS)
        - erosita

    :param list<str> parameters: names of parameters (i.e. column headers) that exist in chosen survey
    :param list<str> errors: names of errors on these parameters (i.e. column headers) that exist in chosen survey
    :param list<str> notes: any notes to include on each of the chosen parameters/errors

    |

    For customised entries outside of supported surveys:

    .. code-block:: python

        survey: {
                "parameters": parameters
                "values": values
                "errors": errors
                "notes": notes
                }

    where:

    :param str survey: name of survey to which the data belongs
    :param list<str> parameters: names of parameters
    :param list<float> values: values of these parameters
    :param list<float> errors: error values on these parameters
    :param list<str> notes: any notes to include on each of the chosen parameters/errors

    :return: populated datatable
    :rtype: bokeh figure

    |

    """

    from .DatapageElements.metadata_table import gettable
    from .Input.input_validation import check_inputs

    inputs = check_inputs(
        {"selection": selection, "source": source, "pos": pos, "radius": radius},
        "datatable",
    )
    selection, source, pos, radius = (
        inputs["selection"],
        inputs["source"],
        inputs["pos"],
        inputs["radius"],
    )

    config.read_config()

    if radius == "config":
        radius = float(config.datapage_datatable_radius)

    return gettable(selection=selection, source=source, pos=pos, radius=radius)


def gridsetup(dimensions, plots, grid_size="config"):
    """gridsetup(dimensions, plots, *)

    Assists in setting up panels for use in datapages.

    :param dict<int> dimensions: dimensions of plot in grid units, in format:

    .. code-block:: python

        'width': width
        'height': height

    :param list<dict> plots: plots for use in datapage, with each entry in format:

    .. code-block:: python

        'name': name
        'figure': figure
        'width': width
        'height': height

    where:

    :param str name: name to assign to the figure
    :param class or bokeh figure or None figure: The figure to give this panel, this can either be an ATK structure that supports plotting (lightcurve, image, etc.) or a Bokeh figure. If filling empty space, None can instead be passed to create a blank panel.
    :param int width: width of panel in grid units
    :param int height: height of panel in grid units

    :return: dictionary of plots with keys = names assigned to each panel
    :rtype: dict

    |

    """
    from .Input.input_validation import check_inputs
    from .Misc.gridsetup import format_grid_plots

    inputs = check_inputs(
        {"dimensions": dimensions, "plots": plots, "grid_size": grid_size}, "gridsetup"
    )
    dimensions, plots, grid_size = (
        inputs["dimensions"],
        inputs["plots"],
        inputs["grid_size"],
    )

    config.read_config()
    if grid_size == "config":
        grid_size = int(config.datapage_grid_size)

    return format_grid_plots(dimensions, plots, grid_size)

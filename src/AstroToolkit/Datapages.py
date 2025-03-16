"""
This module assists in creating datapages, and provides additional panels for this purpose. It is greatly recommended to see the datapage examples alongside this section.
"""

from bokeh.models import Button, DataTable

from .Configuration.baseconfig import ConfigStruct
from .Input.input_validation import check_inputs

config = ConfigStruct()
config.read_config()

newline = "\n"


def buttons(radius: float = "config", grid_size: int = "config", pos: list[float] = None, source: int = None) -> Button:
    """buttons(source/pos, **kwargs)

    Generates Vizier and SIMBAD search buttons for a given target.

    :param source: Target GAIA DR3 Source ID
    :type source: int
    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>
    :param radius: Search radius in arcseconds, default given by :ref:`datapage_search_button_radius <cfg_datapage_search_button_radius>` config key
    :type radius: float, optional
    :param grid_size: datapage grid size used for scaling, default given by :ref:`datapage_grid_size <cfg_datapage_grid_size>` config kkey
    :type grid_size: int, optional

    :return: Vizier and SIMBAD search buttons
    :rtype: bokeh figure

    |

    """

    from .DatapageElements.datapage_buttons import get_search_buttons

    if config.enable_notifications:
        print(f"Generating datapage SIMBAD/Vizier buttons{newline}")

    corrected_inputs = check_inputs(
        {"radius": [radius, float], "pos": [pos, list], "source": [source, int], "grid_size": [grid_size, int]},
        "buttons",
    )
    radius, pos, source, grid_size = corrected_inputs

    config.read_config()
    if radius == "config":
        radius = float(config.datapage_search_button_radius)
    if grid_size == "config":
        grid_size = int(config.datapage_grid_size)

    buttons = get_search_buttons(radius=radius, source=source, pos=pos, grid_size=grid_size)
    return buttons


def datatable(entries: dict, source: int = None, pos: list[float] = None, radius: float = "config") -> DataTable:
    """datatable(entries, source/pos, **kwargs)
    Generates a datatable for a given target.

    :param source: Target GAIA DR3 Source ID
    :type source: int
    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>
    :param radius: search radius to use when fetching data from supported surveys in arcseconds, default given by `datapage_datatable_radius <cfg_datapage_datatable_radius>` config key
    :type radius: float, optional
    :param entries: datatable entry or list of datatable entries
    :type entries: dict/list<dict>

    Each datatable entry in **entries** can take three different forms.

    :return: ATK Datatable
    :rtype: Bokeh figure

    |

    .. rubric:: Default Entries
        :heading-level: 3

    For data surveys that are supported by ATK, default entries are available. These contain key astrometry and photometry. To use these default entries:

    .. code-block:: python

        entry = {
                "kind": "atk_defaults",
                "surveys": <surveys>,
                }

    where:

    :param surveys: list of :ref:`supported ATK data surveys <Data Surveys>` to include in the datatable
    :type surveys: list<str>

    |

    .. rubric:: Vizier Entries
        :heading-level: 3

    To generate a table entry using any catalogue available in Vizier:

    .. code-block:: python

        entry = {
                "kind": "vizier",
                "survey": <survey>,
                "parameters": <parameters>,
                "errors": <errors>,
                "notes": <notes",
                }

    where:

    :param survey: any Vizier catalogue ID, or the name of a user-defined catalogue :ref:`alias <Adding Catalogue Aliases>`
    :type survey: str
    :param parameters: names of parameters (i.e. Vizier column headers) that exist in chosen catalogue
    :type parameters: list<str>
    :param errors: names of errors on these parameters (i.e. Vizier column headers) that exist in chosen catalogue
    :type errors: list<str>
    :param notes: any notes to include on each of the chosen parameters/errors
    :type notes: list<str>

    |

    .. rubric:: Custom Entries
        :heading-level: 3

    To include entries using external data:

    .. code-block:: python

        entry = {
                "kind": "external"
                "survey": <survey>,
                "parameters": <parameters>,
                "values": <values>,
                "errors": <errors>,
                "notes": <notes>,
                }

    where:

    :param survey: name of source from which data originates
    :type survey: str
    :param parameters: names of parameters
    :type parameters: list<str>
    :param values: values of these parameters
    :type values: list<float>
    :param errors: errors on these parameters
    :type errors: list<float>
    :param notes: any notes to include on each of the chosen parameters/errors
    :type notes: list<str>

    |

    """

    if config.enable_notifications:
        print(f"Generating datapage datatable{newline}")

    from .DatapageElements.metadata_table import gettable

    corrected_inputs = check_inputs(
        {"selection": [entries, list], "source": [source, int], "pos": [pos, list], "radius": [radius, float]},
        "datatable",
    )
    entries, source, pos, radius = corrected_inputs

    config.read_config()

    if radius == "config":
        radius = float(config.datapage_datatable_radius)

    return gettable(selection=entries, source=source, pos=pos, radius=radius)


def datapage(dimensions: dict, panels: list[dict], grid_size: int = "config", layout: list = None) -> dict:
    """datapage(dimensions, panels, layout, **kwargs)

    Assists in generating datapages from a set of figures or ATK data structures by adjusting scaling, sizing, style and layout. See the :ref:`Creating a Datapage` tutorial for an example.

    :param dimensions: dimensions of datapage, in the form:
    :type dimensions: dict

    .. code-block:: python

        dimensions = {
                     'width': width (grid units)
                     'height': height (grid units)
                     }

    :param plots: list of plot entries with each entry in format:
    :type plots: list<dict>

    .. code-block:: console

        entry = {
                'name': name to assign to the figure
                'figure': The figure to give this panel, this can either be an ATK structure that supports plotting (lightcurve, image, etc.)
                          or a Bokeh figure. If filling empty space, None can instead be passed to create a blank panel
                'width': width of panel in grid units
                'height': height of panel in grid units
                }

    :param layout: layout of grid, e.g.:
    :type layout: list<list>

    .. code-block:: console

        layouts = [
            [plot0, plot1, [plot2, plot3]],
            [plot4, plot5, plot6],
            [plot7, plot8, plot9],
        ]

    Which would generate the following datapage (assuming all panel sizes were set correctly):

    .. code-block:: console

        +---------+---------+---------+
        |         |         |  plot2  |
        |  plot0  |  plot1  +---------+
        |         |         |  plot3  |
        +---------+---------+---------+
        |         |         |         |
        |  plot4  |  plot5  |  plot6  |
        |         |         |         |
        +---------+---------+---------+
        |         |         |         |
        |  plot7  |  plot8  |  plot9  |
        |         |         |         |
        +---------+---------+---------+

    :return: :class:`Datapage <AstroToolkit.Misc.grid.Datapage>`

    |

    """

    from .Misc.grid import format_grid_plots

    corrected_inputs = check_inputs(
        {"dimensions": [dimensions, dict], "plots": [panels, list], "grid_size": [grid_size, int]}, "datapage"
    )
    dimensions, plots, grid_size = corrected_inputs

    config.read_config()
    if grid_size == "config":
        grid_size = int(config.datapage_grid_size)

    return format_grid_plots(dimensions, plots, grid_size, layout)

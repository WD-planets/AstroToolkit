from typing import Sequence

from bokeh.models import GridBox

from ..plotting.datapage.plot_datapage import get_datapage
from ..structures.DataSet import DataSet


def grid(layout: Sequence[Sequence[DataSet]]) -> GridBox:
    """
    Combines any number of plottable :class:`~ATK.Models.DataSet` objects into a neat grid of figures (see :doc:`here </auto_tutorials/datapages/datapage_plotting>`).

    Parameters
    ----------
    layout : Sequence of Sequence of :class:`~ATK.Models.DataSet`
        Each internal list forms a **row** of the output grid.

    Examples
    --------
    A 5x6 grid layout containing **three rows**:

    - 1st row = 2x2 :class:`~ATK.Models.Image` and 3x2 :class:`~ATK.Models.Lightcurve`
    - 2nd row = 2x2 :class:`~ATK.Models.HRD` and 3x2 :class:`~ATK.Models.Spectrum`
    - 3rd row = 5x2 :class:`~ATK.Models.DataTable`

    |

    Should be passed as:

    .. code-block:: python

        layout = [[image, image, lc    , lc  , lc   ],
                  [image, image, lc    , lc  , lc   ],
                  [hrd  , hrd  , spec , spec , spec ],
                  [hrd  , hrd  , spec , spec , spec ],
                  [table, table, table, table, table],
                  [table, table, table, table, table]]

    where:

    - ``image`` is a :class:`~ATK.Models.DataSet` containing :class:`~ATK.Models.Image` objects
    - ``lc`` is a :class:`~ATK.Models.DataSet` containing :class:`~ATK.Models.Lightcurve` objects
    - ``hrd`` is a :class:`~ATK.Models.DataSet` containing :class:`~ATK.Models.HRD` objects
    - ``spec`` is a :class:`~ATK.Models.DataSet` containing :class:`~ATK.Models.Spectrum` objects
    - ``table`` is a :class:`~ATK.Models.DataSet` containing :class:`~ATK.Models.DataTable` objects
    """

    return get_datapage(layout)

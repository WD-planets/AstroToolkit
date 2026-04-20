import warnings

from astropy.coordinates import SkyCoord

from ..queries.arguments import get_query_arguments
from ..queries.checksum import make_cache_key
from ..queries.query_core import general_query, setup_targeting
from ..structures.DataSet import DataSet
from ..structures.Target import Target
from ..utilities.defaults import RETURNS


def query(kind: str, targets: int | SkyCoord | Target | list[int | SkyCoord | Target], **kwargs) -> DataSet:
    """
    Performs a query to retrieve astronomical data.

    Parameters
    ----------
    kind : {'vizier', 'image', 'lightcurve', 'spectrum', 'sed', 'hrd'}, optional
        Type of query to perform.

    targets : int, :class:`~astropy.coordinates.SkyCoord`, :class:`~ATK.Models.Target`, or iterable of these
        Targets of search (see :doc:`here </auto_tutorials/getting_started/data_query>`).

    path : str or :class:`~pathlib.Path`, optional
        Path of local FITS file to which result should be saved. :func:`~ATK.Tools.query` will attempt to read this file instead of retrieving new data if arguments remain unchanged.

    corrections : bool, optional
        If ``True``, perform automatic proper motion correction (where possible, see :doc:`here </auto_tutorials/getting_started/data_query>`).

        Default is ``True``.

    **kwargs
        Arguments specific to each ``kind``. See below.


    Returns
    -------
    :class:`~ATK.Models.DataSet`
        Dataset containing data matching ``kind``.


    Notes
    -----
    The keyword arguments that are available to :func:`~ATK.Tools.query` depend on the ``kind`` of data that is being requested. These are outlined below:

    |

    ``kind='vizier'``

    survey : str
        `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue or `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue alias (see :doc:`here </auto_tutorials/extension/vizier>`).

    radius : float or :class:`~astropy.units.Quantity`, optional
        Search radius around each target. If float, unit taken from ``query_settings.default_scale`` config key (see :doc:`here </auto_tutorials/configuration/config>`).

        Default taken from ``query_settings.query_radius`` config key
        (see :doc:`here </auto_tutorials/configuration/config>`).

    |

    ``kind='image'``

    survey : str
        Target survey (see :doc:`here </auto_tutorials/images/image_query>`).

    band : str
        Photometric band of image (see :doc:`here </auto_tutorials/images/image_query>`).

    size : float or :class:`~astropy.units.Quantity`, optional
        Width and height of image. If float, unit taken from ``query_settings.default_scale`` config key (see :doc:`here </auto_tutorials/configuration/config>`).

        Default taken from ``query_settings.image_size`` config key
        (see :doc:`here </auto_tutorials/configuration/config>`).

    overlays : list of str, optional
        List of detection overlays to generate (see :doc:`here </auto_tutorials/images/image_query_2>`).

        Default is ``None``

    |

    ``kind='lightcurve'``

    survey : str
        Target survey (see :doc:`here </auto_tutorials/lightcurves/lightcurve_query>`).

    radius : float or :class:`~astropy.units.Quantity`, optional
        Search radius around each target. If float, unit taken from ``query_settings.default_scale`` config key (see :doc:`here </auto_tutorials/configuration/config>`).

        Default taken from ``query_settings.query_radius`` config key
        (see :doc:`here </auto_tutorials/configuration/config>`).

    split : bool, optional
        If ``True``, photometry is separated by a per-survey ID.

        Default is ``False``

    filter : bool, optional
        If ``True``, optional data quality filtering is performed.

        Default is ``True``

    username : str, optional
        ATLAS forced photometry username, only required in ATLAS queries.

    password : str, optional
        ATLAS forced photometry password, only required in ATLAS queries.

    |

    ``kind='spectrum'``

    survey : str
        Target survey (see :doc:`here </auto_tutorials/spectra/spectrum_query>`).

    radius : float or :class:`~astropy.units.Quantity`, optional
        Search radius around each target. If float, unit taken from ``query_settings.default_scale`` config key (see :doc:`here </auto_tutorials/configuration/config>`).

        Default taken from ``query_settings.query_radius`` config key
        (see :doc:`here </auto_tutorials/configuration/config>`).

    |

    ``kind='sed'``

    radius : float or :class:`~astropy.units.Quantity`, optional
        Search radius around each target. If float, unit taken from ``query_settings.default_scale`` config key (see :doc:`here </auto_tutorials/configuration/config>`).

        Default taken from ``query_settings.query_radius`` config key
        (see :doc:`here </auto_tutorials/configuration/config>`).

    |

    ``kind='hrd'``

    colour : str, optional
        Gaia colour in format ``band1-band2``, where each band is a photometry column in `Vizier <https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3>`_.

        Default is ``'BPmag-RPmag'``.

    mag : str, optional
        Gaia absolute magnitude, as listed in `Vizier <https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3>`_.

        Default is ``'Gmag'``.

    |

    ``kind='datatable'``

    columns : dict
        Table columns in format ``survey: cols``, where survey is a `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue or `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue alias (see :doc:`here </auto_tutorials/extension/vizier>`) and cols is a list of columns in that catalogue.

        Example: ``{'gaia': ['Gmag', 'BPmag', 'RPmag'], 'II/328/allwise': ['W1mag', 'W2mag', 'W3mag', 'W4mag']}``

    See Also
    --------
    :class:`~ATK.Models.Record` : Vizier data container, stored in :class:`~ATK.Models.DataSet`\ s returned by ``vizier`` queries.
    :class:`~ATK.Models.Image` : Image data container, stored in :class:`~ATK.Models.DataSet`\ s returned by ``image`` queries.
    :class:`~ATK.Models.Lightcurve` : Light curve data container, stored in :class:`~ATK.Models.DataSet`\ s returned by ``lightcurve`` queries.
    :class:`~ATK.Models.Spectrum` : Spectrum data container, stored in :class:`~ATK.Models.DataSet`\ s returned by ``spectrum`` queries.
    :class:`~ATK.Models.SED` : SED data container, stored in :class:`~ATK.Models.DataSet`\ s returned by ``sed`` queries.
    :class:`~ATK.Models.HRD` : HRD data container, stored in :class:`~ATK.Models.DataSet`\ s returned by ``hrd`` queries.
    :class:`~ATK.Models.DataTable` : Table data container, stored in :class:`~ATK.Models.DataSet`\ s returned by ``datatable`` queries.
    """

    arguments = get_query_arguments(kind, kwargs)

    if targets is None:
        raise ValueError("No targets specified for query.")

    # get flattened list of targets
    targets = setup_targeting(targets)
    if targets is RETURNS.NULL:
        raise ValueError("Query received no targets, likely due to no valid targets being provided.")

    if targets is RETURNS.EXCEPTION:
        raise Exception("Failed to generate requested targets, this is likely due to a Vizier fault.")

    for target in targets:
        target.radius = arguments.get("radius")

    # targets may return exception if Vizier is down
    if any(target is RETURNS.EXCEPTION for target in targets):
        warnings.warn("Failed to generate requested targets, this is likely due to a Vizier fault.")

        structure = DataSet[kind](
            kind=None, survey=arguments.get("survey", None), targets=None, radius=arguments.get("radius", None), exception=True
        )

        return structure

    # disable proper motion correction
    if not arguments.get("corrections", True):
        for target in targets:
            initial_coords = target.initial_coords
            target.initial_coords = SkyCoord(
                ra=initial_coords.ra, dec=initial_coords.dec, frame=initial_coords.frame, obstime=initial_coords.obstime
            )
            target.coords = initial_coords
            target.identifier = None
            target.survey = None
            target.correction = "none"

    structure = general_query(kind, targets, **arguments)

    # attach cache key
    structure._cache_key = make_cache_key(kind, targets, arguments)

    if structure.exception:
        return structure

    # save if path provided
    if arguments.get("path"):
        structure.store(arguments["path"])

    return structure

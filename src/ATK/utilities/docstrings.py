import textwrap

# fmt: off

DOCSTRINGS = {
    "bin": 
        """
        Bins all array-like attributes into a number of bins or a bin size in x (i.e. {x}).
        
        Exactly one of ``bin`` or ``size`` must be provided.

        Parameters
        ----------
        bins : int, optional
            Number of bins in which to bin data.

        size : float or :class:`~astropy.units.Quantity`, optional
            Size of each bin in which to bin data. If a :class:`~astropy.units.Unit` is not provided, ``size`` is assumed to be in the same unit as x.

        inplace : bool, optional
            If ``True``, modify the current :class:`~ATK.Models.{name}` inplace. If ``False``, operate on and return a copy - leaving the original unchanged.

        Returns
        -------
        ``Self``
            The binned :class:`~ATK.Models.{name}`. Returns ``self`` if ``inplace=True``, otherwise returns a new instance.
        """,

    "clip":
        """
        Sigma clips a :class:`~ATK.Models.{name}` to a given sigma (or range in sigma) in y (i.e. {y}).

        This method utilises :func:`astropy.stats.sigma_clip`, see the Astropy documentation for details.

        Parameters
        ----------
        sigma : float
            Number of standard deviations to use for clipping.

        sigma_lower : float or None, optional
            Lower bound for clipping. If ``None``, defaults to ``sigma``.

        sigma_upper : float or None, optional
            Upper bound for clipping. If ``None``, defaults to ``sigma``.

        inplace : bool, optional
            If ``True``, modify the current :class:`~ATK.Models.{name}` inplace. If ``False``, operate on and return a copy - leaving the original unchanged.

        Returns
        -------
        ``Self``
            The clipped :class:`~ATK.Models.{name}`. Returns ``self`` if ``inplace=True``, otherwise returns a new instance.
        """,

    "crop":
        """
        Crops all array-like attributes to a given range in x (i.e. {x}).

        At least one of ``min``, ``max`` must be provided.
        
        Parameters
        ----------
        min : float or :class:`~astropy.units.Quantity`, optional
            Minimum x value. If not provided, bottom range is not clipped.

            If a :class:`~astropy.units.Unit` is not provided, ``min`` is assumed to be in the same unit as x.

        max : float or :class:`~astropy.units.Quantity`, optional. 
            Maximum x value. If not provided, top range is not clipped.

            If a :class:`~astropy.units.Unit` is not provided, ``max`` is assumed to be in the same unit as x.


        inplace : bool, optional
            If ``True``, modify the current :class:`~ATK.Models.{name}` inplace. If ``False``, operate on and return a copy - leaving the original unchanged.

        Returns
        -------
        ``Self``
            The cropped :class:`~ATK.Models.{name}`. Returns ``self`` if ``inplace=True``, otherwise returns a new instance.
        """,
    
    "to_hdu":
        """
        Converts structure into a FITS :class:`~astropy.io.fits.{hdu_type}`.

        Returns
        -------
        :class:`~astropy.io.fits.{hdu_type}`
        """,

    "from_dataframe":
        """
        Construct a :class:`~ATK.Models.{obj}` from a :class:`~pandas.DataFrame`.

        Parameters
        ----------
        target : :class:`~ATK.Models.Target`, int, or :class:`~astropy.coordinates.SkyCoord`
            Astronomical target with which to associate input data. 

            Can be a :class:`~ATK.Models.Target`, a Gaia Source ID (``int``), or a :class:`~astropy.coordinates.SkyCoord`.

        data : :class:`~pandas.DataFrame`
            Tabular data containing the relevant fields (i.e. array-like attributes) required to
            construct a :class:`~ATK.Models.{obj}`. 
            Units are assumed to be as listed above.

        **kwargs
            Additional keyword arguments to be forwarded to the internal parser.

            The following keyword arguments are required: {args}

        Returns
        -------
        :class:`~ATK.Models.{obj}`
        """,

    "from_table":
        """
        Construct a :class:`~ATK.Models.{obj}` from a :class:`~astropy.table.Table`.

        Parameters
        ----------
        target : :class:`~ATK.Models.Target`, int, or :class:`~astropy.coordinates.SkyCoord`
            Astronomical target with which to associate input data. 

            Can be a :class:`~ATK.Models.Target`, a Gaia Source ID (``int``), or a :class:`~astropy.coordinates.SkyCoord`.

        data : :class:`~astropy.table.Table`
            Tabular data containing the relevant fields (i.e. array-like attributes) required to
            construct a :class:`~ATK.Models.{obj}`. 
            Units are taken from the table where available, with missing units assumed to be those listed above.

        **kwargs
            Additional keyword arguments to be forwarded to the internal parser.

            The following keyword arguments are required: {args}

        Returns
        -------
        :class:`~ATK.Models.{obj}`
        """,

    "to_table": 
        """
        Combines all array-like attributes of a structure into a :class:`~astropy.table.Table`, preserving units.

        Returns
        -------
        :class:`~astropy.table.Table`
        """,

    "to_dataframe":
        """
        Combines all array-like attributes of a structure into a :class:`~pandas.DataFrame`.

        Units are not preserved.

        Returns
        -------
        :class:`~pandas.DataFrame`
        """,

    "show":
        """show(self, show_types = False, show_all = False)
        Prints structure to stdout in a human-readable format.

        Parameters
        ----------
        show_types : bool, optional
            If True, print data types of structure attributes.

            Default is ``False``

        show_all : bool, optional
            If True, do not truncate printing of large iterables.

            Default is ``False``.

        Returns
        -------
        ``self``
        """
}

attr_docstrings = {
    "survey": "Survey from which the stored data originates.",
    "band": "Photometric band of stored data.",
    "correction":
        """
        Achieved degree of proper motion correction.
        
        - ``'full'`` = complete 3-dimensional projection on the sky.
        - ``'partial'`` = 2-dimensional plane projection.
        - ``'none'`` = no correction.
        """,
    "search_pos": "Position of search at time of execution (i.e. post-correction).",
    "separation": "Separation between position of the search and the returned data.",
    "overlay":
        """
        Photometric overlay data.

        Each row corresponds to a single detection, and contains:
        
        - ``survey`` : str
            Survey name (e.g. ``'gaia'``).

        - ``ra`` : float
            Right ascension (deg).

        - ``dec`` : float
            Declination (deg).

        - ``pm_ra_cosdec`` : float
            Proper motion in RA direction (mas/yr).

        - ``pm_dec`` : float
            Proper motion in Dec direction (mas/yr).

        - ``dist`` : float
            Distance to source, as derived from parallax (pc).

        - ``correction`` : {'full', 'partial', 'none'}
            Achieved level of proper motion correction.

        - ``mag_name`` : str
            Name of magnitude band.

        - ``mag`` : float
            Apparent magnitude.

        - ``err_name`` : str
            Name of magnitude uncertainty field.

        - ``err`` : float
            Magnitude uncertainty.

        - ``simbad_id`` : str or None
            Primary name in `SIMBAD <https://simbad.cds.unistra.fr/simbad/>`_ (``None`` if no match was found).
        """,
    
    "table":
        """
        Target-matched subset of one or more `Vizier <https://vizier.cds.unistra.fr/>`_ catalogues, combined into a single :class:`~astropy.table.Table`.

        Each row contains:

        - ``survey`` : str
            Alias to `Vizier <https://vizier.cds.unistra.fr/>`_ ``catalogue``.

        - ``catalogue`` : str
            Identifier of the source catalogue providing the measurement, (e.g. ``'I/355/gaiadr3'`` for Gaia DR3, ``'II/349/ps1'`` for Pan-STARRS).

        - ``correction`` : {'full', 'partial', 'none'}
            Achieved level of proper motion correction.

        - ``parameter`` : str
            `Vizier <https://vizier.cds.unistra.fr/>`_ column name, (e.g. ``'Gmag'``, ``'BPmag'``, ``'rmag'``).

        - ``value`` : float
            Value in column given by ``parameter``.

        - ``unit`` : str
            Unit of ``value``.

        - ``separation`` : float
            Separation between position of the search and the returned `Vizier <https://vizier.cds.unistra.fr/>`_ row.
        """,
}

# Remove leading and trailing whitespace (allows indentation for clarity when updating docstrings)
ATTR_DOCSTRINGS = {attr: textwrap.dedent(doc).strip() for attr,doc in attr_docstrings.items()}

"""
Might eventually need to update above to take an attribute name (e.g. #: DOC_OVERRIDE image_overlay), otherwise I am slowly limiting the attribute names I can use (although in this case, only if I want a separate override for 'overlay', which is unlikely)
"""

def get_docstring(func: str, **kwargs):
    base_str = DOCSTRINGS[func]

    # strict check: extract fields
    import string

    fields = [f for _, f, _, _ in string.Formatter().parse(base_str) if f]

    missing = set(fields) - set(kwargs)
    extra = set(kwargs) - set(fields)

    if missing:
        raise ValueError(f"Missing args: {missing}")
    if extra:
        raise ValueError(f"Unexpected args: {extra}")

    return base_str.format(**kwargs)

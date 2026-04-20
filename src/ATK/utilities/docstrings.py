# fmt: off

DOCSTRINGS = {
    "bin": 
        """
        Bins all array-like attributes into a number of bins or a bin size in ``{x}``.
        
        Exactly one of ``bin`` or ``size`` must be provided.

        Parameters
        ----------
        bins : int, optional
            Number of bins in which to bin data.

        size : :class:`~astropy.units.Quantity` or float, optional
            Size of each bin in which to bin data. If a unit is not provided, ``size`` is assumed to be in days.

        inplace : bool, optional
            If ``True``, operate on and return a copy of the :class:`~ATK.Models.{name}` - leaving the original unchanged.
        """,

    "clip":
        """
        Sigma clips a :class:`~ATK.Models.{name}` to a given sigma (or range in sigma) in brightness.

        This method utilises Astropy's :func:`~astropy.stats.sigma_clip`, see the Astropy documentation for details.

        Parameters
        ----------
        sigma : float
            Number of standard deviations to use for clipping.

        sigma_lower : float or None, optional
            Lower bound for clipping. If ``None``, defaults to ``sigma``.

        sigma_upper : float or None, optional
            Upper bound for clipping. If ``None``, defaults to ``sigma``.

        inplace : bool, optional
            If ``True``, operate on and return a copy of the :class:`~ATK.Models.{name}` - leaving the original unchanged.
        """,

    "crop":
        """
        Crops all array-like attributes to a given range in {x}.

        At least one of ``min``, ``max`` must be provided.
        
        Parameters
        ----------
        min : float, optional
            Minimum {x} value. If not provided, bottom range is not clipped.

        max : float, optional. 
            Maximum {x} value. If not provided, top range is not clipped.

        inplace : bool, optional
            If ``True``, operate on and return a copy of the :class:`~ATK.Models.{name}` - leaving the original unchanged.
        """
}


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

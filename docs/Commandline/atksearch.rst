ATKsearch
=========

This tool calls :func:`search() <AstroToolkit.Tools.search>`, searching for a given position or Gaia source in Vizier/SIMBAD.

.. rubric:: Usage 
    :heading-level: 2

.. code-block:: console

    ATKsearch [-h] {simbad,vizier} (--source SOURCE | --pos RA DEC) [-r radius]

**Positional Arguments:**

    .. code-block:: console
    
        {simbad,vizier}  Where to perform search

**Optional Arguments:**

    .. code-block:: console

        -h, --help       Show help message and exit
        --source SOURCE  Gaia DR3 Source ID
        --pos RA DEC     Position in degrees
        -r RADIUS        Radius of search in arcseconds

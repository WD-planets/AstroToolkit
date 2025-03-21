ATKquery
========

This tool runs an ATK query as if :func:`query() <AstroToolkit.Tools.query>` was called. The standard :ref:`data structure <Data Structures>` methods are available as "jobs" for any returned data.

.. rubric:: Usage 1 
    :heading-level: 2
    
.. code-block:: console

    ATKquery {data,reddening,lightcurve,spectrum} [-h] <survey> (--source SOURCE | --pos RA DEC) [-r RADIUS]


Performs a data, reddening, lightcurve or spectrum :func:`query() <AstroToolkit.Tools.query>` for a given position or source in a chosen <survey>.

**Positional Arguments**:

    .. code-block:: console

        {data,reddening,lightcurve,spectrum}  Kind of query to perform
        survey  Target survey

**Optional Arguments:**

    .. code-block:: console

        -h, --help       Show help message and exit
        --source SOURCE  Gaia DR3 Source ID
        --pos RA DEC     Position in degrees
        -r RADIUS        Radius of search in arcseconds

|

.. rubric:: Usage 2 
    :heading-level: 2

.. code-block:: console
        
    ATKquery lightcurve [-h] <survey> (--source SOURCE | --pos RA DEC) [-r RADIUS] [--username USERNAME] [--password PASSWORD]


Performs a lightcurve :func:`query() <AstroToolkit.Tools.query>` for a given position or source in a chosen <survey>.

**Positional Arguments:**

    .. code-block:: console

        survey               Target survey

**Optional Arguments:**

    .. code-block:: console

        -h, --help           Show help message and exit
        --source SOURCE      Gaia DR3 Source ID
        --pos RA DEC         Position in degrees
        -r RADIUS            Radius of search in arcseconds
        --username USERNAME  ATLAS username (only needed in ATLAS queries)
        --password PASSWORD  ATLAS password (only needed in ATLAS queries)

|

.. rubric:: Usage 3 
    :heading-level: 2

.. code-block:: console

    ATKquery image [-h] <survey> (--source SOURCE | --pos RA DEC) [-s SIZE]



Performs an image :func:`query() <AstroToolkit.Tools.query>` for a given position or source in a chosen <survey>.

**Positional Arguments:**
    
    .. code-block:: console

        survey           Target survey

**Optional Arguments:**

    .. code-block:: console

        -h, --help       Show help message and exit
        --source SOURCE  Gaia DR3 Source ID
        --pos RA DEC     Position in degrees
        -s SIZE          Size of image in arcseconds

|

.. rubric:: Usage 4
    :heading-level: 2

.. code-block:: console
    
    ATKquery {bulkdata,sed} [-h] (--source SOURCE | --pos RA DEC) [-r RADIUS]


Performs a bulkdata or SED :func:`query() <AstroToolkit.Tools.query>` for a given position or source.

**Optional Arguments:**

    .. code-block:: console

        -h, --help       Show help message and exit
        --source SOURCE  Gaia DR3 Source ID
        --pos RA DEC     Position in degrees
        -r RADIUS        Radius of search in arcseconds

|

.. rubric:: Usage 5 
    :heading-level: 2

.. code-block:: console

    ATKquery hrd [-h] sources [sources ...]

Performs a hrd :func:`query() <AstroToolkit.Tools.query>` for a given set of sources.

**Positional Arguments:**

    .. code-block:: console

        sources     Sources to overlay

**Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit

|


**Note:**
 - This tool is provided for quick command-line queries, and is therefore more limited in functionality compared to queries performed through either the :ref:`ATK GUI <Gui>` or scripts.

ATKread
=======

This tool calls :func:`readdata() <AstroToolkit.Tools.readdata>`, reading a local ATK data file. The standard :ref:`data structure <Data Structures>` methods are available as "jobs" for any returned data.


.. rubric:: Usage 
    :heading-level: 2
        
.. code-block:: console

    ATkread [-h] [fname]

Reads a file at the location <fname>. If no file name is provided, a file dialogue will open in which a file may be selected.

**Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit
        fname       File path

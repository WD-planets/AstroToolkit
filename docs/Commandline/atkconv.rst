ATKconv
=======

This tool calls :func:`deg2hms() <AstroToolkit.Tools.deg2hms>` or :func:`hms2deg() <AstroToolkit.Tools.hms2deg>` tools, converting coordinates in degrees to HMS±DMS format and vice versa.

.. rubric:: Usage 
    :heading-level: 2

.. code-block:: console

    ATKconv [-h] Position [Position ...]

**Positional Arguments:**

    .. code-block:: console

        Position    2 Arguments:
                      [RA] [DEC] in degrees
                    1 Argument:
                      [POSITION] in HHMMSS.SS±DDMMSS.SS format

**Optional Arguments:**
    
    .. code-block:: console

        -h, --help  Show help message and exit

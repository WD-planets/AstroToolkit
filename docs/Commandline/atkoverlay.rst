ATKoverlay
==========

This tool allows for the various :ref:`Overlays <Overlays>` tools to be utilised from the command line. The first argument sets the job to perform on the overlay definition list. See :ref:`tutorials <Defining New Overlays>` for examples. 

.. rubric:: Usage 1
    :heading-level: 2

.. code-block:: console

    ATKoverlay open [-h]

Calls :func:`openOverlays() <AstroToolkit.Overlays.openOverlays>`, opening the overlay definition list in the default text editor.

**Optional Arguments:**

    .. code-block:: console

        -h, --help Show help message and exit

|

.. rubric:: Usage 2
    :heading-level: 2

.. code-block:: console

    ATKoverlay add [-h] <survey> --ra RA_NAME --dec DEC_NAME --id ID_NAME [--mags MAG_NAMES [MAG_NAMES ...]]

Calls :func:`addOverlay() <AstroToolkit.Overlays.addOverlay>`, adding an overlay definition for a given catalogue alias.

**Positional Arguments:**

    .. code-block:: console

         alias                            Overlay catalogue lias

**Arguments:**
    
    .. code-block:: console

        --ra RA_NAME                      Name of right ascension column in Vizier
        --dec DEC_NAME                    Name of declination column in Vizier
        --id ID_NAME                      Name of ID column in Vizier

**Optional Arguments:**

    .. code-block:: console

        -h, --help                        Show help message and exit
        --mags MAG_NAMES [MAG_NAMES ...]  Names of magnitude columns in Vizier. If provided,
                                          a 'scaled_detection' overlay defintion is generated.
                                          Otherwise, a 'detection' overlay is generated.

|

.. rubric:: Usage 3 
    :heading-level: 2

.. code-block:: console
        
    ATKoverlay del [-h] {scaled_detection,detection} <alias>

Calls :func:`delOverlay() <AstroToolkit.Overlays.delOverlay>`, deleting an overlay definition for catalogue alias <survey>. 

**Positional Arguments:**

    .. code-block:: console

        {scaled_detection,detection}  Kind of overlay definition to delete
        <alias>                       Overlay alias

**Optional Arguments:**

    .. code-block:: console

        -h, --help                    Show help message and exit

|

.. rubric:: Usage 4 
    :heading-level: 2

.. code-block:: console

    ATKoverlay reset [-h]

Calls :func:`resetOverlays() <AstroToolkit.Overlays.resetOverlays>`, resetting the overlay definition list to its default state.

**Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit

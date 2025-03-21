ATKconfig
=========

This tool allows for the various :ref:`Config <Config>` tools to be utilised from the command line. The first argument sets the job to perform on the config. See :ref:`tutorials <The Config>` for examples.

.. rubric:: Usage 1 
    :heading-level: 2

.. code-block:: console
        
    ATKconfig show [-h]

Calls :func:`showConfig() <AstroToolkit.Config.showConfig>`, printing the current state of the config to stdout.

**Optional Arguments:**

    .. code-block:: console

        -h, --help Show help message and exit

|

.. rubric:: Usage 2 
    :heading-level: 2

.. code-block:: console

    ATKconfig open [-h]

Calls :func:`openConfig() <AstroToolkit.Config.openConfig>`, opening the config in the default text editor.

**Optional Arguments:**

    .. code-block:: console

        -h, --help Show help message and exit

|

.. rubric:: Usage 3 
    :heading-level: 2

.. code-block:: console

    ATKconfig edit [-h] <key> <value>

Calls :func:`editConfig() <AstroToolkit.Config.editConfig>`, setting the <value> of a config <key>. See :ref:`Config Keys` for a list of available keys.

**Positional Arguments:**
    
    .. code-block:: console

        key         Name of config key to edit
        value       New value of config key

**Optional Arguments:**
    
    .. code-block:: console

        -h, --help  Show help message and exit

|

.. rubric:: Usage 4
    :heading-level: 2

.. code-block:: console
        
    ATKconfig reset [-h]

Calls :func:`resetConfig() <AstroToolkit.Config.resetConfig>`, resetting the config to its default state. 

**Optional Arguments:**
    
    .. code-block:: console

        -h, --help  Show help message and exit

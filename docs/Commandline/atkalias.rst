ATKalias
========

This tool allows for the various :ref:`Aliases <Aliases>` tools to be utilised from the command line. The first argument sets the job to perform on the alias list. See :ref:`tutorials <Adding Catalogue Aliases>` for examples. 

.. rubric:: Usage 1 
    :heading-level: 2

.. code-block:: console
        
    ATKalias show [-h]

Calls :func:`showAliases() <AstroToolkit.Aliases.showAliases>`, printing the current alias definitions to stdout.

**Optional Arguments:**

    .. code-block:: console

        -h, --help Show help message and exit

|

.. rubric:: Usage 2
    :heading-level: 2

.. code-block:: console

    ATKalias open [-h]

Calls :func:`openAliases() <AstroToolkit.Aliases.openAliases>`, opening the alias list in the default text editor.

**Optional Arguments:**

    .. code-block:: console

        -h, --help Show help message and exit

|

.. rubric:: Usage 3
    :heading-level: 2

.. code-block:: console

    ATKalias add [-h] <alias> <id>

Calls :func:`addAlias() <AstroToolkit.Aliases.addAlias>`, adding an alias definition with name <alias> that refers to a Vizier catalogue <id>.

**Positional Arguments:**
    
    .. code-block:: console

        alias       Alias name
        id          Vizier Catalogue ID

**Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit

|

.. rubric:: Usage 4 
    :heading-level: 2

.. code-block:: console
        
    ATKalias del [-h] <alias> 

Calls :func:`delAlias() <AstroToolkit.Aliases.delAlias>`, deleting an alias definition with name <alias>. 

**Positional Arguments:**

    .. code-block:: console

        alias       Alias name

**Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit

|

.. rubric:: Usage 5 
    :heading-level: 2

.. code-block:: console

    ATKalias reset [-h]

Calls :func:`resetAliases() <AstroToolkit.Aliases.resetAliases>`, resetting the alias list to its default state.

**Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit

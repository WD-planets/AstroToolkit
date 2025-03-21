ATKepoch
========

This tool allows for the various :ref:`Epochs <Epochs>` tools to be utilised from the command line. The first argument sets the job to perform on the epoch list. See :ref:`tutorials <Setting Survey and Alias Epochs>` for examples.

.. rubric:: Usage 1 
    :heading-level: 2

.. code-block:: console
        
    ATKepoch show [-h]

Calls :func:`showEpochs() <AstroToolkit.Epochs.showEpochs>`, printing the current epoch definitions to stdout.

**Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit

|

.. rubric:: Usage 2
    :heading-level: 2

.. code-block:: console

    ATKepoch open [-h]

Calls :func:`openEpochs() <AstroToolkit.Epochs.openEpochs>`, opening the epoch list in the default text editor.

**Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit

|

.. rubric:: Usage 3 
    :heading-level: 2

.. code-block:: console

    ATKepoch set [-h] <section> <survey> <epoch>

Calls :func:`setEpoch() <AstroToolkit.Epochs.setEpoch>`, setting the epoch of a survey (or catalogue alias) <survey> in a given <section> to <epoch>. 


**Positional Arguments:**

    .. code-block:: console

        section     Section of the epochs file in which survey or Vizier
                    catalogue alias is found
        survey      Survey for which epoch should be set
        epoch       Epoch to set (e.g. 2016,0] for Jan 2016)

*Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit

|

.. rubric:: Usage 4
    :heading-level: 2

.. code-block:: console
        
    ATKepoch del [-h] <alias> 

Calls :func:`delEpoch() <AstroToolkit.Epochs.delEpoch>`, deleting an epoch definition for a given <alias>. 

**Positional Arguments:**

    .. code-block:: console

        alias       Name of alias for which an epoch definition should be
                    deleted

**Optional Arguments:**

    .. code-block:: console
    
        -h, --help  Show help message and exit

|

.. rubric:: Usage 5
    :heading-level: 2

.. code-block:: console

    ATKepoch reset [-h]

Calls :func:`resetEpochs() <AstroToolkit.Epochs.resetEpochs>`, resetting the epoch list to its default state.

**Optional Arguments:**

    .. code-block:: console

        -h, --help  Show help message and exit

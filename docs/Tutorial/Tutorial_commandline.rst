The Command Line
================

Now that we have seen how ATK can be used in scripts, we can have a look at its command-line integration. A good first example is the command-line version of query: :ref:`ATKquery`. The simple example of a Gaia data query can be performed via the command line with:

.. code-block:: console

    ATKquery data gaia 587316166180416640

.. code-block:: console

    Running gaia data query
    source = 587316166180416640
    pos = None
    radius = 3.0

    Available Jobs: showdata, savedata <filename: str, optional>, exit

    Job?

Upon the successful retrieval of data, a more limited selection of the relevant :ref:`data structure's <Data Structures>` methods are available as "jobs". In this case, we have access to the :func:`showdata() <AstroToolkit.Data.dataquery.DataStruct.showdata>` and :func:`savedata() <AstroToolkit.Data.dataquery.DataStruct.savedata>` methods.

.. code-block:: console

    Running gaia data query
    source = 587316166180416640
    pos = None
    radius = 3.0

    Available Jobs: showdata, savedata <filename: str, optional>, exit

    Job? showdata

    .kind:       data
    .subkind:    data
    .survey:     gaia
    .catalogue:  I/355/gaiadr3
    .source:     587316166180416640
    .pos:        [141.18533044458, 8.03083432206]
    .identifier: J092444.48+080151.00
    .dataname:   J092444.48+080151.00_587316166180416640_gaia_ATKdata.fits

    .data:
        designation:                     ['Gaia DR3 587316166180416640']
        ra:                              [141.18526028]
        dec:                             [8.0308964]
        solution_id:                     [1636148068921376768]
        source_id:                       [587316166180416640]
        ... 

    Available Jobs: showdata, savedata <filename: str, optional>, exit

    Job? 

Jobs can be continuously entered until "exit" is entered, at which point the process stops. We can also target a system by its position:

.. code-block:: console

    ATKquery lightcurve ztf 141.185 8.031

.. code-block:: console

    Running ztf lightcurve query
    source = None
    pos = [141.185, 8.031]
    radius = 3.0

    Available Jobs: showdata, savedata <filename: str, optional> showplot <filename: str, optional>, saveplot<filename: str, optional>, exit

    Job? showplot
    Plot Type? lightcurve
    Plotting lightcurve data...
    Saving plot to local storage: J092444.40+080150.88_ztf_ATKlightcurve.html

.. raw:: html
    
    <div align="center"><embed src="../_static/commandline_lightcurve.html" width=100% height=500></embed></div>

We can also read files from the command line using the :ref:`ATKread` command, for example if we want to read our data file from the :ref:`previous tutorial <Saving to Local Files>`:

.. code-block:: console

    ATKread test_data.fits

.. code-block:: console

    Recreating data from local storage: test_data.fits
    Available Jobs: showdata, savedata <filename: str, optional>, exit

    Job? showdata

    .kind:       data
    .subkind:    data
    .survey:     galex
    .catalogue:  II/335/galex_ais
    .source:     587316166180416640
    .pos:        [141.18533044458, 8.03083432206]
    .identifier: J092444.48+080151.00
    .dataname:   J092444.48+080151.00_587316166180416640_galex_ATKdata.fits

    .data:
        RAJ2000:  [141.185551]
        DEJ2000:  [8.031037]
        Name:     ['GALEX J092444.5+080151']
        objid:    [6377741628902215075]
        FUVmag:   [19.6878]
        e_FUVmag: [0.113]
        NUVmag:   [19.5523]
        e_NUVmag: [0.0704]
        ...


Many other command line tools are available, some of which will be used in future sections. See :ref:`Command Line` for a full list and description of the available commands.

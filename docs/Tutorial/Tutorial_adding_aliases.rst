Adding Catalogue Aliases
========================

Along with the config, ATK supports `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue aliases. These save you from having to remember the Vizier catalogue IDs in :ref:`data queries <data-query>`, :ref:`bulkdata queries <bulkdata-query>` and :func:`datatable <AstroToolkit.Datapages.datatable>` entries.

As with the config, the easiest way to manages aliases is via the commandline. Upon first using ATK, a default aliases file will be generated which contains only the default data surveys that are :ref:`supported <Supported Surveys>` by the package out of the box. We can have a look at these using the :ref:`ATKalias` command-line tool:

.. code-block:: console
    
    ATKaliases show

.. code-block:: console
    
    Current ATKAliases.ini values:

    [default_catalogues]
    gaia = I/355/gaiadr3
    panstarrs = II/349/ps1
    skymapper = II/379/smssdr4
    galex = II/335/galex_ais
    rosat = IX/11/rosatsrc
    sdss = V/154/sdss16
    wise = II/311/wise
    twomass = II/246/out
    erosita = J/A+A/682/A34/erass1-m
    gaia_lc = I/355/epphot

    [additional_catalogues]

Currently, the only aliases are those that are supported by default. We can add one using the ATKaddalias command. In this case, we will add an alias to the AllWISE catalogue:

.. code-block:: console

    ATKalias add allwise II/328/allwise

.. code-block:: console

    Added alias for II/328/allwise with label allwise.

If we now show our alias file again, we will see that our new alias has been added.

.. code-block:: console

    ATKalias add

.. code-block:: console

    Current ATKAliases.ini values:

    [default_catalogues]
    gaia = I/355/gaiadr3
    panstarrs = II/349/ps1
    skymapper = II/379/smssdr4
    galex = II/335/galex_ais
    rosat = IX/11/rosatsrc
    sdss = V/154/sdss16
    wise = II/311/wise
    twomass = II/246/out
    erosita = J/A+A/682/A34/erass1-m
    gaia_lc = I/355/epphot

    [additional_catalogues]
    allwise = II/328/allwise

And we can now use this alias in a :ref:`data query <data-query>`!

.. code-block:: python

    from AstroToolkit.Tools import query

    query(kind="data",survey="allwise",source=6050296829033196032).showdata()

.. code-block:: console

    Running allwise data query
    source = 6050296829033196032
    pos = None
    radius = 5.0


    .kind:       data
    .subkind:    data
    .survey:     allwise
    .catalogue:  II/328/allwise
    .source:     6050296829033196032
    .pos:        [245.44701332769, -22.88621824697]
    .identifier: J162147.28-225310.39
    .dataname:   J162147.28-225310.39_6050296829033196032_allwise_ATKdata.fits

    .data:
        _r:        [0.485]
        _x:        [-0.01]
        _y:        [0.485]
        _p:        [358.8]
        AllWISE:   ['J162147.29-225310.7']
        RAJ2000:   [245.4470654]
        DEJ2000:   [-22.8863085]
        eeMaj:     [0.0367]
        eeMin:     [0.0343]
        eePA:      [97.]
        Im:        ['Im']
        W1mag:     [10.87300014]
        e_W1mag:   [0.023]
        W2mag:     [10.06299973]
        e_W2mag:   [0.022]
        W3mag:     [7.48899984]
        e_W3mag:   [0.023]
        W4mag:     [5.11399984]
        e_W4mag:   [0.038]
        ...

        Available Methods: .savedata(), .showdata()

Aliases will also be used when searching through all catalogues via :ref:`bulkdata queries <bulkdata-query>`:

.. code-block:: python
    
    query(kind="bulkdata",source=6050296829033196032).showdata()
   
.. code-block:: console

    Running bulkdata query
    source = 6050296829033196032
    pos = None
    radius = 3.0

    .data:
        gaia:
            designation: ['Gaia DR3 6050296829033196032']
            ra:          [245.44706008]
            dec:         [-22.88644709]
            ...

        panstarrs:
            RAJ2000: [245.44705339]
            DEJ2000: [-22.88642441]
            objID:   [80532454470486736]
            ...

        skymapper:
            SMSS:   ['162147.30-225311.1']
            RAICRS: [245.447072]
            DEICRS: [-22.886438]
            ...

        wise:
            WISE:    ['J162147.30-225310.6']
            RAJ2000: [245.447087]
            DEJ2000: [-22.886283]
            ...


        twomass:
            RAJ2000: [245.447007]
            DEJ2000: [-22.886181]
            errMaj:  [0.06]
            ...

        erosita:    
            .recno:  [4107]
            IAUName: ['1eRASS J162147.2-225311']
            DetUId:  ['em01_245114_020_ML00008_002_c010']
            ...

        allwise:
            AllWISE: ['J162147.29-225310.7']
            RAJ2000: [245.4470654]
            DEJ2000: [-22.8863085]
            ...


Any number of aliases can be added, and any that you no longer need can be deleted using the del mode of the :ref:`ATKalias` command:

.. code-block:: console

    ATKalias del allwise

Alternatively, the alias file can be reset entirely:

.. code-block:: console

    ATKalias reset

Just like the config, aliases can be viewed and edited at the same time by opening the alias file in the default text editor:

.. code-block:: console

   ATKalias open


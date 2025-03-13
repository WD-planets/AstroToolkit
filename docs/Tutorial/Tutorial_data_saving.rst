Saving to Local Files
=====================

Now that we have seen some of ATK's data structures, we can look at the package's local storage capabilities.

All ATK :ref:`data structures <Data Structures>` support local file saving. Returning to our :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>` example from :ref:`tutorial two <Fetching Data>`, we can save this to a local file using the :func:`savedata() <AstroToolkit.Data.dataquery.DataStruct.savedata>` method:

.. code-block:: python

    from AstroToolkit.Tools import query

    gaia_data = query(kind="data",source=587316166180416640,survey="galex")
    gaia_data.showdata()
    gaia_data.savedata("test_data.fits")

.. code-block:: console

    .kind:       data
    .subkind:    data
    .survey:     galex
    .catalogue:  II/335/galex_ais
    .source:     587316166180416640
    .pos:        [141.18526027626, 8.03089639753]
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

    Available Methods: .savedata(), .showdata()

This generates a .fits file which contains the data structure. We have overriden the file name here, otherwise it would default to the name stored in the **dataname** attribute. A key feature of file saving in ATK is that it is entirely *lossless* - i.e. the original data structure can be recreated exactly. This is done by reading the file with the :func:`readdata() <AstroToolkit.Tools.readdata>` tool:

.. code-block:: python

    from AstroToolkit.Tools import readdata

    data=readdata("test_data.fits")
    data.showdata()

.. code-block:: console

    .kind:       data
    .subkind:    data
    .survey:     galex
    .catalogue:  II/335/galex_ais
    .source:     587316166180416640
    .pos:        [141.18526027626, 8.03089639753]
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

    Available Methods: .savedata(), .showdata()

The above is equally applicable to any ATK :ref:`data structure <Data Structures>` returned by :func:`query() <AstroToolkit.Tools.query>`, although figures must be saved separately if you wish to avoid having to re-plot the data structure.

The :func:`query() <AstroToolkit.Tools.query>` tool also provides a **check_exists** flag, which searches for an existing file and recreates the data from this if possible. This can save a significant amount of time when running scripts multiple times. If the file is not found the query will run as normal, and the returned data will be saved to the provided file name for next time.

.. code-block:: python

    lightcurve_data = query(kind="lightcurve",source=6050296829033196032,survey="ztf",check_exists="test_lightcurve.fits")

The first time this is run, a query will be performed since the file does not currently exist. The resulting light curve will be saved to local files under the filename "test_lightcurve.fits". Any future executions of this script will then recreate the data structure from local files instead of running a new query.

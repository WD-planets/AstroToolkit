Fetching Data
=============

The cornerstone of ATK is its data structures, which standardise the storage and manipulation of data for use throughout the package. 

These data structures can be accessed via the :func:`query() <AstroToolkit.Tools.query>` tool, which is used to gather any data from :ref:`supported surveys <Supported Surveys>`. To target a system, we need either its J2000 coordinates (called 'pos' throughout ATK) or its Gaia DR3 source ID ('source'). As noted in the :ref:`previous section <Basic Info>`, targeting a system with a Gaia source ID is better where possible as this enables proper motion correction to be used throughout the package. This is the same for most tools in ATK, and is denoted by **pos/source** in a tool's input arguments - meaning that one of these is required (see the :ref:`documentation <Modules>`).

The :func:`query() <AstroToolkit.Tools.query>` tool takes different arguments depending on the kind of query being performed (as described in its documentation), but a data query to Gaia is perhaps the most simple:

.. code-block:: python
    
    from AstroToolkit.Tools import query
    
    gaia_data = query(kind="data",source=587316166180416640,survey="gaia")

where the cataclysmic variable Hu Leo has been targeted using its Gaia DR3 source ID. This returns any Gaia DR3 catalogue data that is found for that source. We can also query the same system in another of the supported surveys, e.g. GALEX:

.. code-block:: python

    galex_data = query(kind="data",source=587316166180416640,survey="galex")

Here, the system's Gaia DR3 coordinates are corrected for proper motion back to GALEX's epoch, and any data is returned. We have not provided a radius, and so one has been taken from the :ref:`config <Config Keys>`.

We have now fetched some data, but what can we do with it? Data, bulkdata, and reddening queries all return a specific ATK data structure: a :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>`. Since these forms of data aren't going to be plotted, we have two methods available: :func:`showdata() <AstroToolkit.Data.dataquery.DataStruct.showdata>` and :func:`savedata() <AstroToolkit.Data.dataquery.DataStruct.savedata>`. The latter will be covered later, so let's focus on the former.

The :func:`showdata() <AstroToolkit.Data.dataquery.DataStruct.showdata>` method is available on all ATK data structures, and prints the data structure to stdout in a readable format. Continuing from the above:

.. code-block:: python

    galex_data.showdata()

.. code-block:: console

    Running galex data query
    source = 587316166180416640
    pos = None
    radius = 3.0


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

**Note:** For readability in this tutorial, the majority of the returned GALEX columns have been omitted.

The first section in the above output notifies us that the query is running, and the rest is the result of :func:`showdata() <AstroToolkit.Data.dataquery.DataStruct.showdata>`. We can now look at the data structure's attributes.

- The **kind** attribute simply describes the type of data being stored (in this case, "data" means catalogue data). **Subkind** is an attribute only found in :class:`DataStructs <AstroToolkit.Data.data.DataStruct>` (as these are shared between data, bulkdata and reddening queries) and denotes which of these the structure is storing. 

- The **survey** attribute describes which survey the data originates from, and **catalogue** holds the `Vizier <https://vizier.cds.unistra.fr/>`_ ID of that survey (for data queries such as those performed above, any Vizier catalogue can be queried).

- The **source** attribute holds the Gaia source to which the data pertains, **pos** holds its J2000 coordinates [right ascension, declination] in degrees, and identifier gives these as a string in HHMMSS.SS±DDMMSS.SS format.

- Finally, the **dataname** attribute gives the default file name to which data will be saved locally using the :func:`savedata() <AstroToolkit.Data.dataquery.DataStruct.savedata>` method (again, this will be covered later). 

In :class:`DataStructs <AstroToolkit.Data.data.DataStruct>`, the resulting data is stored as a dictionary, and hence we can access the value of a certain parameter using its column heading:

.. code-block:: python

    print(galex_data.data["FUVmag"][0],galex_data.data["NUVmag"][0])

.. code-block:: console

    19.6878 19.5523 

|

We have now seen an example of a data query, but this is only one of many kinds of data that can be fetched through ATK. Below is a full list of the various kinds of query and the data structures that they return:

- data query, returns a :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>`
- bulk data ('bulkdata') query, returns a :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>`
- reddening query, returns a :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>`
- image query, returns an :class:`ImageStruct <AstroToolkit.Data.imagequery.ImageStruct>`
- lightcurve query, returns a :class:`LightcurveStruct <AstroToolkit.Data.lightcurvequery.LightcurveStruct>`
- HRD query, returns a :class:`HrdStruct <AstroToolkit.Data.hrdQuery.HrdStruct>`
- SED query, returns an :class:`SedStruct <AstroToolkit.Data.sedquery.SedStruct>`
- spectrum query, returns a :class:`SpectrumStruct <AstroToolkit.Data.spectrumquery.SpectrumStruct>`

While there are some differences between the various :ref:`data structures <Data Structures>` in ATK (more specifically, the format of their .data attribute will of course differ significantly), all share a similar form to the :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>` explored above.

|

Note: The code used in this tutorial can be executed using:

.. code-block:: python

    from AstroToolkit.Examples import data_fetching

Setting Survey and Alias Epochs
===============================

In the :ref:`previous example <Adding Catalogue Aliases>`, you may have noticed some warnings stating that ATK failed to account for proper motion in our queries to AllWISE. This is true, but the reason is simple: we have not yet set an epoch for our new catalogue alias, and so ATK doesn't know the magnitude of any necessary corrections. While this wasn't a problem for the source used in our example due to its small proper motion, this issue is nonetheless worth solving since this will not always be the case.

To highlight this, we return to our high proper motion source from the :ref:`introduction <Basic Info>`, Wolf 28 (Van Maanen 2). With no epoch definition, a query to AllWISE (in this case via the command line) returns no data:

.. code-block:: console

    ATKquery data allwise --source 2552928187080872832

.. code-block:: console

    Running allwise data query
    source = 2552928187080872832
    pos = None
    radius = 3.0

    Note: allwise has no epoch definition. Proper motion has therefore not been corrected.
    Note: allwise data query unsuccessful or returned no data.

However, we can remedy this. The :ref:`Epochs` module allows us to access and edit all epochs that are used by ATK. Much like the :ref:`Alias <Aliases>` and :ref:`Config <Config>` modules, this module is most easily utilised via the command line. In this case, we can use the :ref:`ATKepoch` tool. We can start by showing the current epoch definitions:

.. code-block:: console

    ATKepoch show

.. code-block:: console

    Current ATKEpochs.ini values:

    [default_data_surveys]
    gaia = 2016,0
    panstarrs = 2012,0
    skymapper = 2016,0
    galex = 2006,8
    rosat = 1991,0
    sdss = 2017,0
    wise = 2010,5
    twomass = 1999,0
    erosita = 2022,0

    [additional_data_surveys]

    [lightcurve_surveys]
    atlas = 2021,0
    ztf = 2019,0
    crts = 2000,0
    asassn = 2015,0
    gaia = 2016,0
    tess = 2020,0

    [spectrum_surveys]
    sdss = 2017,0

The above shows the current epoch definitions in our instance of ATK. The **Default_data_surveys** section lists those of the :ref:`supported surveys <Supported Surveys>` in :ref:`data queries <data-query>`, while **lightcurve_surveys** and **spectrum_surveys** list those that are available in :ref:`lightcurve queries <lightcurve-query>` and :ref:`spectrum queries <spectrum-query>`. As expected, the **additional_data_surveys** section is empty. To remedy this, we use the *set* mode of the :ref:`ATKepoch` command:

.. code-block:: console

    ATKepoch set additional_data_surveys allwise 2010,5 

Returning to our previous query, we now see that data is successfully retrieved!

.. code-block:: console

    ATKquery data allwise --source 2552928187080872832

.. code-block:: console

    .kind:       data
    .subkind:    data
    .survey:     allwise
    .catalogue:  ii/328/allwise
    .source:     2552928187080872832
    .pos:        [12.293161429221492, 5.384403492110555]
    .identifier: J004910.36+052303.85
    .dataname:   J004910.36+052303.85_2552928187080872832_allwise_ATKdata.fits
    .trace:      start -> extracted pos from source query, assumed [2016, 0] -> allwise: [2010, 5] -> allwise query performed -> [2000,0] -> end

    .data:
        AllWISE:   ['J004910.78+052249.8']
        RAJ2000:   [12.2949231]
        DEJ2000:   [5.3805217]
        W1mag:     [11.40999985]
        e_W1mag:   [0.023]
        W2mag:     [11.35700035]
        e_W2mag:   [0.02]
        W3mag:     [11.38000011]
        e_W3mag:   [0.197]
        W4mag:     [8.88000011]
        e_W4mag:   [nan]
        ...

For reference: finding this detection via the object's name in Vizier would require a search radius of over 30 arcsec, resulting in the retrieval of three possible detections. Furthermore, this detection is actually listed as the furthest one from the source! To prove that this is in fact an accurately assigned detection, we can use an image with an allwise overlay:

[INSERT AN IMAGE ONCE I SORT OUT CUSTOM OVERLAYS]

Much like the :ref:`ATKalias` command, we can also delete epoch definitions for a specific alias:

.. code-block:: console

    ATKepoch del allwise

We can also reset the list of epoch definitions to its default state:

.. code-block:: console

    ATKepoch reset

and can open it in the default text editor:

.. code-block:: console

    ATKepoch open

**Note:** Epochs do not need to be set for image surveys when used in :ref:`image queries <image-query>`, as the proper motion correction in this case makes use of the specific image epoch rather than a single definition. Similarly, epoch definitions are not required for :ref:`SED <sed-query>` and :ref:`bulkdata <bulkdata-query>` queries, as these simply use the definitions of the relevant surveys above.

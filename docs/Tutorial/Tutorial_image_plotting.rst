Plotting an Image
=================

Besides the :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>` that we have just encountered, all ATK :ref:`data structures <Data Structures>` support plotting functionality. A good first example would be to recreate the image seen in :ref:`tutorial one <Basic Info>`. To fetch the data we need, a similar query can be performed to the one we used in the :ref:`previous tutorial <Fetching Data>`:

.. code-block:: python
    
    from AstroToolkit.Tools import query

    image_data = query(
        kind="image",
        source=2552928187080872832,
        survey="panstarrs",
        overlays=["gaia", "galex"],
        size=120,
    )

where we have used some additional parameters that are specific to :ref:`image queries <image-query>`: **overlays** tells the tool which surveys to fetch overlay data for, and **size** gives the required size of the image in arcseconds. To see the resulting data structure, we can once again make use of the :func:`showdata <AstroToolkit.Data.imagequery.ImageStruct.showdata>` method:

.. code-block:: python

   image_data.showdata()

.. code-block:: console

    .kind:       image
    .survey:     panstarrs
    .source:     2552928187080872832
    .pos:        [12.29124307015, 5.38860939472]
    .identifier: J004909.90+052318.99
    .figure:     None
    .dataname:   J004909.90+052318.99_2552928187080872832_panstarrs_ATKimage.fits
    .plotname:   J004909.90+052318.99_2552928187080872832_panstarrs_ATKimage.html

    .data:
        image_data:   <Image Data>
        image_header: <Image Header>
        size:         120
        image_time:   [2012    1]
        wcs:          <WCS Object>
        image_focus:  [12.29539461  5.37950704]
        overlay:      <Overlay Data>

    Available Methods: .exportplot(), .plot(), .savedata(), .saveplot(), .showdata(), .showplot()

Most of the above should be familiar from the :ref:`previous tutorial <Fetching Data>`, but there are now some additional attributes and methods to have a look at.

- The first additional attribute is **figure**, which stores any plot which is derived from the data structure's data. Since we haven't made one yet, this is currently empty.

- The only other new attribute is **plotname**, which stores the default file name to which this figure will be saved locally using the :func:`showplot() <AstroToolkit.Data.imagequery.ImageStruct.showplot>`, :func:`saveplot() <AstroToolkit.Data.imagequery.ImageStruct.saveplot>` and :func:`exportplot() <AstroToolkit.Data.imagequery.ImageStruct.exportplot>` methods.

|

We can now use one of the new methods to plot our data! The :func:`plot() <AstroToolkit.Data.imagequery.ImageStruct.plot>` method plots the data structure's data and assigns the resulting Bokeh figure to the **figure** attribute:

.. code-block:: python
    
    image_data.plot()

Using :func:`showdata() <AstroToolkit.Data.imagequery.ImageStruct.showdata>` again reveals that the **figure** attribute of our data structure has now been filled:

   .. code-block:: console

    Plotting image data...

    .kind:       image
    .survey:     panstarrs
    .source:     2552928187080872832
    .pos:        [12.29124307015, 5.38860939472]
    .identifier: J004909.90+052318.99
    .figure:     figure(id='p1001', ...)
    .dataname:   J004909.90+052318.99_2552928187080872832_panstarrs_ATKimage.fits
    .plotname:   J004909.90+052318.99_2552928187080872832_panstarrs_ATKimage.html

    .data:
        image_data:   <Image Data>
        image_header: <Image Header>
        size:         120
        image_time:   [2012    1]
        wcs:          <WCS Object>
        image_focus:  [12.29539461  5.37950704]
        overlay:      <Overlay Data>

    Available Methods: .exportplot(), .plot(), .savedata(), .saveplot(), .showdata(), .showplot()

We can open the resulting figure in the web browser using another of the new methods: :func:`showplot() <AstroToolkit.Data.imagequery.ImageStruct.showplot>`.

.. code-block:: python

   image_data.showplot()

.. _example-image:

.. raw:: html
    
    <div align="center"><embed src="../_static/image.html" width=100% height=450></embed></div>

And with that, we have recreated the image from the :ref:`first tutorial <Basic Info>`! Here, we first generated the plot using :func:`plot() <AstroToolkit.Data.imagequery.ImageStruct.plot>` and then opened it with :func:`showplot() <AstroToolkit.Data.imagequery.ImageStruct.showplot>`, but we since we did not customise the plot in any way we also could have just used :func:`showplot() <AstroToolkit.Data.imagequery.ImageStruct.showplot>` and one would have automatically been generated. The use case for :func:`plot() <AstroToolkit.Data.imagequery.ImageStruct.plot>` will be shown in the :ref:`next section <Plotting a Light Curve>`.

|

Note: The code used in this tutorial can be executed using:

.. code-block:: python

    from AstroToolkit.Examples import image_plotting 

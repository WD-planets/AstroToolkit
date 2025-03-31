Defining New Image Overlays
===========================
In the :ref:`previous tutorial <Setting Survey and Alias Epochs>`, we managed to retreive some additional AllWISE data for Wolf 28 using proper motion correction. But how can we be sure that the returned data actually belongs to our target system? The best way to check this is to make use of ATK's image overlays. Since we are using a user-defined alias to AllWISE (see :ref:`previous tutorial <Adding Catalogue Aliases>`), we first need to tell ATK how to generate an AllWISE overlay. We do this using the :ref:`Overlays <Overlays>` module, which we can once again access from the command line via the :ref:`ATKoverlay` tool. 

To generate an overlay definition for a user-defined catalogue alias, we need to provide ATK with the name of the survey's right ascension, declination, ID, and magnitude columns (if available) in Vizier:

.. code-block:: console 
    
    ATKoverlay add allwise --ra RAJ2000 --dec DEJ2000 --id AllWISE --mags W1mag W2mag W3mag W4mag

.. code-block:: console

    Added scaled_detection overlay defintion for alias 'allwise'. 

Since we provided magnitude columns, the overlay will be scaled with brightness (see :ref:`here <Plotting an Image>`) - hence a 'scaled_detection' overlay definition has been generated. We can now generate images with AllWISE overlays just like we would with any survey that ATK supports out of the box!

.. code-block:: python

    from AstroToolkit.Tools import query

    image = query(kind="image", survey="panstarrs", source=2552928187080872832, overlays=["allwise"], size=120).showplot()

.. code-block:: console

    Running panstarrs image query
    source = 2552928187080872832
    pos = None
    size = 120

    Plotting image data...

.. _user-defined_overlay:

.. raw:: html
    
    <div align="center"><embed src="../_static/user_defined_overlay.html" width=100% height=450></embed></div>

While we can already see that the detection lines up almost perfectly with our target source, we can make sure by hovering over the central detection (i.e. the one that pertains to Wolf 28) and comparing the object ID (J004910.78+052249.8) with that of the data that we retrieved in the :ref:`previous tutorial <Setting Survey and Alias Epochs>`. As expected, they match!

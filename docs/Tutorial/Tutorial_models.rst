Using External Data
===================

While support for further surveys is being gradually added, it is also possible to use data from external sources in ATK through the :ref:`Models` module. This module provides empty :ref:`data structures <Data Structures>` for use throughout the package, which can be filled with any relevant data. Some attributes can be filled automatically if the model is provided with a **pos** or **source**.

As an example, we will use some photometry from the Thai National Telescope (TNT) for a system that we have :ref:`already seen <Plotting a Light Curve>`: the white dwarf pulsar, AR Sco. For reference, the data takes the form:

.. code-block:: console
    
    mjd            flux     error
    2457517.209956 0.454216 0.0185251
    2457517.210012 0.123707 0.0158678
    2457517.210074 -0.0441742 0.015052
    2457517.210136 0.160731 0.0161515
    2457517.210198 0.398715 0.0172122
    ...            ...      ...

and is saved as AR_Sco_TNT.txt in this example. We must first read the data into Python, in this case we will use Pandas:

.. code-block:: python

    import pandas as pd

    data = pd.read_csv("AR_Sco_TNT.txt", delimiter="\s+")

We then import and create a custom :class:`ATK Lightcurve Structure <AstroToolkit.Models.CustomLightcurveStruct>`:

.. code-block:: python

    from AstroToolkit.Models import CustomLightcurveStruct
    
    lightcurve = CustomLightcurveStruct().showdata()

.. code-block:: console
    
    .kind:       lightcurve
    .survey:     None
    .source:     None
    .pos:        None
    .identifier: None
    .figure:     None
    .dataname:   None
    .plotname:   None
    .data:       None
    Available Methods: .bin(), .crop(), .exportplot(), .plot(), .savedata(), .saveplot(), .showdata(), .showplot(), .sigmaclip()

As we can see, we now have an empty :class:`LightcurveStruct <AstroToolkit.Data.lightcurvequery.LightcurveStruct>` which we can fill with data. Since we know the source (or position) of our data, we can get a head start by providing this to :class:`CustomLightcurveStruct <AstroToolkit.Models.CustomLightcurveStruct>`:

.. code-block:: python

    lightcurve = CustomLightcurveStruct(source=6050296829033196032).showdata()

.. code-block:: console 

    Running gaia data query
    source = 6050296829033196032
    pos = None
    radius = 3.0


    .kind:       lightcurve
    .survey:     None
    .source:     6050296829033196032
    .pos:        [245.44701332769, -22.88621824697]
    .identifier: J162147.28-225310.39
    .figure:     None
    .dataname:   J162147.28-225310.39_6050296829033196032_ATKlightcurve.fits
    .plotname:   J162147.28-225310.39_6050296829033196032_ATKlightcurve.html
    .data:       None
    Available Methods: .bin(), .crop(), .exportplot(), .plot(), .savedata(), .saveplot(), .showdata(), .showplot(), .sigmaclip()

With a number of the structure's fields being filled out for us, we simply fill those that remain with the name of the survey and our custom data (in the same format as a normal :class:`LightcurveStruct <AstroToolkit.Data.lightcurvequery.LightcurveStruct>`):

.. code-block:: python 
    
    lightcurve.survey = "TNT" 
    lightcurve.data = [
        {
            "band": "g",
            "hjd": data["mjd"].tolist(),
            "mag": data["flux"].tolist(),
            "mag_err": data["error"].tolist(),
        }
    ]

And that is it! We can now use our external data throughout ATK as if it were officially supported:

.. code-block:: python

    lightcurve.plot(colours=["green"]).showplot()

.. raw:: html
    
    <div align="center"><embed src="../_static/AR_Sco_TNT_Lightcurve.html" width=100% height=500></embed></div>

.. code-block:: python

    lightcurve.plot(kind="powspec", start_freq=650, stop_freq=800).showplot()

.. raw:: html
    
    <div align="center"><embed src="../_static/AR_Sco_TNT_Powspec.html" width=100% height=425></embed></div>

.. code-block:: python

    lightcurve.plot(kind="phasefold", bins=300, foverlay=False, freq=6.74157303371).showplot()

.. raw:: html
    
    <div align="center"><embed src="../_static/AR_Sco_TNT_Phasefold.html" width=100% height=500></embed></div>

|

Note: The code used in this tutorial can be executed using:

.. code-block:: python

    from AstroToolkit.Examples import external_data

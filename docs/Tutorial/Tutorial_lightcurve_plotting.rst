Plotting a Light Curve
======================

We can now look at plotting another kind of data: light curves. Through this, we will see some of the customisation options that are available through :func:`plot() <AstroToolkit.Data.lightcurvequery.LightcurveStruct.plot>`.

The first step is to fetch some light curve data. For this, we introduce another system: the white dwarf pulsar, AR Sco. 

.. code-block:: python

    from AstroToolkit.Tools import query

    lightcurve_data = query(kind="lightcurve",source=6050296829033196032,survey="ztf")
    lightcurve_data.showdata()

.. code-block:: console

    .kind:       lightcurve
    .survey:     ztf
    .source:     6050296829033196032
    .pos:        [245.44701332846378, -22.886218247040002]
    .identifier: J162147.28-225310.39
    .figure:     None
    .dataname:   J162147.28-225310.39_6050296829033196032_ztf_ATKlightcurve.fits
    .plotname:   J162147.28-225310.39_6050296829033196032_ztf_ATKlightcurve.html
    .trace:      start -> extracted pos from source query, assumed [2016, 0] -> ztf: [2019, 0] -> ztf query performed -> [2000,0] -> end

    .data:
        band: g
            ra:      [245.447071  245.4470765 245.447069  ... 245.4471056 245.4470825 245.4470859]
            dec:     [-22.8864202 -22.8864231 -22.8864182 ... -22.8865103 -22.8865078 -22.8865203]
            hjd:     [2458235.88088966 2458235.94319608 2458246.86172527 ... 2460524.68224955 2460527.69099994 2460529.69086604]
            mag:     [16.7012329 14.7080278 16.4311085 ... 15.7799644 14.800931  16.0791607]
            mag_err: [0.02124788 0.019023   0.02014201 ... 0.01882575 0.01893662 0.01923533]

        band: r
            ra:      [245.4470594 245.4470709 245.4470644 ... 245.4471027 245.4470833 245.4470933]
            dec:     [-22.8864405 -22.8864381 -22.8864126 ... -22.8865397 -22.8865631 -22.8865343]
            hjd:     [2458218.90986354 2458218.97188203 2458246.9025947  ... 2460504.77009979 2460510.75659636 2460527.65778532]
            mag:     [15.6002302 15.0777416 14.619998  ... 15.7173185 14.579071  15.1474562]
            mag_err: [0.01373806 0.01352318 0.01368265 ... 0.01386647 0.01370837 0.01352345]

        band: i
            ra:      None
            dec:     None
            hjd:     None
            mag:     None
            mag_err: None

    Available Methods: .bin(), .crop(), .exportplot(), .plot(), .savedata(), .saveplot(), .showdata(), .showplot(), .sigmaclip()

The attributes of the resulting :class:`LightcurveStruct <AstroToolkit.Data.lightcurvequery.LightcurveStruct>` are all familiar. We do have some new methods (:func:`bin() <AstroToolkit.Data.lightcurvequery.LightcurveStruct.bin>`, :func:`crop() <AstroToolkit.Data.lightcurvequery.LightcurveStruct.crop>` and :func:`sigmaclip() <AstroToolkit.Data.lightcurvequery.LightcurveStruct.sigmaclip>`), but we won't look at those here.

Let's now plot our light curve, skipping :func:`plot() <AstroToolkit.Data.lightcurvequery.LightcurveStruct.plot>`:

.. code-block:: python

    lightcurve_data.showplot()

.. code-block:: console

    Note: No plot was found. One will be generated with a default configuration.
    Plotting lightcurve data...
    Saving plot to local storage: J092444.48+080151.00_587316166180416640_ztf_ATKlightcurve.html

.. raw:: html
    
    <div align="center"><embed src="../_static/lightcurve1.html" width=100% height=500></embed></div>

|

This looks a little bland, and the bands all merge together. We can pass some additional arguments to :func:`plot() <AstroToolkit.Data.lightcurvequery.LightcurveStruct.plot>` to add some colour:

.. code-block:: python
    
    lightcurve_data.plot(bands=["g","r","i"],colours=["green","red","blue"])
    lightcurve_data.showplot()

.. raw:: html
    
    <div align="center"><embed src="../_static/lightcurve2.html" width=100% height=500></embed></div>

|

**Note:** i-band data is missing since AR Sco has no i-band photometry in ZTF.

That looks better. However, the options that we can pass to :func:`plot() <AstroToolkit.Data.lightcurvequery.LightcurveStruct.plot>` are not purely stylistic. :class:`LightcurveStructs <AstroToolkit.Data.lightcurvequery.LightcurveStruct>` support three kinds of plotting: "lightcurve" (as seen above), "powspec", and "phasefold". The latter two utilise the Lomb-Scargle periodogram to perform time series analysis.

Let's try generating a power spectrum:

.. _power-spectrum:

.. code-block:: python

    lightcurve_data.plot(kind="powspec").showplot()

.. raw:: html
    
    <div align="center";><embed src="../_static/powspec.html" width=100% height=425></embed></div>

and a phase folded light curve:

.. code-block:: python

    lightcurve_data.plot(kind="phasefold").showplot()

.. raw:: html
    
    <div align="center";><embed src="../_static/phasefold1.html" width=100% height=425></embed></div>

Additional plotting options are available depending on the kind of data being plotted, and the way in which it is being plotted (e.g. the lightcurve/powspec/phasefold plots above). As an example, we could bin the phase-folded light curve and shift it slightly to the right to better align it with our phase axis:

.. code-block:: python

    lightcurve_data.plot(kind="phasefold",bins=100,shift=0.115)

.. raw:: html
    
    <div align="center";><embed src="../_static/phasefold2.html" width=100% height=425></embed></div>

For a full description of available parameters, see each :ref:`data structure's <Data Structures>` documentation.

|

Note: The code used in this tutorial can be executed using:

.. code-block:: python

    from AstroToolkit.Examples import lightcurve_plotting 

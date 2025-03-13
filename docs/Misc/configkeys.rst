.. _config-keys:

Config Keys
===========

Below is a list and description of all available config keys along with their default values. See the :ref:`Config` module for information on how to customise them.

|

Global Settings
---------------
    .. _cfg_enable_notifications:

    enable_notifications: *bool*, defaults to True
        toggles the output of notifications to stdout when a tool is running 

    .. _cfg_unit_size:

    unit_size: *int*, defaults to 500
        sets the scaling size of all ATK figures. 

    .. _cfg_output_backend:

    output_backend: *str*, defaults to canvas
        sets the output backend used to generate bokeh figures, from: svg, canvas, webgl. The former is the slowest, and should only be used when printing to PDF. The latter is the fastest, but support is limited (see the `Bokeh documentation <https://docs.bokeh.org/en/latest/docs/user_guide/output/webgl.html>`_ for details)

    .. _cfg_show_toolbars:    

    show_toolbars: *str*, defaults to True
        enables/disables the toolbar in all plots

    .. _cfg_show_grids:

    show_grids: *str*, defaults to True
      enables/disables the background grid in all ATK figures 

    .. _cfg_show_titles:

    show_titles: *str*, defaults to True
        enables/disables titles in all ATK figures

    .. _cfg_font_size:

    font_size: int, defaults to 14
        sets the font size of all text in all ATK figures

    .. _cfg_font:

    font: *str*, defaults to Helvetica
        sets the font used in any text in ATK figures 
    
    |

Query Settings
--------------
    .. _cfg_query_data_radius:

    query_data_radius: *float*, defaults to 3
        sets the default search radius used in :ref:`data queries <data-query>` in arcseconds

    .. _cfg_query_phot_radius:

    query_phot_radius: *float*, defaults to 3
        sets the default search radius for :ref:`photometry queries <phot-query>` in arcseconds

    .. _cfg_query_bulkphot_radius:

    query_bulkphot_radius: *float*, defaults to 3
        sets the default search radius for :ref:`bulkphot queries <bulkphot-query>` in arcseconds

    .. _cfg_query_lightcurve_radius:

    query_lightcurve_radius: *float*, defaults to 3
        sets the default search radius for :ref:`lightcurve queries <lightcurve-query>` in arcseconds

    .. _cfg_query_spectrum_radius:

    query_spectrum_radius: *float*, defaults to 3
        sets the default search radius for :ref:`spectrum queries <spectrum-query>` in arcseconds

    .. _cfg_query_sed_radius:

    query_sed_radius: *float*, defaults to 3
        sets the default search radius for :ref:`sed queries <sed-query>` in arcseconds

    .. _cfg_query_reddening_radius:

    query_reddening_radius: *float*, defaults to 3
        sets the default search radius for :ref:`reddening queries <reddening-query>` in arcseconds

    .. _cfg_query_image_size:

    query_image_size: *int*, defaults to 30
        sets the default image size to be used in :ref:`image queries <image-query>` in arcseconds

    .. _cfg_query_image_overlays:

    query_image_overlays: *str*, defaults to gaia
        sets the default detection overlays to use in :ref:`image queries <image-query>` 

    .. _cfg_query_image_band:

    query_image_band: *str*, defaults to g
        sets the default band to use in :ref:`image queries <image-query>`

    .. _cfg_query_lightcurve_atlas_username:

    query_lightcurve_atlas_username: *str*, defaults to None
        sets the default ATLAS username to use in ATLAS :ref:`lightcurve queries <lightcurve-query>`

    .. _cfg_query_lightcurve_atlas_password:

    query_lightcurve_atlas_password: *str*, defaults to None
        sets the default ATLAS password to use in ATLAS :ref:`lightcurve queries <lightcurve-query>`

    |

Image Overlay Settings
----------------------
    .. _cfg_gaia_overlay_mag:

    gaia_overlay_mag: *str*, defaults to phot_g_mean_mag
        sets the default band to use in Gaia detection overlays

    .. _cfg_galex_overlay_mag:

    galex_overlay_mag: *str*, defaults to NUVmag
        sets the default band to use in GALEX detection overlays

    .. _cfg_wise_overlay_mag:

    wise_overlay_mag: *str*, defaults to W1mag
        sets the default band to use in WISE detection overlays

    .. _cfg_sdss_overlay_mag:

    sdss_overlay_mag: *str*, defaults to gPmag
        sets thedefault band to use in SDSS detection overlays

    .. _cfg_2mass_overlay_mag:

    twomass_overlay_mag: *str*, defaults to Jmag
        sets the default band to use in 2MASS detection overlays

    .. _cfg_panstarrs_overlay_mag:

    panstarrs_overlay_mag: *str*, defaults to gmag
        sets the default band to use in Pan-STARRS detection overlays

    .. _cfg_skymapper_overlay_mag:

    skymapper_overlay_mag: *str*, defaults to gPSF
        sets the default band to use in SkyMapper detection overlays

    .. _cfg_overlay_piggyback_radius:

    overlay_piggyback_radius: *float*, defaults to 5
        sets the default 'piggyback' radius to use in arcseconds when proper motion correcting non-Gaia detections in imaging overlays. Larger radius will result in more detections being corrected, but a larger chance of this correction being erroneous. 

    |

Search Settings
---------------
    .. _cfg_search_radius:

    search_radius: *float*, defaults to 3
        sets the default search radius used in Vizier/SIMBAD searches using the :func:`search() <AstroToolkit.Tools.search>` tool in arcseconds

    |

Datapage Settings
-----------------
    .. _cfg_datapage_search_button_radius:

    datapage_search_button_radius: *float*, defaults to 3
        sets the default search radius used in Vizier/SIMBAD :func:`buttons <AstroToolkit.Datapages.buttons>` in arcseconds

    .. _cfg_datapage_datatable_radius:

    datapage_datatable_radius: *float*, defaults to 3
        sets the default search radius used for populating :func:`datatables <AstroToolkit.Datapages.datatable>` in arcseconds

    .. _cfg_datapage_font_size:

    datapage_font_size: *int*, defaults to 12
        sets the default font size when formatting plots with :func:`grid() <AstroToolkit.Datapages.grid>`

    .. _cfg_datapage_grid_size:

    datapage_grid_size: *int*, defaults to 250
        sets the default size of datapage :func:`grids <AstroToolkit.Datapages.grid>`


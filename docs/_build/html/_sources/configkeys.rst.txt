Config Keys
===========

Here, all config keys and there functions will be outlined. See the config module for information on how to customise them.

|

Global Settings
---------------
   **enable_notifications:** bool = True
      if True, outputs notifications of queries to stdout when executed

   **unit_size:** int = 400
      scales all ATK figures

   **output_backend:** str = canvas
      sets the output backend used to generate bokeh figures, from svg, canvas, webgl

   **font_size:** int = 18
      sets the font size of all text in ATK figures

   **font:** str = Times New Roman
      sets the font used in all text in ATK figures, from: any font supported by bokeh

|

Query Settings
--------------
   **query_data_radius:** int = 3
      sets the default search radius for data queries in arcseconds

   **query_phot_radius:** int = 3
      sets the default search radius for phot queries in arcseconds

   **query_bulkphot_radius:** int = 3
      sets the default search radius for bulkphot queries in arcseconds

   **query_reddening_radius:** int = 5
      sets the default search radius for reddening queries in arcseconds
   
   **query_lightcurve_radius:** int = 3
      sets the default search radius for lightcurve queries in arcseconds

   **query_spectrum_radius:** int = 3
      sets the default search radius for spectrum queries in arcseconds

   **query_sed_radius:** int = 3
      sets the default search radius for sed queries in arcseconds

   **query_image_size:** int = 30
      sets the default image size in arcseconds 

   **query_image_overlays:** str = gaia
      sets the default detection overlays to use in image queries 

   **query_image_band:** str = g
      sets the default band to use in image queries

   **query_lightcurve_atlas_username:** str = None
      sets the default ATLAS username to use in ATLAS lightcurve queries

   **query_lightcurve_atlas_password:** str = None
      sets the default ATLAS password to use in ATLAS lightcurve queries

|

Image Overlay Settings
----------------------
   **gaia_overlay_mag:** str = phot_g_mean_mag
      sets the default band to use in gaia detection overlays
   
   **galex_overlay_mag:** str = NUVmag
      sets the default band to use in gaia detection overlays

   **wise_overlay_mag:** str = W1mag
      sets the default band to use in gaia detection overlays
    
   **sdss_overlay_mag:** str = gPmag
      sets the default band to use in gaia detection overlays

   **twomass_overlay_mag:** str = Jmag
      sets the default band to use in gaia detection overlays

   **panstarrs_overlay_mag:** str = gmag
      sets the default band to use in gaia detection overlays

   **skymapper_overlay_mag:** str = gPSF
      sets the default band to use in gaia detection overlays

   **overlay_piggyback_radius:** int = 5
      sets the default 'piggyback' radius to use in arcseconds when proper motion correcting non-Gaia surveys in imaging overlays. Larger radius will result in more detections being corrected, but a larger chance of this correction being erroneous. 

|

Search Settings
---------------
   **seach_radius:** int = 3
      sets the default search radius in arcseconds for Vizier/SIMBAD searches using the search() tool

|

Datapage Settings
-----------------
   **datapage_search_button_radius:** int = 3
      sets the default search radius in arcseconds for Vizier/SIMBAD buttons from the Datapages module

   **datapage_datatable_radius:** int = 3
      sets the default search radius used for populating a datatable from the Datapages module

   **datapage_grid_size:** int = 250
      sets the default size of the grid in datapages

The Config
==========

You may have noticed that many specifics in previous tutorials have been chosen for us (e.g. the magnitudes which were used in :ref:`image overlays <example-image>` or the search radius of any of our queries). Because we have not explicitly provided these, they have been taken from the :ref:`config <Config>`. Upon first using ATK, a config file is generated with default values. This config contains a list of keys which set the defaults used throughout the package.

To have a look at the current state of the config (which is also the default state assuming you haven't edited it), we can print it to stdout via the :ref:`ATKshowconfig` command in the command line:

.. code-block:: console

    ATKshowconfig

.. code-block:: console
    
    Current ATKConfig.ini values:

    enable_notifications = True
    unit_size = 500
    output_backend = canvas
    show_toolbars = True
    show_grids = True
    show_titles = True
    font_size = 14
    font = Helvetica
    query_data_radius = 3
    query_phot_radius = 3
    query_bulkphot_radius = 3
    query_lightcurve_radius = 3
    query_spectrum_radius = 3
    query_sed_radius = 3
    query_reddening_radius = 3
    query_image_size = 30
    query_image_overlays = gaia
    query_image_band = g
    query_lightcurve_atlas_username = None
    query_lightcurve_atlas_password = None
    gaia_overlay_mag = phot_g_mean_mag
    galex_overlay_mag = NUVmag
    wise_overlay_mag = W1mag
    sdss_overlay_mag = gPmag
    twomass_overlay_mag = jmag
    skymapper_overlay_mag = gPSF
    panstarrs_overlay_mag = gmag
    overlay_piggyback_radius = 5
    overlay_simbad_search_radius = 3
    search_radius = 3
    datapage_search_button_radius = 3
    datapage_datatable_radius = 3
    datapage_font_size = 12
    datapage_grid_size = 250

A full description of each of these can be found :ref:`here <Config Keys>`, but hopefully some of them are self-explanatory.

Config keys can be edited from the command line via the :ref:`ATKeditconfig` command:

.. code-block:: console

    ATKeditconfig query_data_radius 5

.. code-block:: console

    Written change to ATKConfig.ini. New Values:

    enable_notifications = True
    unit_size = 500
    output_backend = canvas
    show_toolbars = True
    show_grids = True
    show_titles = True
    font_size = 14
    font = Helvetica
    query_data_radius = 5
    ...

Or the config can be viewed and edited at the same time by opening it in the default text editor:

.. code-block:: console

   ATKopenconfig

The config can be reset to its default state at any time via the :ref:`ATKresetconfig` command:

.. code-block:: console

    ATKresetconfig

.. code-block:: console

    Resetting ATKConfig.ini to default values...

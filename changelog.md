=======
Changes
=======
- Rewrote entire package from the ground up, should be a lot easier to develop in the future

- All coordinate/source handling now done via a central Target class
    - A Target is automatically generated if a source_id or SkyCoord is entered, or one can be manually created via Target.from_coord() or Target.from_id()
    - Reduces explicit dependency on Gaia, so future astrometric surveys like LSST will be a lot easier to implement, just need to provide an equivalent to get_gaia_target (i.e. converts LSST id -> ATK Target with LSST astrometry)
    - Added an "astrometric_backend" config key, which will (in future) select the default astrometry survey (currently only Gaia is supported)
    - Coordinates can be provided in any frame, automatically transformed to ICRS
- Configuration setup now far more robust
    - Improved cross-platform file opening
    - Config files now stored in hidden HOME/.AstroToolkit directory
- Improved .show() (previously .showdata()), now recursively handles the printing of arbitrarily complex structures in an improved format + can optionally show types via show_all_types = True
    - .show() is now supported on all ATK objects (e.g. DataSets, Containers, and Targets)
- Improved file saving to adaptively store structures as fits files
- Improved file reading to adaptively generate ATK structures from fits files
- Opened figures now save to HOME/.AstroToolkit/cached_figures/ temporarily, with a duration given by the config
- Image bands now correctly supported, and added 2 new colour maps - viridis (default) and false colour. Latter converts filter wavelength to an approximate real colour. Can choose colour map by passing cmap = 'viridis' | 'false_colour' | 'grey' to plot()
- Image plotting can now utilise relatives axes (i.e. +- arcsec from the centre, now default)
- Unified all structures into a single class DataSet
- Unified .data attribute - all query types now stored data as a list of ATK data containers
- SkyMapper image queries updated to SkyMapper DR4, and now sorted by exposure time (desc) and air mass (asc) to return best image
- Image queries to DSS1/DSS2 now properly implemented
- Added support for WISE, 2MASS and SDSS image queries 
- Improved image plotting
    - images should now always be near-perfect squares, regardless of plotting options
- Added an optional parameter "disable_corrections" to query(), which disable all astrometric corrections if True (defaults to False)
- SEDs now retain all detections in specified radius
- Added support for DESI DR1 spectral queries
- Improved filtering of ATLAS photometry
- Light curve queries now accept a kwarg "split", which will assign each returned light curve a per-survey object ID (where possible)
    - This defaults to True. Disabling it essentially makes every light curve survey a forced photometry (e.g. ATLAS) query
    - Light curve plotting now automatically splits light curves into per-survey and per-id figures
- Added kwarg 'cmap' to light curve plotting, which may be "mean" or "flat". The former uses a gradient colour map with a band of darker colour at the mean magnitude, while the latter uses a single flat colour per band. Defaults to "mean"
- User is now warned if ATK version does not match at time of fits file reading/writing
- Added integration with astropy units. Where relevant, parameters are treated and stored as astropy Quantities.
    - Input parameters may also use astropy quantities, e.g. a radius of 2 * u.arcmin may be requested for a 2 arcmin search
    - Added global config option 'unit_format' = 'text'/'symbol' to print units as text or symbol representations (defaults to 'symbol')
- Added query settings config option 'default_scale' = 'arcsec'/'arcmin'/'deg' to choose the default unit for query radius/image sizes (defaults to 'arcsec')
- Significantly reduced the number of dependencies
- added ability to query multiple targets at once
    - DataSets now store targets and a mapping between the query target and the returned data containers, added .fetch_by_id(), .fetch_by_coord() and .fetch_by_target() methods to extract data per-source
- added kwarg "background" to HRD plotting, which limits the background sample to a given fraction of the full sample (e.g. 0.5 would half the number of background points to reduced file size and lag)
- added kwarg 'combine' to HRD plotting, if False one HRD will be plotted for each source (primarily for datapage creation where each datapage will want its own hrd)
- added tab titles to plots when opened in browser
- added targeting information to figure titles (i.e. identifier or coordinates)
- improved SED/spectra overlay functionality
    - uses target matching to automatically match SEDs and spectra
    - plots with overlays now keep all interactive elements from both plot types
    - can now overlay both ways, i.e. can overlay spectra on an SED as before, but can also overlay an SED over spectra
    - overlaying multiple spectra for a single source (i.e. if a survey returned multiple spectra for the same source) now produces duplicated SEDs with a separate spectrum overlay for each one
- 'check_exists' kwarg in queries replaced by 'path' kwarg with same functionality
    - DataSets now use a checksum system -> changes in any query parameters will now automatically re-run the query and overwrite the local file
    - local file will also be overwritten if working ATK version doesn't match that at time of local file creation
- light curve binning now much faster
- usability of all data methods (e.g. lightcurve .bin(),.crop() etc.) significantly improved:
    - added method .apply() to DataSet, applies a given method to all stored containers
    - argument 'inplace' can be used to modify the structure/container in-place or perform modifications on a copy which is then returned
- powspec ('pspec') and phasefolding ('fold') now done via the .apply() method, allowing this data to be stored rather than being performed at plot-time
    - added 'multiband' kwarg to both of the above. If True, analyses all bands simultaneously to increase SNR and produce a single freq/period. If False, each band is treated entirely separately
    - added 'subtract' kwarg to phase folding. If 'median', all bands are median-subtracted. If 'mean', all bands are mean-subtracted. If None, no subtraction is performed.
- crop data method now supported by SEDs, spectra and power spectra
- bin data method now supported by spectra
- improved structure of user-accessible imports
    - main functions now available via: from ATK import ...
    - all ATK models available via: from ATK.Models import ...
- added spectral feature fitting functionality via plot method 'fit'
- added spectral radial velocity functionality via plot method 'rv_fit'
    - can detect multi-component radial velocities
- simplified datapage creation, layout now given as a list (columns) of lists (rows) of datasets
- optimized datapage plotting to minimise whitespace, align axes and force square images as best as possible
- greatly simplified DataTable creation
    - now just takes a dict in form survey: <list of cols> to include, automatically provides survey/catalogue/correction/separation columns and fetches parameter/value columns using Vizier queries
- DataTables now shrink in height to match table, and will grow in height up to requested size
- DataTables store data as an astropy Table under the .table attribute
    - when plotted, DataTables now have a unit column
- Datapages now returned as a DataPages object, with show(), show_by_target(), show_by_id(), and show_by_coords() methods to show all datapages or single out those of individual targets
- data is now saved via the .store() method, while plots are saved with the .save() method (the latter is also used for DataPages)
- Made spectral lines hidden by default
- Rewrote + significantly improved all documentation
    - Full tutorials added
    - Tutorials written with sphinx-gallery, now automatically run code and embed bokeh plots/terminal output, etc.
    - Tutorials now downloadable
- figure legends now dynamically shrink font size to remain inside bounds
- local files can now be compressed by providing a compressed file extension (e.g. .fits.gz)
- phase folding now uses phase dispersion minimisation to discern between the real frequency and its harmonics
- folded light curves can now be automatically aligned by passing align = "mean", "median", "max", "min" to align to the mean/median magnitude or a maxima/minima
- Vizier data containers (Records) now store their data as an astropy Table under the .table attribute
    - Returned Vizier data now maintains units
- old 'raw' argument in light curve queries is now 'filter', disables all non-required filtering if False (default=True)

+ other stuff that I forgot to write down



=========
To-Do Now
=========

Plotting
--------

Data
----
- talk to boris about my rv_fit process
- ztf light curve API not working
    - get better light curve sigma clipping example

- test all types of custom data set
- default units for Quantity arrays
    - needed to make sure .to() etc. doesn't fail

- sort container manipulations (to/from dataframe/table/etc), this can maybe wait but needs to at least work internally
- open by id, save by id etc.

Docs
====
- add a general note to the docs about how ATK implicitly converts SkyCoords and IDs to Target objects, and uses these to link data
- check phase folding in docs and finish this section, hopefully once ztf is actually working
- add note about TESS and ATLAS light curve filtering toggle
- rename tutorial .py files
- figures in docs need to scale with screen resolution, some html scaling thing and keep rest the same, maybe in plot_formatting under _utilities.py?
- check docstrings / comments
- check type hints
- finish docs

Other
=====
- check dependencies -> scipy?



===========
To-Do Later
===========
- clean up pm correction
    - main issue is e.g. overlays.py cannot be vectorised -> astropy doesn't let you define non-scalar SkyCoords with varying validity of distance/pm information, so you cannot reliably correct a grouped SkyCoord. Luckily, not used much as usually things are handled as single coordinates anyway
    - distance and pm ONLY matter for corrections, in searches its fine to lose this information and hence stay vectorised
- add a way to tell if any significant periodicity was detected in powspec
- pdm implementation
    - improve true frequency detection across range of orbital morphologies
- dataset operations, e.g. merge etc.
- gui/website (?)
- providing a SkyCoord with proper motion measurements
- clean up and improve generalisation of data methods (pass struct instead of arrays)
- look into using container methods on the containers themselves rather than via .apply on a DataSet
- rich text output option (https://realpython.com/python-rich-package/)
- calibrating + combining multiple light curves to make one massive light curve
- decouple from Gaia with a properly implemented astrometric backend system
- allow user to use config from within scripts, e.g. ATK.CONFIG[...][...] = ...
- see if I can get crts working, possibly a temporary outage
- add best-epoch separation to light curves?
- add matplotlib as an optional plotting backend to avoid issues with many data points (e.g. hrd/tess/power spectra)?
- don't change coordinates to icrs/celestial straight away, keep these and just change in query() without affecting Target.initial_coords?
- possibly store metadata, e.g. object IDs from light curve surveys in a .meta attribute
- add search_pos to SEDs?
- allow any magnitudes to be used in hrd query - i.e. user searches with Gaia ID to make sure there is a valid distance (if not return no data) -> supplied two bands are fetched from a Vizier query to get the colour?
- Let Vizier catalogue names be used in overlays (as a fallback if alias not in alias file)
- move most globals (or things that need to be edited on occasion) to one place (?) + include global prefix to make sure that these aren't edited (?)
- sort defaults for kwargs (should be in function definitions/config - or somewhere else, not as default arg in kwargs.get()) (?)
- used ValueError too much, should only be used for argument errors
- struct vs ctnr vs etc.
- option to only return closest photometry from each survey in SED queries
- recursive show_types=True in .show()
- light curve overlays in images as way to show how to recombine data
    - more generally, be able to pass any structure to image plotting -> extracts positional data + overlay
- improve warnings/logging (no print statements?)
- add annotations to ATK keywords in fits headers
- let spectral peak fitting work in velocity-space (?)

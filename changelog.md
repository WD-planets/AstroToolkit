Changes
-------
- Rewrote entire package, should be a lot easier to develop in the future

- All coordinate/source handling now done via a central Target class
    - A Target is automatically generated if a source_id or SkyCoord is entered, or one can be manually created via Target.from_pos() or Target.from_id()
    - Reduces explicit dependency on Gaia, so future astrometric surveys like LSST will be a lot easier to implement, just need to provide an equivalent to get_gaia_target (i.e. converts LSST id -> ATK Target with LSST astrometry)
    - Added an "astrometric_backend" config key, which will (in future) select the default astrometry survey (currently only Gaia)
    - Coordinates can be provided in any frame, automatically transformed to ICRS
- Configuration setup now far more robust
    - Added ATKoverlay show
    - Improved cross-platform file opening
    - Config files now stored in hidden HOME/.AstroToolkit directory
- Improved .show() (previously .showdata()), now recursively handles the printing of arbitrarily complex structures in an improved format + can optionally show types via show_all_types = True
    - .show() is now supported on all ATK objects (e.g. QueryResults, Containers, and Targets)
- Improved file saving to adaptively store structures as fits files
- Improved file reading to adaptively generate ATK structures from fits files
- Opened figures now save to HOME/.AstroToolkit/cached_figures/ temporarily, with a duration given by the config
- Images bands now correctly supported, and added 2 new colour maps - viridis (default) and false colour. Latter converts filter wavelength to an approximate real colour. Can choose colour map by passing cmap = 'viridis' | 'false_colour' | 'grey' to plot()
- Image plotting can now use relatives axes (i.e. +- arcsec from the centre)
- Unified all structures into a single class BaseQueryResult, from which QueryResult and PlottableQueryResult inherit
- Unified .data attribute - all query types now stored data as a list of pandas DataFrames (vizier queries) or new ATK data objects (basically everything else).
- SkyMapper image queries updated to SkyMapper DR4, and now sorted by exposure time (desc) and air mass (asc) to return best image
- Image queries to DSS1/DSS2 now properly implemented
- Added support for WISE, 2MASS and SDSS image queries 
- Improved image plotting
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
- Added query settings config option 'default_unit' = 'arcsec'/'arcmin'/'deg' to choose the default unit for query radius/image sizes (defaults to 'arcsec')
- Significantly reduced the number of dependencies
- added ability to query multiple targets at once
    - QueryResults now store targets and a mapping between the query target and the returned data containers, added .fetch_by_id() and .fetch_by_coord() methods to extract data per-source

To-Do Now
---------
- Let Vizier catalogue names be used in overlays (as a fallback if alias not in alias file)
- move most globals (or things that need to be edited on occasion) to one place (?) + include global prefix to make sure that these aren't edited (?)
- try to remove unnecessary dependencies
    - reproject
- add unit tests for each survey (known working examples to check if survey is not working, can auto run thes on exception optionally) these should also save and read a file to test this
    - maybe add this as an option for the user - i.e. if an exception occurs and the unit test then fails, retry every ~ 5 mins (not too much traffic, only intended for large studies)
- include distances in overlay corrections
- add survey ID to SED hovertool
- add spectral line fitting tool
- possibly store metadata, e.g. object IDs from light curve surveys in a .meta attribute
- properly sort warnings/logging (no print statements?)
- add filter kwarg to light curve queries to disable all unrequired filtering
- add annotation to ATK keywords in fits headers
- sort defaults for kwargs (should be in function definitions/config - or somewhere else, not as default arg in kwargs.get())
- check docstrings / comments
- check type hints, especially for astropy quantities after change was made
- consider turning off split = True for ZTF lightcurves as default, too many "objects"
- default units for Quantity arrays?


QueryResults should store a target, not positions and identifiers separately -> print as e.g. "123.456 12.345 (identifier if present)"
    - now these are coupled together, removes complexity of some targets having an identifier while others do not

identifier/position should go into container names if one is being stored, e.g. <1234567 gaia G vs bp-rp HRD>
store position/identifier of each container in header + reconstruct individually, don't need query position/identifier here


To-Do Later
-----------
- calibrating + combining multiple light curves to make one massive light curve
- decouple from Gaia with a properly implemented astrometric backend system
- allow user to use config from within scripts, e.g. ATK.CONFIG[...][...] = ...
- see if I can get crts working, possibly a temporary outage
- add best-epoch separation to light curves
- add matplotlib as an optional plotting backend to avoid issues with many data points (e.g. hrd/tess/power spectra)

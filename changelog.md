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

To-Do 
-----
- Let Vizier catalogue names be used in overlays (as a fallback if alias not in alias file)
- use astropy units to define axes
- move most globals to one place + include global prefix to make sure that these aren't edited
- try to remove unnecessary dependencies
    - reproject
- Make an actual_bounds() function for Image
- test file saving on vizier and image data
- add colours to image plots
- add unit tests for each survey (known working examples to check if survey is not working, can auto run thes on exception optionally)
- use astropy units in query radius/size, assume arcsec if no unit given but accept other units + convert
- decouple from Gaia with a properly implemented astrometric backend system
- include distances in overlay corrections
- add survey ID to SED hovertool
- add spectral line fitting tool
- possibly store metadata, e.g. object IDs from light curve surveys in a .meta attribute
- what happens if plotting no data


- calibrating + combining multiple light curves to make one massive light curve

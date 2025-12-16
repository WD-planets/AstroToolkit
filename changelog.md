Changes
-------
- Rewrote entire package, should be a lot easier to develop in the future

- Configuration setup now far more robust
    - Added ATKoverlay show
    - Improved cross-platform file opening
    - Config files now stored in hidden HOME/.AstroToolkit directory
- Improved .show() (previously .showdata()), now recursively handles the printing of arbitrarily complex structures in an improved format + can optionally show types via show_all_types = True
- Improved file saving to adaptively store structures as fits files
- Improved file reading to adaptively generate ATK structures from fits files
- Opened figures now save to HOME/.AstroToolkit/cached_figures/ temporarily, with a duration given by the config
- Images bands now correctly supported, and added 2 new colour maps - viridis (default) and false colour. Latter converts filter wavelength to an approximate real colour. Can choose colour map by passing cmap = 'viridis' | 'false_colour' | 'grey' to plot()
- Unified all structures into a single class BaseQueryResult, from which QueryResult and PlottableQueryResult inherit
- Unified .data attribute - all query types now stored data as a list of pandas DataFrames (vizier queries) or new ATK data objects (basically everything else).
- Image queries to DSS1/DSS2 now properly implemented
- Added support for WISE, 2MASS and SDSS image queries 
- Improved image plotting

To-Do 
-----
- Let Vizier catalogue names be used in overlays (as a fallback if alias not in alias file)
- use astropy units to define axes
- move most globals to one place + include global prefix to make sure that these aren't edited

Changes
-------
- Rewrote entire package, should be a lot easier to develop in the future
- Configuration setup now far more robust
    - Added ATKoverlay show
    - Improved cross-platform file opening
- Improved .show() (previously .showdata()), now recursively handles the printing of arbitrarily complex structures

To-Do 
-----
- Let Vizier catalogue names be used in overlays (as a fallback if alias not in alias file)
- use astropy units to define axes

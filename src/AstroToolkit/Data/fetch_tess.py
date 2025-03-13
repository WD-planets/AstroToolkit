# Adapted from code given by Keith
def get_TESS_lightcurve(pos, radius):
    """
    Note on converting to magnitudes:

    https://heasarc.gsfc.nasa.gov/docs/tess/observing-technical.html
    says that 15000 e/s gives a magnitude of 10

    Tess returns data in astropy quantities - flux in e/s

    u.Magnitude(lc['flux']) is the same as -2.5*log10(lc['flux'])

    so plugging in 15000 and m=10 gives zp=20.44*u.ABmag
    """

    import warnings

    import astropy.units as u
    import lightkurve as lk
    import numpy as np
    from astropy.table import Table, vstack

    ra, dec = pos

    Search = lk.search_lightcurve(
        str(ra) + " " + str(dec), mission="TESS", radius=radius
    )

    if len(Search) > 0:
        lcurve = []
        for s in Search:
            try:
                lqq = s.download_all()
                lcurve += lqq
            except:
                print("Note: experiencing issues with TESS.")
        sectors = []
        LC = Table([])
        for lc in lcurve:
            if lc.meta["SECTOR"] not in sectors:
                sectors += [lc.meta["SECTOR"]]
                warnings.simplefilter("ignore")
                lc = lc[lc["flux"] > 0]
                lc["mag"] = u.Magnitude(lc["flux"]) + 20.44 * u.ABmag
                lc["mjd"] = lc["time"].mjd.astype("float")
                try:
                    tab = lc.to_table()[
                        [
                            "time",
                            "flux",
                            "flux_err",
                            "sap_flux",
                            "quality",
                            "mag",
                            "mjd",
                        ]
                    ]
                except:
                    tab = lc.to_table()[
                        ["time", "flux", "flux_err", "quality", "mag", "mjd"]
                    ]
                    tab["sap_flux"] = 0
                has_nan = np.zeros(len(tab), dtype=bool)
                for col in tab.itercols():
                    if col.info.dtype.kind in ["f", "i"]:
                        has_nan |= np.isnan(col)
                tab_nonan = tab[~has_nan]
                if len(LC) == 0:
                    LC = tab_nonan
                else:
                    LC = vstack([LC, tab_nonan])
        return LC.to_pandas()
    else:
        print("Note: TESS query returned no data.")
        return None

class LightcurveStruct(object):
    """LightcurveStruct()
    This structure is returned from lightcurve queries, when read from a data file that was originally created by a lightcurve query, or through the Models module (in which case all attributes are set to None).

    Attributes
    ----------

    kind : str
        "lightcurve"

    survey: str
        survey to which the data pertains, see lightcurve query description for details.

    source: int
        Gaia source ID (if one was used to target the system, else None)

    pos: list<int>
        Target position [right ascension, declination] in degrees

    identifier: str
        Target position in JHHMMSS.SS±DDMMSS.SS format

    dataname: str
        Default file name that data will be saved to when using the savedata() method

    plotname: str
        Default file name that plot will be saved to when using the showplot() and saveplot() methods

    figure: None
        Stores any figure resulting from the plot() method

    data: list<dict>
        Returned data as a list of bands, with each band in format:

        .. code-block:: python

            "band": band
            "ra": ra
            "dec": dec
            "hjd"/"mjd": hjd
            "hjd_ori"/"hjd_ori": hjd_ori
            "mag": mag
            "mag_err": mag_err

    where:

        :**band**: str, band of data entry

        :**ra**: list<float>, right ascension in degrees of each detection

        :**dec**: list<float>, declination in degrees of each detection

        :**hjd or mjd**: list<float>, hjd or mjd of detection, reduced to start at zero

        :**hjd_ori or mjd_ori**: list<float>, non-reduced hjd or mjd of detection

        :**mag**: list<float>, magnitude of detection

        :**mag_err**: list<float>, error on the magnitude

    and the bands for each supported light curve survey are:

        - ztf: g, r, i
        - atlas: o, c, i
        - gaia: g, bp, rp
        - asassn: g, v
        - crts: v
        - tess: i

    |

    """

    def __init__(self, survey, source, pos, data, identifier=None):
        self.kind = "lightcurve"
        self.survey = survey
        self.source = source
        self.pos = pos
        self.identifier = identifier
        self.data = data
        self.figure = None
        self.dataname = None
        self.figure = None

    def __str__(self):
        return "<ATK Lightcurve Structure>"

    def plot(
        self,
        kind="lightcurve",
        colours=None,
        bands=None,
        freq=None,
        bins=None,
        timeformat="reduced",
        method="ls",
        foverlay=True,
        repeat=1,
        shift=0,
    ):
        """
        plot(kind, *)
        Plots light curve data and assigns the resulting figure to the **figure** attribute.

        :param str kind: type of figure to create, from: lightcurve, powspec (power spectrum), phasefold (phase-folded light curve)

        There are then additional optional parameters depending on the kind of figure requested:

        For lightcurve plotting:

        :param list<str> bands: list of bands to plot, see listed bands for each supported survey above, defaults to all supported bands
        :param list<str> colours: list of colours to apply to selected bands. Supported colours: green, red, blue, black, orange, purple. Default =  black for all
        :param str timeformat: time format from: reduced, original. Default = reduced

        |

        For powspec plotting:

        :param str method: time series analysis routine to use, from:

            - ls (Lomb-Scargle)
            - amhw (Multiharmonic analysis of variance)
            - aovw (Analysis of variance for phase bins)
            - atrw (Analysis of variance for planetary transits/eclipses)
            - pspw (Discrete power spectrum)
            - f_mw (Generalised Lomb_Scargle)
            - lomw (Analysis of variance)

        **Note:** all time series analysis routines except for ls require PyAOV to be setup. See setup module for more information.

        |

        For phasefold plotting:

        :param float freq: frequency on which to fold the data in :math:`\\text{days}^{-1}`, defaults to calcualted Lomb-Scargle frequency
        :param int bins: number of equally spaced bins in time to bin phase folded data into, default = no binning
        :param int repeat: number of repetitions to plot on the folded period. Default = 2
        :param float shift: shift to apply to phase folded data in units of phase, default = 0
        :param bool foverlay: overlays a sine wave of the fold frequency, default = True

        |

        """

        from ..Plotting.plotmap import map_to_plot

        return map_to_plot(
            self,
            kind=kind,
            colours=colours,
            bands=bands,
            freq=freq,
            bins=bins,
            timeformat=timeformat,
            method=method,
            foverlay=foverlay,
            repeat=repeat,
            shift=shift,
        )

    def sigmaclip(self, sigma=3):
        """
        Sigma clips light curve data.

        :param float sigma: number of standard deviations beyond which to clip data
        :return: self
        :rtype: class

        |

        """

        from .sigmaclip import sigma_clip

        if self.data:
            self.data = [
                sigma_clip(band, sigma) if band["mag"] else band for band in self.data
            ]
        return self

    def bin(self, bins=None, binsize=None):
        """
        Bins light curve data into a given number of equally spaced bins in time or a given bin size in days, hours or minutes.

        :param int bins: number of bins in which to bin data

        or

        :param str binsize: bin size ending with the unit, e.g. "10d", "10h" or "10m" for 10 days, hours or minutes, respectively

        :return: self
        :rtype: class

        |

        """

        from .bin_lightcurves import binning

        if not bins and not binsize:
            raise Exception("Must provide bins or binsize.")
        elif bins and binsize:
            raise Exception("Both bins and binsize provided.")

        if self.data:
            self.data = [
                binning(band, bins, binsize) if band["mag"] else band
                for band in self.data
            ]
        return self

    def crop(
        self,
        start=None,
        stop=None,
        start_percent=None,
        stop_percent=None,
        timeformat="reduced",
    ):
        """crop(*)
        Crops light curve data between a given start and end time in reduced or original time format, or between two percentages of the total coverage.

        :param str timeformat: time format, from: reduced, original
        :param float start: start time in format given by **timeformat**
        :param float stop: end time in format given by **timeformat**

        or

        :param float start_percent: percentage of the total time covered by the light curve at which to start crop
        :param float stop_percent: percentage of the total time covered by the light curve at which to end crop

        :return: self
        :rtype: class

        |

        """

        from .crop_lightcurves import crop_lightcurve

        self.data = crop_lightcurve(
            self, start, stop, start_percent, stop_percent, timeformat
        )
        return self

    def savedata(self, fname=None):
        """
        Saves data to local files. Default file name = dataname attribute.

        :param str fname: overrides file name
        :return: name of file to which data was saved
        :rtype: str

        |

        """

        from .data_saving import savedata

        fname = savedata(self, fname)
        return fname

    def showdata(self, raw=False):
        """
        Prints data structure to stdout in a readable format.

        :param bool raw: if True, doesn't collapse arrays to improve readability, default = False
        :return: self
        :rtype: class

        |

        """

        from .data_printing import print_data

        print_data(self, raw)
        return self

    def showplot(self, fname=None):
        """
        Opens the figure stored in the **figure** attribute in the default web browser, and saves it to local files.

        :param str fname: file name to save the figure to, default file name = plotname attribute
        :return: file name to which the figure was saved
        :rtype: str

        |

        """

        from ..Plotting.plotmap import showplot

        fname = showplot(self, fname=fname)
        return fname

    def saveplot(self, fname=None):
        """
        Saves the figure stored in the **figure** attribute to local files without opening it in the web browser.

        :param str fname: file name to save the figure to, default file name = plotname attribute
        :return: file name to which the figure was saved
        :rtype: str

        |

        """

        from ..Plotting.plotmap import saveplot

        fname = saveplot(self, fname=fname)
        return fname


class GeneralQuery(object):
    def __init__(self, survey, radius, pos, username, password, pmra, pmdec, raw):
        self.survey = survey
        self.radius = radius
        self.pos = pos
        self.url = ""
        self.username = username
        self.password = password
        self.pmra = pmra
        self.pmdec = pmdec
        self.raw = raw

    def get_response(self):
        import requests
        from requests.adapters import HTTPAdapter, Retry

        s = requests.Session()
        retries = Retry(
            total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504]
        )
        s.mount("http://", HTTPAdapter(max_retries=retries))

        try:
            response = s.get(self.url, timeout=180)
        except:
            print(f"Note: experiencing issues with {self.survey}")
            return None

        if response.status_code != 200:
            print(f"Note: experiencing issues with {self.survey}")
            return None

        return response


class ZtfQuery(GeneralQuery):
    def set_url(self):
        url = f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE {self.pos[0]} {self.pos[1]} {self.radius/3600}&BANDNAME=g,r,i&FORMAT=CSV"
        self.url = url

    def generate_data(self, response):
        from io import BytesIO

        import pandas as pd

        data = pd.read_csv(BytesIO(response.content))
        if len(data) == 0:
            print(f"Note: {self.survey} query returned no data.")
            return None

        return data

    def format_data(self, data):
        data_arr = []
        for filter_code in ["zg", "zr", "zi"]:
            current_band_data = data.loc[data["filtercode"] == filter_code].reset_index(
                drop=True
            )
            if current_band_data.empty:
                data_arr.append(
                    {
                        "band": filter_code[1:],
                        "ra": None,
                        "dec": None,
                        "mjd": None,
                        "mjd_ori": None,
                        "mag": None,
                        "mag_err": None,
                    }
                )
                continue

            mag, mag_err, ra, dec, hjd_ori = (
                current_band_data["mag"].tolist(),
                current_band_data["magerr"].tolist(),
                current_band_data["ra"].tolist(),
                current_band_data["dec"].tolist(),
                current_band_data["hjd"].tolist(),
            )
            hjd = [x - min(hjd_ori) for x in hjd_ori]
            data_arr.append(
                {
                    "band": filter_code[1:],
                    "ra": ra,
                    "dec": dec,
                    "hjd": hjd,
                    "hjd_ori": hjd_ori,
                    "mag": mag,
                    "mag_err": mag_err,
                }
            )

        return data_arr


class AtlasQuery(GeneralQuery):
    def set_url(self):
        base_url = "https://fallingstar-data.com/forcedphot"
        self.url = base_url

    def get_response(self):
        import re
        import time

        import requests

        r = requests.post(
            url=f"{self.url}/api-token-auth/",
            data={"username": self.username, "password": self.password},
        )
        if r.status_code == 200:
            token = r.json()["token"]
            headers = {"Authorization": f"Token {token}", "Accept": "application/json"}
            self.headers = headers
        else:
            print(
                f"Note: experiencing issues with {self.survey}. Please ensure that you have provided valid login credentials."
            )
            return None

        do_correction = False
        if self.pmra and self.pmdec:
            do_correction = True

        task_url = None
        while not task_url:
            with requests.Session() as s:
                if do_correction:
                    r = s.post(
                        f"{self.url}/queue/",
                        headers=headers,
                        data={
                            "ra": self.pos[0],
                            "dec": self.pos[1],
                            "mjd_min": 50000.0,
                            "propermotion_ra": self.pmra,
                            "propermotion_dec": self.pmdec,
                            "radec_epoch_year": 2000,
                            "use_reduced": True,
                        },
                    )
                else:
                    r = s.post(
                        f"{self.url}/queue/",
                        headers=headers,
                        data={
                            "ra": self.pos[0],
                            "dec": self.pos[1],
                            "mjd_min": 50000.0,
                            "radec_epoch_year": 2000,
                            "use_reduced": True,
                        },
                    )

                if r.status_code == 201:
                    task_url = r.json()["url"]
                elif r.status_code == 429:
                    message = r.json()["detail"]
                    print(f"{r.status_code} {message}")
                    t_sec = re.findall(r"Available in (\d+) seconds", message)
                    t_min = re.findall(r"Available in (\d+) minutes", message)

                    if t_sec:
                        waittime = int(t_sec[0])
                    elif t_min:
                        waittime = int(t_min[0])
                    else:
                        waittime = 10

                        print(f"Waiting {waittime} seconds.")
                        time.sleep(waittime)
                else:
                    print(
                        f"{self.survey} login failed. Please ensure that you have provided valid login credentials."
                    )
                    return None

        result_url = None
        while not result_url:
            with requests.Session() as s:
                r = s.get(task_url, headers=headers)
                if r.status_code == 200:
                    if r.json()["finishtimestamp"]:
                        return r
                    elif r.json()["starttimestamp"]:
                        print(
                            f"Task is running (started at {r.json()['starttimestamp']})"
                        )
                    else:
                        print("Waiting for job to start. Checking again in 10 seconds.")
                    time.sleep(10)
                else:
                    print(f"Note: experiencing issues with {self.survey}")
                    return None

    def generate_data(self, response):
        import io

        import pandas as pd
        import requests

        result_url = response.json()["result_url"]

        with requests.Session() as s:
            textdata = s.get(result_url, headers=self.headers).text

        df = pd.read_csv(io.StringIO(textdata.replace("###", "")), sep="\\s+")
        if df.empty:
            print("Note: ATLAS query returned no data.")
            return None
        else:
            return df

    def format_data(self, data):
        import itertools

        mjd_ori = data["MJD"].tolist()
        mjd = [x - min(mjd_ori) for x in mjd_ori]
        mag = data["m"].tolist()
        mag_err, ra_arr, dec_arr, obs_arr, chi_arr, flux, flux_err = (
            data["dm"].tolist(),
            data["RA"].tolist(),
            data["Dec"].tolist(),
            data["Obs"].tolist(),
            data["chi/N"].tolist(),
            data["uJy"].tolist(),
            data["duJy"].tolist(),
        )

        if not self.raw:
            chi_arr = [float(x) for x in chi_arr]
            flux = [float(x) for x in flux]
            flux_err = [float(x) for x in flux_err]

            # abs(flux) < 3
            bad_indices_flux = set(
                [i for i, element in enumerate(flux) if abs(element) < 3]
            )
            # flux_err > 4000
            bad_indices_flux_err = set(
                [i for i, element in enumerate(flux_err) if element > 4000]
            )

            """
            # mag_err > 0.5 * mag_range
            bad_indices_mag = set(
                [
                    i
                    for i, element in enumerate(mag_err)
                    if element > 0.5 * (max(mag) - min(mag))
                ]
            )
            """

            bad_indices_mag_hard_lim = set(
                [i for i, element in enumerate(mag_err) if element > 0.5]
            )

            bad_indices = list(
                set(
                    itertools.chain(
                        *[
                            bad_indices_flux,
                            bad_indices_flux_err,
                            bad_indices_mag_hard_lim,
                        ]
                    )
                )
            )

            # apply the above filters
            filtered_list = []
            for _, val in enumerate(
                [obs_arr, mag, mag_err, mjd, mjd_ori, ra_arr, dec_arr]
            ):
                filtered_list.append(
                    [element for i, element in enumerate(val) if i not in bad_indices]
                )
            obs_arr, mag, mag_err, mjd, mjd_ori, ra_arr, dec_arr = filtered_list

            filtered_list = []
            mask = [i for i, val in enumerate(mag) if val < 0]
            for _, val in enumerate(
                [obs_arr, mag, mag_err, mjd, mjd_ori, ra_arr, dec_arr]
            ):
                filtered_list.append(
                    [element for i, element in enumerate(val) if i not in mask]
                )
            obs_arr, mag, mag_err, mjd, mjd_ori, ra_arr, dec_arr = filtered_list

            filtered_list = []
            mask = [i for i, val in enumerate(mag_err) if val == 0]
            for _, val in enumerate(
                [obs_arr, mag, mag_err, mjd, mjd_ori, ra_arr, dec_arr]
            ):
                filtered_list.append(
                    [element for i, element in enumerate(val) if i not in mask]
                )
            obs_arr, mag, mag_err, mjd, mjd_ori, ra_arr, dec_arr = filtered_list

        # split data into bands
        data_split = []
        for band in ["o", "c", "i"]:
            split_arr = []
            for _, val in enumerate([ra_arr, dec_arr, mjd, mjd_ori, mag, mag_err]):
                split_arr.append(
                    [element for i, element in enumerate(val) if obs_arr[i][-1] == band]
                )
            if len(split_arr[4]) > 0:
                data_split.append(
                    {
                        "band": band,
                        "ra": split_arr[0],
                        "dec": split_arr[1],
                        "mjd": split_arr[2],
                        "mjd_ori": split_arr[3],
                        "mag": split_arr[4],
                        "mag_err": split_arr[5],
                    }
                )
            else:
                data_split.append(
                    {
                        "band": band,
                        "ra": None,
                        "dec": None,
                        "mjd": None,
                        "mjd_ori": None,
                        "mag": None,
                        "mag_err": None,
                    }
                )

        return data_split


class GaiaQuery(GeneralQuery):
    def get_data(self, source=None, pos=None):
        from ..Tools import query

        data = query(
            survey="gaia_lc",
            kind="data",
            pos=pos,
            source=source,
            radius=self.radius,
            level="internal",
        )
        if data:
            return data.data

    def format_data(self, data):
        import math

        import numpy as np
        from astropy.time import Time

        # magnitude error formula taken from https://astronomy.stackexchange.com/questions/38371/how-can-i-calculate-the-uncertainties-in-magnitude-like-the-cds-does
        bands = ["g", "bp", "rp"]
        cols = []

        data_arr = []
        for band in bands:
            if band == "g":
                # time, flux, flux_err, mag
                cols = [
                    "_tab8_5",
                    "_tab8_6",
                    "_tab8_7",
                    "_tab8_9",
                    "RA_ICRS",
                    "DE_ICRS",
                ]
            elif band == "bp":
                # time, flux, flux_err, mag
                cols = [
                    "_tab8_11",
                    "_tab8_12",
                    "_tab8_13",
                    "_tab8_15",
                    "RA_ICRS",
                    "DE_ICRS",
                ]
            elif band == "rp":
                # time, flux, flux_err, mag
                cols = [
                    "_tab8_16",
                    "_tab8_17",
                    "_tab8_18",
                    "_tab8_20",
                    "RA_ICRS",
                    "DE_ICRS",
                ]

            time, flux, flux_err, mag, ra, dec = (
                data[cols[0]],
                data[cols[1]],
                data[cols[2]],
                data[cols[3]],
                data[cols[4]],
                data[cols[5]],
            )

            mag_err = [(2.5 / np.log(10)) * (f_e / f) for f, f_e in zip(flux, flux_err)]

            bad_indices = [i for i, val in enumerate(time) if math.isnan(val)]

            # apply the above filter
            filtered_list = []
            for _, val in enumerate([mag, mag_err, time, ra, dec]):
                filtered_list.append(
                    [element for i, element in enumerate(val) if i not in bad_indices]
                )
            mag, mag_err, time, ra, dec = filtered_list

            mjd_ori = [Time(t + 2455197.5, format="jd").value for t in time]
            mjd = [t - min(mjd_ori) for t in mjd_ori]

            if len(mag) > 0:
                data_arr.append(
                    {
                        "band": band,
                        "ra": ra,
                        "dec": dec,
                        "mjd": mjd,
                        "mjd_ori": mjd_ori,
                        "mag": mag,
                        "mag_err": mag_err,
                    }
                )
            else:
                data_arr.append(
                    {
                        "band": band,
                        "ra": None,
                        "dec": None,
                        "mjd": None,
                        "mjd_ori": None,
                        "mag": None,
                        "mag_err": None,
                    }
                )

        return data_arr


class AsassnQuery(GeneralQuery):
    def set_url(self):
        url = f"https://asas-sn.osu.edu/photometry.json?action=index&controller=photometry&dec={self.pos[1]}&epochs_max=&epochs_min=&ra={self.pos[0]}&radius={self.radius/60}&rms_max=&rms_min=&sort_by=raj2000&utf8=%E2%9C%93&vmag_max=&vmag_min="
        self.url = url

    def generate_data(self, response):
        import json

        data = json.loads(response.content)
        if data["count"] > 1:
            print(f"Note: {self.survey} query returned data for multiple objects")
        elif data["count"] == 0:
            print("Note: asassn query returned no data.")
            return None

        return data

    def format_data(self, data):
        import json

        import requests

        v_ra_arr, v_dec_arr, v_mag, v_mag_err, v_hjd_ori = [], [], [], [], []
        g_ra_arr, g_dec_arr, g_mag, g_mag_err, g_hjd_ori = [], [], [], [], []

        for i in range(0, len(data["results"])):
            link = data["results"][i]["link"]

            ra = data["results"][i]["raj2000"]
            dec = data["results"][i]["dej2000"]

            v_cameras = ["a", "b", "c", "d", "e", "f", "g", "h"]
            g_cameras = ["i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t"]

            complete = False
            while not complete:
                r = requests.get(link, timeout=60)

                lc_data = json.loads(r.content)

                for j in range(0, len(lc_data["results"])):
                    current_data = lc_data["results"][j]
                    camera = current_data["camera"][-1:]
                    if camera in v_cameras:
                        v_mag.append(current_data["mag"])
                        v_mag_err.append(current_data["mag_err"])
                        v_hjd_ori.append(current_data["hjd"])
                        v_ra_arr.append(ra)
                        v_dec_arr.append(dec)
                    elif camera in g_cameras:
                        g_mag.append(current_data["mag"])
                        g_mag_err.append(current_data["mag_err"])
                        g_hjd_ori.append(current_data["hjd"])
                        g_ra_arr.append(ra)
                        g_dec_arr.append(dec)

                if lc_data["next"] is not None:
                    link = lc_data["next"]
                else:
                    complete = True

        v_hjd = [t - min(v_hjd_ori) for t in v_hjd_ori]
        g_hjd = [t - min(g_hjd_ori) for t in g_hjd_ori]

        if not self.raw:
            v_bad_indices = set(
                [
                    i
                    for i, element in enumerate(v_mag_err)
                    if element > 0.5 * (max(v_mag) - min(v_mag))
                ]
            )
            g_bad_indices = set(
                [
                    i
                    for i, element in enumerate(g_mag_err)
                    if element > 0.5 * (max(g_mag) - min(g_mag))
                ]
            )

            # apply the above filters
            filtered_list = []
            for _, val in enumerate(
                [v_mag, v_mag_err, v_hjd, v_hjd_ori, v_ra_arr, v_dec_arr]
            ):
                filtered_list.append(
                    [element for i, element in enumerate(val) if i not in v_bad_indices]
                )
            v_mag, v_mag_err, v_hjd, v_hjd_ori, v_ra_arr, v_dec_arr = filtered_list

            # apply the above filters
            filtered_list = []
            for _, val in enumerate(
                [g_mag, g_mag_err, g_hjd, g_hjd_ori, g_ra_arr, g_dec_arr]
            ):
                filtered_list.append(
                    [element for i, element in enumerate(val) if i not in g_bad_indices]
                )
            g_mag, g_mag_err, g_hjd, g_hjd_ori, g_ra_arr, g_dec_arr = filtered_list

        data_arr = []
        if len(v_mag) > 0:
            data_arr.append(
                {
                    "band": "v",
                    "ra": v_ra_arr,
                    "dec": v_dec_arr,
                    "hjd": v_hjd,
                    "hjd_ori": v_hjd_ori,
                    "mag": v_mag,
                    "mag_err": v_mag_err,
                }
            )
        else:
            data_arr.append(
                {
                    "band": "v",
                    "ra": None,
                    "dec": None,
                    "mjd": None,
                    "mjd_ori": None,
                    "mag": None,
                    "mag_err": None,
                }
            )

        if len(g_mag) > 0:
            data_arr.append(
                {
                    "band": "g",
                    "ra": g_ra_arr,
                    "dec": g_dec_arr,
                    "hjd": g_hjd,
                    "hjd_ori": g_hjd_ori,
                    "mag": g_mag,
                    "mag_err": g_mag_err,
                }
            )
        else:
            data_arr.append(
                {
                    "band": "g",
                    "ra": None,
                    "dec": None,
                    "mjd": None,
                    "mjd_ori": None,
                    "mag": None,
                    "mag_err": None,
                }
            )

        return data_arr


class CrtsQuery(GeneralQuery):
    def get_data(self):
        from .fetch_crts import get_CRTS_lightcurve

        data = get_CRTS_lightcurve(self.pos, self.radius)
        return data

    def format_data(self, data):
        ra_arr = data["RA"].tolist()
        dec_arr = data["Dec"].tolist()
        mjd_ori = data["MJD"].tolist()
        mjd = [t - min(mjd_ori) for t in mjd_ori]
        mag = data["Mag"].tolist()
        mag_err = data["Magerr"].tolist()

        data_dict = {
            "band": "v",
            "ra": ra_arr,
            "dec": dec_arr,
            "mjd": mjd,
            "mjd_ori": mjd_ori,
            "mag": mag,
            "mag_err": mag_err,
        }

        return [data_dict]


class TessQuery(GeneralQuery):
    def get_data(self):
        from .fetch_tess import get_TESS_lightcurve

        data = get_TESS_lightcurve(self.pos, self.radius)
        return data

    def format_data(self, data):
        import math

        import numpy as np

        mjd_ori = data["mjd"].tolist()
        mjd = [t - min(mjd_ori) for t in mjd_ori]
        mag = data["mag"].tolist()
        ra_arr = [self.pos[0]] * len(mag)
        dec_arr = [self.pos[1]] * len(mag)
        flux = data["flux"].tolist()
        flux_err = data["flux_err"].tolist()

        mag_err = (
            (2.5 / np.log(10)) * (np.asarray(flux_err) / np.asarray(flux))
        ) ** 2 + 0.05**2
        mag_err = [math.sqrt(x) for x in mag_err]

        if not self.raw:
            mask = [i for i, val in enumerate(mag_err) if val > 1]

            def filter_data(arr, mask):
                filtered_arr = [val for i, val in enumerate(arr) if i not in mask]
                return filtered_arr

            mjd_ori = filter_data(mjd_ori, mask)
            mjd = filter_data(mjd, mask)
            mag = filter_data(mag, mask)
            mag_err = filter_data(mag_err, mask)
            ra_arr = filter_data(ra_arr, mask)
            dec_arr = filter_data(dec_arr, mask)

        data_dict = {
            "band": "i",
            "ra": ra_arr,
            "dec": dec_arr,
            "mjd": mjd,
            "mjd_ori": mjd_ori,
            "mag": mag,
            "mag_err": mag_err,
        }

        return [data_dict]


def get_f_return(survey):
    f_return = []

    from .dataquery import SurveyInfo

    bands = SurveyInfo().lightcurve_bands[survey]
    for band in bands:
        f_return.append(
            {
                "band": band,
                "ra": None,
                "dec": None,
                "mjd": None,
                "mjd_ori": None,
                "mag": None,
                "mag_err": None,
            }
        )

    return f_return


def query(survey, source, pos, radius, raw, username=None, password=None):
    f_return = LightcurveStruct(survey=survey, source=source, pos=pos, data=None)

    def get_lightcurve():
        query_object = globals()[f"{survey.capitalize()}Query"](
            pos=pos,
            radius=radius,
            survey=survey,
            username=username,
            password=password,
            pmra=pmra,
            pmdec=pmdec,
            raw=raw,
        )

        if survey == "gaia":
            if source:
                data = query_object.get_data(source=source)
            else:
                data = query_object.get_data(pos=pos)
        elif survey in ["crts", "tess"]:
            data = query_object.get_data()
        else:
            query_object.set_url()
            response = query_object.get_response()
            if response is None:
                return get_f_return(survey)
            data = query_object.generate_data(response)

        if data is None:
            return get_f_return(survey)
        data = query_object.format_data(data)

        return data

    pmra, pmdec = None, None
    if source and survey != "atlas":
        from ..Tools import correctpm

        pos = correctpm(source=source, target_survey=survey)
        if not pos:
            return f_return

    elif source and survey == "atlas":
        from ..Tools import query

        gaia_data = query(kind="data", source=source, survey="gaia", level="internal")
        if gaia_data:
            import math

            ra2000, dec2000, pmra, pmdec = (
                gaia_data.data["ra2000"][0],
                gaia_data.data["dec2000"][0],
                gaia_data.data["pmra"][0],
                gaia_data.data["pmdec"][0],
            )
            pos = [ra2000, dec2000]
            if math.isnan(pmra) or math.isnan(pmdec):
                print("Note: could not correct coordinates due to missing pmra/pmdec")
                pmra, pmdec = None, None

    lightcurve = LightcurveStruct(
        survey=survey, source=source, pos=pos, data=get_lightcurve()
    )

    return lightcurve

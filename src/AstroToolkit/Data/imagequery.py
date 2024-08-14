class ImageStruct(object):
    """ImageStruct()
    This structure is returned from image queries, when read from a data file that was originally created by an image query, or through the Models module (in which case all attributes are set to None).

    Attributes
    ----------

    kind : str
        "image"

    survey: str
        survey to which the data pertains, see image query description for details.

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

    data: dict<list>
        Returned data in format:

        .. code-block:: python

            "image_data": image_data
            "image_header": image_header
            "size": size
            "image_time": image_time
            "wcs": wcs
            "image_focus": image_focus
            "overlay": overlay

        where:

        :**image_data**: list<float>, raw image data

        :**image_header**: astropy image header

        :**size**: int, image size in arcseconds

        :**image_time**: list<int>, [year,month] time at which the image was taken

        :**wcs**: astropy wcs object

        :**image_focus**: list<float>, [right ascension, declination] of image focus in degrees

        :**overlay**: list<dict>, detection information to overlay when plotting images with the plot() method

    |

    """

    def __init__(self, survey, source, pos, data, identifier=None):
        self.kind = "image"

        self.survey = survey
        self.source = source
        self.pos = pos
        self.identifier = identifier
        self.data = data
        self.figure = None
        self.dataname = None

    def __str__(self):
        return "<ATK Image Structure>"

    def plot(self):
        """
        Plots image data and assigns the resulting figure to the **figure** attribute.

        |

        """

        from ..Plotting.plotmap import map_to_plot

        return map_to_plot(self)

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
    """Base class for image queries"""

    def __init__(self, survey, size, band, pos):
        self.survey = survey

        self.size = size

        self.band = band
        self.pos = pos
        self.url = ""

    def get_image_data(self):
        from io import BytesIO

        import numpy as np
        import requests
        from astropy.io import fits
        from astropy.visualization import AsinhStretch, PercentileInterval

        # send request using the URL returned for the selected survey

        try:
            r = requests.get(self.url, timeout=15)

        except:
            print(f"Note: experiencing issues with {self.survey}.")

            return None, None

        if r.status_code != 200:
            print(f"Note: experiencing issues with {self.survey}.")

            return None, None

        try:
            fh = fits.open(BytesIO(r.content))

        except:
            print(f"Note: {self.survey} image query returned no data.")

            return None, None

        # read image fits file, apply a contrast filter
        (fh[0].data)[np.isnan(fh[0].data)] = 0.0

        transform = AsinhStretch() + PercentileInterval(95)

        self.image_header = fh[0].header

        self.image_data = transform(fh[0].data)

        return self.image_data, self.image_header

    @property
    def image_time(self):
        from astropy.time import Time

        mjd = self.image_header["MJD-OBS"]

        imageTime = Time(mjd, format="mjd").to_datetime()

        imageTime = [imageTime.year, imageTime.month]

        return imageTime

    @property
    def get_wcs(self):
        from astropy.wcs import WCS

        return WCS(self.image_header)


class PanstarrsQuery(GeneralQuery):
    def set_url(self):
        import numpy as np
        from astropy.table import Table

        url_size = self.size * 4

        url = f"https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={self.pos[0]}&dec={self.pos[1]}&band={self.band}"

        try:
            table = Table.read(url, format="ascii")

            if not len(table) > 0:
                print(f"Note: {self.survey} image query returned no data.")

                return None

            sub_url = f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra={self.pos[0]}&dec={self.pos[1]}&size={url_size}&format=fits"

            filter_list = ["yzirg".find(x) for x in table["filter"]]

            sub_table = table[np.argsort(filter_list)]

            url_base = f"{sub_url}&red="

            url_main = []

            for filename in sub_table["filename"]:
                url_main.append(url_base + filename)

            url_main = url_main[0]
            self.url = url_main

        except:
            print(f"Note: experiencing issues with {self.survey}.")

            return None


class SkymapperQuery(GeneralQuery):
    def set_url(self):
        import pandas as pd
        from astropy.table import Table

        url_size = self.size / 3600

        url = f"https://api.skymapper.nci.org.au/public/siap/dr2/query?POS={self.pos[0]},{self.pos[1]}&SIZE={url_size}&BAND={self.band}&FORMAT=image/fits&VERB=3&INTERSECT=covers&RESPONSEFORMAT=CSV"

        try:
            table = pd.read_csv(url)

            table = Table.from_pandas(table)

        except:
            print(f"Note: experiencing issues with {self.survey}.")

            return None

        if len(table) > 0:
            url_main = table["get_image"][0]
            self.url = url_main

        else:
            print(f"Note: {self.survey} image query returned no data.")

            return None


class DssQuery(GeneralQuery):
    def set_url(self):
        url_size = self.size / 60

        url_main = f"http://archive.stsci.edu/cgi-bin/dss_search?ra={self.pos[0]}&d={self.pos[1]}&v=3&e=J2000&f=fits&h={url_size}&w={url_size}"
        self.url = url_main

    @property
    def image_time(self):
        from astropy.time import Time

        time = self.image_header["DATE-OBS"]

        mins = int(time[14:16])

        hours = int(time[11:13])

        if mins >= 60:
            mins = mins - 60

            hours += 1

        mins = str(mins).zfill(2)

        hours = str(hours).zfill(2)

        time = time[0:10] + "T" + str(hours) + ":" + str(mins) + ":" + time[17:20]

        imageTime = Time(time, format="fits")

        mjd = imageTime.mjd

        imageTime = Time(mjd, format="mjd").to_datetime()

        imageTime = [imageTime.year, imageTime.month]

        return imageTime


def query(survey, size, band, pos=None, source=None, overlays=None):
    f_return = ImageStruct(survey=survey, source=source, pos=pos, data=None)

    def getimage():
        query_object = globals()[f"{survey.capitalize()}Query"](
            pos=pos, size=size, band=band, survey=survey
        )

        query_object.set_url()

        if query_object.url is None:
            return f_return

        image_data, image_header = query_object.get_image_data()

        if image_data is None or image_header is None:
            return f_return

        image_time = query_object.image_time

        wcs = query_object.get_wcs

        data_dict = {
            "image_data": image_data,
            "image_header": image_header,
            "size": size,
            "image_time": image_time,
            "wcs": wcs,
            "image_focus": pos,
        }

        data = ImageStruct(survey=survey, source=source, pos=pos, data=data_dict)

        return data

    if source:
        from ..Data.dataquery import SurveyInfo
        from ..Tools import correctpm
        from ..Tools import query as data_query

        survey_times = SurveyInfo().times

        gaia_data = data_query(
            kind="data", survey="gaia", source=source, level="internal"
        ).data
        if gaia_data:
            pos = [gaia_data["ra"][0], gaia_data["dec"][0]]

        else:
            print("Note: No gaia object found with given Source.")
            return f_return

        # get an initial image to get image_time

        image = getimage()

        if image.data:
            image_time = image.data["image_time"]

        else:
            return f_return

        # correct coords of source to image_time

        pos = correctpm(
            input_time=survey_times["gaia"], target_time=image_time, source=source
        )

    image = getimage()

    if overlays:
        from ..Data.imageoverlay import overlay_query

        overlay_data = overlay_query(image, overlays)
        image.data["overlay"] = overlay_data

    return image

from ..Data.dataquery import DataStruct


class GeneralQuery(object):
    def __init__(self, survey, source=None, pos=None, distance=None, radius=None):
        self.radius = None
        self.survey = survey
        self.source = source
        self.pos = pos
        self.data = None
        self.url = ""
        self.distance = distance
        self.radius = radius

        self.f_return = DataStruct(
            catalogue=None,
            survey=self.survey,
            source=self.source,
            pos=self.pos,
            data=None,
        )

    def send_request(self):
        import requests
        from requests.adapters import HTTPAdapter, Retry

        s = requests.Session()
        retries = Retry(total=5, backoff_factor=1)
        s.mount("http://", HTTPAdapter(max_retries=retries))
        r = requests.get(self.url, allow_redirects=True, timeout=30)

        if r.ok:
            return r
        else:
            return None


class StilismQuery(GeneralQuery):
    def set_url(self, distance, gal_lon, gal_lat):
        self.url = f"http://stilism.obspm.fr/reddening?frame=galactic&vlong={gal_lon}&ulong=deg&vlat={gal_lat}&ulat=deg&distance={distance}"

    def generate_data(self, response):
        from io import StringIO

        import pandas as pd

        file = StringIO(response.content.decode("utf-8"))
        df = pd.read_csv(file)
        return df

    def read_response(self):
        import pandas as pd

        from ..Tools import query as data_query

        gaia_data = data_query(
            kind="data", survey="gaia", source=self.source, level="internal"
        ).data
        if gaia_data:
            ra, dec, gal_lon, gal_lat, parallax, parallax_error = (
                gaia_data["ra"][0],
                gaia_data["dec"][0],
                gaia_data["l"][0],
                gaia_data["b"][0],
                gaia_data["parallax"][0],
                gaia_data["parallax_error"][0],
            )
        else:
            return None

        self.pos = [ra, dec]

        distance = 1 / (parallax / 1000)
        self.distance = distance

        upper_distance = 1 / ((parallax + parallax_error) / 1000)
        lower_distance = 1 / ((parallax - parallax_error) / 1000)

        data = {
            "distance": pd.DataFrame(),
            "upper_distance": pd.DataFrame(),
            "lower_distance": pd.DataFrame(),
        }
        for dist, label in zip(
            [distance, upper_distance, lower_distance], list(data.keys())
        ):
            self.set_url(distance=dist, gal_lon=gal_lon, gal_lat=gal_lat)
            response = self.send_request()
            if not response:
                return GeneralQuery(
                    survey=self.survey,
                    source=self.source,
                    pos=self.pos,
                    radius=self.radius,
                )
            gen_data = self.generate_data(response)
            if response is None:
                print("Note: Reddening query returned no data.")
                return GeneralQuery(
                    survey=self.survey,
                    source=self.source,
                    pos=self.pos,
                    radius=self.radius,
                ).f_return

            data[label] = gen_data

        return data

    def format_data(self, data):
        reddening_upper_error = (
            data["upper_distance"]["reddening[mag]"].tolist()[0]
            - data["distance"]["reddening[mag]"].tolist()[0]
        )
        reddening_lower_error = (
            data["lower_distance"]["reddening[mag]"].tolist()[0]
            - data["distance"]["reddening[mag]"].tolist()[0]
        )

        dist, dist_err, reddening, reddening_min_err, reddening_max_err = (
            data["distance"]["distance[pc]"].tolist()[0],
            data["distance"]["distance_uncertainty[pc]"].tolist()[0],
            data["distance"]["reddening[mag]"].tolist()[0],
            data["distance"]["reddening_uncertainty_min[mag]"].tolist()[0],
            data["distance"]["reddening_uncertainty_max[mag]"].tolist()[0],
        )

        import math

        red_err_lower = math.sqrt(reddening_min_err**2 + reddening_lower_error**2)
        red_err_upper = math.sqrt(reddening_max_err**2 + reddening_upper_error**2)

        return DataStruct(
            survey=self.survey,
            catalogue=None,
            pos=self.pos,
            source=self.source,
            data={
                "dist": [self.distance],
                "reddening_distance": [dist],
                "reddening_distance_err": [dist_err],
                "reddening": [reddening],
                "reddening_upper_limit": [round(reddening + red_err_upper, 3)],
                "reddening_lower_limit": [round(reddening - red_err_lower, 3)],
            },
            sub_kind="reddening",
        )


class GdreQuery(GeneralQuery):
    def set_url(self):
        if self.source:
            from ..Tools import correctpm

            self.pos = correctpm(source=self.source, target_time=[2000, 0])
            if not self.pos:
                return None

        self.url = f"https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={self.pos[0]}+{self.pos[1]}+equ+j2000&regSize={float(self.radius)}"

    def read_response(self):
        import xml.etree.ElementTree as ET

        self.set_url()
        response = self.send_request()
        if not response:
            return GeneralQuery(
                survey=self.survey, source=self.source, pos=self.pos, radius=self.radius
            ).f_return

        tree = ET.ElementTree(ET.fromstring(response.text))
        root = tree.getroot()

        status = root.attrib["status"]
        if status != "ok":
            return GeneralQuery(
                survey=self.survey, source=self.source, pos=self.pos, radius=self.radius
            ).f_return

        def get_value(element, param, kind):
            if kind == "ext":
                return float(
                    element.find("statistics")
                    .find(param)
                    .text.lstrip()
                    .rstrip()
                    .removesuffix(" (mag)")
                )
            elif kind == "mic":
                return float(
                    element.find("statistics")
                    .find(param)
                    .text.lstrip()
                    .rstrip()
                    .removesuffix(" (MJy/sr)")
                )
            elif kind == "dust":
                return float(
                    element.find("statistics")
                    .find(param)
                    .text.lstrip()
                    .rstrip()
                    .removesuffix(" (K)")
                )

        try:
            for element in root:
                if element.find("desc") is not None:
                    if (
                        element.find("desc").text.lstrip().rstrip()
                        == "E(B-V) Reddening"
                    ):
                        extinction_table_url = (
                            element.find("data").find("table").text.lstrip().rstrip()
                        )
                        extinction_image_url = (
                            element.find("data").find("image").text.lstrip().rstrip()
                        )
                        extinction_pixel_value_2011 = get_value(
                            element, "refPixelValueSandF", "ext"
                        )
                        extinction_pixel_value_1998 = get_value(
                            element, "refPixelValueSFD", "ext"
                        )
                        extinction_mean_value_2011 = get_value(
                            element, "meanValueSandF", "ext"
                        )
                        extinction_mean_value_1998 = get_value(
                            element, "meanValueSFD", "ext"
                        )
                        extinction_stdev_2011 = get_value(element, "stdSandF", "ext")
                        extinction_stdev_1998 = get_value(element, "stdSFD", "ext")
                        extinction_max_2011 = get_value(element, "maxValueSandF", "ext")
                        extinction_max_1998 = get_value(element, "maxValueSFD", "ext")
                        extinction_min_2011 = get_value(element, "minValueSandF", "ext")
                        extinction_min_1998 = get_value(element, "minValueSFD", "ext")
                    elif (
                        element.find("desc").text.lstrip().rstrip()
                        == "100 Micron Emission"
                    ):
                        micron_image_url = (
                            element.find("data").find("image").text.lstrip().rstrip()
                        )
                        micron_pixel_value = get_value(element, "refPixelValue", "mic")
                        micron_mean_value = get_value(element, "meanValue", "mic")
                        micron_stdev = get_value(element, "std", "mic")
                        micron_max = get_value(element, "maxValue", "mic")
                        micron_min = get_value(element, "minValue", "mic")
                    elif (
                        element.find("desc").text.lstrip().rstrip()
                        == "Dust Temperature"
                    ):
                        dust_image_url = (
                            element.find("data").find("image").text.lstrip().rstrip()
                        )
                        dust_pixel_value = get_value(element, "refPixelValue", "dust")
                        dust_mean_value = get_value(element, "meanValue", "dust")
                        dust_stdev = get_value(element, "std", "dust")
                        dust_max = get_value(element, "maxValue", "dust")
                        dust_min = get_value(element, "minValue", "dust")
        except:
            print("Note: reddening query returned data, but failed to read response.")
            return GeneralQuery(
                survey=self.survey, source=self.source, pos=self.pos, radius=self.radius
            ).f_return

        data_dict = {
            "extinction_2011": extinction_pixel_value_2011,
            "extinction_mean_2011": extinction_mean_value_2011,
            "extinction_stdev_2011": extinction_stdev_2011,
            "extinction_max_2011": extinction_max_2011,
            "extinction_min_2011": extinction_min_2011,
            "extinction_1998": extinction_pixel_value_1998,
            "extinction_mean_1998": extinction_mean_value_1998,
            "extinction_stdev_1998": extinction_stdev_1998,
            "extinction_max_1998": extinction_max_1998,
            "extinction_min_1998": extinction_min_1998,
            "extinction_table_url": extinction_table_url,
            "extinction_image_url": extinction_image_url,
            "100micron_emission": micron_pixel_value,
            "100micron_emission_mean": micron_mean_value,
            "100micron_emission_stdev": micron_stdev,
            "100micron_emission_max": micron_max,
            "100micron_emission_min": micron_min,
            "100micron_image": micron_image_url,
            "dust_temp": dust_pixel_value,
            "dust_temp_mean": dust_mean_value,
            "dust_temp_stdev": dust_stdev,
            "dust_temp_max": dust_max,
            "dust_temp_min": dust_min,
            "dust_image": dust_image_url,
        }

        return data_dict

    def format_data(self, data):
        return DataStruct(
            survey=self.survey,
            catalogue=None,
            pos=None,
            source=self.source,
            data=data,
            sub_kind="reddening",
        )


def query(survey, source=None, pos=None, radius=None):
    query_object = globals()[f"{survey.capitalize()}Query"](
        survey=survey, source=source, pos=pos, radius=radius
    )
    data = query_object.read_response()
    if data:
        data = query_object.format_data(data)
    if not data:
        return GeneralQuery(
            survey=survey, source=source, pos=pos, radius=radius
        ).f_return

    return data

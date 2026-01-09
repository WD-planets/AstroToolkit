import xml.etree.ElementTree as ET

import astropy.units as u
from requests import Response

from ...structures.definitions import Target
from ...utilities.defaults import RETURNS
from ...utilities.requests import send_request
from .image_core import get_image_data

# maps of ATK band names to actual IRSA band names
BAND_MAP = {
    "dss1": {"blue": "DSS1 Blue", "red": "DSS1 Red"},
    "dss2": {"blue": "DSS2 Blue", "red": "DSS2 Red", "ir": "DSS2 IR"},
    "wise": {"w1": "w1", "w2": "w2", "w3": "w3", "w4": "w4"},
    "2mass": {"j": "J", "h": "H", "k": "K"},
}


def parse_xml(survey: str, band: str, xml_string: str) -> Response | RETURNS:
    """
    Parses an IRSA XML tree to get a url to a fits image
    """

    root = ET.fromstring(xml_string)

    # Check status
    status = root.get("status")
    if status != "ok":
        return RETURNS.EXCEPTION

    # Navigate into <result>
    result = root.find("result")
    if result is None:
        return RETURNS.EXCEPTION

    # For example, get total images
    total_images = result.findtext("totalimages")
    if not int(total_images):
        return RETURNS.NULL

    # Loop through image entries
    urls = []
    for image_elem in result.findall("image"):
        if image_elem.findtext("band") == band:
            urls.append(image_elem.findtext("fitsurl"))

    return urls[0]


def check_inputs(survey: str, band: str, size: int):
    """
    Check if requested band and size are acceptable
    """

    if band not in BAND_MAP[survey]:
        raise ValueError(f"Invalid {survey} band. Supported bands are {list(BAND_MAP[survey].keys())}.")
    if not 6 * u.arcsec < size < 3600 * u.arcsec:
        raise ValueError(f"Size too large. Size of {survey} images must be between 6 and 3600 arcsec.")


def irsa_query(survey: str, target: Target, size: int, band: str, **kwargs: dict):
    """
    Perform an IRSA image query (DSS1/2, WISE, 2MASS)
    """

    check_inputs(survey, band, size)

    # get rid of dss generation
    query_band = BAND_MAP[survey][band]
    if survey in ["dss1", "dss2"]:
        query_survey = survey[:-1]
    else:
        query_survey = survey

    # get size in arcmin
    url_size = size.to(u.arcmin).value
    url = f"https://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?locstr={target.coords.ra.value}%20{target.coords.dec.value}&subsetsize={url_size}&survey={query_survey.capitalize()}&mode=prog&"

    # get response from url
    response = send_request(survey, url)
    if response is RETURNS.EXCEPTION:
        return response

    main_url = parse_xml(query_survey, query_band, response.content)
    if main_url in (RETURNS.NULL, RETURNS.EXCEPTION):
        return main_url

    image_hdu = get_image_data(main_url, survey, band, size, kwargs.get("epoch_fetcher"), kwargs.get("epoch_key"))

    return image_hdu

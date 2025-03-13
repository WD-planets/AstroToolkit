import re
import warnings

from astroquery.simbad import Simbad

warnings.simplefilter("ignore", category=Warning)


def identifier_query(identifier):
    url = f"https://simbad.cds.unistra.fr/simbad/sim-id?output.format=ASCII&Ident={identifier}"

    table = Simbad.query_object(identifier)
    if not table:
        return None
    df = table.to_pandas()
    response = df["MAIN_ID"].tolist()[0]
    response = re.sub(r"\s+", " ", response)

    return response


def pos_query(pos):
    ra, dec = pos[0], pos[1]

    from astropy.coordinates import SkyCoord

    simbad = Simbad()
    simbad.Row_LIMIT = 10
    coordinates = SkyCoord(ra, dec, unit=("deg", "deg"))
    response = simbad.query_region(coordinates, radius="3arcsec")
    if not response:
        return None
    response = response[0]["MAIN_ID"]
    response = re.sub(r"\s+", " ", response)

    return response

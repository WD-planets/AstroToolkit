from astropy.io.fits import Header
from astropy.time import Time

from ...structures.definitions import Target
from .irsa_queries import irsa_query


def parse_2mass_epoch(hdr: Header) -> Time:
    # Take UT_DATE and UT
    date_str = hdr.get("UT_DATE")
    time_str = hdr.get("UT")

    if date_str is None or time_str is None:
        raise ValueError("Failed to parse date/time keywords in header")

    year = int(date_str[:2])
    month = int(date_str[2:4])
    day = int(date_str[4:6])

    # 2MASS took place in ~1997–2001/2002 -> if year < 50, assume 2000+, otherwise 1900+
    year = 2000 + year if year < 50 else 1900 + year

    iso = f"{year:04d}-{month:02d}-{day:02d}T{time_str}"

    return Time(iso, scale="utc")


def query(target: Target, **kwargs: dict):
    band, size = kwargs["band"], kwargs["size"]

    return irsa_query("2mass", target, size, band, epoch_fetcher=parse_2mass_epoch, epoch_key=None)

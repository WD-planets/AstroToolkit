from astropy.coordinates import SkyCoord

from .image_core import fits_to_epoch
from .irsa_queries import irsa_query


def query(search_pos: SkyCoord, **kwargs: any):
    band, size = kwargs["band"], kwargs["size"]

    return irsa_query("dss2", search_pos, size, band, epoch_fetcher=fits_to_epoch, epoch_key="DATE-OBS")

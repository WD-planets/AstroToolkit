from astropy.coordinates import SkyCoord

from .image_core import fits_to_epoch
from .irsa_queries import irsa_query


def query(search_pos: SkyCoord, **kwargs: any):
    """
    Performs a WISE image query
    """

    band, size = kwargs["band"], kwargs["size"]

    return irsa_query("wise", search_pos, size, band, epoch_fetcher=fits_to_epoch, epoch_key="MIDOBS")

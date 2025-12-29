from ...structures.definitions import Target
from .image_core import fits_to_epoch
from .irsa_queries import irsa_query


def query(target: Target, **kwargs: any):
    """
    Perform a DSS1 image query
    """

    band, size = kwargs["band"], kwargs["size"]

    return irsa_query("dss1", target, size, band, epoch_fetcher=fits_to_epoch, epoch_key="DATE-OBS")

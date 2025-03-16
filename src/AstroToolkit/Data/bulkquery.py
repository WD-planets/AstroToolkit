from ..Configuration.catalogue_setup import CatalogueStruct

aliases = CatalogueStruct()
bulkdata_surveys = aliases.get_catalogue_list()


def bulkdata_query(radius, pos=None, source=None):
    from ..Tools import query

    bulk_data = {}
    for survey in bulkdata_surveys:
        data = query(kind="data", pos=pos, source=source, radius=radius, survey=survey, level="internal").data
        bulk_data[survey] = data

    if source:
        gaia_data = query(kind="data", source=source, survey="gaia", level="internal").data
        if gaia_data:
            pos = [gaia_data["ra"][0], gaia_data["dec"][0]]

    from ..Data.dataquery import DataStruct

    dataStruct = DataStruct(survey="all", catalogue=None, pos=pos, source=source, data=bulk_data, sub_kind="bulkdata")

    return dataStruct

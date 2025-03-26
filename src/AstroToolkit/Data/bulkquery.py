from ..Configuration.catalogue_setup import CatalogueStruct
from ..Configuration.epochs import EpochStruct

epochs = EpochStruct().epoch_list

aliases = CatalogueStruct()
bulkdata_surveys = aliases.get_catalogue_list()

break_str = "|"


def bulkdata_query(radius, pos=None, source=None):
    import math

    from ..Tools import query

    if source:
        gaia_data = query(kind="data", survey="gaia", source=source).data
        if math.isnan(gaia_data["pmra"][0]) or math.isnan(gaia_data["pmdec"][0]):
            success = False
        else:
            success = True

        trace = f"start -> extracted pos from source query, assumed {epochs['gaia']}{break_str}"
        if success:
            for survey in bulkdata_surveys:
                if survey != "gaia" and survey in epochs:
                    trace += f" -> {survey}: {epochs[survey]} -> {survey} query performed{break_str}"
                elif survey != "gaia" and survey not in epochs:
                    trace += f" -> {survey} lacking epoch definition, proper motion correction failed{break_str}"
            trace += " -> [2000,0] -> end"
        else:
            for survey in bulkdata_surveys:
                if survey != "gaia" and survey in epochs:
                    trace += f" -> proper motion correction failed -> {survey} query performed{break_str}"
                elif survey != "gaia" and survey not in epochs:
                    trace += f" -> {survey} lacking epoch definition, proper motion correction failed{break_str}"
            trace += " -> proper motion correction failed -> end"
    elif pos:
        trace = None

    bulk_data = {}
    for survey in bulkdata_surveys:
        data = query(kind="data", pos=pos, source=source, radius=radius, survey=survey, level="internal")
        data = data.data
        bulk_data[survey] = data

    if source:
        from ..Tools import correctpm

        gaia_data = query(kind="data", source=source, survey="gaia", level="internal").data
        ra, dec, pmra, pmdec = (gaia_data["ra"][0], gaia_data["dec"][0], gaia_data["pmra"][0], gaia_data["pmdec"][0])
        pos = [ra, dec]
        final_pos = correctpm(pos=pos, input_time=epochs["gaia"], target_time=[2000, 0], pmra=pmra, pmdec=pmdec)
    else:
        final_pos = pos

    from ..Data.dataquery import DataStruct

    dataStruct = DataStruct(survey="all", catalogue=None, pos=pos, source=source, data=bulk_data, sub_kind="bulkdata")
    dataStruct.trace = trace
    dataStruct.pos = final_pos

    return dataStruct

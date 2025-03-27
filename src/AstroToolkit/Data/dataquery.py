import warnings
from functools import wraps

import astropy.coordinates as coord
import astropy.units as u
import pandas as pd
from astroquery.vizier import Vizier

from ..Configuration.epochs import EpochStruct
from ..Misc.pmcorrection import correctpm
from ..Misc.rename_headers import renameHeadersDR3
from ..PackageInfo import SurveyInfo
from ..StructureMethods.method_definitions import savedata
from ..StructureMethods.method_definitions import showdata as showdata

warnings.simplefilter(action="ignore", category=UserWarning)

surveyInfo = SurveyInfo()
epochs = EpochStruct().epoch_list

# ensure that the row limit in returned data is infinite
row_limit = -1
Vizier.ROW_LIMIT = -1


class DataStruct(object):
    """DataStruct()
    This structure is returned from data, bulkdata and reddening queries, when read from a data file that was originally created by a data, bulkdata or reddening query, or through the Models module (in which case all attributes are set to None).

    .. rubric:: Attributes
        :heading-level: 1

    kind : *str*
        "data"

    subkind: *str*
        data query kind, from: "data", "bulkdata", "reddening"

    survey: *str*
        survey from which data originates. For bulkdata queries, defaults to None.

    catalogue: *str*
        Vizier catalogue ID from which data originates

    source : *int*
        Gaia source ID of target system (if given, else None)

    pos : *list<float>*
        Position of target system [right ascension, declination] in degrees

    identifier : *str*
        Position of target system in JHHMMSS.SS±DDMMSS.SS format

    dataname : *str*
        Default file for the :func:`savedata` method

    data: *dict*
        returned data in the form:

        .. code-block:: python

            {
            <Vizier column header>: <value>
            ...
            }

    |

    """

    def __init__(self, survey, catalogue, source, pos, data, identifier=None, sub_kind="data", trace=None):
        self.kind = "data"
        self.subkind = sub_kind
        self.survey = survey
        self.catalogue = catalogue
        self.source = source
        self.pos = pos
        self.identifier = identifier
        self.data = data
        self.dataname = None
        self.trace = trace

    def __str__(self):
        return "<ATK Data Structure>"

    @wraps(showdata)
    def showdata(self, pprint=True, print_methods=True):
        showdata(self, pprint, print_methods)
        return self

    @wraps(savedata)
    def savedata(self, fname=None):
        fname = savedata(self, fname)
        return fname


class VizierQuery(object):
    """Performs Vizier queries"""

    def __init__(self, catalogue, radius, survey=None, pos=None, source=None):
        self.catalogue = catalogue
        self.pos = pos
        self.radius = radius
        self.survey = survey
        self.source = source
        self.data = None

        self.f_return = DataStruct(
            survey=self.survey, catalogue=self.catalogue, source=self.source, pos=self.pos, data=None
        )

    # performs queries by coordinates
    def pos_query(self):
        data = []
        v = Vizier(columns=["**"], row_limit=row_limit)
        data.append(
            v.query_region(
                coord.SkyCoord(ra=self.pos[0], dec=self.pos[1], unit=(u.deg, u.deg), frame="icrs"),
                width=self.radius * u.arcsec,
                catalog=self.catalogue,
            )
        )
        try:
            data = data[0][0].to_pandas().reset_index(drop=True)
            return self.check_data(data)
        except:
            return self.f_return

    # performs queries by Gaia source
    def source_query(self):
        v = Vizier(columns=["**"], column_filters={"Source": "==" + str(self.source)}, row_limit=row_limit)
        data = v.query_object(f"GAIA DR3 {self.source}", catalog=self.catalogue)
        try:
            data = data[0].to_pandas().reset_index(drop=True)
            return self.check_data(data)
        except:
            return self.f_return

    # checks if any data was returned
    def check_data(self, data):
        try:
            if not data.empty:
                data = pd.DataFrame.to_dict(data, orient="list")
                return DataStruct(
                    survey=self.survey, catalogue=self.catalogue, source=self.source, pos=self.pos, data=data
                )
            else:
                return DataStruct(
                    survey=self.survey, catalogue=self.catalogue, source=self.source, pos=self.pos, data=None
                )
        except:
            return DataStruct(survey=self.survey, catalogue=self.catalogue, source=self.source, pos=self.pos, data=None)


# maps coordinates to vizier surveys, performing proper motion correction for source queries.
def query(survey, radius, pos=None, source=None):
    from ..Configuration.catalogue_setup import CatalogueStruct

    aliases = CatalogueStruct()
    all_surveys = aliases.get_catalogue_list()

    # get the necessary basic survey info
    supported_catalogues = surveyInfo.getCatalogueDict

    # if survey isn't a supported survey, take the 'survey' to be a vizier catalogue ID
    if survey not in all_surveys and survey != "gaia_lc":
        catalogue = survey
    elif survey == "gaia_lc":
        catalogue = supported_catalogues[survey]
    else:
        catalogue = all_surveys[survey]

    # perform coordinate Vizier query
    if pos:
        data = VizierQuery(survey=survey, catalogue=catalogue, radius=radius, pos=pos).pos_query()
        trace = None
        final_pos = pos

    # perform source Vizier query
    elif source:
        gaia_data = (
            VizierQuery(survey=survey, catalogue=supported_catalogues["gaia"], radius=radius, source=source)
            .source_query()
            .data
        )
        if not gaia_data:
            raise ValueError(f"Gaia source {source} not found.")
        ra, dec, pmra, pmdec = (
            gaia_data["RA_ICRS"][0],
            gaia_data["DE_ICRS"][0],
            gaia_data["pmRA"][0],
            gaia_data["pmDE"][0],
        )

        if catalogue == "I/355/gaiadr3":
            data = VizierQuery(survey=survey, catalogue=catalogue, radius=radius, source=source).source_query()
            trace = f"start -> extracted pos from source query, assumed {epochs['gaia']} -> [2000,0] -> end"
            final_pos = correctpm(
                input_time=epochs["gaia"], target_time=[2000, 0], ra=ra, dec=dec, pmra=pmra, pmdec=pmdec
            )
        elif catalogue == "I/355/epphot":
            data = VizierQuery(survey=survey, catalogue=catalogue, radius=radius, source=source).source_query()
            trace = f"start -> extracted pos from source query, assumed {epochs['gaia']} -> [2000,0] -> end"
            final_pos = correctpm(
                input_time=epochs["gaia"], target_time=[2000, 0], ra=ra, dec=dec, pmra=pmra, pmdec=pmdec
            )
        else:
            if survey in epochs and survey != "gaia":
                pos, success = correctpm(epochs["gaia"], epochs[survey], ra, dec, pmra, pmdec, check_success=True)
                if success:
                    trace = f"start -> extracted pos from source query, assumed {epochs['gaia']} -> {survey}: {epochs[survey]} -> {survey} query performed -> [2000,0] -> end"
                    final_pos = correctpm(
                        input_time=epochs[survey], target_time=[2000, 0], ra=ra, dec=dec, pmra=pmra, pmdec=pmdec
                    )
                else:
                    trace = f"start -> extracted pos from source query, assumed {epochs['gaia']} -> proper motion correction failed -> {survey} query performed -> end"
                    final_pos = pos
            else:
                print(f"Note: {survey} has no epoch definition. Proper motion has therefore not been corrected.")
                trace = None
                pos = [ra, dec]
                final_pos = pos

            data = VizierQuery(survey=survey, catalogue=catalogue, radius=radius, pos=pos, source=source).pos_query()
    else:
        raise Exception("No source or coordinates provided.")

    data.pos = final_pos
    data.trace = trace

    if catalogue == "I/355/gaiadr3" and data.data:
        data = renameHeadersDR3(data)

    if not data.data:
        print(f"Note: {survey} data query unsuccessful or returned no data.")

    return data

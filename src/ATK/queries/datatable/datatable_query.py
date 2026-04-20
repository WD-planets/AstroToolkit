import warnings

import numpy as np
from astropy.table import Table, vstack
from astropy.utils.metadata import MergeConflictWarning
from numpy.ma import is_masked

from ...configuration.survey_config import SURVEY_CONFIG
from ...structures.DataTable import DataTable
from ...structures.Target import Target
from ...Tools.query import query as general_query
from ...utilities.defaults import RETURNS

warnings.simplefilter("ignore", MergeConflictWarning)


def query(target: Target, **kwargs):
    radius = kwargs["radius"]
    cols = kwargs["columns"]

    aliases = SURVEY_CONFIG._get_aliases()

    data = []
    for survey, cols in cols.items():
        survey_data = general_query(kind="vizier", targets=target, survey=survey, radius=radius)
        if survey in aliases:
            catalogue = aliases[survey]
        else:
            catalogue = survey

        if survey_data.exception:
            return RETURNS.EXCEPTION

        if not survey_data.data:
            continue

        for col in cols:
            if col not in survey_data.data[0].table.colnames:
                raise ValueError(f"Couldn't find column {col} in {survey} catalogue.")

        base_tbl = survey_data.data[0].table[cols]

        correction = survey_data.data[0].correction

        if "_r" in survey_data.data[0].table.colnames:
            sep_base = np.asarray(survey_data.data[0].table["_r"])
        else:
            sep_base = None

        sep_str = f"separation ({radius.unit.to_string('fits')})"

        tables = []

        for colname in cols:
            col = base_tbl[colname]
            t = Table()

            t["survey"] = np.full(len(col), survey, dtype=object)
            t["catalogue"] = np.full(len(col), catalogue, dtype=object)
            t["correction"] = np.full(len(col), correction, dtype=object)
            t["parameter"] = np.full(len(col), colname)

            unit_str = f" {col.unit}" if col.unit is not None else ""
            value_col, unit_col = [], []
            for val in col:
                if is_masked(val) or (isinstance(val, float) and np.isnan(val)):
                    value_col.append("")
                    unit_col.append("")
                else:
                    value_col.append(f"{val:.3g}")
                    unit_col.append(unit_str)
            t["value"] = value_col
            t["unit"] = unit_col

            if sep_base is not None:
                t[sep_str] = sep_base
            else:
                t[sep_str] = np.nan

            tables.append(t)

        tbl_long = vstack(tables)

        tbl_long.sort([sep_str, "parameter"])

        param_order = {c: i for i, c in enumerate(cols)}
        tbl_long["param_order"] = [param_order[p] for p in tbl_long["parameter"]]
        tbl_long.sort([sep_str, "param_order"])
        tbl_long.remove_column("param_order")

        data.append(tbl_long)

    final_tbl = vstack(data)

    if not len(final_tbl):
        return RETURNS.NULL

    dt = DataTable(_target_key=target._key, table=final_tbl)

    return [dt]

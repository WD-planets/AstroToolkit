import numpy as np
import pandas as pd

from ...configuration.survey_config import SURVEY_CONFIG
from ...structures.DataTable import DataTable
from ...structures.Target import Target
from ...Tools.query import query as general_query
from ...utilities.defaults import RETURNS


def query(target: Target, **kwargs):
    radius = kwargs["radius"]
    rows = kwargs["rows"]

    aliases = SURVEY_CONFIG._get_aliases()

    data = []
    for survey, cols in rows.items():
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
            if col not in survey_data.data[0].data:
                raise ValueError(f"Couldn't find column {col} in {survey} catalogue.")

        df = survey_data.data[0].data[cols]
        df = df.round(3)
        df = df.astype(str)
        df = df.melt(var_name="parameter", value_name="value")
        df["survey"] = survey
        df["catalogue"] = catalogue
        df["correction"] = survey_data.data[0].correction

        sep_str = f"separation ({radius.unit.to_string('fits')})"
        if "_r" in survey_data.data[0].data:
            sep = np.tile(np.asarray(survey_data.data[0].data["_r"]), len(cols))
        else:
            sep = np.nan
        df[sep_str] = sep

        df["parameter"] = pd.Categorical(df["parameter"], categories=cols, ordered=True)
        df = df.sort_values([sep_str, "parameter"])

        front_cols = ["survey", "catalogue", "correction", sep_str]
        df = df[front_cols + [c for c in df.columns if c not in front_cols]]

        data.append(df)

    final_df = pd.concat(data)

    if final_df.empty:
        return RETURNS.NULL

    dt = DataTable(_target_key=target._key, data=final_df)

    return [dt]

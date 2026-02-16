import pandas as pd

from ...configuration.alias_config import ALIAS_CONFIG
from ...structures.DataTable import DataTable
from ...structures.Target import Target
from ...Tools.query import query as general_query
from ...utilities.defaults import RETURNS


def query(target: Target, **kwargs):
    radius = kwargs["radius"]
    rows = kwargs["rows"]

    aliases = ALIAS_CONFIG.as_dict()["vizier_aliases"]

    data = []
    for survey, cols in rows.items():
        if survey in aliases:
            survey_data = general_query(kind="vizier", target=target, survey=survey, radius=radius)

            catalogue = aliases[survey]
        else:
            survey_data = general_query(kind="vizier", target=target, catalogue=survey, radius=radius)

            catalogue = survey
            survey = None

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
        df["separation"] = survey_data.data[0].separation

        front_cols = ["survey", "catalogue", "correction", "separation"]
        df = df[front_cols + [c for c in df.columns if c not in front_cols]]

        data.append(df)

    final_df = pd.concat(data)

    if final_df.empty:
        return RETURNS.NULL

    dt = DataTable(_target_key=target._key, data=final_df)

    return [dt]

def gettable(selection, source, pos, radius):
    from bokeh.models import ColumnDataSource, DataTable, TableColumn

    from ..Data.dataquery import SurveyInfo
    from ..Misc.identifier_generation import identifier_from_pos
    from ..Tools import query

    metadata_defaults = SurveyInfo().metadata_defaults
    supported_surveys = SurveyInfo().list

    survey_col = []
    parameters_col = []
    values_col = []
    errors_col = []
    notes_col = []

    if source:
        gaia_data = query(
            kind="data", survey="gaia", pos=pos, source=source, level="internal"
        )
        if gaia_data.data:
            ra, dec = gaia_data.data["ra2000"][0], gaia_data.data["dec2000"][0]
            identifier = identifier_from_pos([ra, dec])
    elif pos:
        identifier = identifier_from_pos(pos)

    survey_col += ["ATK", "ATK", "ATK"]
    parameters_col += ["ATK pos", "ATK source", "ATK identifier"]
    values_col += [str(source), str(pos), str(identifier)]
    errors_col += [None, None, None]
    notes_col += ["ATK input source", "ATK input pos", None]

    for survey, columns in selection.items():
        if columns != "default":
            if not isinstance(selection[survey]["parameters"], list):
                selection[survey]["parameters"] = [selection[survey]["parameters"]]
            if not isinstance(selection[survey]["errors"], list):
                selection[survey]["errors"] = [selection[survey]["errors"]]
            if not isinstance(selection[survey]["notes"], list):
                selection[survey]["notes"] = [selection[survey]["notes"]]

        if survey in supported_surveys:
            data = query(
                kind="data",
                survey=survey,
                source=source,
                pos=pos,
                radius=radius,
                level="internal",
            ).data
            if data:
                if columns == "default":
                    selection[survey] = metadata_defaults[survey]
                for parameter, error, note in zip(
                    selection[survey]["parameters"],
                    selection[survey]["errors"],
                    selection[survey]["notes"],
                ):
                    survey_col.append(survey)
                    parameters_col.append(parameter)
                    values_col.append(str(data[parameter][0]))
                    errors_col.append(str(data[error][0]) if error else "---")
                    notes_col.append(note if note else "---")
        else:
            if not isinstance(selection[survey]["values"], list):
                selection[survey]["values"] = [selection[survey]["values"]]

            for parameter, value, error, note in zip(
                selection[survey]["parameters"],
                selection[survey]["values"],
                selection[survey]["errors"],
                selection[survey]["notes"],
            ):
                survey_col.append(survey)
                parameters_col.append(parameter)
                values_col.append(str(value))
                errors_col.append(str(error))
                notes_col.append(str(note))

    survey_col = [str(x) for x in survey_col]
    parameters_col = [str(x) for x in parameters_col]
    values_col = [str(x) for x in values_col]
    errors_col = [str(x) for x in errors_col]
    notes_col = [str(x) for x in notes_col]

    survey_col = ["---" if x == "None" else x for x in survey_col]
    parameters_col = ["---" if x == "None" else x for x in parameters_col]
    values_col = ["---" if x == "None" else x for x in values_col]
    errors_col = ["---" if x == "None" else x for x in errors_col]
    notes_col = ["---" if x == "None" else x for x in notes_col]

    data_structure = dict(
        survey=survey_col,
        parameter=parameters_col,
        value=values_col,
        error=errors_col,
        notes=notes_col,
    )

    data_source = ColumnDataSource(data_structure)
    table_columns = [
        TableColumn(field="survey", title="Survey"),
        TableColumn(field="parameter", title="Parameter"),
        TableColumn(field="value", title="Value"),
        TableColumn(field="error", title="Error"),
        TableColumn(field="notes", title="Notes"),
    ]

    from ..Configuration.baseconfig import ConfigStruct

    config = ConfigStruct()
    config.read_config()

    table = DataTable(
        source=data_source,
        columns=table_columns,
        width=int(3 * int(config.unit_size)),
        height=int(config.unit_size),
    )

    return table

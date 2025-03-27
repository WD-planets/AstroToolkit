from bokeh.models import InlineStyleSheet

from ..Configuration.baseconfig import ConfigStruct
from ..PackageInfo import metadataInfo

metadataInfo = metadataInfo()

config = ConfigStruct()
config.read_config()


class DataTable_Style:
    def __init__(self):
        temp_config = ConfigStruct()
        temp_config.read_config()

        font = temp_config.font
        font_size = temp_config.font_size
        if not font_size.endswith("pt"):
            font_size += "pt"

        style_sheet = InlineStyleSheet(
            css=f".slick-header-columns {{background-color: #e0e0e0 !important;font-family: {font.lower()};font-size: {int(font_size[:-2])}pt; font-weight: normal}}.slick-row {{font-size: {int(font_size[:-2]) - 1}pt; font-weight: normal}}"
        )

        self.style_sheet = style_sheet

        datapage_font_size = temp_config.datapage_font_size
        if not datapage_font_size.endswith("pt"):
            datapage_font_size += "pt"

        self.datapage_style_sheet = InlineStyleSheet(
            css=f".slick-header-columns {{background-color: #e0e0e0 !important;font-family: {font.lower()};font-size: {int(datapage_font_size[:-2])}pt; font-weight: normal}}.slick-row {{font-size: {int(datapage_font_size[:-2]) - 1}pt; font-weight: normal}}"
        )


def gettable(selection, source, pos, radius):
    from bokeh.models import ColumnDataSource, DataTable, TableColumn

    from ..Misc.identifier_generation import identifier_from_pos
    from ..Tools import query

    metadata_defaults = metadataInfo.metadataDefaults

    survey_col = []
    parameters_col = []
    values_col = []
    errors_col = []
    notes_col = []

    if source:
        gaia_data = query(kind="data", survey="gaia", pos=pos, source=source, level="internal")
        if gaia_data.data:
            ra, dec = gaia_data.data["ra2000"][0], gaia_data.data["dec2000"][0]
            identifier = identifier_from_pos([ra, dec])
    elif pos:
        identifier = identifier_from_pos(pos)

    survey_col += ["ATK", "ATK", "ATK"]
    parameters_col += ["ATK source", "ATK pos", "ATK identifier"]
    values_col += [str(source), str(pos), str(identifier)]
    errors_col += [None, None, None]
    notes_col += ["ATK input source", "ATK input pos", None]

    for entry in selection:
        kind = entry["kind"]
        if kind == "vizier":
            surveys = [entry["survey"]]
        elif kind == "atk_defaults":
            surveys = entry["surveys"]
        if kind in ["vizier", "atk_defaults"]:
            for survey in surveys:
                data = query(kind="data", survey=survey, source=source, pos=pos, radius=radius, level="internal").data
                if data:
                    if kind == "atk_defaults":
                        entry["parameters"] = metadata_defaults[survey]["parameters"]
                        entry["errors"] = metadata_defaults[survey]["errors"]
                        entry["notes"] = metadata_defaults[survey]["notes"]
                    for parameter, error, note in zip(entry["parameters"], entry["errors"], entry["notes"]):
                        survey_col.append(survey)
                        parameters_col.append(parameter)
                        values_col.append(str(data[parameter][0]))
                        errors_col.append(str(data[error][0]) if error else "---")
                        notes_col.append(note if note else "---")
        elif kind == "external":
            survey = entry["survey"]
            for parameter, value, error, note in zip(
                entry["parameters"], entry["values"], entry["errors"], entry["notes"]
            ):
                survey_col.append(survey)
                parameters_col.append(parameter)
                values_col.append(value)
                errors_col.append(error)
                notes_col.append(note)

    survey_col = [str(x) if x else None for x in survey_col]
    parameters_col = [str(x) if x else None for x in parameters_col]
    values_col = [str(x) if x else None for x in values_col]
    errors_col = [str(x) if x else None for x in errors_col]
    notes_col = [str(x) if x else None for x in notes_col]

    survey_col = ["---" if not x else x for x in survey_col]
    parameters_col = ["---" if not x else x for x in parameters_col]
    values_col = ["---" if not x else x for x in values_col]
    errors_col = ["---" if not x else x for x in errors_col]
    notes_col = ["---" if not x else x for x in notes_col]

    data_structure = dict(
        survey=survey_col, parameter=parameters_col, value=values_col, error=errors_col, notes=notes_col
    )

    data_source = ColumnDataSource(data_structure)
    table_columns = [
        TableColumn(field="survey", title="Survey"),
        TableColumn(field="parameter", title="Parameter"),
        TableColumn(field="value", title="Value"),
        TableColumn(field="error", title="Error"),
        TableColumn(field="notes", title="Notes"),
    ]

    table = DataTable(
        source=data_source, columns=table_columns, width=int(3 * int(config.unit_size)), height=int(config.unit_size)
    )

    table.stylesheets = [DataTable_Style().style_sheet]

    return table

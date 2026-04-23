from ..structures.Base import Container


def get_axis_label(ctnr: Container, attr: str):
    """
    Returns a Bokeh-suitable axis label for an astropy quantity
    """

    unit = ctnr._get_attr_unit(attr)
    if unit is None:
        return rf"{attr.capitalize()}"

    label = f"{attr.capitalize()} / {unit.to_string('unicode')}"

    return label

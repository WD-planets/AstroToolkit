from ..structures.definitions import BaseContainer


def get_axis_label(ctnr: BaseContainer, attr: str):
    """
    Returns a Bokeh-suitable axis label for a quantity with units using MathJax and latex math
    """

    unit = ctnr._get_attr_unit(attr)
    if unit is None:
        return rf"{attr.capitalize()}"

    label = f"{attr.capitalize()} / {unit.to_string('unicode')}"

    return label

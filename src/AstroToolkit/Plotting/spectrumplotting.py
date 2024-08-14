from bokeh import events
from bokeh.models import CustomJS, Label, Range1d
from bokeh.plotting import figure


def get_overlay():
    wavelengths = []
    wavelengths.append(
        {
            "label": "Hydrogen",
            "wavelengths": [8503, 8454, 8598, 8665, 6562, 4862, 4340, 4101.734],
            "line_labels": [
                None,
                None,
                None,
                None,
                r"\[H\alpha\]",
                r"\[H\beta\]",
                r"\[H\gamma\]",
                r"\[H\delta\]",
            ],
        }
    )
    wavelengths.append(
        {
            "label": "Helium",
            "wavelengths": [4472, 4686, 4713, 4921, 5016, 5876, 6678],
            "line_labels": [None, None, None, None, None, None, None],
        }
    )
    wavelengths.append(
        {
            "label": "Sodium",
            "wavelengths": [8183, 8195],
            "line_labels": [None, None],
        }
    )
    wavelengths.append(
        {
            "label": "Calcium",
            "wavelengths": [3934, 3967, 8498, 8542, 8662],
            "line_labels": [None, None, None, None, None],
        }
    )
    wavelengths.append(
        {
            "label": "Calcium II",
            "wavelengths": [3608, 3854, 4109, 4383, 4737, 5165, 5636, 6191],
            "line_labels": [None, None, None, None, None, None, None, None],
        }
    )
    colours = ["red", "blue", "purple", "orange", "green", "lime", "brown"]

    return wavelengths, colours


def plot_spectrum(struct):
    x, y = struct.data["wavelength"], struct.data["flux"]

    y = [i * pow(10, 16) for i in y]

    plot = figure(
        width=400,
        height=400,
        title=f"{struct.survey} Spectrum",
        x_axis_label=r"\[\lambda\text{ }[\text{AA}]\]",
        y_axis_label=r"\[\text{flux [}10^{-16}\text{ erg}\text{ cm }^{-2}\text{ s }^{-1}\text{AA}^{-1}]\]",
    )

    plot.line(x, y, color="black", line_width=1)

    wavelengths, colours = get_overlay()

    bottom_margin = 0.05
    top_margin = 0.4
    y_min = min(y) - bottom_margin * (max(y) - min(y))
    y_max = max(y) + top_margin * (max(y) - min(y))
    plot.y_range = Range1d(y_min, y_max)
    plot.y_range.max_interval = y_max
    plot.y_range.min_interval = y_min

    for index, element in enumerate(wavelengths):
        annotation_count = 0
        annotation_group_size = len(
            [x for x in element["line_labels"] if x is not None]
        )
        for wavelength, annotation in zip(
            element["wavelengths"], element["line_labels"]
        ):
            line_renderer = plot.line(
                x=[wavelength, wavelength],
                y=[y_min, y_max],
                color=colours[index],
                legend_label=element["label"],
                level="underlay",
            )
            if annotation:
                annotation_height = max(y) + (
                    annotation_count / annotation_group_size
                ) * 0.3 * max(y)
                label = Label(
                    x=wavelength,
                    y=annotation_height,
                    x_offset=2,
                    text=annotation,
                    text_font_size="12pt",
                )

                plot.add_layout(label)
                annotation_count += 1

                line_renderer.js_on_change(
                    "visible",
                    CustomJS(args=dict(ls=label), code="ls.visible = cb_obj.visible;"),
                )

    toggle_legend_js = CustomJS(
        args=dict(leg=plot.legend[0]),
        code="""
            if (leg.visible) {
                leg.visible = false
                }
            else {
                leg.visible = true
            }
    """,
    )

    plot.js_on_event(events.DoubleTap, toggle_legend_js)
    plot.legend.click_policy = "hide"

    return plot

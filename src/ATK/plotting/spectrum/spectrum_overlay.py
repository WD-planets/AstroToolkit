import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource, CustomJS, Label, Span
from bokeh.plotting import figure

from ...configuration.base_config import BASE_CONFIG
from ...structures.methods.spectrum.radial_velocities import get_velocities
from ...structures.Spectrum import Spectrum
from ..colours import get_palette

# wavelengths in Angstroms
OVERLAY_LINES = [
    {"label": "Hydrogen", "wavelength": 8503},
    {"label": "Hydrogen", "wavelength": 8454},
    {"label": "Hydrogen", "wavelength": 8598},
    {"label": "Hydrogen", "wavelength": 8665},
    {"label": "Hydrogen", "wavelength": 6562, "line_label": r"\[\text{H}\alpha\]"},
    {"label": "Hydrogen", "wavelength": 4862, "line_label": r"\[\text{H}\beta\]"},
    {"label": "Hydrogen", "wavelength": 4340, "line_label": r"\[\text{H}\gamma\]"},
    {"label": "Hydrogen", "wavelength": 4101.734, "line_label": r"\[\text{H}\delta\]"},
    *[{"label": "Helium", "wavelength": w} for w in [4472, 4686, 4713, 4921, 5016, 5876, 6678]],
    *[{"label": "Sodium", "wavelength": w} for w in [8183, 8195]],
    *[{"label": "Calcium", "wavelength": w} for w in [3934, 3967, 8498, 8542, 8662]],
    *[{"label": "Calcium II", "wavelength": w} for w in [3608, 3854, 4109, 4383, 4737, 5165, 5636, 6191]],
]


def plot_overlay(plot: figure, spectrum: Spectrum):
    flux = spectrum._get_attr_value("flux")
    x = spectrum._get_attr_value(spectrum.x_type)

    overlay_lines = OVERLAY_LINES

    if spectrum.x_type == "velocity":
        for line in overlay_lines:
            line["velocity"] = get_velocities(line["wavelength"], spectrum.wav_ref.to(u.angstrom).value).value

    elements = list(set([line["label"] for line in overlay_lines]))
    colours = get_palette(len(elements))

    label_counters = {}
    element_renderers = {element: {"spans": [], "labels": []} for element in elements}

    text_size = str(BASE_CONFIG._get("plot_settings", "font_size"))
    if not text_size.endswith("pt"):
        text_size += "pt"
    text_font = str(BASE_CONFIG._get("plot_settings", "font"))

    y_min = float(np.min(flux))
    y_max = float(np.max(flux))

    for line in overlay_lines:
        xpos = line[spectrum.x_type]

        if not (np.min(x) < xpos < np.max(x)):
            continue

        label_name = line["label"]
        colour = colours[elements.index(label_name)]

        span = Span(location=xpos, dimension="height", line_color=colour, line_width=1)
        span.visible = False
        plot.add_layout(span)
        element_renderers[label_name]["spans"].append(span)

        if "line_label" in line:
            n_labels = label_counters.get(label_name, 0)
            total_labels = sum(1 for ln in overlay_lines if ln.get("label") == label_name and "line_label" in ln)

            y_pos = y_max + (n_labels / max(1, total_labels)) * 0.3 * y_max
            lbl = Label(x=xpos, y=y_pos, x_offset=2, text=line["line_label"], text_font_size=text_size, text_font=text_font)
            lbl.visible = False
            plot.add_layout(lbl)
            element_renderers[label_name]["labels"].append(lbl)
            label_counters[label_name] = n_labels + 1

    for idx, element in enumerate(elements):
        colour = colours[idx]

        y_min = float(np.min(flux))
        y_max = float(np.max(flux))

        if y_max == y_min:
            y_inside = y_min
        else:
            y_inside = y_min + 0.01 * (y_max - y_min)

        x_min = float(np.min(x))
        x_span = (float(np.max(x)) - x_min) * 1e-6  # extremely small

        source = ColumnDataSource(dict(x=[x_min, x_min + x_span], y=[y_inside, y_inside]))

        dummy_renderer = plot.line("x", "y", source=source, line_color=colour, line_width=1, legend_label=element)
        dummy_renderer.visible = False

        # CustomJS to toggle spans + labels
        callback = CustomJS(
            args=dict(spans=element_renderers[element]["spans"], labels=element_renderers[element]["labels"]),
            code="""
                const visible = cb_obj.visible;
                for (let i = 0; i < spans.length; i++) {
                    spans[i].visible = visible;
                }
                for (let j = 0; j < labels.length; j++) {
                    labels[j].visible = visible;
                }
            """,
        )
        dummy_renderer.js_on_change("visible", callback)

    plot.legend.click_policy = "hide"
    return plot

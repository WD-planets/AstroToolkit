import astropy.units as u
import numpy as np
from bokeh.models import CustomJS, Label, Range1d
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

    elements = list(set([line["label"] for line in OVERLAY_LINES]))
    colours = get_palette(len(elements))

    y_min, y_max = 0, max(flux) * 1.4
    plot.y_range = Range1d(y_min, y_max)
    plot.y_range.min_interval = y_min
    plot.y_range.max_interval = y_max

    # Track how many labels per element for stacked annotations
    label_counters = {}

    text_size = str(BASE_CONFIG.get("plot_settings", "font_size"))
    if not text_size.endswith("pt"):
        text_size += "pt"
    text_font = str(BASE_CONFIG.get("plot_settings", "font"))

    for idx, line in enumerate(OVERLAY_LINES):
        if not (line[spectrum.x_type] > np.min(x)) and (line[spectrum.x_type] < np.max(x)):
            continue

        label_name = line["label"]

        colour = colours[elements.index(label_name)]

        # Create vertical line
        line_renderer = plot.line(
            x=[line[spectrum.x_type], line[spectrum.x_type]],
            y=[1.5 * y_min, 1.5 * y_max],
            color=colour,
            legend_label=label_name,
            level="underlay",
        )

        # Handle annotation if it exists
        if "line_label" in line:
            # count how many annotations for this element so far
            n_labels = label_counters.get(label_name, 0)
            total_labels = sum(1 for line in OVERLAY_LINES if line.get("label") == label_name and "line_label" in line)
            y_pos = max(flux) + (n_labels / max(1, total_labels)) * 0.3 * max(flux)

            label = Label(
                x=line[spectrum.x_type], y=y_pos, x_offset=2, text=line["line_label"], text_font_size=text_size, text_font=text_font
            )
            plot.add_layout(label)

            # Sync label visibility with line visibility
            line_renderer.js_on_change("visible", CustomJS(args=dict(lbl=label), code="lbl.visible = cb_obj.visible;"))

            label_counters[label_name] = n_labels + 1

    return plot

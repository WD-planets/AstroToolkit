import astropy.units as u
import numpy as np
from astropy.timeseries import LombScargle
from bokeh import events
from bokeh.models import CustomJS
from bokeh.plotting import figure

newline = "\n"


def lomb_scargle(
    data, freq=None, bins=None, foverlay=True, repeat=1, shift=0, survey=None
):
    class timeseries_data(object):
        def __init__(
            self, time, mag, mag_err, power, frequency, freq, foverlay, repeat, shift
        ):
            self.time = time
            self.mag = mag
            self.mag_err = mag_err
            self.power = power
            self.frequency = frequency
            self.phase_freq = freq
            self.phase_bins = bins
            self.foverlay = foverlay
            self.repeat = repeat
            self.shift = shift

        @property
        def powspec_plot(self):
            plot = figure(
                width=400,
                height=400,
                x_axis_label=r"Frequency / \[\text{days}^{-1}\]",
                y_axis_label="Lomb-Scargle Power",
                title=f"{survey} L-S Power Spectrum",
            )

            max_freq = self.frequency[np.nanargmax(self.power)]
            period = 24 / max_freq

            plot.line(
                x=self.frequency,
                y=self.power,
                legend_label=f"Max Frequency: {round(max_freq.value,2)} 1/d {newline}Period: {round(period.value,2)} h",
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

            return plot, max_freq

        @property
        def phasefold_plot(self):
            # code for this function adapted from Keith
            def do_binning(bins, phase, mag, mag_err):
                bin_phase = []
                bin_y = []
                bin_y_err = []

                mag = np.asarray(mag)
                mag_err = np.asarray(mag_err)

                for bphase in np.linspace(0, 1 - 1 / bins, num=bins):
                    mask = (phase < bphase + 1 / bins) & (phase >= bphase)
                    if len(mag[mask]) > 0:
                        bin_phase += [bphase + 0.5 / bins]
                        weights = 1 / (mag_err[mask] ** 2)
                        baverage, norm = np.average(
                            mag[mask], weights=weights, returned=True
                        )
                        bin_y += [baverage]
                        bin_y_err += [1 / np.sqrt(norm)]

                return bin_phase, bin_y, bin_y_err

            if not self.phase_freq:
                best_frequency = self.frequency[np.nanargmax(self.power)]
                freq_multiplier = 1.0
            else:
                best_frequency = self.phase_freq / u.day
                freq_multiplier = (
                    self.frequency[np.nanargmax(self.power)].value / best_frequency
                ) / u.day

            t_fit = np.linspace(0, 1 / best_frequency.value, 1000) * u.day
            ls = LombScargle(self.time, self.mag, self.mag_err)

            # this uses the frequency determined by L-S
            y_fit = ls.model(t=t_fit, frequency=freq_multiplier * best_frequency)

            phase = (self.time % (1 / best_frequency)) * best_frequency

            t_fit = t_fit * best_frequency

            cut_indices = [i for i, val in enumerate(t_fit) if val > 1]

            t_fit_formatted = [
                val for i, val in enumerate(t_fit) if i not in cut_indices
            ]
            y_fit_formatted = [
                val for i, val in enumerate(y_fit) if i not in cut_indices
            ]

            median_mag = np.median(self.mag).value

            phase = [x.value for x in phase]
            self.mag = [x.value - median_mag for x in self.mag]
            t_fit_formatted = [x.value for x in t_fit_formatted]
            y_fit_formatted = [x.value - median_mag for x in y_fit_formatted]
            self.mag_err = [x.value for x in self.mag_err]

            mag = self.mag
            mag_err = self.mag_err

            if self.phase_bins:
                phase, mag, mag_err = do_binning(self.phase_bins, phase, mag, mag_err)

            if self.repeat > 1:
                base_mag = mag.copy()
                base_mag_err = mag_err.copy()
                base_phase = phase.copy()
                t_fit_base = t_fit_formatted.copy()
                y_fit_base = y_fit_formatted.copy()

                for i in range(1, self.repeat):
                    mag += base_mag
                    mag_err += base_mag_err
                    phase += [x + i for x in base_phase]
                    t_fit_formatted += [x + i for x in t_fit_base]
                    y_fit_formatted += y_fit_base

            if self.shift != 0:
                self.shift = self.shift - int(self.shift)
                phase = [
                    (
                        (x + self.shift)
                        if (x + self.shift <= self.repeat and x + self.shift >= 0)
                        else (
                            (x + self.shift + self.repeat)
                            if (x + self.shift < 0)
                            else (x + self.shift - self.repeat)
                        )
                    )
                    for x in phase
                ]
                t_fit_formatted = [
                    (
                        (x + self.shift)
                        if (x + self.shift <= self.repeat and x + self.shift >= 0)
                        else (
                            (x + self.shift + self.repeat)
                            if (x + self.shift < 0)
                            else (x + self.shift - self.repeat)
                        )
                    )
                    for x in t_fit_formatted
                ]

                t_fit_formatted, y_fit_formatted = zip(
                    *sorted(zip(t_fit_formatted, y_fit_formatted))
                )

            err_xs = [[x, x] for x in phase]
            err_ys = [[y - y_err, y + y_err] for y, y_err in zip(mag, mag_err)]

            plot = figure(
                width=400,
                height=400,
                x_axis_label="Phase",
                y_axis_label="Magnitude (relative to median)",
            )
            plot.scatter(x=phase, y=mag, level="guide")

            if self.foverlay:
                plot.line(
                    x=t_fit_formatted,
                    y=y_fit_formatted,
                    line_width=3,
                    line_color="black",
                    level="overlay",
                    legend_label=f"Period: {round(24/(best_frequency.value),2)} h",
                    alpha=0.75,
                )

            plot.multi_line(
                xs=err_xs,
                ys=err_ys,
                line_width=0.5,
                level="glyph",
                line_cap="square",
            )

            plot.y_range.flipped = True

            if self.foverlay:
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

            return plot

    from .timeseries_formatting import format_data

    try:
        time, mag, mag_err = format_data(data)
    except:
        print(
            "Note: No data passed to timeseries tool, suggests no light curve data was found."
        )
        return None

    time, mag, mag_err = time * u.day, mag * u.mag, mag_err * u.mag

    freqs = np.linspace(0, 60, 1500000) / u.day
    power = LombScargle(time, mag, mag_err, fit_mean=False).power(freqs)

    data_class = timeseries_data(
        time, mag, mag_err, power, freqs, freq, foverlay, repeat, shift
    )

    return data_class

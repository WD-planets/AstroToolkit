from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from bokeh.models import Arrow, ColumnDataSource, CustomJS, HoverTool, OpenHead
from bokeh.plotting import figure
from scipy.signal import find_peaks

from .fitting import do_fitting

if TYPE_CHECKING:
    from ....structures.Spectrum import Spectrum

C_KMS = 299792.458


def get_velocities(wav: Quantity, wav_ref: Quantity) -> np.ndarray:
    """
    Calculates radial velocity from a reference wavelength
    """

    return (wav - wav_ref) / wav_ref * C_KMS * u.Unit("km s-1")


def refine_rv(rvs: np.ndarray, scores: np.ndarray, idx: int, window: int = 2) -> float:
    """
    Refine RV peak using a quadratic fit around a given index
    """

    # window mask
    idx_start = max(idx - window, 0)
    idx_end = min(idx + window + 1, len(rvs))

    # mask arrays
    x = rvs[idx_start:idx_end]
    y = scores[idx_start:idx_end]

    mask = np.isfinite(y)
    if np.sum(mask) < 3:
        return rvs[idx]

    a, b, c = np.polyfit(x[mask], y[mask], 2)

    # vertex of quadratic fit
    rv_refined = -b / (2 * a)

    return rv_refined


def match_lines(wav_obs_val: np.ndarray, wav_lab_val: np.ndarray, rv: float, dv_tol_val: float, used_obs: np.ndarray):
    # shift wavelengths by rv increment
    shifted = wav_lab_val * (1.0 + rv / C_KMS)
    used = used_obs.copy()
    residuals = []

    # for each wavelength in shifted lab wavelengths
    for w in shifted:
        # get separation
        diffs = np.abs(wav_obs_val - w)
        # don't re-use the same features
        diffs[used] = np.inf
        # get best-matched feature
        j = np.argmin(diffs)

        # get separation in velocity
        dv = diffs[j] / w * C_KMS
        if dv < dv_tol_val:
            used[j] = True
            residuals.append(dv)

    return used, residuals


def find_rv_components(
    wav_obs: u.Quantity,
    wav_lab: u.Quantity,
    rv_min: float = -500.0,
    rv_max: float = 500.0,
    rv_step: float = 0.01,
    dv_tol: u.Quantity = 30 * u.km / u.s,
    min_matches: int = 2,
    min_rv_sep: u.Quantity = 30 * u.km / u.s,
):
    # convert Quantities to scalar arrays
    wav_obs_val = wav_obs.to_value(u.AA)
    wav_lab_val = wav_lab.to_value(u.AA)
    dv_tol_val = dv_tol.to_value(u.km / u.s)
    min_rv_sep_val = min_rv_sep.to_value(u.km / u.s)

    # set up radial velocity array at a given resolution
    rvs = np.arange(rv_min, rv_max + rv_step, rv_step)
    scores = np.full(len(rvs), -np.inf, dtype=float)

    # score radial velocity features
    for i, rv in enumerate(rvs):
        used, residuals = match_lines(wav_obs_val, wav_lab_val, rv, dv_tol_val, used_obs=np.zeros(len(wav_obs_val), dtype=bool))

        # append average separation to scores array
        if len(residuals) >= min_matches:
            scores[i] = -np.median(np.abs(residuals))

    # find peaks in scores array - minimum distance in pixels, prominence given by std dev
    peak_indices, _ = find_peaks(scores, distance=int(min_rv_sep_val / rv_step), prominence=np.std(scores[np.isfinite(scores)]) * 0.5)

    # sort peaks by score
    peak_indices = sorted(peak_indices, key=lambda i: scores[i], reverse=True)

    rv_components, matched_lines = [], []
    globally_used_obs = np.zeros(len(wav_obs_val), dtype=bool)

    # iterate through peaks
    for idx in peak_indices:
        # improve rv estimate
        rv = refine_rv(rvs, scores, idx)

        # ignore components that are too close to those that have already been found
        if any(abs(rv - rvi) < min_rv_sep_val for rvi in rv_components):
            continue

        used, _ = match_lines(wav_obs_val, wav_lab_val, rv, dv_tol_val, used_obs=globally_used_obs)

        newly_matched = used & ~globally_used_obs
        if np.sum(newly_matched) < min_matches:
            continue

        rv_components.append(rv)
        matched_lines.append(wav_obs[newly_matched])
        globally_used_obs |= newly_matched

    rv_components = [rv * u.km / u.s for rv in rv_components]
    rvs = rvs * u.km / u.s

    return rv_components, matched_lines, rvs, scores


def plot_rv(
    plot: figure,
    spectrum: Spectrum,
    wav_lab: u.Quantity,
    rv_components: list[u.Quantity],
    matched_lines: list[u.Quantity],
    features: u.Quantity,
    peak_vals: np.ndarray,
):
    colors = itertools.cycle(["red", "blue", "green", "orange", "purple"])
    lab_vals = wav_lab.to_value(u.AA)
    features_val = features.to_value(u.AA)

    flux_min = np.min(spectrum.flux.value) * 0.95
    flux_max = np.max(spectrum.flux.value) * 1.05

    for rv, matched_obs in zip(rv_components, matched_lines):
        color = next(colors)
        rv_val = rv.to_value(u.km / u.s)

        xs, ys = [], []
        lab_rest_list, obs_shift_list, rv_list, y_end_list = [], [], [], []

        for obs in matched_obs.to_value(u.AA):
            # Find the lab line that was matched to this observed feature
            lab_idx = np.argmin(np.abs(lab_vals - obs / (1 + rv_val / C_KMS)))
            lab_rest = lab_vals[lab_idx]
            lab_shift = lab_rest * (1 + rv_val / C_KMS)

            # Find y-value from closest feature
            feat_idx = np.argmin(np.abs(features_val - obs))
            y_peak = peak_vals[feat_idx]

            xs.append([lab_rest, lab_shift])
            ys.append([y_peak, y_peak])
            lab_rest_list.append(lab_rest)
            obs_shift_list.append(lab_shift)
            rv_list.append(rv_val)
            y_end_list.append(y_peak)

        source = ColumnDataSource(dict(xs=xs, ys=ys, lab_rest=lab_rest_list, obs_shift=obs_shift_list, rv=rv_list, y_end=y_end_list))

        # Invisible hover lines
        hover_line = plot.multi_line(xs="xs", ys="ys", source=source, line_width=20, alpha=0, legend_label=f"RV = {rv_val:.1f} km/s")
        plot.add_tools(HoverTool(renderers=[hover_line], tooltips=[("rest λ", "@lab_rest Å"), ("shifted λ", "@obs_shift Å"), ("RV", "@rv km/s")]))

        head = OpenHead(line_color=color, line_width=2, size=6)

        # Draw arrows and vertical lines
        for lab_start, lab_end, y_end in zip(lab_rest_list, obs_shift_list, y_end_list):
            arrow = Arrow(
                end=head,
                line_color=color,
                line_width=2,
                x_start=lab_start,  # start at lab wavelength
                y_start=y_end,
                x_end=lab_end,  # end at RV-shifted lab wavelength
                y_end=y_end,
            )

            vline = plot.line(
                x=[lab_end, lab_end],
                y=[1.5 * flux_min, 1.5 * flux_max],
                line_dash="dashed",
                line_color=color,
                line_width=2,
                alpha=0.6,
                legend_label=f"RV = {rv_val:.1f} km/s",
            )

            vline.js_on_change("visible", CustomJS(args=dict(arw=arrow), code="arw.visible = cb_obj.visible;"))

            plot.add_layout(arrow)


def get_rvs(plot: figure, spectrum: Spectrum, prominence: float = 2, smoothing: int = 3, snr: float = 3, **kwargs):
    """
    Calculates and plots multi-component radial velocities of spectral features
    """

    from ....plotting.spectrum.spectrum_overlay import OVERLAY_LINES

    # get wavelengths and peak values of spectral features
    features, peak_vals = do_fitting(plot, spectrum, prominence, smoothing, snr, get_features=True, **kwargs)

    # extract wavelengths
    wav_obs = features * u.AA
    wav_lab = np.array([line["wavelength"] for line in OVERLAY_LINES if min(spectrum.wavelength.value) <= line["wavelength"] <= max(spectrum.wavelength.value)]) * u.AA

    # find radial velocities
    rv_components, matched_lines, rvs, scores = find_rv_components(wav_obs, wav_lab)

    # plot radial velocity overlay
    plot_rv(
        plot=plot,
        spectrum=spectrum,
        wav_lab=wav_lab,
        rv_components=rv_components,
        matched_lines=matched_lines,
        features=wav_obs,
        peak_vals=peak_vals,
    )

    return plot

"""
fig_utils.py — Shared plotting utilities for the Pericyte Paper 2026.

All figure notebooks import from this module. To use:
    from fig_utils import *

Functions are grouped into:
    1. Style configuration
    2. Scalebar utilities
    3. ABF trace plotting
    4. Statistical plots (boxplots, paired dots, heatmaps, correlations)
    5. Fluorescence / calcium imaging plots

Dependencies:
    matplotlib, numpy, pandas, pyabf, seaborn, scipy
"""

import json
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import pyabf
import seaborn as sns
from pathlib import Path
from scipy import stats
from scipy.stats import wilcoxon
from matplotlib.ticker import MultipleLocator

# =============================================================================
# 0 — ABF / CSV LOADING
# =============================================================================

class TraceData:
    """Lightweight stand-in for pyabf.ABF that loads from .npz or .csv + .json.

    Provides the same attributes used by the plotting functions:
    sweepX, sweepY, sweepC, sweepLabelX, sweepLabelY, sweepUnitsX,
    sweepUnitsY, sweepCount, channelCount, sweepList, sampleRate,
    and setSweep().
    """

    def __init__(self, data_path, json_path=None):
        data_path = Path(data_path)
        json_path = Path(json_path) if json_path else data_path.with_suffix('.json')

        # Load metadata
        with open(json_path) as f:
            meta = json.load(f)

        self.sweepCount = meta['sweep_count']
        self.channelCount = meta['channel_count']
        self.sampleRate = meta['sample_rate']
        self.sweepLabelX = meta['sweep_label_x']
        self.sweepLabelY = meta['sweep_label_y']
        self.sweepUnitsX = meta['sweep_units_x']
        self.sweepUnitsY = meta['sweep_units_y']
        self.sweepList = list(range(self.sweepCount))

        # Load arrays from .npz or .csv
        if data_path.suffix == '.npz':
            self._data = dict(np.load(str(data_path)))
        else:
            df = pd.read_csv(data_path)
            self._data = {col: df[col].values for col in df.columns}

        # Initialise to sweep 0, channel 0
        self.sweepX = self._data['time']
        self.sweepY = self._data['ch0_s0']
        self.sweepC = self._data.get('cmd_s0', np.zeros(len(self.sweepX)))

    def setSweep(self, sweepNumber=0, channel=0):
        """Mimic pyabf.ABF.setSweep()."""
        self.sweepY = self._data[f'ch{channel}_s{sweepNumber}']
        self.sweepC = self._data.get(f'cmd_s{sweepNumber}',
                                     np.zeros(len(self.sweepX)))


def load_abf_or_csv(file_path):
    """
    Load a trace from .abf, .npz, or .csv (each with companion .json).

    Auto-detects format by extension. If no extension, tries .abf → .npz → .csv.

    Returns an object with the pyabf.ABF-compatible interface.
    """
    p = Path(file_path)

    if p.suffix == '.npz':
        return TraceData(p)

    if p.suffix == '.csv':
        return TraceData(p)

    if p.suffix == '.abf':
        return pyabf.ABF(str(p))

    # No extension: try all formats
    for ext, loader in [('.abf', lambda x: pyabf.ABF(str(x))),
                        ('.npz', TraceData),
                        ('.csv', TraceData)]:
        candidate = p.with_suffix(ext)
        if candidate.exists():
            return loader(candidate)

    raise FileNotFoundError(f"No .abf, .npz, or .csv found for: {file_path}")


# =============================================================================
# 1 — STYLE CONFIGURATION & SAVE CONTROL
# =============================================================================

# Defaults (override per-notebook by passing kwargs or calling set_style again)
DEFAULT_FONT = 'arial'
DEFAULT_FONTSIZE = 15
DEFAULT_DPI = 600
SAVE_FIGURES = False  # Set to True to enable saving; all save_path is then active


def enable_saving(on=True):
    """
    Toggle figure saving globally.

    Usage (top of notebook):
        enable_saving()       # turn on — all figures with save_path will be saved
        enable_saving(False)  # turn off (default)
    """
    global SAVE_FIGURES
    SAVE_FIGURES = on


def save_if_enabled(save_path):
    """Internal helper: save current figure if SAVE_FIGURES is True and path given."""
    if SAVE_FIGURES and save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=DEFAULT_DPI)

def set_style(font=DEFAULT_FONT, fontsize=DEFAULT_FONTSIZE):
    """
    Apply consistent publication-quality style defaults.
    Call once at the top of each notebook, or per-cell if needed.
    """
    plt.rc('font', family=font)
    plt.rcParams.update({
        'font.size': fontsize,
        'axes.labelsize': fontsize,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'axes.titlesize': fontsize + 1,
        'figure.dpi': DEFAULT_DPI,
    })
    sns.set_style('ticks')


# =============================================================================
# 2 — SCALEBAR UTILITIES
# =============================================================================

def scalebar(abf=None, hideTicks=True, hideFrame=True, fontSize=8,
             scaleXsize=None, scaleYsize=None,
             scaleXunits="", scaleYunits="", lineWidth=2,
             padX=0.10, padY=0.05, textOffset=None):
    """
    Add an L-shaped scalebar to the current axes, removing ticks and frame.

    Parameters
    ----------
    padX : float
        Horizontal padding from the right edge (fraction of x-span). Default 0.10.
    padY : float
        Vertical padding from the bottom edge (fraction of y-span). Default 0.05.

    Source: Scott Harden, pyABF
    https://github.com/swharden/pyABF/blob/main/src/pyabf/plot.py
    License: MIT  |  Retrieved: 2026-03-19
    """
    if abf:
        scaleXunits = abf.sweepUnitsX
        scaleYunits = abf.sweepUnitsY

    x1, x2, y1, y2 = plt.axis()
    xs, ys = abs(x2 - x1), abs(y2 - y1)

    if not scaleXsize:
        scaleXsize = abs(plt.xticks()[0][1] - plt.xticks()[0][0]) / 2
    if not scaleYsize:
        scaleYsize = abs(plt.yticks()[0][1] - plt.yticks()[0][0]) / 2

    lblX = _format_scalebar_label(scaleXsize, scaleXunits)
    lblY = _format_scalebar_label(scaleYsize, scaleYunits)

    _draw_scalebar(x1, x2, y1, y2, xs, ys,
                   scaleXsize, scaleYsize, lblX, lblY,
                   hideTicks, hideTicks, hideFrame, fontSize, lineWidth,
                   padX, padY, textOffset=textOffset)

def scalebar_concat(abfs=None, hideXTicks=True, hideYTicks=True, hideFrame=True,
                    fontSize=8, scaleXsize=None, scaleYsize=None,
                    scaleXunits="", scaleYunits="", lineWidth=2,
                    padX=0.10, padY=0.05, textOffset=None):
    """
    Scalebar variant for concatenated ABF traces where tick-based auto-sizing
    is unreliable. Computes round scalebar sizes from the visible axis range
    (~20 % of span).

    Uses the first ABF's units. Otherwise identical layout to scalebar().

    Source: adapted from pyABF (Scott W Harden, MIT license)
    """
    if abfs and len(abfs) > 0:
        scaleXunits = getattr(abfs[0], 'sweepUnitsX',
                              abfs[0].sweepLabelX.split('(')[-1].split(')')[0])
        scaleYunits = getattr(abfs[0], 'sweepUnitsY',
                              abfs[0].sweepLabelY.split('(')[-1].split(')')[0])

    x1, x2, y1, y2 = plt.axis()
    xs, ys = abs(x2 - x1), abs(y2 - y1)

    if not scaleXsize:
        if xs < 1:
            scaleXsize = round(xs * 0.2 * 1000) / 1000
        elif xs < 10:
            scaleXsize = round(xs * 0.2)
        elif xs < 60:
            scaleXsize = round(xs * 0.2 / 5) * 5
        else:
            scaleXsize = round(xs * 0.2 / 10) * 10

    if not scaleYsize:
        if ys < 10:
            scaleYsize = round(ys * 0.2 * 10) / 10
        elif ys < 50:
            scaleYsize = round(ys * 0.2)
        else:
            scaleYsize = round(ys * 0.2 / 5) * 5

    lblX = _format_scalebar_label(scaleXsize, scaleXunits)
    lblY = _format_scalebar_label(scaleYsize, scaleYunits)

    _draw_scalebar(x1, x2, y1, y2, xs, ys,
                   scaleXsize, scaleYsize, lblX, lblY,
                   hideXTicks, hideYTicks, hideFrame, fontSize, lineWidth,
                   padX, padY, textOffset=textOffset)
    


# ── Scalebar internals ──────────────────────────────────────────────────────

def _format_scalebar_label(size, units):
    """Format a scalebar value + units string, converting sec→ms if needed."""
    lbl = str(size)
    if lbl.endswith(".0"):
        lbl = lbl[:-2]
    if units in ("sec", "s", "seconds") and "." in lbl:
        lbl = str(int(float(lbl) * 1000))
        units = "ms"
    return f"{lbl} {units}".strip()


def _draw_scalebar(x1, x2, y1, y2, xs, ys,
                   scaleXsize, scaleYsize, lblX, lblY,
                   hideXTicks, hideYTicks, hideFrame,
                   fontSize, lineWidth, padX=0.10, padY=0.05,
                   textOffset=None):
    """Shared drawing logic for both scalebar variants."""
    bx = x2 - padX * xs
    by = y1 + padY * ys
    bar_xs = [bx - scaleXsize, bx, bx]
    bar_ys = [by, by, by + scaleYsize]
    if textOffset is None:
        lbl_pad = xs * (0.005 + 0.002 * lineWidth)
    else:
        lbl_pad = textOffset

    ax = plt.gca()
    if hideYTicks:
        ax.get_yaxis().set_visible(False)
    if hideXTicks:
        ax.get_xaxis().set_visible(False)
    if hideFrame:
        for spine in ax.spines.values():
            spine.set_visible(False)

    plt.plot(bar_xs, bar_ys, 'k-', lw=lineWidth)
    plt.text((bx - scaleXsize + bx) / 2, by - lbl_pad, lblX,
             ha='center', va='top', fontsize=fontSize)
    plt.text(bx + lbl_pad, (by + by + scaleYsize) / 2, lblY,
             ha='left', va='center', fontsize=fontSize)


# =============================================================================
# 3 — ABF TRACE PLOTTING
# =============================================================================

def despine(ax, top=True, right=True, bottom=False, left=False):
    """Hide selected spines on an axis."""
    for side, hide in [('top', top), ('right', right),
                       ('bottom', bottom), ('left', left)]:
        ax.spines[side].set_visible(not hide)


def plot_single_trace(abf_file, xlim=None, ylim=None,
                      color='black', figsize=(15, 8), lw=0.5,
                      fontsize=DEFAULT_FONTSIZE,
                      scalebar_fontsize=20,
                      scaleXsize=None, scaleYsize=None,
                      padX=0.10, padY=0.05,
                      use_scalebar=True, save_path=None):
    """
    Plot a single-sweep, single-channel ABF recording with optional scalebar.

    Parameters
    ----------
    abf_file : str          Path to ABF file.
    xlim, ylim : tuple      Axis limits.
    color : str             Trace colour.
    figsize : tuple         Figure size in inches.
    lw : float              Line width.
    fontsize : int          Font size for axis labels and ticks.
    scalebar_fontsize : int Font size for scalebar text.
    scaleXsize, scaleYsize  Custom scalebar dimensions.
    use_scalebar : bool     Scalebar (True) or standard axes (False).
    save_path : str         File path for saving. Only saves if enable_saving() was called.
    """
    set_style(fontsize=fontsize)
    abf = load_abf_or_csv(abf_file)
    fig, ax = plt.subplots(figsize=figsize, dpi=DEFAULT_DPI)

    ax.plot(abf.sweepX, abf.sweepY, color=color, lw=lw)
    ax.set_xlabel(abf.sweepLabelX, fontsize=fontsize)
    ax.set_ylabel(abf.sweepLabelY, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    despine(ax)

    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    if use_scalebar:
        scalebar(abf, fontSize=scalebar_fontsize, hideFrame=True,
                 scaleXsize=scaleXsize, scaleYsize=scaleYsize,
                 padX=padX, padY=padY)

    plt.tight_layout()
    save_if_enabled(save_path)
    plt.show()


def plot_two_channels(abf_file, xlim_0=None, ylim_0=None,
                      xlim_1=None, ylim_1=None,
                      figsize=(8, 8), fontsize=DEFAULT_FONTSIZE,
                      scalebar_fontsize=20,
                      scaleXsize=None, scaleYsize=None,
                      padX=0.10, padY=0.05, textOffset=None,
                      save_path=None):
    """
    Plot channels 0 and 1 from a single ABF file as vertically
    stacked subplots with scalebars.

    Parameters
    ----------
    scaleXsize, scaleYsize : float, optional
        Custom scalebar dimensions. Auto-calculated if None.
    """
    set_style(fontsize=fontsize)
    abf = load_abf_or_csv(abf_file)
    fig = plt.figure(figsize=figsize, dpi=DEFAULT_DPI)

    for idx, (ch, xlim, ylim) in enumerate([
        (0, xlim_0, ylim_0),
        (1, xlim_1, ylim_1),
    ]):
        plt.subplot(2, 1, idx + 1)
        abf.setSweep(sweepNumber=0, channel=ch)
        plt.plot(abf.sweepX, abf.sweepY, color='black')
        plt.xlabel(abf.sweepLabelX, fontsize=fontsize)
        plt.ylabel(abf.sweepLabelY, fontsize=fontsize)
        plt.tick_params(labelsize=fontsize)
        if xlim and ylim:
            plt.axis([*xlim, *ylim])
        scalebar(abf, fontSize=scalebar_fontsize, hideFrame=True,
                 scaleXsize=scaleXsize, scaleYsize=scaleYsize,
                 padX=padX, padY=padY, textOffset=textOffset)

    plt.tight_layout()
    save_if_enabled(save_path)
    plt.show()


def plot_average_sweeps(abf_file, x_limits=None, y_limits=None,
                        figsize=None, add_scalebar=False, dpi=DEFAULT_DPI,
                        fontsize=DEFAULT_FONTSIZE, scalebar_fontsize=20,
                        padX=0.10, padY=0.05,
                        save_path=None):
    """
    Average all sweeps per channel and plot as two vertically stacked subplots.
    Used for dual-patch current-clamp recordings.
    """
    x_limits = x_limits or (0, 2)
    y_limits = y_limits or (-90, -65)
    figsize = figsize or (6, 6)
    set_style(fontsize=fontsize)

    abf = load_abf_or_csv(abf_file)
    num_sweeps = abf.sweepCount

    avg_ch0 = np.zeros_like(abf.sweepY)
    avg_ch1 = np.zeros_like(abf.sweepY)
    for sweep in range(num_sweeps):
        abf.setSweep(sweepNumber=sweep, channel=0)
        avg_ch0 += abf.sweepY
        abf.setSweep(sweepNumber=sweep, channel=1)
        avg_ch1 += abf.sweepY
    avg_ch0 /= num_sweeps
    avg_ch1 /= num_sweeps

    fig = plt.figure(figsize=figsize, dpi=dpi)
    for idx, (avg, color, label) in enumerate([
        (avg_ch0, '#DE2514', "Average sweep — Channel 0"),
        (avg_ch1, '#11A238', "Average sweep — Channel 1"),
    ]):
        plt.subplot(2, 1, idx + 1)
        plt.plot(abf.sweepX, avg, label=label, color=color)
        plt.xlabel(abf.sweepLabelX, fontsize=fontsize)
        plt.ylabel(abf.sweepLabelY, fontsize=fontsize)
        plt.tick_params(labelsize=fontsize)
        plt.axis([*x_limits, *y_limits])
        if add_scalebar:
            scalebar(abf=abf, fontSize=scalebar_fontsize, hideFrame=True,
                     padX=padX, padY=padY)

    plt.tight_layout()
    save_if_enabled(save_path)
    plt.show()


def plot_concatenated_traces(abf_files, color='black',
                             xlim=None, ylim=None,
                             figsize=(4, 4), lw=0.5,
                             fontsize=DEFAULT_FONTSIZE,
                             scalebar_fontsize=20,
                             scaleXsize=None, scaleYsize=None,
                             padX=0.10, padY=0.05,
                             use_concat_scalebar=True,
                             save_path=None):
    """
    Load multiple ABF files and plot them end-to-end as a single
    continuous trace. Uses scalebar_concat by default for proper
    auto-sizing on concatenated time axes.

    Parameters
    ----------
    use_concat_scalebar : bool
        If True (default), use scalebar_concat. If False, use standard scalebar.
    """
    set_style(fontsize=fontsize)
    abfs = [load_abf_or_csv(f) for f in abf_files]
    fig, ax = plt.subplots(figsize=figsize, dpi=DEFAULT_DPI)

    t_offset = 0.0
    for abf in abfs:
        t_shifted = abf.sweepX - abf.sweepX[0] + t_offset
        ax.plot(t_shifted, abf.sweepY, color, lw=lw)
        t_offset = t_shifted[-1]

    ax.set_xlabel(abfs[0].sweepLabelX, fontsize=fontsize)
    ax.set_ylabel(abfs[0].sweepLabelY, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    despine(ax)

    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)

    sb_func = scalebar_concat if use_concat_scalebar else scalebar
    sb_kwargs = dict(fontSize=scalebar_fontsize, hideFrame=True,
                     scaleXsize=scaleXsize, scaleYsize=scaleYsize)
    if use_concat_scalebar:
        sb_func(abfs, padX=padX, padY=padY, **sb_kwargs)
    else:
        sb_func(abf=abfs[0], padX=padX, padY=padY, **sb_kwargs)

    plt.tight_layout()
    save_if_enabled(save_path)
    plt.show()


def plot_stacked_sweeps(abf_file, y_offset=15, duration_s=1.0,
                        figsize=(10, 10), fontsize=DEFAULT_FONTSIZE,
                        scalebar_fontsize=20,
                        scaleXsize=None, scaleYsize=None,
                        padX=0.10, padY=0.05,
                        save_path=None):
    """
    Plot all sweeps stacked vertically with a y-offset.
    Typically used for IV step protocols.
    """
    set_style(fontsize=fontsize)
    abf = load_abf_or_csv(abf_file)
    fig = plt.figure(figsize=figsize, dpi=DEFAULT_DPI)

    for sweep_num in abf.sweepList:
        abf.setSweep(sweep_num)
        i2 = int(abf.sampleRate * duration_s)
        dataX = abf.sweepX[:i2]
        dataY = abf.sweepY[:i2] + y_offset * sweep_num
        plt.plot(dataX, dataY, color='black', alpha=0.5)

    plt.xlabel(abf.sweepLabelX, fontsize=fontsize)
    plt.ylabel(abf.sweepLabelY, fontsize=fontsize)
    plt.tick_params(labelsize=fontsize)
    scalebar(abf, fontSize=scalebar_fontsize, hideFrame=True,
             scaleXsize=scaleXsize, scaleYsize=scaleYsize,
             padX=padX, padY=padY)

    plt.tight_layout()
    save_if_enabled(save_path)
    plt.show()


def visualize_abf_data(abf_files, xlim=None, ylim=None,
                       ylabel=None,
                       figsize=(15, 4), fontsize=DEFAULT_FONTSIZE,
                       scalebar_fontsize=20,
                       padX=0.10, padY=0.05,
                       scalebar_x=None, scalebar_y=None, 
                       yticks=None,
                       add_scalebar=False, show_yaxis=True,
                       hlines=None, save_path=None):
    set_style(fontsize=fontsize)
    abfs = [load_abf_or_csv(f) for f in abf_files]
    fig, ax = plt.subplots(figsize=figsize, dpi=DEFAULT_DPI)
    t_offset = 0.0
    for abf in abfs:
        ax.plot(t_offset + abf.sweepX, abf.sweepY, 'black', lw=0.5)
        t_offset += abf.sweepX[-1] - abf.sweepX[0]
    ax.set_xlabel(abfs[0].sweepLabelX, fontsize=fontsize + 5)
    ax.set_ylabel(ylabel, fontsize=fontsize + 5)
    ax.tick_params(labelsize=fontsize + 5)
    despine(ax)
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    if yticks:
        ax.set_yticks(yticks)
        ax.spines['left'].set_bounds(min(yticks), max(yticks))
    if add_scalebar:
        _x1, _x2, _y1, _y2 = plt.axis()
        _xs, _ys = abs(_x2 - _x1), abs(_y2 - _y1)
        _padX = (_x2 - scalebar_x) / _xs if scalebar_x is not None else padX
        _padY = (scalebar_y - _y1) / _ys if scalebar_y is not None else padY
        scalebar(abf=abfs[0], fontSize=scalebar_fontsize,
                 hideFrame=False, hideTicks=False,
                 padX=_padX, padY=_padY, textOffset=1)
        ax.get_xaxis().set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if not show_yaxis:
            ax.get_yaxis().set_visible(False)
            ax.spines['left'].set_visible(False)
    if hlines:
        for y in hlines:
            ax.axhline(y=y, color='black', lw=1, dashes=(20, 20))
    plt.tight_layout()
    save_if_enabled(save_path)
    plt.show()


# =============================================================================
# 4 — STATISTICAL / SUMMARY PLOTS
# =============================================================================

def plot_boxplot_scatter(data_list, tick_labels, box_colors=None,
                         ylabel=r'$\Delta$ Membrane potential (mV)',
                         ylim=None, figsize=(3.5, 5), box_width=1.0,
                         fontsize=DEFAULT_FONTSIZE, tick_fontsize=None,
                         ytick_interval=None, rotation=60, spines=False,
                         xlim_right_pad=0, hline_y=None, text_annotations=None, 
                         save_path=None):
    """
    Boxplot with individual data points overlaid and optional
    horizontal reference line.
 
    Parameters
    ----------
    data_list : list of array-like   One array per group.
    tick_labels : list of str        X-axis labels.
    box_colors : list/str/None       Box fill colours. Defaults to grey.
    ylabel : str                     Y-axis label.
    ylim : tuple                     Y-axis limits.
    figsize : tuple                  Figure size.
    box_width : float                Width of each box. Default 1.0.
    fontsize : int                   Label/tick font size.
    tick_fontsize : int              X-tick font size (defaults to fontsize).
    ytick_interval : float           Major tick interval for y-axis.
    rotation : int                   X-tick label rotation in degrees. Default 45.
    hline_y : float                  Dashed reference line.
    save_path : str                  If given, save figure.
    """
    set_style(fontsize=fontsize)
    tick_fontsize = tick_fontsize or fontsize
 
    if box_colors is None:
        box_colors = ['#BAB9BD'] * len(data_list)
    elif isinstance(box_colors, str):
        box_colors = [box_colors] * len(data_list)
 
    fig, ax = plt.subplots(figsize=figsize, dpi=DEFAULT_DPI)
    n = len(data_list)
    positions = list(range(0, n * 2, 2)) 
        
    bpl = ax.boxplot(data_list, positions=positions, sym='',
                     widths=box_width, patch_artist=True)
 
    for box, color in zip(bpl['boxes'], box_colors):
        box.set_facecolor(color)
        box.set_linewidth(1)
        box.set(edgecolor='black')
    for med in bpl['medians']:
        med.set_color('black')
        med.set_linewidth(1.8)
    for wh in bpl['whiskers']:
        wh.set_linewidth(2)
 
    for pos, vals in zip(positions, data_list):
        ax.scatter(np.ones(len(vals)) * pos, vals,
                   color='black', edgecolor='black', zorder=10, alpha=0.5)
 
    if hline_y is not None:
        ax.axhline(hline_y, color='black', linestyle='--', zorder=2, alpha=0.2)
    
    ax.set_xticks(positions)
    ax.set_xticklabels(tick_labels, fontsize=tick_fontsize,
                       rotation=rotation, va='top')
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_xlim(-2, n * 2 + xlim_right_pad) 
    if ylim:
        ax.set_ylim(*ylim)
    if ytick_interval:
        ax.yaxis.set_major_locator(MultipleLocator(ytick_interval))
    if not spines:
        despine(ax)
 
    if text_annotations:
        for annotation in text_annotations:
            x = annotation.pop('x')
            y = annotation.pop('y')
            text = annotation.pop('text')
            ax.text(x, y, text, **annotation)
    else:
        plt.tight_layout()

    save_if_enabled(save_path)
    plt.show()


def plot_paired_dots(csv_file, col_a, col_b,
                     label_a="Group A", label_b="Group B",
                     color_lines='black', figsize=(2.2, 1.6),
                     fontsize=10, save_path=None):
    """
    Paired dot-plot showing coupling coefficients in both directions,
    connected by lines for each recording pair.

    Set color_lines='direction' to colour by increase (grey) vs. decrease
    (darkblue).
    """
    set_style(fontsize=fontsize)
    df = pd.read_csv(csv_file)
    a = df[col_a].values
    b = df[col_b].values

    fig = plt.figure(figsize=figsize, dpi=DEFAULT_DPI)
    plt.scatter(np.zeros(len(a)), a, color='black', marker='o', s=20)
    plt.scatter(np.ones(len(b)), b, color='black', marker='o', s=20)

    for i in range(len(a)):
        c = ('grey' if b[i] > a[i] else 'darkblue') \
            if color_lines == 'direction' else color_lines
        plt.plot([0, 1], [a[i], b[i]], c=c, linewidth=0.5)

    plt.xticks([0, 1], [label_a, label_b], rotation=60, fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.ylabel('Coupling Coeffi- \ncient: Creci / Cstim', fontsize=fontsize)

    ax = plt.gca()
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    ax.set_yticklabels([f'{y:.1f}' for y in np.arange(0, 1.1, 0.2)])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='y', which='major', length=4, width=0.8, 
               left=True, right=False, top=False, bottom=False)
    ax.tick_params(axis='x', which='major', length=2, width=0.8,
               left=False, right=False, top=False, bottom=True)
    plt.ylim(0, 1)
    plt.xlim(-0.5, 1.5)
    plt.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.1)

    save_if_enabled(save_path)
    plt.show()


def plot_coupling_heatmap(csv_file, columns, labels,
                          figsize=(12, 8), square=True, annot_size=8,
                          fontsize=DEFAULT_FONTSIZE,
                          title=None, save_path=None, cbar_ticksize=None):
    """
    Coupling-strength heatmap from an Excel matrix file.
    """
    set_style(fontsize=fontsize)
    df = pd.read_csv(csv_file)[columns]
    fig, ax = plt.subplots(figsize=figsize, dpi=DEFAULT_DPI)
    sns.set_context('paper', font_scale=1.4)
    h = sns.heatmap(
        df.T, annot=True, fmt=".2f", cmap='YlGnBu',
        vmin=0, vmax=1.0, linewidths=0.5, linecolor='white',
        square=square,
        cbar_kws={'label': 'Coupling Strength', 'orientation': 'vertical',
                'shrink': 0.3, 'ticks': [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]},
        yticklabels=labels,
        annot_kws={'size': annot_size},
        ax=ax,
    )
    if cbar_ticksize:
        h.collections[0].colorbar.ax.tick_params(labelsize=cbar_ticksize)
    ax.set_xticks([])
    ax.yaxis.set_visible(True)
    ax.tick_params(axis='y', left=True)
    ax.set_yticklabels(labels, rotation=0)
    if title:
        plt.title(title, pad=20, y=1.1)
    plt.tight_layout()
    save_if_enabled(save_path)
    plt.show()


# =============================================================================
# 5 — FLUORESCENCE / CALCIUM IMAGING PLOTS
# =============================================================================

def plot_fluorescence_with_command(fluor_csv, fluor_x_col, fluor_y_col,
                                   abf_file, fluor_color='#076c89',
                                   title="", fluor_ylim=None,
                                   figsize=(6.5,5.5), fontsize=DEFAULT_FONTSIZE,
                                   save_path=None):
    """
    Dual-axis plot: fluorescence ΔF/F₀ (left y-axis) overlaid with the
    ABF command signal (right y-axis, dashed grey).

    The right axis is automatically aligned so that −100 mV on the
    command axis corresponds to 0.0 on the fluorescence axis.
    """
    df = pd.read_csv(fluor_csv)
    x_fluor = df[fluor_x_col]
    y_fluor = df[fluor_y_col]

    abf = load_abf_or_csv(abf_file)
    cmd_x_all, cmd_y_all = [], []
    for sweep_num in abf.sweepList:
        abf.setSweep(sweep_num)
        t_offset = sweep_num * abf.sweepX[-1]
        cmd_x_all.append(abf.sweepX + t_offset)
        cmd_y_all.append(abf.sweepC - np.min(abf.sweepC))
    cmd_x = np.concatenate(cmd_x_all)
    cmd_y = np.concatenate(cmd_y_all) - 100

    set_style(fontsize=fontsize)
    fig, ax1 = plt.subplots(figsize=figsize, dpi=DEFAULT_DPI)

    ax1.plot(x_fluor, y_fluor, color=fluor_color, zorder=1)
    ax1.set_ylabel(r"$\Delta F / F_0$", color=fluor_color, fontsize=fontsize)
    ax1.set_xlabel("Time (s)", fontsize=fontsize)
    ax1.tick_params(axis='y', colors='black', labelsize=fontsize)
    ax1.tick_params(axis='x', labelsize=fontsize)
    ax1.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
    ax1.set_yticks([0.25, 0.75, 1.25, 1.75], minor=True)
    despine(ax1)
    if title:
        ax1.text(0.05, 0.95, title, transform=ax1.transAxes,
             fontsize=fontsize + 1, va='top')
    if fluor_ylim:
        ax1.set_ylim(*fluor_ylim)

    ax2 = ax1.twinx()
    ax2.plot(cmd_x, cmd_y, color='darkgray', lw=1, ls='-', alpha=0.7, zorder=-2)
    ax2.set_ylabel("Command\nVoltage (mV)", color='darkgray',
                    rotation=90, labelpad=10, fontsize=fontsize)
    ax2.tick_params(axis='y', colors='black', labelsize=fontsize)
    ax2.spines['right'].set_color('black')
    ax2.spines['top'].set_visible(False)
    ax2.set_yticks([-100, -60, -20, 20, 60])
    ax2.set_yticks([-80, -40, 0, 40], minor=True)

    # Align −100 mV (right) with 0.0 ΔF/F₀ (left)
    l_min, l_max = ax1.get_ylim()
    r_max = 60
    frac = (0 - l_min) / (l_max - l_min)
    r_min = (-100 - r_max * frac) / (1 - frac)
    ax2.set_ylim(r_min, r_max)

    plt.tight_layout()
    save_if_enabled(save_path)
    plt.show()
"""Shared typography, input loading and export for publication panels."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#000000", "#F0E442"]
PROGRAMMES = ["Replication", "Redox", "Matrix"]
PROGRAMME_COLORS = dict(zip(PROGRAMMES, [PALETTE[0], PALETTE[2], PALETTE[1]]))
SEQUENTIAL = colors.LinearSegmentedColormap.from_list("OkabeBlue", ["#FFFFFF", "#56B4E9", "#0072B2"])
DIVERGING = colors.LinearSegmentedColormap.from_list("OkabeDiverging", ["#0072B2", "#FFFFFF", "#D55E00"])
METHODS = ["No change", "Target mean", "Paired mean", "Unpaired mean", "Nearest pair", "Global anchor", "Paired midpoint", "PCHT"]


def configure_style(font_dir=None):
    """Select Liberation Sans, the Okabe–Ito palette and embedded vector fonts."""
    if font_dir:
        for path in sorted(Path(font_dir).glob("*.ttf")):
            font_manager.fontManager.addfont(str(path))
    available = {font.name for font in font_manager.fontManager.ttflist}
    if "Liberation Sans" not in available:
        raise RuntimeError("Liberation Sans is required. Supply its font directory with --font-dir.")
    plt.rcParams.update({
        "font.family": "Liberation Sans", "font.size": 7,
        "axes.labelsize": 7, "axes.titlesize": 7.5, "axes.titleweight": "bold",
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": .65, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "lines.linewidth": .8, "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none", "savefig.facecolor": "white",
        "mathtext.fontset": "dejavusans",
    })


def wrap_label(value, width=16):
    return "\n".join(textwrap.wrap(str(value), width, break_long_words=False, break_on_hyphens=False))


def load_panel(data_dir, panel):
    """Load panel-ready tables and their external display metadata."""
    directory = Path(data_dir)
    metadata = json.loads((directory / f"{panel}.json").read_text(encoding="utf-8"))
    if metadata.get("panel") != panel:
        raise ValueError(f"Panel identity mismatch for {panel}")
    inputs = {"metadata": metadata}
    for key, descriptor in metadata["datasets"].items():
        source = directory / descriptor["file"]
        inputs[key] = pd.read_csv(source, index_col=0 if descriptor.get("indexed") else None)
    return inputs


def scatter(ax, x, y, color=PALETTE[0], label=None, size=12, alpha=.8, **kwargs):
    return ax.scatter(x, y, s=size, c=color, alpha=alpha, linewidths=0, label=label, **kwargs)


def legend(ax, **kwargs):
    item = ax.legend(frameon=False, handlelength=1.2, handletextpad=.45, borderpad=.1, **kwargs)
    if item:
        for handle in item.legend_handles:
            handle.set_alpha(1)
    return item


def colorbar(ax, artist, label):
    cb = ax.figure.colorbar(artist, cax=ax.inset_axes([1.035, .05, .035, .9]))
    cb.set_label(wrap_label(label, 22), fontsize=7, labelpad=3)
    cb.ax.tick_params(labelsize=7, length=2)
    from matplotlib.ticker import MaxNLocator
    cb.locator = MaxNLocator(4)
    cb.update_ticks()
    return cb


def label_indices(labels, max_labels, preferred=()):
    if max_labels is None or len(labels) <= max_labels:
        return list(range(len(labels)))
    step = max(1, int(np.ceil(len(labels) / max_labels)))
    keep = list(range(0, len(labels), step))
    for j, label in enumerate(labels):
        if any(word in str(label) for word in preferred):
            keep = [k for k in keep if abs(k - j) >= step]
            keep.append(j)
    return sorted(set(keep))


def heatmap(ax, data, label, cmap=SEQUENTIAL, vmin=None, vmax=None, labels=True, rotation=0, max_y_labels=None, max_x_labels=None):
    image = ax.imshow(data.to_numpy(dtype=float), cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax, interpolation="nearest")
    if labels:
        rows, columns = [str(x) for x in data.index], [str(x) for x in data.columns]
        max_y_labels = max_y_labels if max_y_labels is not None else getattr(ax, "_max_heatmap_y_labels", None)
        max_x_labels = max_x_labels if max_x_labels is not None else getattr(ax, "_max_heatmap_x_labels", None)
        yi = label_indices(rows, max_y_labels, ["TGFB1"])
        xi = label_indices(columns, max_x_labels)
        ax.set_yticks(yi, [wrap_label(rows[j], 24) for j in yi])
        ax.set_xticks(xi, [wrap_label(columns[j], 16) for j in xi], rotation=rotation, ha="right" if rotation else "center")
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar(ax, image, label)
    return image


def spatial_map(ax, data, label, metadata):
    image = scatter(ax, data.x, data.y, data.value, size=3.5, alpha=1, cmap=SEQUENTIAL,
                    vmin=metadata["vmin"], vmax=metadata["vmax"], rasterized=True)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar(ax, image, label)
    return image


def annotate_genes(ax, data, offsets):
    indexed = data.set_index("gene")
    for gene, offset in offsets.items():
        if gene in indexed.index:
            row = indexed.loc[gene]
            ax.annotate(gene, (row.x, row.y), xytext=offset, textcoords="offset points", fontsize=7,
                        arrowprops={"arrowstyle": "-", "color": "black", "lw": .4})


def export_figure(fig, output_dir, stem):
    """Export vector masters, a 350 dpi preview and an automatic layout report."""
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    axes = list(fig.axes)
    for ax in fig.axes:
        axes.extend(ax.child_axes)
    for ax in axes:
        ax.tick_params(axis="both", which="both", labelsize=7, length=2.5, pad=2)
        for axis in (ax.xaxis, ax.yaxis):
            if axis.get_scale() == "log":
                axis.set_major_formatter(FuncFormatter(lambda x, position: f"{x:.3g}"))
                axis.set_minor_formatter(FuncFormatter(lambda x, position: ""))
        for text in ax.findobj(matplotlib.text.Text):
            if text.get_text():
                text.set_fontfamily("Liberation Sans")
                text.set_fontsize(max(7, text.get_fontsize()))
    for line in fig.findobj(matplotlib.lines.Line2D):
        if line.get_linewidth() > 0:
            line.set_linewidth(float(np.clip(line.get_linewidth(), .25, 1)))
    for collection in fig.findobj(matplotlib.collections.Collection):
        widths = collection.get_linewidths()
        if len(widths):
            collection.set_linewidths(np.where(np.asarray(widths) > 0, np.clip(widths, .25, 1), 0))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    visible_text = list(fig.texts)
    for ax in axes:
        visible_text.extend(list(ax.texts) + [ax.title, ax._left_title, ax._right_title, ax.xaxis.label, ax.yaxis.label])
        if ax.get_legend():
            visible_text.extend(ax.get_legend().get_texts())
        for axis in (ax.xaxis, ax.yaxis):
            lower, upper = sorted(axis.get_view_interval())
            tolerance = (upper - lower) * 1e-8
            for tick in list(axis.get_major_ticks()) + list(axis.get_minor_ticks()):
                if lower - tolerance <= tick.get_loc() <= upper + tolerance:
                    visible_text.extend([tick.label1, tick.label2])
    for text in visible_text:
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
            outside.append(text.get_text())
    for extension in ("pdf", "svg", "png"):
        fig.savefig(destination / f"{stem}.{extension}", dpi=350)
    report = {"figure": stem, "dimensions_mm": [round(v * 25.4, 2) for v in fig.get_size_inches()],
              "font_family": "Liberation Sans", "minimum_font_pt": 7, "text_outside_canvas": outside}
    (destination / f"{stem}_layout.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def run_panel(draw, panel, title):
    parser = argparse.ArgumentParser(description=f"Render {panel}: {title}.")
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--font-dir", type=Path)
    args = parser.parse_args()
    configure_style(args.font_dir)
    fig = plt.figure(figsize=(110 / 25.4, 95 / 25.4))
    ax = fig.add_axes([.34, .25, .48, .58])
    draw(ax, load_panel(args.data_dir, panel))
    ax.set_title(wrap_label(title, 42), loc="left", pad=8, fontsize=7.5)
    report = export_figure(fig, args.output_dir, panel)
    plt.close(fig)
    print(json.dumps(report))

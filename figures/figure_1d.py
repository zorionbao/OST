"""Render pure-state matrix calibration from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    for method, color, marker in [('CRMP fixed', PALETTE[1], 's'), ('Pure-state CRMP', PALETTE[0], 'o')]:
        points = data[data.method == method]
        scatter(ax, points.truth, points.prediction, color, method, size=20, marker=marker)
    ax.plot([meta['identity_min'], meta['identity_max']], [meta['identity_min'], meta['identity_max']], ':', color='black', lw=.7)
    ax.set(xlabel='Pure osteogenic matrix rank', ylabel='Predicted matrix rank')
    legend(ax, loc='lower right')


if __name__ == "__main__":
    run_panel(draw, "figure_1d", 'Pure-state matrix calibration')

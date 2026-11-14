"""Render patient and model state space from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    for category, color, marker, size, alpha in [('Patient tumours', PALETTE[0], 'o', 10, .5), ('Cell models', PALETTE[1], '^', 22, .8)]:
        points = data[data.category == category]
        scatter(ax, points.x, points.y, color, category, size=size, alpha=alpha, marker=marker)
    ax.set(xlabel='Matrix score', ylabel='Replication score')
    legend(ax)


if __name__ == "__main__":
    run_panel(draw, "figure_2d", 'Patient and model state space')

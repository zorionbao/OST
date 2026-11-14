"""Render matched context depletion from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    for category, color in [('Both contexts', PALETTE[0]), ('Lung threshold only', PALETTE[1]), ('Other genes', PALETTE[2])]:
        points = data[data.category == category]
        scatter(ax, points.x, points.y, color, category, size=13, alpha=.7)
    annotate_genes(ax, data, {'BCL6': (4, -18), 'KLF4': (-35, 10), 'STAT3': (8, 4)})
    ax.set(xlabel='Culture log2 fold change', ylabel='Lung log2 fold change')
    legend(ax, loc='lower right')


if __name__ == "__main__":
    run_panel(draw, "figure_2a", 'Matched context depletion')

"""Render independent model dependencies from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    order = inputs['metadata']['gene_order']
    for i, gene in enumerate(order):
        points = data[data.gene == gene]
        scatter(ax, points.value, i + np.linspace(-.16, .16, len(points)), PALETTE[0], size=8, alpha=.7)
    ax.set_yticks(range(len(order)), order)
    ax.set_xlabel('Chronos gene effect')
    ax.invert_yaxis()


if __name__ == "__main__":
    run_panel(draw, "figure_2c", 'Independent model dependencies')

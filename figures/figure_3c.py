"""Render distinct environmental response axes from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    scatter(ax, data.x, data.y, PALETTE[0], size=8, alpha=.6)
    annotate_genes(ax, data, {gene: (4, 4) for gene in ['ITGB1', 'PTK2', 'GPX4', 'CDK6']})
    ax.set(xlabel='Dome − plastic rank', ylabel='Serum-free − serum rank')


if __name__ == "__main__":
    run_panel(draw, "figure_3c", 'Distinct environmental response axes')

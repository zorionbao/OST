"""Render response projection in osteosarcoma from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    spatial_map(ax, inputs['data'], 'TGFB1 response', inputs['metadata'])


if __name__ == "__main__":
    run_panel(draw, "figure_4c", 'Response projection in osteosarcoma')

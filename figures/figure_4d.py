"""Render matrix programme in the same section from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    spatial_map(ax, inputs['data'], 'Matrix programme', inputs['metadata'])


if __name__ == "__main__":
    run_panel(draw, "figure_4d", 'Matrix programme in the same section')

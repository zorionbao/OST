"""Render gene-disjoint spatial associations from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    heatmap(ax, data, 'Median partial ρ', cmap=DIVERGING, vmin=-meta['color_limit'], vmax=meta['color_limit'])


if __name__ == "__main__":
    run_panel(draw, "figure_4e", 'Gene-disjoint spatial associations')

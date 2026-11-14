"""Render paired redistribution of tissue states from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    heatmap(ax, data, 'Relapse − primary', cmap=DIVERGING,
            vmin=-meta['color_limit'], vmax=meta['color_limit'], rotation=45)


if __name__ == "__main__":
    run_panel(draw, "figure_4h", 'Paired redistribution of tissue states')

"""Render the measured dependency landscape from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    heatmap(ax, data, 'Fitness rank', vmin=0, vmax=1, labels=False)
    ax.set(xlabel=f'{len(data.columns)} genes', ylabel=f'{len(data)} screens')


if __name__ == "__main__":
    run_panel(draw, "figure_3a", 'The measured dependency landscape')

"""Render gene-level held-out predictions from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    image = ax.hexbin(data.observed, data.predicted, gridsize=inputs['metadata']['hexbin_gridsize'],
                      mincnt=1, cmap=SEQUENTIAL, linewidths=0, rasterized=True)
    ax.set(xlabel='Observed fitness rank', ylabel='Predicted fitness rank')
    colorbar(ax, image, 'Gene–context pairs')


if __name__ == "__main__":
    run_panel(draw, "figure_3f", 'Gene-level held-out predictions')

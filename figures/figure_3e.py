"""Render error structure across models from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data'].copy()
    data.index = [str(identifier).replace('ACH-', '') for identifier in data.index]
    image = heatmap(ax, data, 'RMSE', vmin=min(.1, float(data.min().min())), vmax=max(.4, float(data.max().max())))
    ax.set_xticks(range(8), [f'M{i}' for i in range(1, 9)])


if __name__ == "__main__":
    run_panel(draw, "figure_3e", 'Error structure across models')

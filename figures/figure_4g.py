"""Render primary–relapse state coupling from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    heatmap(ax, data, 'Coupled spot mass', vmin=0, vmax=meta['vmax'])
    ax.set(xlabel='Relapse metacell', ylabel='Primary metacell')


if __name__ == "__main__":
    run_panel(draw, "figure_4g", 'Primary–relapse state coupling')

"""Render patient representation response from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data'].sort_values('rejection_cost')
    ax.plot(data.rejection_cost, data.represented_percent, 'o-', color=PALETTE[0], ms=3)
    ax.set(xlabel='Unmatched-mass cost', ylabel='Represented patient mass (%)')


if __name__ == "__main__":
    run_panel(draw, "figure_2f", 'Patient representation response')

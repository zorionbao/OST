"""Render matrix mixing response from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    for method, color in zip(inputs['metadata']['method_order'], [PALETTE[6], PALETTE[1], PALETTE[0]]):
        curve = data[data.method == method].sort_values('rna_fraction')
        ax.plot(curve.rna_fraction, curve.value, 'o-', color=color, label=method, ms=3)
    ax.set(xlabel='Osteogenic RNA fraction', ylabel='Matrix absolute error')
    legend(ax, loc='upper right')


if __name__ == "__main__":
    run_panel(draw, "figure_1b", 'Matrix mixing response')

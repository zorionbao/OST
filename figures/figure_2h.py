"""Render duplicate-independent solver scaling from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    for method, color, marker in [('Full sparse', PALETTE[1], 's'), ('Ancestry quotient', PALETTE[0], 'o')]:
        curve = data[data.method == method].sort_values('models')
        ax.plot(curve.models, curve.milliseconds, marker=marker, color=color, label=method, ms=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set(xlabel='Total model records', ylabel='Median solve time (ms)')
    legend(ax, loc='upper left')


if __name__ == "__main__":
    run_panel(draw, "figure_2h", 'Duplicate-independent solver scaling')

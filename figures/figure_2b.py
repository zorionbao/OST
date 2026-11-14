"""Render screen-wide depletion distributions from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    for category, color in [('Lung', PALETTE[1]), ('Culture', PALETTE[0])]:
        curve = data[data.category == category]
        ax.plot(curve.x, curve.y, color=color, label=category)
    ax.set(xlabel='Gene depletion log2 fold change', ylabel='Cumulative gene fraction')
    legend(ax)


if __name__ == "__main__":
    run_panel(draw, "figure_2b", 'Screen-wide depletion distributions')

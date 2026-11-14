"""Render culture geometry in dependency space from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    for category, color, marker in [('Plastic', PALETTE[0], 'o'), ('Dome', PALETTE[1], '^')]:
        points = data[data['format'] == category]
        scatter(ax, points.x, points.y, color, category, marker=marker, size=25)
    ax.set(xlabel='Dependency PC1', ylabel='Dependency PC2')
    legend(ax)


if __name__ == "__main__":
    run_panel(draw, "figure_3b", 'Culture geometry in dependency space')

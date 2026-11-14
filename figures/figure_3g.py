"""Render pathway-specific context shifts from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data'].copy()
    data.index = ['PI3K–MTOR Signaling' if str(label) == 'PI3K MTOR Signaling' else str(label) for label in data.index]
    limit = max(.12, float(data.abs().max().max()))
    heatmap(ax, data, 'Mean rank shift', cmap=DIVERGING, vmin=-limit, vmax=limit)
    ax.set_yticks(range(len(data)), data.index)


if __name__ == "__main__":
    run_panel(draw, "figure_3g", 'Pathway-specific context shifts')

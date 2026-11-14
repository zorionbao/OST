"""Render osteosarcoma-linked context atlas from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    image = scatter(ax, data.x, data.y, data.color, size=14, cmap=DIVERGING,
                    vmin=-meta['color_limit'], vmax=meta['color_limit'])
    annotate_genes(ax, data, {gene: (4, 4) for gene in ['ITGB1', 'PTK2', 'GPX4', 'CDK6']})
    ax.set(xlabel='Osteosarcoma mean Chronos', ylabel='Organoid − adherent Chronos')
    colorbar(ax, image, 'Paired format shift')


if __name__ == "__main__":
    run_panel(draw, "figure_3h", 'Osteosarcoma-linked context atlas')

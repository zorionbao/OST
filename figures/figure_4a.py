"""Render spatially annotated tumour clusters from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    for category, color in zip(meta['category_order'], PALETTE[:5]):
        points = data[data.category == category]
        scatter(ax, points.x, points.y, color, category, size=3, alpha=.8, rasterized=True)
    ax.set_aspect('equal')
    ax.set_axis_off()
    legend(ax, loc='upper left', bbox_to_anchor=(1.05, 1.02))


if __name__ == "__main__":
    run_panel(draw, "figure_4a", 'Spatially annotated tumour clusters')

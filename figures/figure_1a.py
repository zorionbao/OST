"""Render measured compartment programmes from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    x = {name: i for i, name in enumerate(meta['programmes'])}
    y = {name: i for i, name in enumerate(meta['compartments'])}
    image = ax.scatter(data.programme.map(x), data.compartment.map(y), s=data.area,
                       c=data.value, cmap=SEQUENTIAL, vmin=meta['vmin'], vmax=meta['vmax'],
                       linewidths=.3, edgecolors='white')
    ax.set_xticks(range(len(x)), list(x), rotation=25, ha='right')
    ax.set_yticks(range(len(y)), list(y))
    ax.invert_yaxis()
    colorbar(ax, image, 'Mean programme rank')


if __name__ == "__main__":
    run_panel(draw, "figure_1a", 'Measured compartment programmes')

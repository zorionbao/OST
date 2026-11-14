"""Render family mass under duplication from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data'].sort_values('copies')
    ax.plot(data.copies, data.duplicated_family_mass, 'o-', color=PALETTE[0], ms=3, label='Ancestry quotient')
    ax.plot(data.copies, data.kernel_family_mass, 's-', color=PALETTE[1], ms=3, label='Row kernel')
    ax.set_xscale('log')
    ax.set(xlabel='Identical records in one family', ylabel='Duplicated-family mass')
    legend(ax, loc='upper left')


if __name__ == "__main__":
    run_panel(draw, "figure_2g", 'Family mass under duplication')

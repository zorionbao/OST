"""Render exact ancestry-quotient allocation from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    heatmap(ax, data, 'Fraction of patient mass', vmin=0, vmax=1, labels=False)
    ax.set_xticks(range(len(data.columns)), meta['family_display_labels'], rotation=0, ha='center')
    ax.set_ylabel(f'{len(data)} represented tumours')


if __name__ == "__main__":
    run_panel(draw, "figure_2e", 'Exact ancestry-quotient allocation')

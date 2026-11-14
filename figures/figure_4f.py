"""Render section-level perturbation contrasts from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    colors = {'Tgfbr2': PALETTE[1], 'Jak2': PALETTE[2], 'Ifngr2': PALETTE[3]}
    for i, perturbation in enumerate(meta['perturbation_order']):
        points = data[data.perturbation == perturbation]
        scatter(ax, np.full(len(points), i), points.difference, colors[perturbation], perturbation, size=25)
    ax.set_xticks(range(len(meta['perturbation_order'])), meta['perturbation_order'])
    ax.set_ylabel('TGFB1 score: cluster − KP')


if __name__ == "__main__":
    run_panel(draw, "figure_4f", 'Section-level perturbation contrasts')

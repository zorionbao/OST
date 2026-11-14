"""Render independent donor uncertainty from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data'].set_index('programme')
    order = inputs['metadata']['programme_order']
    for i, programme in enumerate(order):
        row = data.loc[programme]
        ax.errorbar(row.relative_reduction, i,
                    xerr=[[row.relative_reduction - row.ci_low], [row.ci_high - row.relative_reduction]],
                    fmt='o', color=PROGRAMME_COLORS[programme], capsize=2, ms=4)
    ax.set_yticks(range(3), order)
    ax.set_xlabel('Donor-mean RMSE reduction (%)')


if __name__ == "__main__":
    run_panel(draw, "figure_1g", 'Independent donor uncertainty')

"""Render mixing and state-retention balance from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    for method, color in zip(meta['method_order'], [PALETTE[6], PALETTE[1], PALETTE[0]]):
        values = data[data.method == method].set_index('metric').loc[meta['metric_order'], 'value']
        ax.plot(range(3), values, 'o-', color=color, label=method, ms=3)
    ax.set_xticks(range(3), ['Mixture\nRMSE', 'Pure\nRMSE', 'Score\nSD'])
    ax.set(ylabel='Within-metric relative value', ylim=(-.05, 1.85))
    legend(ax, loc='upper left')


if __name__ == "__main__":
    run_panel(draw, "figure_1f", 'Mixing and state-retention balance')

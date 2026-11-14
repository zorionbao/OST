"""Render whole-model prediction performance from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, meta = inputs['data'], inputs['metadata']
    values = data.pivot(index='model', columns='method', values='rmse').loc[meta['model_order'], meta['method_order']]
    for model, row in values.iterrows():
        ax.plot(range(len(row)), row, color=PALETTE[0], marker='o', ms=2.8, alpha=.45)
    for method, color in [('PCHT', PALETTE[1]), ('Global anchor', PALETTE[2])]:
        scatter(ax, np.full(len(values), meta['method_order'].index(method)), values[method], color, method, size=14, alpha=1, zorder=3)
    ax.set_xticks(range(8), [f'M{i}' for i in range(1, 9)])
    ax.set_ylabel('Held-out model RMSE')
    ax.set_ylim(values.to_numpy().min() - .015, values.to_numpy().max() + .055)
    legend(ax, loc='upper left', ncol=2)


if __name__ == "__main__":
    run_panel(draw, "figure_3d", 'Whole-model prediction performance')

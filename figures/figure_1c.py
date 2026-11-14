"""Render independent donor improvement from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data, density = inputs['data'], inputs['density']
    summaries = inputs['quantiles'].set_index('programme')
    for i, programme in enumerate(inputs['metadata']['programme_order']):
        color = PROGRAMME_COLORS[programme]
        curve = density[density.programme == programme]
        if len(curve):
            ax.fill_between(curve.x, i + .10, i + .10 + curve.height, color=color, alpha=.30, lw=0)
            ax.plot(curve.x, i + .10 + curve.height, color=color, lw=.8)
        points = data[data.programme == programme]
        scatter(ax, points.value, i - .14 + np.linspace(-.06, .06, len(points)), color, size=20)
        row = summaries.loc[programme]
        ax.plot([row.q1, row.q3], [i - .29] * 2, color=color, lw=1)
        scatter(ax, [row['median']], [i - .29], color, size=23)
    ax.set_yticks(range(3), inputs['metadata']['programme_order'])
    ax.set(xlabel='Donor RMSE reduction (%)', ylim=(-.55, 2.65))


if __name__ == "__main__":
    run_panel(draw, "figure_1c", 'Independent donor improvement')

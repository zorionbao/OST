"""Render preserved programme identity from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data']
    order = inputs['metadata']['programme_order']
    statistics = inputs['boxes'].set_index('programme')
    box_input = [{'label': programme, **statistics.loc[programme].to_dict(), 'fliers': []} for programme in order]
    boxes = ax.bxp(box_input, patch_artist=True, widths=.42, showfliers=False,
                   medianprops={'color': 'black', 'linewidth': 1},
                   whiskerprops={'linewidth': .8}, capprops={'linewidth': .8})
    for i, (programme, patch) in enumerate(zip(order, boxes['boxes'])):
        color = PROGRAMME_COLORS[programme]
        patch.set_facecolor(color)
        patch.set_alpha(.22)
        points = data[data.programme == programme]
        scatter(ax, i + 1 + np.linspace(-.13, .13, len(points)), points.value, color, size=15)
    ax.set_ylabel('Pure-profile signed error')
    ax.tick_params(axis='x', rotation=25)


if __name__ == "__main__":
    run_panel(draw, "figure_1e", 'Preserved programme identity')

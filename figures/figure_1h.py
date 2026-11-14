"""Render discovery parameter surface from prepared figure inputs."""
from osteosarcoma_core.plotting import (PALETTE, PROGRAMME_COLORS, SEQUENTIAL, DIVERGING, np, scatter, legend, colorbar, heatmap, spatial_map, annotate_genes, run_panel)


def draw(ax, inputs):
    data = inputs['data'].pivot(index='lambda', columns='x', values='value')
    x, y = np.meshgrid(data.columns.to_numpy(dtype=float), data.index.to_numpy(dtype=float))
    image = ax.contourf(x, y, data.to_numpy(), levels=inputs['metadata']['contour_levels'], cmap=SEQUENTIAL)
    ax.scatter(x, y, s=8, color='black', linewidths=0)
    ax.set(xlabel='log10 coefficient penalty', ylabel='Drift penalty')
    colorbar(ax, image, 'MSE (1e-3)')


if __name__ == "__main__":
    run_panel(draw, "figure_1h", 'Discovery parameter surface')

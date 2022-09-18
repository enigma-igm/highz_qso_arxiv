import numpy as np

from IPython import embed

def plot_contour(x, y, ax, levels=None, range=None, bins=200, color='grey', zorder=1, label=None):
    """Plot a density contour of the points x, y
    modified from https://github.com/dfm/corner.py

    Parameters
    ----------
    x, y : array_like
        The data to plot
    """
    if range is None:
        range = [[x.min(), x.max()], [y.min(), y.max()]]
    
    if levels is None:
        levels = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)

    from matplotlib.colors import LinearSegmentedColormap, colorConverter
    base_color = ax.get_facecolor()
    base_cmap = LinearSegmentedColormap.from_list(
        "base_cmap", [base_color, base_color], N=2
    )
    
    H, X, Y = np.histogram2d(
        x.flatten(),
        y.flatten(),
        bins=bins,
        range=list(map(np.sort, range))
    )

    # Compute the density levels.
    Hflat = H.flatten()
    inds = np.argsort(Hflat)[::-1]
    Hflat = Hflat[inds]
    sm = np.cumsum(Hflat)
    sm = sm / sm[-1]
    V = np.empty(len(levels))
    for i, v0 in enumerate(levels):
        try:
            V[i] = Hflat[sm <= v0][-1]
        except IndexError:
            V[i] = Hflat[0]
    V.sort()
    m = np.diff(V) == 0

    while np.any(m):
        V[np.where(m)[0][0]] *= 1.0 - 1e-4
        m = np.diff(V) == 0
    V.sort()

    # Compute the bin centers.
    X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])

    # Extend the array for the sake of the contours at the plot edges.
    H2 = H.min() + np.zeros((H.shape[0] + 4, H.shape[1] + 4))
    H2[2:-2, 2:-2] = H
    H2[2:-2, 1] = H[:, 0]
    H2[2:-2, -2] = H[:, -1]
    H2[1, 2:-2] = H[0]
    H2[-2, 2:-2] = H[-1]
    H2[1, 1] = H[0, 0]
    H2[1, -2] = H[0, -1]
    H2[-2, 1] = H[-1, 0]
    H2[-2, -2] = H[-1, -1]
    X2 = np.concatenate(
        [
            X1[0] + np.array([-2, -1]) * np.diff(X1[:2]),
            X1,
            X1[-1] + np.array([1, 2]) * np.diff(X1[-2:]),
        ]
    )
    Y2 = np.concatenate(
        [
            Y1[0] + np.array([-2, -1]) * np.diff(Y1[:2]),
            Y1,
            Y1[-1] + np.array([1, 2]) * np.diff(Y1[-2:]),
        ]
    )

    ax.scatter(x, y, alpha=0.1, zorder=zorder, s=0.1, color=color)
    ax.contourf(
        X2,
        Y2,
        H2.T,
        [V.min(), H.max()],
        cmap=base_cmap,
        antialiased=False, zorder=zorder
    )
    ax.contour(X2, Y2, H2.T, V, colors=color, alpha=0.5, zorder=zorder, label=label)

    return ax

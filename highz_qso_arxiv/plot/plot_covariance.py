import numpy as np
from matplotlib.patches import Ellipse

from IPython import embed

def error_cov(ax, xc, yc, cov, sigma=1, color='b', **kwargs):
    w, v = np.linalg.eigh(cov) # assumes symmetric matrix
    theta = np.degrees(np.arctan2(w[0]-cov[0,0],cov[0,1]))
    ellipse = Ellipse(xy=(xc,yc),width=2.*sigma*np.sqrt(w[0]),height=2.*sigma*np.sqrt(w[1]),angle=theta,color=color, zorder=2, **kwargs)
    ellipse.set_facecolor('none')
    ax.add_artist(ellipse)
    return ax

def plot_cov(x, y, cov, ax=None, sigma=1, color='b', **kwargs):
    for i in range(len(x)):
        ax = error_cov(ax, x[i], y[i], cov[i], sigma=sigma, color=color, **kwargs)
    return ax
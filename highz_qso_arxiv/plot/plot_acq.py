import matplotlib.pyplot as plt
from photutils.aperture import CircularAperture

import matplotlib as mpl
mpl.rcParams['image.interpolation'] = 'nearest'
mpl.rcParams['image.origin'] = 'lower'

def plot_acq(image, peak_pos=None, fig=None, ax=None, display=False, save_file=""):
    if ax == None:
        fig, ax = plt.subplots(figsize=(6,6))
    im = ax.imshow(image)
    fig.colorbar(im, ax=ax)

    if peak_pos is not None:
        apertures = CircularAperture(peak_pos, r=5.)
        apertures.plot(color='#0547f9', axes=ax, lw=1.5)
    ax.set_xlabel("x [pix]", fontsize=15)
    ax.set_ylabel("y [pix]", fontsize=15)
    if display:
        plt.show()
    if save_file:
        fig.savefig(save_file)
    return fig, ax

def plot_hist(image, bins=20, fig=None, ax=None, display=False, save_file=""):
    """plot 1d hist of an image
       help us understand the background

    Args:
        image (_type_): _description_
        bins (int, optional): _description_. Defaults to 20.
    """
    if ax == None:
        fig, ax = plt.subplots(figsize=(8,6))
    ax.hist(image.flatten(), bins=bins)
    ax.set_xlabel("counts", fontsize=15)
    if display:
        plt.show()
    if save_file:
        fig.savefig(save_file)
    return fig, ax

def plot_acq_and_hist(image, peak_pos=None, bins=20, title="", display=False, save_file=""):
    fig, axs = plt.subplots(1, 2, figsize=(20,6))
    plot_acq(image, peak_pos, fig, axs[0])
    plot_hist(image, bins, fig, axs[1])
    if title:
        plt.suptitle(title, fontsize=15)
    if display:
        plt.show()
    if save_file:
        fig.savefig(save_file)
    return fig, axs
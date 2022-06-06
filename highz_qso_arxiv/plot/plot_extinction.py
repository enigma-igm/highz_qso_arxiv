from ..util import unix_to_iso
import matplotlib.pyplot as plt

def plot_extinction(extinction_data, star_unix=None, qso_unix=None, offset=0, display=False, save_file=""):
    """plot extinction curve

    Args:
        extinction_data (_type_): _description_
        star_unix (int): star frame observing time, in UNIX
        qso_unix (int): qso frame observing time, in UNIX
        offset (float): constant offset, CFHT staff usually assume an average offset of +0.03
        display (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """
    fig, ax = plt.subplots(figsize=(6,12))
    mask = extinction_data["extinction"] > 0
    ax.scatter(extinction_data["extinction"][mask]-offset, extinction_data["time"][mask], color="black")
    ax.errorbar(extinction_data["extinction"][mask]-offset, extinction_data["time"][mask], 
                xerr=extinction_data["scatter"][mask], linestyle='None', 
                ecolor="black", elinewidth=1, capsize=1, capthick=1)
    ylabel = ax.get_yticks()
    ylabel = [unix_to_iso(y) for y in ylabel]
    ax.set_yticklabels(ylabel, rotation=0)

    xmin, xmax = ax.get_xlim()
    if star_unix is not None:
        ax.hlines(star_unix, xmin, xmax, ls="dashed", label="star", color="red")
    if qso_unix is not None:
        ax.hlines(qso_unix, xmin, xmax, ls="dashed", label="qso", color="blue")
    ax.set_xlim(xmin, xmax)
    ax.legend()
    ax.set_xlabel("extinction [mag]", fontsize=15)
    if display:
        plt.show()
    if save_file:
        fig.savefig(save_file)
    return fig, ax
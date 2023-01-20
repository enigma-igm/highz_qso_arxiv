from re import template
from highz_qso_arxiv.plot import plot_series

name_all = []
fits_all = []

# NIRES-1903

name_list = ["J0756+5744", "J0800+3034",
            "J0803+0030", "J0816+1622",
            "J0953-0853", "J1217+3136", "J1312+5707", "J1423+2901",
            "J1546+2945", "J1550+2558",
            "J1635+5940", "J1638+5412"]

fits_list = [f"../NIRES_1903/reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]

name_all += name_list
fits_all += fits_list

# NIRES-1905

name_list = ["J1432+5746",
            "J1450+3302", "J1458+6300",
            "J1635+1758", "J2212+2040"]

fits_list = [f"../NIRES_1905/reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]

import pandas as pd
cat = pd.read_csv('../hizqa_catalog.csv')
mask_dwarf = cat['label'] == 'star'
mask_QSO = cat['label'] == 'QSO'
mask_uNQ = cat['label'] == 'uNQ'
mask_inconclusive = cat['label'] == 'inconclusive'

idx_dwarf, idx_QSO, idx_uNQ, idx_inconclusive = [], [], [], []
for nm in cat['name'][mask_dwarf]:
    if nm in name_all:
        idx_dwarf.append(name_all.index(nm))
for nm in cat['name'][mask_QSO]:
    if nm in name_all:
        idx_QSO.append(name_all.index(nm))
for nm in cat['name'][mask_uNQ]:
    if nm in name_all:
        idx_uNQ.append(name_all.index(nm))
for nm in cat['name'][mask_inconclusive]:
    if nm in name_all:
        idx_inconclusive.append(name_all.index(nm))

xrange = (9600, 24400)
if len(idx_dwarf) > 0:
    name_dwarf = [name_all[i] for i in idx_dwarf]
    fits_dwarf = [fits_all[i] for i in idx_dwarf]
    idx_all = [1 for i in range(len(name_dwarf))]
    template_all = [True for i in range(len(name_dwarf))]
    qso_all = [False for i in range(len(name_dwarf))]
    fig, axs= plot_series(name_dwarf, fits_dwarf, idx_all, telluric=False, template_list=template_all,
                          xrange=xrange, smooth_window=9, qso_list=qso_all, display=False)
    for ax in axs:
        # mask out the telluric region
        ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
        ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
    fig.savefig('nires_dwarf.pdf')

    for seq, i in enumerate(range(0, len(name_dwarf), 10)):
        fig, axs = plot_series(name_dwarf[i:i+10], fits_dwarf[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, smooth_window=9)
        for ax in axs:
            # mask out the telluric region
            ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
            ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
        fig.savefig(f'nires_dwarf_{seq}.pdf')
    fits_dwarf = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_dwarf]
    fig, axs = plot_series(name_dwarf, fits_dwarf, idx_all, telluric=True, template_list=template_all, 
                           xrange=xrange, smooth_window=9, qso_list=qso_all, display=False)
    for ax in axs:
        # mask out the telluric region
        ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
        ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
    fig.savefig('nires_dwarf_tellcorr.pdf')
    
    for seq, i in enumerate(range(0, len(name_dwarf), 10)):
        fig, axs = plot_series(name_dwarf[i:i+10], fits_dwarf[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, smooth_window=9)
        for ax in axs[:-1]:
            # mask out the telluric region
            ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
            ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
            # ymax * 1.5
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.5)
        fig.savefig(f'nires_dwarf_tellcorr_{seq}.pdf')

if len(idx_QSO) > 0:
    name_QSO = [name_all[i] for i in idx_QSO]
    fits_QSO = [fits_all[i] for i in idx_QSO]
    idx_all = [1 for i in range(len(name_QSO))]
    template_all = [False for i in range(len(name_QSO))]
    qso_all = [True for i in range(len(name_QSO))]
    fig, axs = plot_series(name_QSO, fits_QSO, idx_all, telluric=False, template_list=template_all, 
                            xrange=xrange, smooth_window=9, qso_list=qso_all, display=False)
    for ax in axs:
        # mask out the telluric region
        ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
        ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
    fig.savefig('nires_QSO.pdf')
    for seq, i in enumerate(range(0, len(name_QSO), 10)):
        fig, axs = plot_series(name_QSO[i:i+10], fits_QSO[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, smooth_window=9)
        for ax in axs:
            # mask out the telluric region
            ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
            ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
        fig.savefig(f'nires_QSO_{seq}.pdf')
    fits_QSO = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_QSO]
    fig, axs = plot_series(name_QSO, fits_QSO, idx_all, telluric=True, template_list=template_all, 
                            xrange=xrange, qso_list=qso_all, display=False)
    for ax in axs:
        # mask out the telluric region
        ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
        ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
    fig.savefig('nires_QSO_tellcorr.pdf')

    for seq, i in enumerate(range(0, len(name_QSO), 10)):
        fig, axs= plot_series(name_QSO[i:i+10], fits_QSO[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, smooth_window=9)
        for ax in axs[:-1]:
            # mask out the telluric region
            ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
            ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
            # ymax * 1.5
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.5)
        fig.savefig(f'nires_QSO_tellcorr_{seq}.pdf')

if len(idx_uNQ) > 0:
    name_uNQ = [name_all[i] for i in idx_uNQ]
    fits_uNQ = [fits_all[i] for i in idx_uNQ]
    idx_all = [1 for i in range(len(name_uNQ))]
    template_all = [False for i in range(len(name_uNQ))]
    qso_all = [True for i in range(len(name_uNQ))]
    fig, axs = plot_series(name_uNQ, fits_uNQ, idx_all, telluric=False, template_list=template_all, 
                            xrange=xrange, smooth_window=9, qso_list=qso_all, display=False)
    for ax in axs:
        # mask out the telluric region
        ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
        ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
    fig.savefig('nires_uNQ.pdf')

    for seq, i in enumerate(range(0, len(name_uNQ), 10)):
        fig, axs = plot_series(name_uNQ[i:i+10], fits_uNQ[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, smooth_window=9)
        for ax in axs:
            # mask out the telluric region
            ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
            ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
        fig.savefig(f'nires_uNQ_{seq}.pdf')

    fits_uNQ = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_uNQ]
    telluric_all = [True for i in range(len(name_uNQ))]
    fig, axs = plot_series(name_uNQ, fits_uNQ, idx_all, telluric=True, template_list=template_all, 
                            xrange=xrange, smooth_window=9, qso_list=qso_all, display=False)
    for ax in axs:
        # mask out the telluric region
        ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
        ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
    fig.savefig('nires_uNQ_tellcorr.pdf')

    for seq, i in enumerate(range(0, len(name_uNQ), 10)):
        fig, axs= plot_series(name_uNQ[i:i+10], fits_uNQ[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, smooth_window=9)
        for ax in axs[:-1]:
            # mask out the telluric region
            ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=10)
            ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=10)
            # ymax * 1.5
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.5)
        fig.savefig(f'nires_uNQ_tellcorr_{seq}.pdf')    

# if len(idx_inconclusive) > 0:
#     name_inconclusive = [name_all[i] for i in idx_inconclusive]
#     fits_inconclusive = [fits_all[i] for i in idx_inconclusive]
#     idx_all = [1 for i in range(len(name_inconclusive))]
#     telluric_all = [False for i in range(len(name_inconclusive))]
#     template_all = [False for i in range(len(name_inconclusive))]
#     qso_all = [True for i in range(len(name_inconclusive))]
#     plot_series(name_inconclusive, fits_inconclusive, idx_all, telluric_list=telluric_all, template_list=template_all, 
#                 smooth_window=9, qso_list=qso_all, display=False, save_file='nires_inconclusive.pdf')
#     for seq, i in enumerate(range(0, len(name_inconclusive), 10)):
#         plot_series(name_inconclusive[i:i+10], fits_inconclusive[i:i+10], idx_all[i:i+10], telluric_list=telluric_all[i:i+10], 
#                     template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
#                     smooth_window=9, save_file=f'nires_inconclusive_{seq}.pdf')
#     fits_inconclusive = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_inconclusive]
#     telluric_all = [True for i in range(len(name_inconclusive))]
#     plot_series(name_inconclusive, fits_inconclusive, idx_all, telluric_list=telluric_all, template_list=template_all, 
#                 smooth_window=9, qso_list=qso_all, display=False, save_file='nires_inconclusive_tellcorr.pdf')
#     for seq, i in enumerate(range(0, len(name_inconclusive), 10)):
#         plot_series(name_inconclusive[i:i+10], fits_inconclusive[i:i+10], idx_all[i:i+10], telluric_list=telluric_all[i:i+10], 
#                     template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
#                     smooth_window=9, save_file=f'nires_inconclusive_tellcorr_{seq}.pdf')
from re import template
from highz_qso_arxiv.plot import plot_series

name_all = []
fits_all = []

# MOSFIRE-1911
name_list = ["J0058+1715", "J0208+1749", "J0426+0221",
             "J2132-1434", "J2258-1540", "J2312+3215"]

fits_list = [f"../MOSFIRE_1911/reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]

name_all += name_list
fits_all += fits_list

# MOSFIRE-2005
name_list = ["J1311+5932", "J1600+3127", "J1626+1337"]

fits_list = [f"../MOSFIRE_2005/reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]
fits_list[name_list.index("J1626+1337")] = "../MOSFIRE_2005/reduced/all_part2/coadd2d/J1626+1337_coadd.fits"

name_all += name_list
fits_all += fits_list

# MOSFIRE-2010
name_list = ["J0020+1956", "J0102-3327",
            "J0159-2017", "J0220+3458",
            "J0251+1724", "J0339-1610", "J0705+2850",
            "J0758+0040", "J0759+0735",
            "J0812+3201", "J0813+2219",
            "J0822+5045", "J0826+1217",
            "J0845+3033", "J2033-0428",
            "J2105-0049", "J2119-1335",
            "J2129+0601", "J2212+2040",
            "J2244+2907", "J2247+3247",
            "J2312+3215", "J2328+2755",
            "J2346+2545", "J2352+2655",
            "J2355+2048",]

fits_list = [f"../MOSFIRE_2010/reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]
fits_list[name_list.index('J0020+1956')] = f"../MOSFIRE_2010/reduced/all_day4/coadd2d/J0020+1956_coadd.fits"
fits_list[name_list.index('J0339-1610')] = f"../MOSFIRE_2010/reduced/J0339/coadd2d/J0339-1610_coadd.fits"
fits_list[name_list.index('J0759+0735')] = f"../MOSFIRE_2010/reduced/all_day4/coadd2d/J0759+0735_coadd.fits"
fits_list[name_list.index('J2328+2755')] = f"../MOSFIRE_2010/reduced/all_day4/coadd2d/J2328+2755_coadd.fits"
fits_list[name_list.index('J2352+2655')] = f"../MOSFIRE_2010/reduced/all_day4/coadd2d/J2352+2655_coadd.fits"

name_all += name_list
fits_all += fits_list

# MOSFIRE-2201

name_list = ["J0637+3812", "J0838-0112",
            "J0847+0139", "J0854+2908",
            "J0938+1341"]
fits_list = [f"../MOSFIRE_2201/reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]

name_all += name_list
fits_all += fits_list

# MOSFIRE-2204

name_list = ["J1150-0143",
            "J1202+0129",
            "J1318+2932",
            "J1332+0150",
            "J1333+0919",
            "J1511+0344"]
fits_list = [f"../MOSFIRE_2204/reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]

name_all += name_list
fits_all += fits_list

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

xrange = (9650, 11200)

if len(idx_dwarf) > 0:

    name_dwarf = [name_all[i] for i in idx_dwarf]
    fits_dwarf = [fits_all[i] for i in idx_dwarf]
    idx_all = [1 for i in range(len(name_dwarf))]
    template_all = [True for i in range(len(name_dwarf))]
    qso_all = [False for i in range(len(name_dwarf))]
    plot_series(name_dwarf, fits_dwarf, idx_all, telluric=False, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='mosfire_dwarf.pdf')
    for seq, i in enumerate(range(0, len(name_dwarf), 10)):
        plot_series(name_dwarf[i:i+10], fits_dwarf[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'mosfire_dwarf_{seq}.pdf')
    fits_dwarf = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_dwarf]
    plot_series(name_dwarf, fits_dwarf, idx_all, telluric=True, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='mosfire_dwarf_tellcorr.pdf')
    for seq, i in enumerate(range(0, len(name_dwarf), 10)):
        plot_series(name_dwarf[i:i+10], fits_dwarf[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'mosfire_dwarf_tellcorr_{seq}.pdf')

if len(idx_QSO) > 0:
    name_QSO = [name_all[i] for i in idx_QSO]
    fits_QSO = [fits_all[i] for i in idx_QSO]
    idx_all = [1 for i in range(len(name_QSO))]
    template_all = [False for i in range(len(name_QSO))]
    qso_all = [True for i in range(len(name_QSO))]
    plot_series(name_QSO, fits_QSO, idx_all, telluric=False, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='mosfire_QSO.pdf')
    for seq, i in enumerate(range(0, len(name_QSO), 10)):
        plot_series(name_QSO[i:i+10], fits_QSO[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'mosfire_QSO_{seq}.pdf')
    fits_QSO = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_QSO]
    plot_series(name_QSO, fits_QSO, idx_all, telluric=True, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='mosfire_QSO_tellcorr.pdf')
    for seq, i in enumerate(range(0, len(name_QSO), 10)):
        plot_series(name_QSO[i:i+10], fits_QSO[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'mosfire_QSO_tellcorr_{seq}.pdf')

if len(idx_uNQ) > 0:
    name_uNQ = [name_all[i] for i in idx_uNQ]
    fits_uNQ = [fits_all[i] for i in idx_uNQ]
    idx_all = [1 for i in range(len(name_uNQ))]
    template_all = [False for i in range(len(name_uNQ))]
    qso_all = [True for i in range(len(name_uNQ))]
    plot_series(name_uNQ, fits_uNQ, idx_all, telluric=False, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='mosfire_uNQ.pdf')
    for seq, i in enumerate(range(0, len(name_uNQ), 10)):
        plot_series(name_uNQ[i:i+10], fits_uNQ[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'mosfire_uNQ_{seq}.pdf')
    fits_uNQ = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_uNQ]
    plot_series(name_uNQ, fits_uNQ, idx_all, telluric=True, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='mosfire_uNQ_tellcorr.pdf')
    for seq, i in enumerate(range(0, len(name_uNQ), 10)):
        plot_series(name_uNQ[i:i+10], fits_uNQ[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'mosfire_uNQ_tellcorr_{seq}.pdf')

if len(idx_inconclusive) > 0:
    name_inconclusive = [name_all[i] for i in idx_inconclusive]
    fits_inconclusive = [fits_all[i] for i in idx_inconclusive]
    idx_all = [1 for i in range(len(name_inconclusive))]
    template_all = [False for i in range(len(name_inconclusive))]
    qso_all = [True for i in range(len(name_inconclusive))]
    plot_series(name_inconclusive, fits_inconclusive, idx_all, telluric=False, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='mosfire_inconclusive.pdf')
    for seq, i in enumerate(range(0, len(name_inconclusive), 10)):
        plot_series(name_inconclusive[i:i+10], fits_inconclusive[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'mosfire_inconclusive_{seq}.pdf')
    fits_inconclusive = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_inconclusive]
    plot_series(name_inconclusive, fits_inconclusive, idx_all, telluric=True, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='mosfire_inconclusive_tellcorr.pdf')
    for seq, i in enumerate(range(0, len(name_inconclusive), 10)):
        plot_series(name_inconclusive[i:i+10], fits_inconclusive[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'mosfire_inconclusive_tellcorr_{seq}.pdf')
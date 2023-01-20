from re import template
from highz_qso_arxiv.plot import plot_series

name_all = []
fits_all = []

# LRIS-2201
name_list = ['J0739+2328', 'J0759+2811', 'J0831+2558',
             'J0849+2601', 'J0936+3244', 'J0936+3346', 'J0938+3332',
             'J1003-2610', 'J1030-0031', 'J1145+0933', 'J1146+0039',
             'J1206-0051', 'J1209+0135', 'J1241-0134', 'J1246-3045',
             'J1256-0306', 'J1305-1549', 'J1326+0927']
fits_list = [f"../LRIS_2201/reduced/all/{nm}/{nm}_coadd.fits" for nm in name_list]
fits_list[9] = f"../LRIS_2201/reduced/all/coadd2d/{name_list[9]}_coadd.fits"
fits_list[11] = f"../LRIS_2201/reduced/all/coadd2d/{name_list[11]}_coadd.fits"
fits_list[15] = f"../LRIS_2201/reduced/all/coadd2d/{name_list[15]}_coadd.fits"

name_all += name_list
fits_all += fits_list

# LRIS-2203
targets_0305 = ["J1100+0203", "J1143-0248", "J1200+0112", "J1223+0114",
                "J1250+0347", "J1319+0101", "J1327-0207", "J1335+0103",
                "J1355-0044", "J1355-0111", "J1356+0050", "J1358-0118",
                "J1422+0307"]
targets_0306 = ["J0810+2352", "J0850+0146", "J0901+2906", "J0911+0022",
                "J0947+0111", "J1153-2239", "J1233-2807", "J1412-2757",
                "J1425+0004", "J1427+0202", "J1458+1012", "J1510-0144", 
                "J1514-0121", "J1522+0051", "J1534+0046", "J1655-0051", "J1724+3718"]
name_list = targets_0305 + targets_0306
name_list.sort()

fits_list = []
for nm in name_list:
    if nm in targets_0305:
        fits_list.append(f"../LRIS_2203//LRIS_220305/reduced/all/{nm}/{nm}_coadd.fits")
    else:
        fits_list.append(f"../LRIS_2203//LRIS_220306/reduced/all/{nm}/{nm}_coadd.fits")
fits_list[name_list.index('J1100+0203')] = "../LRIS_2203//LRIS_220305/reduced/all/coadd2d/J1100+0203_coadd.fits"
fits_list[name_list.index('J1200+0112')] = "../LRIS_2203/LRIS_220305/reduced/all/coadd2d/J1200+0112_coadd.fits"
fits_list[name_list.index('J1427+0202')] = "../LRIS_2203/LRIS_220306/reduced/all/coadd2d/J1427+0202_coadd.fits"

name_all += name_list
fits_all += fits_list

# LRIS-2204
name_list = ["J0854+2908", "J0954-0117", "J1100+0237", "J1111+0640",
             "J1238+0234", "J1253-0200", "J1317-0302",
             "J1340+1422", "J1349+0805", "J1351+0128", "J1352+0002",
             "J1401+4542", "J1426+0148", "J1429+0219", "J1434+0240", "J1438-0103", 
             "J1441+0149", "J1450+0228", "J1508-0109", "J1513-0059", 
             "J1523+2935", "J1532-0155", "J1535+6146", 
             "J1625+2414", "J1635+5940", "J2212+2040"]
fits_list = [f"../LRIS_2204/reduced/all/{nm}/{nm}_coadd.fits" for nm in name_list]
fits_list[0] = "../LRIS_2204/reduced/all/coadd2d/J0854+2908_coadd.fits"
fits_list[2] = "../LRIS_2204/reduced/J1100/coadd2d/J1100+0237_coadd.fits"
fits_list[-1] = "../LRIS_2204/reduced/J2212/J2212+2040_coadd.fits"

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

xrange = (7300, 10600)

if len(idx_dwarf) > 0:
    name_dwarf = [name_all[i] for i in idx_dwarf]
    fits_dwarf = [fits_all[i] for i in idx_dwarf]
    idx_all = [1 for i in range(len(name_dwarf))]
    telluric_all = [False for i in range(len(name_dwarf))]
    template_all = [True for i in range(len(name_dwarf))]
    qso_all = [False for i in range(len(name_dwarf))]
    plot_series(name_dwarf, fits_dwarf, idx_all, telluric=False, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='lris_dwarf.pdf')
    for seq, i in enumerate(range(0, len(name_dwarf), 10)):
        plot_series(name_dwarf[i:i+10], fits_dwarf[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'lris_dwarf_{seq}.pdf')
    fits_dwarf = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_dwarf]
    telluric_all = [True for i in range(len(name_dwarf))]
    plot_series(name_dwarf, fits_dwarf, idx_all, telluric=True, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='lris_dwarf_tellcorr.pdf')
    for seq, i in enumerate(range(0, len(name_dwarf), 10)):
        plot_series(name_dwarf[i:i+10], fits_dwarf[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'lris_dwarf_tellcorr_{seq}.pdf')

if len(idx_QSO) > 0:
    name_QSO = [name_all[i] for i in idx_QSO]
    fits_QSO = [fits_all[i] for i in idx_QSO]
    idx_all = [1 for i in range(len(name_QSO))]
    telluric_all = [False for i in range(len(name_QSO))]
    template_all = [False for i in range(len(name_QSO))]
    qso_all = [True for i in range(len(name_QSO))]
    plot_series(name_QSO, fits_QSO, idx_all, telluric=False, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='lris_QSO.pdf')
    for seq, i in enumerate(range(0, len(name_QSO), 10)):
        plot_series(name_QSO[i:i+10], fits_QSO[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'lris_QSO_{seq}.pdf')
    fits_QSO = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_QSO]
    telluric_all = [True for i in range(len(name_QSO))]
    plot_series(name_QSO, fits_QSO, idx_all, telluric=True, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='lris_QSO_tellcorr.pdf')
    for seq, i in enumerate(range(0, len(name_QSO), 10)):
        plot_series(name_QSO[i:i+10], fits_QSO[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'lris_QSO_tellcorr_{seq}.pdf')

if len(idx_uNQ) > 0:
    name_uNQ = [name_all[i] for i in idx_uNQ]
    fits_uNQ = [fits_all[i] for i in idx_uNQ]
    idx_all = [1 for i in range(len(name_uNQ))]
    telluric_all = [False for i in range(len(name_uNQ))]
    template_all = [False for i in range(len(name_uNQ))]
    qso_all = [True for i in range(len(name_uNQ))]
    plot_series(name_uNQ, fits_uNQ, idx_all, telluric=False, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='lris_uNQ.pdf')
    for seq, i in enumerate(range(0, len(name_uNQ), 10)):
        plot_series(name_uNQ[i:i+10], fits_uNQ[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'lris_uNQ_{seq}.pdf')
    fits_uNQ = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_uNQ]
    telluric_all = [True for i in range(len(name_uNQ))]
    plot_series(name_uNQ, fits_uNQ, idx_all, telluric=True, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='lris_uNQ_tellcorr.pdf')
    for seq, i in enumerate(range(0, len(name_uNQ), 10)):
        plot_series(name_uNQ[i:i+10], fits_uNQ[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'lris_uNQ_tellcorr_{seq}.pdf')

if len(idx_inconclusive) > 0:
    name_inconclusive = [name_all[i] for i in idx_inconclusive]
    fits_inconclusive = [fits_all[i] for i in idx_inconclusive]
    idx_all = [1 for i in range(len(name_inconclusive))]
    telluric_all = [False for i in range(len(name_inconclusive))]
    template_all = [False for i in range(len(name_inconclusive))]
    qso_all = [True for i in range(len(name_inconclusive))]
    plot_series(name_inconclusive, fits_inconclusive, idx_all, telluric=False, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='lris_inconclusive.pdf')
    for seq, i in enumerate(range(0, len(name_inconclusive), 10)):
        plot_series(name_inconclusive[i:i+10], fits_inconclusive[i:i+10], idx_all[i:i+10], telluric=False, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'lris_inconclusive_{seq}.pdf')
    fits_inconclusive = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_inconclusive]
    telluric_all = [True for i in range(len(name_inconclusive))]
    plot_series(name_inconclusive, fits_inconclusive, idx_all, telluric=True, template_list=template_all, 
                xrange=xrange, qso_list=qso_all, display=False, save_file='lris_inconclusive_tellcorr.pdf')
    for seq, i in enumerate(range(0, len(name_inconclusive), 10)):
        plot_series(name_inconclusive[i:i+10], fits_inconclusive[i:i+10], idx_all[i:i+10], telluric=True, 
                    template_list=template_all[i:i+10], qso_list=qso_all[i:i+10], display=False, 
                    xrange=xrange, save_file=f'lris_inconclusive_tellcorr_{seq}.pdf')
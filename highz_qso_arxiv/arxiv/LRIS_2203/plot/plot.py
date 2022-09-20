from re import template
from highz_qso_arxiv.plot import plot_series

targets_0305 = ["J1100+0203", "J1143-0248", "J1200+0112", "J1223+0114",
                "J1250+0347", "J1319+0101", "J1327-0207", "J1335+0103",
                "J1355-0044", "J1355-0111", "J1356+0050", "J1358-0118",
                "J1422+0307"]
targets_0306 = ["J0810+2352", "J0850+0146", "J0901+2906", "J0911+0022",
                "J0947+0111", "J1153-2239", "J1233-2807", "J1412-2757",
                "J1425+0004", "J1427+0202", "J1458+1012", "J1510-0144", 
                "J1514-0121", "J1522+0051", "J1534+0046", "J1655-0051", "J1724+3718"]
targets_list = targets_0305 + targets_0306
targets_list.sort()

fits_list = []
for nm in targets_list:
    if nm in targets_0305:
        fits_list.append(f"../LRIS_220305/reduced/all/{nm}/{nm}_coadd.fits")
    else:
        fits_list.append(f"../LRIS_220306/reduced/all/{nm}/{nm}_coadd.fits")
fits_list[targets_list.index('J1100+0203')] = "../LRIS_220305/reduced/all/coadd2d/J1100+0203_coadd.fits"
fits_list[targets_list.index('J1200+0112')] = "../LRIS_220305/reduced/all/coadd2d/J1200+0112_coadd.fits"
fits_list[targets_list.index('J1427+0202')] = "../LRIS_220306/reduced/all/coadd2d/J1427+0202_coadd.fits"
idx_list = [1 for i in range(len(targets_list))]
template_list = [True for i in range(len(fits_list))]
template_list[targets_list.index('J0901+2906')] = False
template_list[targets_list.index('J1319+0101')] = False
template_list[targets_list.index('J1458+1012')] = False
template_list[targets_list.index('J1724+3718')] = False

telluric_list = [False for i in range(len(fits_list))]
qso_list = [False for i in range(len(fits_list))]
qso_list[targets_list.index('J0901+2906')] = True
qso_list[targets_list.index('J1319+0101')] = True
qso_list[targets_list.index('J1458+1012')] = True
qso_list[targets_list.index('J1724+3718')] = True

# divide the list into three parts and plot

plot_series(targets_list[:10], fits_list[:10], idx_list[:10], telluric_list=telluric_list[:10], 
            template_list=template_list[:10], qso_list=qso_list[:10], display=False, save_file="LRIS_2203_part1.pdf")
plot_series(targets_list[10:20], fits_list[10:20], idx_list[10:20], telluric_list=telluric_list[10:20], 
            template_list=template_list[10:20], qso_list=qso_list[10:20], display=False, save_file="LRIS_2203_part2.pdf")
plot_series(targets_list[20:], fits_list[20:], idx_list[20:], telluric_list=telluric_list[20:], 
            template_list=template_list[20:], qso_list=qso_list[20:], display=False, save_file="LRIS_2203_part3.pdf")

plot_series(targets_list, fits_list, idx_list, smooth_window=3, 
            template_list=template_list, telluric_list=telluric_list, qso_list=qso_list, display=False, save_file="LRIS_2203.pdf")

fits_list = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_list]
telluric_list = [True for i in range(len(fits_list))]
plot_series(targets_list[:10], fits_list[:10], idx_list[:10], telluric_list=telluric_list[:10], 
            template_list=template_list[:10], qso_list=qso_list[:10], display=False, save_file="LRIS_2203_tellcorr_part1.pdf")
plot_series(targets_list[10:20], fits_list[10:20], idx_list[10:20], telluric_list=telluric_list[10:20], 
            template_list=template_list[10:20], qso_list=qso_list[10:20], display=False, save_file="LRIS_2203_tellcorr_part2.pdf")
plot_series(targets_list[20:], fits_list[20:], idx_list[20:], telluric_list=telluric_list[20:], 
            template_list=template_list[20:], qso_list=qso_list[20:], display=False, save_file="LRIS_2203_tellcorr_part3.pdf")
plot_series(targets_list, fits_list, idx_list, smooth_window=3,
            template_list=template_list, telluric_list=telluric_list, qso_list=qso_list, display=False, save_file="LRIS_2203_tellcorr.pdf")
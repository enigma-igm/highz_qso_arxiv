from re import template
from highz_qso_arxiv.plot import plot_series

targets_0305 = ["J1100+0203", "J1143-0248", "J1200+0112", "J1223+0114",
                "J1250+0347", "J1319+0101", "J1327-0207", "J1335+0103",
                "J1355-0044", "J1355-0111", "J1356+0050", "J1358-0118",
                "J1422+0307"]
fits_list = [f"../LRIS_220305/reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in targets_0305]
fits_list[0] = "../LRIS_220305/reduced/all/coadd2d/J1100+0203_coadd_tellcorr.fits"
fits_list[2] = "../LRIS_220305/reduced/all/coadd2d/J1200+0112_coadd_tellcorr.fits"
idx_list = [1 for i in range(len(targets_0305))]

# a list of True with the same length as fits_list
template_list = [True for i in range(len(fits_list))]
template_list[5] = False
telluric_list = [True for i in range(len(fits_list))]
plot_series(targets_0305, fits_list, idx_list, smooth_window=3, 
            template_list=template_list, telluric_list=telluric_list, display=False, save_file="LRIS_220305.pdf")

targets_0306 = ["J0810+2352", "J0850+0146", "J0901+2906", "J0911+0022",
                "J0947+0111", "J1153-2239", "J1233-2807", "J1412-2757",
                "J1425+0004", "J1427+0202", "J1458+1012", "J1510-0144", 
                "J1514-0121", "J1522+0051", "J1534+0046", "J1655-0051", "J1724+3718"]
fits_list = [f"../LRIS_220306/reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in targets_0306]
fits_list[9] = "../LRIS_220306/reduced/all/coadd2d/J1427+0202_coadd_tellcorr.fits"
idx_list = [1 for i in range(len(targets_0306))]

template_list = [True for i in range(len(fits_list))]
template_list[2] = False
template_list[10] = False
template_list[-1] = False
telluric_list = [True for i in range(len(fits_list))]
plot_series(targets_0306, fits_list, idx_list, smooth_window=3,
            template_list=template_list, telluric_list=telluric_list, display=False, save_file="LRIS_220306.pdf")

qso_targets = ["J1458+1012", "J1319+0101", "J1724+3718", "J0901+2906"]
fits_list = [f"../LRIS_220306/reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in qso_targets]
fits_list[1] = f"../LRIS_220305/reduced/all/J1319+0101/J1319+0101_coadd_tellcorr.fits"
idx_list = [1 for i in range(len(qso_targets))]

template_list = [False for i in range(len(fits_list))]
telluric_list = [True for i in range(len(fits_list))]
plot_series(qso_targets, fits_list, idx_list, template_list=template_list, telluric_list=telluric_list, display=False, save_file="LRIS_2203_qso.pdf")

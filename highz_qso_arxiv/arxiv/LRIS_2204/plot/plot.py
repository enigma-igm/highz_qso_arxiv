from highz_qso_arxiv.plot import plot_series

"""
    plot all without tellruic
"""
name_list = ["J0739+2645", "J0954-0117", "J1111+0640",
             "J1238+0234", "J1253-0200", "J1317-0302",
             "J1340+1422", "J1349+0805", "J1351+0128", "J1352+0002",
             "J1438-0103", "J1426+0148", "J1429+0219",
             "J1434+0240", "J1523+2935", "J1625+2414", "J1635+5940",
             "J1441+0149", "J1437-0151", "J1450+0228", "J1401+4542",
             "J1535+6146", "J1508-0109", "J1513-0059", "J1532-0155",
             "J2212+2040"]
fits_list = [f"../reduced/all/{nm}/{nm}_coadd.fits" for nm in name_list]
fits_list[-1] = "../reduced/J2212/J2212+2040_coadd.fits"
idx_list = [1 for i in range(len(name_list))]
template_list = [True for i in range(len(fits_list))]
telluric_list = [False for i in range(len(fits_list))]
plot_series(name_list, fits_list, idx_list, template_list=template_list, telluric_list=telluric_list, 
            display=False, save_file="LRIS_2204.pdf")

"""
    plot all with telluric
"""
fits_list = [f"../reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in name_list]
fits_list[-1] = "../reduced/J2212/J2212+2040_coadd_tellcorr.fits"
idx_list = [1 for i in range(len(name_list))]
template_list = [True for i in range(len(fits_list))]
telluric_list = [True for i in range(len(fits_list))]
plot_series(name_list, fits_list, idx_list, template_list=template_list, telluric_list=telluric_list, 
            display=False, save_file="LRIS_2204_tellcorr.pdf")
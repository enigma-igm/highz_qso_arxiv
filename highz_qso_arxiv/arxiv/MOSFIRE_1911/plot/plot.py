from highz_qso_arxiv.plot import plot_series

"""
    plot all without tellruic
"""
name_list = ["J0038-0653", "J0058+1715",
             "J0208+1749", "J0313-1806", "J0426+0221",
             "J2132-1434", "J2258-1540", "J2312+3215"]

fits_list = [f"../reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
template_list = [False for i in range(len(fits_list))]
telluric_list = [False for i in range(len(fits_list))]
plot_series(name_list, fits_list, idx_list, template_list=template_list, telluric_list=telluric_list, 
            display=False, save_file="MOSFIRE_1911.pdf")

"""
    plot all with telluric
"""
fits_list = [f"../reduced/all/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
template_list = [False for i in range(len(fits_list))]
telluric_list = [True for i in range(len(fits_list))]
plot_series(name_list, fits_list, idx_list, template_list=template_list, telluric_list=telluric_list, 
            display=False, save_file="MOSFIRE_1911_tellcorr.pdf")
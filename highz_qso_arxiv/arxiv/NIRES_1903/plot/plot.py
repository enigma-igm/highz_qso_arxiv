from highz_qso_arxiv.plot import plot_series

"""
    plot all without telluric
"""

name_list = ["J0739+4843", "J1217+3136",
             "J0756+5744", "J1423+2901",
             "J0800+3034", "J1546+2945",
             "J0800+5128", "J1550+2558",
             "J0803+0030", "J1635+5940",
             "J0816+1622", "J1638+5412"]
fits_list = [f"../reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
telluric_list = [False for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, display=False, save_file="NIRES_1903.pdf")

"""
    plot all with telluric
"""

name_list = ["J0739+4843", "J1217+3136",
             "J0756+5744", "J1423+2901",
             "J0800+3034", "J1546+2945",
             "J0800+5128", "J1550+2558",
             "J0803+0030", "J1635+5940",
             "J0816+1622", "J1638+5412"]
fits_list = [f"../reduced/all/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
telluric_list = [True for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, display=False, save_file="NIRES_1903_tellcorr.pdf")

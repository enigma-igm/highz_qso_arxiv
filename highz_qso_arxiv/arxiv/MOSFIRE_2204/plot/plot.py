from highz_qso_arxiv.plot import plot_series

"""
    plot all with telluric
"""
name_list = ["J1150-0143", "J1332+0150",
             "J1202+0129", "J1333+0919",
             "J1318+2932", "J1511+0344"]

fits_list = [f"../reduced/all_redo/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
telluric_list = [True for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, display=False, save_file="MOSFIRE_2204_tellcorr.pdf")
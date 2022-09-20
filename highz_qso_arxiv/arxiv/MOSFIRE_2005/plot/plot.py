from highz_qso_arxiv.plot import plot_series

# name_list = ["J1311+5932", "J1342+0928", "J1600+3127"]
name_list = ["J1311+5932", "J1600+3127"]

fits_list = [f"../reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
telluric_list = [False for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
qso_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, 
            qso_list=qso_list, display=False, save_file="MOSFIRE_2005.pdf")

fits_list = [f"../reduced/all/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
telluric_list = [True for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, 
            qso_list=qso_list, display=False, save_file="MOSFIRE_2005_tellcorr.pdf")

from highz_qso_arxiv.plot import plot_series

# name_list = ["J1007+2115", "J1432+5746",
#             "J1450+3302", "J1458+6300",
#             "J1635+1758", "J2212+2040"]

name_list = ["J1432+5746",
            "J1450+3302", "J1458+6300",
            "J1635+1758", "J2212+2040"]

fits_list = [f"../reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
telluric_list = [False for i in range(len(name_list))]
template_list = [True for i in range(len(name_list))]
qso_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, 
            qso_list=qso_list, display=False, save_file="NIRES_1905.pdf")

fits_list = [f"../reduced/all/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
telluric_list = [True for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, 
            qso_list=qso_list, display=False, save_file="NIRES_1905_tellcorr.pdf")
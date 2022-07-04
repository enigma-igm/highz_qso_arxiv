from highz_qso_arxiv.plot import plot_series

"""
    plot all with telluric
"""
name_list = ['J0637+3812', 'J0854+2908',
             'J0838-0112', 'J0938+1341',
             'J0847+0139', 'J1024-0037']
fits_list = [f"../reduced/all/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
telluric_list = [True for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, display=False, save_file="MOSFIRE_2201_tellcorr.pdf")
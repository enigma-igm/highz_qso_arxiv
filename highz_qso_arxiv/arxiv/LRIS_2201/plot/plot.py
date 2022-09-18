from highz_qso_arxiv.plot import plot_series

"""
    plot all without tellruic
"""
name_list = ['J0739+2328', 'J0759+2811', 'J0820+2412', 'J0831+2558',
             'J0849+2601', 'J0936+3244', 'J0936+3346', 'J0938+3332',
             'J1003-2610', 'J1030-0031', 'J1145+0933', 'J1146+0039',
             'J1206-0051', 'J1209+0135', 'J1241-0134', 'J1246-3045',
             'J1256-0306', 'J1305-1549', 'J1326+0927']
fits_list = [f"../reduced/all/{nm}/{nm}_coadd.fits" for nm in name_list]
fits_list[10] = f"../reduced/all/coadd2d/{name_list[10]}_coadd.fits"
fits_list[12] = f"../reduced/all/coadd2d/{name_list[12]}_coadd.fits"
fits_list[16] = f"../reduced/all/coadd2d/{name_list[16]}_coadd.fits"
idx_list = [1 for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, display=False, save_file="LRIS_2201.pdf")

"""
    plot all with telluric
"""
fits_list = [f"../reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in name_list]
fits_list[10] = f"../reduced/all/coadd2d/{name_list[10]}_coadd_tellcorr.fits"
fits_list[12] = f"../reduced/all/coadd2d/{name_list[12]}_coadd_tellcorr.fits"
fits_list[16] = f"../reduced/all/coadd2d/{name_list[16]}_coadd_tellcorr.fits"
idx_list = [1 for i in range(len(name_list))]
telluric_list = [True for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, display=False, save_file="LRIS_2201_tellcorr.pdf")

"""
    plot J1206 and J1256
"""
name_list = ['J1206-0051', 'J1256-0306']
fits_list = [f"../reduced/all/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
telluric_list = [True for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, display=False, save_file="J1206_J1256.pdf")

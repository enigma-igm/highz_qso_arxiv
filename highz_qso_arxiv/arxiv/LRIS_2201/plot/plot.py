from re import template
from highz_qso_arxiv.plot import plot_series

name_list = ['J0739+2328', 'J0759+2811', 'J0831+2558',
             'J0849+2601', 'J0936+3244', 'J0936+3346', 'J0938+3332',
             'J1003-2610', 'J1030-0031', 'J1145+0933', 'J1146+0039',
             'J1206-0051', 'J1209+0135', 'J1241-0134', 'J1246-3045',
             'J1256-0306', 'J1305-1549', 'J1326+0927']
fits_list = [f"../reduced/all/{nm}/{nm}_coadd.fits" for nm in name_list]
fits_list[9] = f"../reduced/all/coadd2d/{name_list[9]}_coadd.fits"
fits_list[11] = f"../reduced/all/coadd2d/{name_list[11]}_coadd.fits"
fits_list[15] = f"../reduced/all/coadd2d/{name_list[15]}_coadd.fits"
idx_list = [1 for i in range(len(name_list))]
template_list = [True for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, template_list=template_list, display=False, save_file="LRIS_2201.pdf")

fits_list = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_list]
telluric_list = [True for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, template_list=template_list, telluric_list=telluric_list, display=False, save_file="LRIS_2201_tellcorr.pdf")
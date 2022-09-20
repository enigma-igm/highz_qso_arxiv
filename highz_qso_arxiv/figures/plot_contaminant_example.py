from re import template
from highz_qso_arxiv.plot import plot_series

name_list = ['J0739+2328',
             'J1206-0051', 'J1319+0101']
fits_list = [f"../arxiv/LRIS_2201/reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in name_list]
fits_list[1] = f"../arxiv/LRIS_2201/reduced/all/coadd2d/{name_list[1]}_coadd_tellcorr.fits"
fits_list[2] = f"../arxiv/LRIS_2203/LRIS_220305/reduced/all/{name_list[2]}/{name_list[2]}_coadd_tellcorr.fits"

idx_list = [1 for i in range(len(name_list))]
template_list = [True, False, False]
telluric_list = [True for i in range(len(name_list))]
qso_list = [False, False, True]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, qso_list=qso_list, display=False, save_file="contaminant_example.pdf")
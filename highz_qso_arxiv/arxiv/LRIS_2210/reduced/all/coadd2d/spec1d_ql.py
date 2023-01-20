from highz_qso_arxiv.plot import plot_series

name_all = []
fits_all = []

name_list = ['J0703+5159', 'J0714+4653']
fits_list = ['J0703+5159_coadd.fits', 'J0714+4653_coadd.fits']
idx_list = [1, 1]
template_list = [True, True]
qso_list = [False, False]

xrange = (7300, 10600)
plot_series(name_list, fits_list, idx_list, telluric=False, template_list=template_list, 
            xrange=xrange, qso_list=qso_list, display=True)

from highz_qso_arxiv.plot import plot_series

name_list = ['J0759+2811', 'J0739+2328', 'J0820+2412', 'J0831+2558',
             'J0849+2601', 'J0936+3244', 'J0936+3346', 'J0938+3332',
             'J1030-0031', 'J1003-2610', 'J1145+0933', 'J1146+0039',
             'J1206-0051', 'J1241-0134', 'J1246-3045', 'J1256-0306',
             'J1305-1549', 'J1326+0927', 'J1209+0135']
fits_list = [f"./{nm}/{nm}_coadd.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, display=False, save_file="LRIS_2201.pdf")
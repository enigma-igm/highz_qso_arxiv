from highz_qso_arxiv.plot import plot_series

nonqso_targets_0305 = ["J1250+0347", "J1200+0112", "J1335+0103", "J1355-0111", \
                       "J1422+0307", "J1327-0207", "J1355-0044", "J1143-0248", \
                       "J1356+0050", "J1100+0203", "J1223+0114", "J1358-0118"]
fits_list = [f"../LRIS_220305/reduced/all/{nm}/{nm}_coadd.fits" for nm in nonqso_targets_0305]
idx_list = [1 for i in range(len(nonqso_targets_0305))]
plot_series(nonqso_targets_0305, fits_list, idx_list, display=False, save_file="LRIS_220305_noqso.pdf")

nonqso_targets_0306 = ["J0810+2352", "J0850+0146", "J0911+0022", "J0947+0111", \
                       "J1153-2239", "J1233-2807", "J1412-2757", "J1425+0004", \
                       "J1427+0202", "J1510-0144", "J1514-0121", "J1522+0051", \
                       "J1534+0046", "J1655-0051"]
fits_list = [f"../LRIS_220306/reduced/all/{nm}/{nm}_coadd.fits" for nm in nonqso_targets_0306]
idx_list = [1 for i in range(len(nonqso_targets_0306))]
plot_series(nonqso_targets_0306, fits_list, idx_list, display=False, save_file="LRIS_220306_noqso.pdf")

qso_targets = ["J1458+1012", "J1319+0101", "J1724+3718", "J0901+2906"]
fits_list = [f"../LRIS_220306/reduced/quasar/{nm}_coadd_tellcorr.fits" for nm in qso_targets]
fits_list[1] = f"../LRIS_220305/reduced/quasar/J1319+0101_coadd_tellcorr.fits"
idx_list = [1 for i in range(len(qso_targets))]
plot_series(qso_targets, fits_list, idx_list, telluric=True, display=False, save_file="LRIS_2203_qso.pdf")

from highz_qso_arxiv.plot import plot_series

"""
    plot all without tellruic
"""
name_list = ["J0854+2908", "J0954-0117", "J1100+0237", "J1111+0640",
             "J1238+0234", "J1253-0200", "J1317-0302",
             "J1340+1422", "J1349+0805", "J1351+0128", "J1352+0002",
             "J1401+4542", "J1426+0148", "J1429+0219", "J1434+0240", "J1438-0103", 
             "J1441+0149", "J1450+0228", "J1508-0109", "J1513-0059", 
             "J1523+2935", "J1532-0155", "J1535+6146", 
             "J1625+2414", "J1635+5940", "J2212+2040"]
fits_list = [f"../reduced/all/{nm}/{nm}_coadd.fits" for nm in name_list]
fits_list[0] = "../reduced/all/coadd2d/J0854+2908_coadd.fits"
fits_list[2] = "../reduced/J1100/coadd2d/J1100+0237_coadd.fits"
fits_list[-1] = "../reduced/J2212/J2212+2040_coadd.fits"
idx_list = [1 for i in range(len(name_list))]
template_list = [True for i in range(len(fits_list))]
telluric_list = [False for i in range(len(fits_list))]

template_list = [True for i in range(len(fits_list))]
template_list[name_list.index('J1111+0640')] = False
template_list[name_list.index('J1401+4542')] = False
template_list[name_list.index('J1523+2935')] = False

qso_list = [False for i in range(len(fits_list))]
qso_list[name_list.index('J1111+0640')] = True
qso_list[name_list.index('J1401+4542')] = True
qso_list[name_list.index('J1523+2935')] = True

# divide list into three parts and plot

plot_series(name_list[:8], fits_list[:8], idx_list[:8], telluric_list=telluric_list[:8], template_list=template_list[:8], 
            qso_list=qso_list[:8], display=False, save_file="LRIS_2204_part1.pdf")
plot_series(name_list[8:16], fits_list[8:16], idx_list[8:16], telluric_list=telluric_list[8:16], template_list=template_list[8:16], 
            qso_list=qso_list[8:16], display=False, save_file="LRIS_2204_part2.pdf")
plot_series(name_list[16:], fits_list[16:], idx_list[16:], telluric_list=telluric_list[16:], template_list=template_list[16:], 
            qso_list=qso_list[16:], display=False, save_file="LRIS_2204_part3.pdf")

plot_series(name_list, fits_list, idx_list, template_list=template_list, telluric_list=telluric_list, 
            qso_list=qso_list, display=False, save_file="LRIS_2204.pdf")

"""
    plot all with telluric
"""
fits_list = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_list]
telluric_list = [True for i in range(len(fits_list))]

plot_series(name_list[:8], fits_list[:8], idx_list[:8], telluric_list=telluric_list[:8], template_list=template_list[:8], 
            qso_list=qso_list[:8], display=False, save_file="LRIS_2204_tellcorr_part1.pdf")
plot_series(name_list[8:16], fits_list[8:16], idx_list[8:16], telluric_list=telluric_list[8:16], template_list=template_list[8:16], 
            qso_list=qso_list[8:16], display=False, save_file="LRIS_2204_tellcorr_part2.pdf")
plot_series(name_list[16:], fits_list[16:], idx_list[16:], telluric_list=telluric_list[16:], template_list=template_list[16:], 
            qso_list=qso_list[16:], display=False, save_file="LRIS_2204_tellcorr_part3.pdf")

plot_series(name_list, fits_list, idx_list, template_list=template_list, telluric_list=telluric_list, 
            qso_list=qso_list, display=False, save_file="LRIS_2204_tellcorr.pdf")
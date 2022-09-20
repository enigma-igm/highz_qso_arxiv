from highz_qso_arxiv.plot import plot_series

# name_list = ["J0020+1956",
#             "J0038-0653", "J0102-3327",
#             "J0159-2017", "J0220+3458",
#             "J0251+1724", "J0339-1610", "J0705+2850",
#             "J0758+0040", "J0759+0735",
#             "J0812+3201", "J0813+2219",
#             "J0822+5045", "J0826+1217",
#             "J0845+3033", "J2033-0428",
#             "J2105-0049", "J2119-1335",
#             "J2129+0601", "J2212+2040",
#             "J2244+2907", "J2247+3247",
#             "J2312+3215", "J2328+2755",
#             "J2346+2545", "J2352+2655",
#             "J2355+2048",]

name_list = ["J0020+1956", "J0102-3327",
            "J0159-2017", "J0220+3458",
            "J0251+1724", "J0339-1610", "J0705+2850",
            "J0758+0040", "J0759+0735",
            "J0812+3201", "J0813+2219",
            "J0822+5045", "J0826+1217",
            "J0845+3033", "J2033-0428",
            "J2105-0049", "J2119-1335",
            "J2129+0601", "J2212+2040",
            "J2244+2907", "J2247+3247",
            "J2312+3215", "J2328+2755",
            "J2346+2545", "J2352+2655",
            "J2355+2048",]

fits_list = [f"../reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]
fits_list[name_list.index('J0020+1956')] = f"../reduced/all_day4/coadd2d/J0020+1956_coadd.fits"
fits_list[name_list.index('J0339-1610')] = f"../reduced/J0339/coadd2d/J0339-1610_coadd.fits"
fits_list[name_list.index('J0759+0735')] = f"../reduced/all_day4/coadd2d/J0759+0735_coadd.fits"
fits_list[name_list.index('J2328+2755')] = f"../reduced/all_day4/coadd2d/J2328+2755_coadd.fits"
fits_list[name_list.index('J2352+2655')] = f"../reduced/all_day4/coadd2d/J2352+2655_coadd.fits"
idx_list = [1 for i in range(len(name_list))]
telluric_list = [False for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
qso_list = [False for i in range(len(name_list))]

# divide lists into three parts and plot
plot_series(name_list[0:9], fits_list[0:9], idx_list[0:9], template_list=template_list[0:9], telluric_list=telluric_list[0:9],
            qso_list=qso_list[0:9], display=False, save_file="MOSFIRE_2010_part1.pdf")
plot_series(name_list[9:18], fits_list[9:18], idx_list[9:18], template_list=template_list[9:18], telluric_list=telluric_list[9:18],
            qso_list=qso_list[9:18], display=False, save_file="MOSFIRE_2010_part2.pdf")
plot_series(name_list[18:], fits_list[18:], idx_list[18:], template_list=template_list[18:], telluric_list=telluric_list[18:],
            qso_list=qso_list[18:], display=False, save_file="MOSFIRE_2010_part3.pdf")

plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, qso_list=qso_list, display=False, save_file="MOSFIRE_2010.pdf")

fits_list = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_list]
telluric_list = [True for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
qso_list = [False for i in range(len(name_list))]

plot_series(name_list[0:9], fits_list[0:9], idx_list[0:9], template_list=template_list[0:9], telluric_list=telluric_list[0:9],
            qso_list=qso_list[0:9], display=False, save_file="MOSFIRE_2010_tellcorr_part1.pdf")
plot_series(name_list[9:18], fits_list[9:18], idx_list[9:18], template_list=template_list[9:18], telluric_list=telluric_list[9:18],
            qso_list=qso_list[9:18], display=False, save_file="MOSFIRE_2010_tellcorr_part2.pdf")
plot_series(name_list[18:], fits_list[18:], idx_list[18:], template_list=template_list[18:], telluric_list=telluric_list[18:],
            qso_list=qso_list[18:], display=False, save_file="MOSFIRE_2010_tellcorr_part3.pdf")
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, display=False, save_file="MOSFIRE_2010_tellcorr.pdf")

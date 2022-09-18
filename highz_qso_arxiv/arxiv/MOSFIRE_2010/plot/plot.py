from highz_qso_arxiv.plot import plot_series

"""
    plot all with telluric
"""
name_list = ["J0020+1956",
            "J0038-0653", "J0102-3327",
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
fits_list[0] = f"../reduced/all_day4/coadd2d/J0020+1956_coadd.fits"
fits_list[6] = f"../reduced/J0339/coadd2d/J0339-1610_coadd.fits"
fits_list[9] = f"../reduced/all_day4/coadd2d/J0759+0735_coadd.fits"
fits_list[23] = f"../reduced/all_day4/coadd2d/J2328+2755_coadd.fits"
fits_list[25] = f"../reduced/all_day4/coadd2d/J2352+2655_coadd.fits"
idx_list = [1 for i in range(len(name_list))]
telluric_list = [False for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, display=False, save_file="MOSFIRE_2010.pdf")

fits_list = [f"../reduced/all/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
fits_list[0] = f"../reduced/all_day4/coadd2d/J0020+1956_coadd_tellcorr.fits"
fits_list[6] = f"../reduced/J0339/coadd2d/J0339-1610_coadd_tellcorr.fits"
fits_list[9] = f"../reduced/all_day4/coadd2d/J0759+0735_coadd_tellcorr.fits"
fits_list[23] = f"../reduced/all_day4/coadd2d/J2328+2755_coadd_tellcorr.fits"
fits_list[25] = f"../reduced/all_day4/coadd2d/J2352+2655_coadd_tellcorr.fits"
idx_list = [1 for i in range(len(name_list))]
telluric_list = [True for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, display=False, save_file="MOSFIRE_2010_tellcorr.pdf")

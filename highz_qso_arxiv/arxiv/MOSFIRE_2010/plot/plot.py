from highz_qso_arxiv.plot import plot_series

"""
    plot all with telluric
"""
name_list = ["J0038-0653", "J0813+2219", "J2129+0601",
            "J0102-3327", "J0822+5045", "J2212+2040",
            "J0159-2017", "J0826+1217", "J2244+2907",
            "J0220+3458", "J0845+3033", "J2247+3247",
            "J0251+1724", "J2033-0428", "J2312+3215",
            "J0705+2850", "J2105-0049", "J2346+2545",
            "J0758+0040", "J2119-1335", "J2355+2048"]
fits_list = [f"../reduced/all/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
telluric_list = [True for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, display=False, save_file="MOSFIRE_2010_tellcorr.pdf")

name_list = ["J0020+1956", "J2328+2755", "J0759+0735", "J2352+2655"]
fits_list = [f"../reduced/all_day4/coadd2d/{nm}_coadd_tellcorr.fits" for nm in name_list]
idx_list = [1 for i in range(len(name_list))]
telluric_list = [True for i in range(len(name_list))]
template_list = [False for i in range(len(name_list))]
plot_series(name_list, fits_list, idx_list, telluric_list=telluric_list, template_list=template_list, display=False, save_file="MOSFIRE_2010_day4_tellcorr.pdf")
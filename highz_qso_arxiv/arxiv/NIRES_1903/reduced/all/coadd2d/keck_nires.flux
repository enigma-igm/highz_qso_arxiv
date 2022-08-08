# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-07-28

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True
  
# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_NR.20190314.21243-NR.20190314.21243-J0739+4843.fits ../GD71_nires_sens.fits
 Science_coadd/spec1d_NR.20190314.22418-NR.20190314.22418-J0800+5128.fits
 Science_coadd/spec1d_NR.20190314.23543-NR.20190314.25436-J0756+5744.fits
 Science_coadd/spec1d_NR.20190314.52459-NR.20190314.52813-J1312+5707.fits
 Science_coadd/spec1d_NR.20190314.54277-NR.20190314.54631-J1546+2945.fits
 Science_coadd/spec1d_NR.20190315.24463-NR.20190315.24463-J0800+3034.fits
 Science_coadd/spec1d_NR.20190315.31485-NR.20190315.31838-J0953.fits
 Science_coadd/spec1d_NR.20190315.33275-NR.20190315.33275-J0803+0030.fits
 Science_coadd/spec1d_NR.20190315.34007-NR.20190315.34361-J0816+1622.fits
 Science_coadd/spec1d_NR.20190315.36689-NR.20190315.37042-J1217+3136.fits
 Science_coadd/spec1d_NR.20190315.45339-NR.20190315.53538-J1635+5940.fits
 Science_coadd/spec1d_NR.20190315.47367-NR.20190315.47721-J1550+2558.fits
 Science_coadd/spec1d_NR.20190315.49129-NR.20190315.49483-1638+5412.fits
 Science_coadd/spec1d_NR.20190315.55055-NR.20190315.55409-J1423+2901.fits
flux end


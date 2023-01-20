# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-09-25

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_MF.20200529.42212-MF.20200529.42694-J1517+2551.fits ../Feige110_mosfire_sens.fits
 Science_coadd/spec1d_MF.20200529.43510-MF.20200529.44669-J1626+1337.fits
 Science_coadd/spec1d_MF.20200529.45548-MF.20200529.45708-J1526+2848.fits
 Science_coadd/spec1d_MF.20200529.47793-MF.20200529.47953-J1642+2731.fits
 Science_coadd/spec1d_MF.20200529.50195-MF.20200529.52961-J2212+2040.fits
flux end


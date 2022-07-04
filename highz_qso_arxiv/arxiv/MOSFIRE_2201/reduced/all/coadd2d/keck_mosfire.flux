# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-07-03

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_MF.20220111.33636-MF.20220111.34558-J0637+3812.fits ../GD153_mosfire_sens.fits
 Science_coadd/spec1d_MF.20220111.35725-MF.20220111.36647-J0730+5949.fits
 Science_coadd/spec1d_MF.20220111.39267-MF.20220111.40200-J0838.fits
 Science_coadd/spec1d_MF.20220111.42075-MF.20220111.42993-J0847+0139.fits
 Science_coadd/spec1d_MF.20220111.44492-MF.20220111.44674-J0854+2908.fits
 Science_coadd/spec1d_MF.20220111.50331-MF.20220111.50514-J1024.fits
 Science_coadd/spec1d_MF.20220111.56422-MF.20220111.56605-J1506+3316.fits
 Science_coadd/spec1d_MF.20220111.46711-MF.20220111.48430-J0938+1341.fits
flux end


# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2023-01-06

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_r221027-r221027-J0104+1533.fits ../GD71_lris_sens.fits
 Science_coadd/spec1d_r221027-r221027-J0633+5524.fits
 Science_coadd/spec1d_r221027-r221027-J0703+5159.fits
 Science_coadd/spec1d_r221027-r221027-J0708+4058.fits
 Science_coadd/spec1d_r221027-r221027-J0714+4653.fits
 Science_coadd/spec1d_r221027-r221027-J0736+2743.fits
 Science_coadd/spec1d_r221027-r221027-J0738+2749.fits
 Science_coadd/spec1d_r221028-r221028-J0126.fits
 Science_coadd/spec1d_r221028-r221028-J0201.fits
 Science_coadd/spec1d_r221028-r221028-J0216.fits
 Science_coadd/spec1d_r221028-r221028-J0249.fits
 Science_coadd/spec1d_r221028-r221028-J0740+2637.fits
 Science_coadd/spec1d_r221028-r221028-J0741+2615.fits
 Science_coadd/spec1d_r221028-r221028-J0748+2331.fits
 Science_coadd/spec1d_r221028-r221028-J0748+2716.fits
 Science_coadd/spec1d_r221028-r221028-J0749+2005.fits
 Science_coadd/spec1d_r221028-r221028-J0749+2356.fits
 Science_coadd/spec1d_r221028-r221028-J0750+2036.fits
 Science_coadd/spec1d_r221028-r221028-J0751+2653.fits
flux end


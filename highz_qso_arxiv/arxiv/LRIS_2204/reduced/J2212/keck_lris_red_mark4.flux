# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-08-08

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True
  
# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science/spec1d_r220423_00111-J2212+2040_OFF_LRISr_20220423T143645.158.fits GD153_lris_sens.fits
 Science/spec1d_r220423_00112-J2212+2040_OFF_LRISr_20220423T144221.773.fits
 Science/spec1d_r220423_00115-J2212+2040_OFF_LRISr_20220423T145416.819.fits
 Science/spec1d_r220423_00117-lds749_LRISr_20220423T150657.744.fits
flux end


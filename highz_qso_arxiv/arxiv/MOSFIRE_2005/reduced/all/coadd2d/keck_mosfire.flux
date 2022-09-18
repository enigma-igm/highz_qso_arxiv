# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-09-19

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_MF.20200528.36909-MF.20200529.29112-J1410+3927.fits ../Feige110_mosfire_sens.fits
 Science_coadd/spec1d_MF.20200528.41023-MF.20200528.41184-J1600+3127.fits
 Science_coadd/spec1d_MF.20200528.48923-MF.20200528.49083-J2133+0605.fits
 Science_coadd/spec1d_MF.20200528.50202-MF.20200528.50362-J2129+0601.fits
 Science_coadd/spec1d_MF.20200529.21866-MF.20200529.22026-J1056+4613.fits
 Science_coadd/spec1d_MF.20200529.24519-MF.20200529.24680-J1100+5215.fits
 Science_coadd/spec1d_MF.20200529.27374-MF.20200529.27535-J1141+2822.fits
 Science_coadd/spec1d_MF.20200529.30553-MF.20200529.30713-J1342+0928.fits
 Science_coadd/spec1d_MF.20200529.34819-MF.20200529.34980-J1311+5932.fits
 Science_coadd/spec1d_MF.20200529.40064-MF.20200529.40225-J1517+2846.fits
flux end


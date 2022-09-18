# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-09-16

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_MF.20191118.30271-MF.20191118.30943-J0007+2517.fits ../GD71_mosfire_sens.fits
 Science_coadd/spec1d_MF.20191118.33822-MF.20191118.34492-J0038.fits
 Science_coadd/spec1d_MF.20191118.39696-MF.20191118.40352-J0058+1715.fits
 Science_coadd/spec1d_MF.20191118.42213-MF.20191118.42866-J0208+1749.fits
 Science_coadd/spec1d_MF.20191118.44344-MF.20191118.45000-J0313.fits
 Science_coadd/spec1d_MF.20191118.45957-MF.20191118.46603-J0426+0221.fits
 Science_coadd/spec1d_MF.20191118.47544-MF.20191118.47544-J0427+0327.fits
 Science_coadd/spec1d_MF.20191119.26551-MF.20191120.20095-J2140+1608.fits
 Science_coadd/spec1d_MF.20191119.28208-MF.20191119.28840-J2245+2709.fits
 Science_coadd/spec1d_MF.20191119.29889-MF.20191119.30521-J2247+3247.fits
 Science_coadd/spec1d_MF.20191119.31698-MF.20191119.32330-J2250.fits
 Science_coadd/spec1d_MF.20191120.18157-MF.20191120.18490-J2128+2231.fits
 Science_coadd/spec1d_MF.20191120.22373-MF.20191120.22556-J2132.fits
 Science_coadd/spec1d_MF.20191120.23456-MF.20191120.24921-J2258.fits
 Science_coadd/spec1d_MF.20191120.25568-MF.20191120.25901-J2333.fits
 Science_coadd/spec1d_MF.20191120.27025-MF.20191120.27359-J0022.fits
 Science_coadd/spec1d_MF.20191120.28915-MF.20191120.29247-J2312+3215.fits
 Science_coadd/spec1d_MF.20191120.29875-MF.20191120.30207-J2328+2755.fits
flux end


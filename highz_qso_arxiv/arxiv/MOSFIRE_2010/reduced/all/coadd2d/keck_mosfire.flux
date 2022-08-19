# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-08-18

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_MF.20201022.20660-MF.20201022.22319-z7.fits ../GD71_mosfire_sens.fits
 Science_coadd/spec1d_MF.20201022.23624-MF.20201022.26033-z7.fits
 Science_coadd/spec1d_MF.20201022.28846-MF.20201022.29769-z7.fits
 Science_coadd/spec1d_MF.20201022.31023-MF.20201022.31943-z7.fits
 Science_coadd/spec1d_MF.20201022.33150-MF.20201022.34069-z7.fits
 Science_coadd/spec1d_MF.20201022.35878-MF.20201022.36061-z7.fits
 Science_coadd/spec1d_MF.20201022.41614-MF.20201022.43274-z7.fits
 Science_coadd/spec1d_MF.20201022.46190-MF.20201022.47109-z7.fits
 Science_coadd/spec1d_MF.20201022.48080-MF.20201022.49739-z7.fits
 Science_coadd/spec1d_MF.20201022.50645-MF.20201022.51567-lo.fits
 Science_coadd/spec1d_MF.20201022.52787-MF.20201022.52970-z7.fits
 Science_coadd/spec1d_MF.20201022.55450-MF.20201022.55633-z7.fits
 Science_coadd/spec1d_MF.20201023.18335-MF.20201023.22464-z7.fits
 Science_coadd/spec1d_MF.20201023.20421-MF.20201023.20603-z7.fits
 Science_coadd/spec1d_MF.20201023.23643-MF.20201023.23826-z7.fits
 Science_coadd/spec1d_MF.20201023.25205-MF.20201023.25387-z7.fits
 Science_coadd/spec1d_MF.20201023.28048-MF.20201023.28966-z7.fits
 Science_coadd/spec1d_MF.20201023.30213-MF.20201023.30396-z7.fits
 Science_coadd/spec1d_MF.20201023.32388-MF.20201023.34788-z7.fits
 Science_coadd/spec1d_MF.20201023.51390-MF.20201023.51572-z7.fits
 Science_coadd/spec1d_MF.20201023.53311-MF.20201023.53494-z7.fits
 Science_coadd/spec1d_MF.20201023.54617-MF.20201023.54800-z7.fits
 Science_coadd/spec1d_MF.20201023.55556-MF.20201023.55738-z7.fits
 Science_coadd/spec1d_MF.20201024.17223-MF.20201024.17405-z7.fits
 Science_coadd/spec1d_MF.20201024.18434-MF.20201024.21229-z7.fits
 Science_coadd/spec1d_MF.20201024.22629-MF.20201024.25022-J0038.fits
 Science_coadd/spec1d_MF.20201024.50126-MF.20201024.51778-z7.fits
 Science_coadd/spec1d_MF.20201024.53028-MF.20201024.53953-z7.fits
 Science_coadd/spec1d_MF.20201024.55149-MF.20201024.55332-z7.fits
flux end


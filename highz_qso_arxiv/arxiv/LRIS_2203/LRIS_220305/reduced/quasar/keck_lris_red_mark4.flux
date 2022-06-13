# Auto-generated PypeIt file using PypeIt version: 1.7.1.dev544+gbb7e776b9
# 2022-03-10

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science/spec1d_r220305_00190-J1319+0101_OFF_LRISr_20220305T143022.320.fits GD153_IR_sens.fits
 Science/spec1d_r220305_00191-J1319+0101_OFF_LRISr_20220305T143559.885.fits
 Science/spec1d_r220305_00206-GD153_LRISr_20220305T155330.106.fits
flux end


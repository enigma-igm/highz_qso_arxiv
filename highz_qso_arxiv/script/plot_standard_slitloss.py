import os
import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.interpolate import interp1d

from pypeit import io
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit.sensfunc import IRSensFunc
from pypeit.images.detector_container import DetectorContainer

narrow_spec1df = os.path.join('../arxiv/LRIS_2203/LRIS_220305/reduced/all/Science/', 
                              'spec1d_r220305_00206-GD153_LRISr_20220305T155330.106.fits')
narrow_spec2df = os.path.join('../arxiv/LRIS_2203/LRIS_220305/reduced/all/Science/',
                                'spec2d_r220305_00206-GD153_LRISr_20220305T155330.106.fits')
wide_spec1df = os.path.join('../arxiv/LRIS_2203/LRIS_220305/reduced/long_8.7/Science',
                            'spec1d_r220305_00207-GD153_LRISr_20220305T155549.814.fits')
wide_spec2df = os.path.join('../arxiv/LRIS_2203/LRIS_220305/reduced/long_8.7/Science',
                                'spec2d_r220305_00207-GD153_LRISr_20220305T155549.814.fits')
detname = DetectorContainer.get_name(det=1)
narrow_sobjs = specobjs.SpecObjs.from_fitsfile(narrow_spec1df, chk_version=False)
narrow_spec2DObj = spec2dobj.Spec2DObj.from_file(narrow_spec2df, detname, chk_version=False)

wide_sobjs = specobjs.SpecObjs.from_fitsfile(wide_spec1df, chk_version=False)
wide_spec2DObj = spec2dobj.Spec2DObj.from_file(wide_spec2df, detname, chk_version=False)

pixel_scale = 0.123 # arcsec/pixel
binning = 2 # in spatial direction
wide_fwhm_func = interp1d(wide_sobjs[0].OPT_WAVE, wide_sobjs[0].FWHMFIT, fill_value="extrapolate")
wide_fwhmfit = wide_fwhm_func(narrow_sobjs[0].OPT_WAVE)
narrow_fwhmfit = narrow_sobjs[0].FWHMFIT

ratio = erf((8.7/2)/(np.sqrt(2)*wide_fwhmfit*pixel_scale*binning/2.355)) / \
        erf((1/2)/(np.sqrt(2)*narrow_fwhmfit*pixel_scale*binning/2.355))

# ratio = erf((8.7/2)/(np.sqrt(2)*wide_sobjs[0].FWHM*pixel_scale*binning/2.355)) / \
#         erf((1/2)/(np.sqrt(2)*narrow_sobjs[0].FWHM*pixel_scale*binning/2.355))
print(ratio)
plt.plot(narrow_sobjs[0].OPT_WAVE, narrow_sobjs[0].OPT_COUNTS*ratio, label='Narrow (corrected)')
plt.plot(narrow_sobjs[0].OPT_WAVE, narrow_sobjs[0].OPT_COUNTS, label='Narrow', alpha=0.8)
plt.plot(wide_sobjs[0].OPT_WAVE, wide_sobjs[0].OPT_COUNTS, label='Wide')
plt.legend()
plt.show()

plt.plot(narrow_sobjs[0].OPT_WAVE, ratio)
plt.show()
# wide_func = interp1d(wide_sobjs[0].OPT_WAVE, wide_sobjs[0].OPT_COUNTS, fill_value="extrapolate")
# plt.plot(narrow_sobjs[0].OPT_WAVE, narrow_sobjs[0].OPT_COUNTS/wide_func(narrow_sobjs[0].OPT_WAVE))
# plt.show()
import numpy as np
from scipy.special import erf

from IPython import embed

__all__ = ["gauss_comb"]

def gauss_comb(shape, center, sig=5.):
    """
    Creat gassian profiles along x-axis (spatial) in an image
    """
    img = np.zeros(shape, dtype=float)
    x = np.arange(shape[1])
    for i,_c in enumerate(center):
        img[i,:] = (erf((x-center[i]+0.5)/np.sqrt(2)/sig) 
                   - erf((x-center[i]-0.5)/np.sqrt(2)/sig))/2.
    return img
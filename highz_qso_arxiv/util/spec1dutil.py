import numpy as np

from IPython import embed

__all__ = ["gp_trough", "damping_wing", "extend_to_lower"]

def gp_trough(wave, flux, redshift):
    """
    add gunn-peterson trough on the spectrum
    TODO: we need real absorption
    """
    wl_lya = 1215.67 * (1 + redshift)
    trough = wave < wl_lya
    flux[trough] = 0
    return flux

def damping_wing(wave, flux, redshift):
    """
    add damping wing on the spectrum
    TODO: add this
    """
    pass

def extend_to_lower(wave, counts, wave_min):
    """
    expand wave to wave_min, fill in zero
    """
    if wave[0] < wave_min:
        return wave, counts
    try:
        unit = counts.unit
    except AttributeError:
        unit = 1

    wave_min_old = np.min(wave)
    wave_new = np.arange(wave_min, wave_min_old, 1)
    counts_new = np.zeros_like(wave_new) * unit
    wave_new = np.append(wave_new, wave)
    counts_new = np.append(counts_new, counts)
    return wave_new, counts_new
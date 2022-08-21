import numpy as np
from scipy import interpolate

from IPython import embed

__all__ = ['add_gp_trough', 'add_damping_wing', 'add_telluric', 'extend_to_lower', 'rescale']

def add_gp_trough(wave, flux, redshift):
    """
    add gunn-peterson trough on the spectrum
    TODO: we need real absorption
    """
    wl_lya = 1215.67 * (1 + redshift)
    trough = wave < wl_lya
    flux[trough] = 0
    return flux

def add_damping_wing(wave, flux, redshift):
    """
    add damping wing on the spectrum
    TODO: add this
    """
    pass

def add_telluric(wave, flux, telluric):
    wave_tell, tell = telluric.WAVE[0], telluric.TELLURIC[0]
    tell_intep = interpolate.interp1d(wave_tell, tell)
    mask = (wave > wave_tell[0]) & (wave < wave_tell[-1])
    flux[mask] = flux[mask] * tell_intep(wave[mask])
    return wave, flux

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

def rescale(wl, flux, wl_base, flux_base, err_base=None):
    """
    rescale the spectrum to the base spectrum
    """
    if err_base is None:
        err_base = np.ones_like(flux_base)
    template_interp_func = interpolate.interp1d(wl, flux, kind="cubic")
    mask_for_scale = (wl_base < min(np.max(wl), np.max(wl_base))) & (wl_base > max(np.min(wl), np.min(wl_base)))

    # interpolate the template to the same wavelength as the data
    flux_template_interp = template_interp_func(wl_base[mask_for_scale])

    # calculate the scale factor
    # d(chi_sq)/d(scale_factor) = 0
    scale = np.sum(flux_base[mask_for_scale]*flux_template_interp/err_base[mask_for_scale]**2) / \
            np.sum(flux_template_interp**2/err_base[mask_for_scale]**2)
    flux_template_interp = flux_template_interp * scale
    flux = scale * flux
    chi_sq = np.sum((flux_base[mask_for_scale]-flux_template_interp)**2/err_base[mask_for_scale]**2)
    chi_sq = chi_sq / (len(flux_base[mask_for_scale]) - 1)
    return flux, chi_sq
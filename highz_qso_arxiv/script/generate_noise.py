import numpy as np
from astropy.io import fits
from scipy import interpolate
from astropy.table import Table

from highz_qso_arxiv.util import get_project_root

from IPython import embed

def interp(err_f):
    # Interpolates and generates the error cumulatives from different flux bins
    n = len(err_f)
    C = np.arange(n) / (n - 1)
    sig_sort = err_f[err_f.argsort()]
    return interpolate.interp1d(C, sig_sort)

def add_noise(syn_flux, bin_edges_ref, **kwargs):
    # This function convolves the simulated / deconvolved fluxes 
    # with errors taken from the cumulative error distributions coming from real flux
    # Here it opens the real noisy data

    root = get_project_root()

    XD_params = {
        'table_name': root / 'resource/catalog/VIKING_catalog_clean_nobright.fits',
        'ref_mag':'J_mag_aper_3p0',
        'fluxes':['J_flux_aper_3p0', 'flux_z','flux_w1'],
        'fluxes_err':['J_flux_aper_err_3p0', 'flux_z_err','flux_w1_err'],
    }

    hdu_list = fits.open(XD_params['table_name'], memmap=True)
    output = Table(hdu_list[1].data)

    flux = np.array([output[XD_params['fluxes'][i]] for i in range(len(XD_params['fluxes']))])
    flux_err = np.array([output[XD_params['fluxes_err'][i]] for i in range(len(XD_params['fluxes_err']))])
    mag_ref=output[XD_params['ref_mag']]

    flux = flux[:, (mag_ref <= bin_edges_ref[1]) & (mag_ref > bin_edges_ref[0])]
    flux_err = flux_err[:, (mag_ref <= bin_edges_ref[1]) & (mag_ref > bin_edges_ref[0])]
    mag_ref = mag_ref[(mag_ref <= bin_edges_ref[1]) & (mag_ref > bin_edges_ref[0])]

    # Interpolate the noise from real data
    bin = kwargs.get('bin', 150)
    #bin_edges = np.zeros((len(flux), bin + 1))
    data_points_per_bin = len(mag_ref) // bin
    interp_err = []
    syn_flux_err = np.zeros(syn_flux.shape)
    for i in range(len(flux)):
        data = np.sort(flux[i])
        bins = [data[_ * data_points_per_bin: (_ + 1) * data_points_per_bin] for _ in range(bin)]
        #bin_edges[i] = np.append([min(bins[j]) for j in range(bin)], max(data))
        bin_edges = np.array([min(bins[j]) for j in range(bin)], max(data))
        ind = np.digitize(flux[i], bin_edges)
        interp_err += [[interp(flux_err[i, ind == j]) for j in range(1, bin + 1)]]

        index_sim = np.digitize(syn_flux[i], bin_edges)
        for j in range(len(syn_flux[0, :])):
            if index_sim[j]>=bin+1: index_sim[j]=bin
            if index_sim[j] == 0: index_sim[j] = 1
            syn_flux_err[i, j] = interp_err[i][index_sim[j] - 1](np.random.uniform())
            syn_flux[i, j] = np.random.normal(syn_flux[i, j], syn_flux_err[i, j])

    return syn_flux, syn_flux_err
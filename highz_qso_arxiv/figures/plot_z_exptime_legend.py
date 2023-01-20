import numpy as np
import matplotlib.pyplot as plt

import highz_qso_arxiv.plot as hzp

from IPython import embed

import matplotlib as mpl
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
CB_color_cycle = ['#b2182b', '#ef8a62', '#67a9cf', '#2166ac']

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--show', action='store_true')
args = parser.parse_args()

snr_thresh = 6
snr = np.array([15.2255819 , 18.4389655 , 20.86029965, 11.36261624, 13.44583139,
       20.16298293, 13.65737015, 13.01396451, 12.36416336, 19.4181692 ,
       21.5769301 ,  8.739841  ,  8.90679418,  9.64735524, 12.57169699,
       10.80387453,  7.02037961,  6.35158582,  4.63115404,  3.63421048,
        2.47087022,  2.33758011])

star_type = ['M9', 'L1']
snr_star = np.array([3.460383286699083, 2.6086089194558952])

exptime_mag_22 = 300*(snr_thresh/snr)**2
exptime_mag_22_star = 300*(snr_thresh/snr_star)**2

redshift = np.array([5.5, 5.6, 5.7, 5.8, 5.9, 6. , 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7,
                     6.8, 6.9, 7. , 7.1, 7.2, 7.3, 7.4, 7.5, 7.6])

# redshift = np.array([6. , 6.5, 7. , 7.5, 7.6])
# exptime_mag_22 = np.array([11.80678548, 10.31007746, 41.12275227, 786.2150739, 878.43216636])
exptime_mag_21_5 = exptime_mag_22 * 10**(-0.5*2/2.5)
exptime_mag_21_5_star = exptime_mag_22_star * 10**(-0.5*2/2.5)

exptime_mag_22_5 = exptime_mag_22 * 10**(0.5*2/2.5)
exptime_mag_22_5_star = exptime_mag_22_star * 10**(0.5*2/2.5)

exptime_mag_23 = exptime_mag_22 * 10**(1*2/2.5)
exptime_mag_23_star = exptime_mag_22_star * 10**(1*2/2.5)

# fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 8), gridspec_kw={'width_ratios': [5, 1]}, sharey=True)
fig, ax = plt.subplots(figsize=(30, 2))
ax.scatter([], [], s=200, label=r'$m_{J}=21.5$', color=CB_color_cycle[0])
ax.scatter([], [], marker='s', s=200, label=r'$m_{J}=22$', color=CB_color_cycle[1])
ax.scatter([], [], s=200, label=r'$m_{J}=22.5$', color=CB_color_cycle[2])
ax.scatter([], [], s=200, label=r'$m_{J}=23$', color=CB_color_cycle[3])
ax.plot([], [], linestyle='--', lw=2, color='k', label='observing efficiency=0.5')

ax.legend(ncol=5, loc='center', fontsize=40)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# plt.show()
plt.savefig('z_exptime_legend.pdf')
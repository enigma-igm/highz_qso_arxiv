import numpy as np
import matplotlib.pyplot as plt

import highz_qso_arxiv.plot as hzp

from IPython import embed

import matplotlib as mpl
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
CB_color_cycle = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3']

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--show', action='store_true')
args = parser.parse_args()

fig, ax = plt.subplots(figsize=(12, 8))

snr_thresh = 6

snr = np.array([15.2255819 , 18.4389655 , 20.86029965, 11.36261624, 13.44583139,
       20.16298293, 13.65737015, 13.01396451, 12.36416336, 19.4181692 ,
       21.5769301 ,  8.739841  ,  8.90679418,  9.64735524, 12.57169699,
       10.80387453,  7.02037961,  6.35158582,  4.63115404,  3.63421048,
        2.47087022,  2.33758011])
exptime_mag_22 = 300*(snr_thresh/snr)**2
redshift = np.arange(5.5, 7.7, 0.1)
ax.plot(redshift, exptime_mag_22, color=CB_color_cycle[0])
ax.scatter(redshift, exptime_mag_22, marker='s', s=50, label=r'LRIS', color=CB_color_cycle[0])

snr = np.array([3.9751916 , 3.55878821, 4.2817164 , 4.54080957, 4.160312  ,
       4.00182748, 7.99286595, 9.43244198, 6.40594656, 5.06752347])
exptime_mag_22 = (150*4)*(snr_thresh/snr)**2
redshift = np.arange(7, 8, 0.1)
ax.plot(redshift, exptime_mag_22, color=CB_color_cycle[1])
ax.scatter(redshift, exptime_mag_22, marker='s', s=50, label=r'MOSFIRE', color=CB_color_cycle[1])

snr = np.array([0.85749994, 1.54652816, 2.48852709, 3.37261222, 3.86033387,
       4.27733366, 6.82454749, 8.08224094, 7.30044491, 6.81444445,
       7.73099821, 9.73429894, 8.06743271, 5.98758221, 6.45939754])
exptime_mag_22 = (360*2)*(snr_thresh/snr)**2
redshift = np.arange(7.0, 8.5, 0.1)
ax.plot(redshift, exptime_mag_22, color=CB_color_cycle[2])
ax.scatter(redshift, exptime_mag_22, marker='s', s=50, label=r'NIRES', color=CB_color_cycle[2])

ax.set_xlabel('Redshift', fontsize=30)
ax.set_ylabel('Exposure time (s)', fontsize=30)
ax.set_yscale('log')
ax.text(0.8, 0.95, r'$m_{J}=22$', transform=ax.transAxes, fontsize=30, verticalalignment='top')

xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
ax.hlines(300, xmin, xmax, linestyle='--', color='k', label='observing efficiency=0.5')
ax.fill_between([5.,10], 0, 300, color='grey', alpha=0.3)
ax.set_xlim(xmin, xmax)
ax.tick_params(axis='both', which='major', labelsize=30, width=1, size=6)
ax.tick_params(axis='both', which='minor', labelsize=30, width=1, size=3)
ax.legend(fontsize=25)


fig.tight_layout()
fig.subplots_adjust(wspace=0)

if args.show:
    plt.show()
else:
    plt.savefig('z_exptime_all.pdf')
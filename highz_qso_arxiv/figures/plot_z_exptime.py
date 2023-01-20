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
fig, (ax2, ax) = plt.subplots(1, 2, figsize=(10, 8), gridspec_kw={'width_ratios': [1, 5]}, sharey=True)
ax.plot(redshift, exptime_mag_21_5, color=CB_color_cycle[0])
ax.scatter(redshift, exptime_mag_21_5, label=r'$m_{J}=21.5$', color=CB_color_cycle[0])

ax.plot(redshift, exptime_mag_22, color=CB_color_cycle[1])
ax.scatter(redshift, exptime_mag_22, marker='s', s=50, label=r'$m_{J}=22$', color=CB_color_cycle[1])

ax.plot(redshift, exptime_mag_22_5, color=CB_color_cycle[2])
ax.scatter(redshift, exptime_mag_22_5, label=r'$m_{J}=22.5$', color=CB_color_cycle[2])

ax.plot(redshift, exptime_mag_23, color=CB_color_cycle[3])
ax.scatter(redshift, exptime_mag_23, label=r'$m_{J}=23$', color=CB_color_cycle[3])

ax.set_xlabel('Redshift', fontsize=35)
# ax.set_ylabel('Exposure time (s)', fontsize=20)
ax.set_yscale('log')
ax.set_xticks([5.5, 6., 6.5, 7., 7.5])
# ax.set_title(f'LRIS: S/N={snr_thresh}', fontsize=25)
ax.set_ylim(5e0, 3e5)

xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
ax.hlines(300, xmin, xmax, linestyle='--', color='k', label='observing efficiency=0.5')
ax.fill_between([5.,8], 0, 300, color='grey', alpha=0.3)
ax.set_xlim(xmin, xmax)
ax.tick_params(axis='both', which='major', labelsize=30, width=1, size=8)
ax.tick_params(axis='both', which='minor', labelsize=30, width=1, size=4)
# ax.legend(fontsize=20)
ax.set_title('LRIS: SNR=6', fontsize=35)

x = np.arange(len(star_type))
ax2.scatter(x, exptime_mag_21_5_star, color=CB_color_cycle[0])
ax2.scatter(x, exptime_mag_22_star, marker='s', s=50, color=CB_color_cycle[1])
ax2.scatter(x, exptime_mag_22_5_star, color=CB_color_cycle[2])
ax2.scatter(x, exptime_mag_23_star, color=CB_color_cycle[3])
ax2.set_xticks([])
ax2.set_xlim(-0.5, 1.5)
ax2.set_ylim(ymin, ymax)
ax2.set_yscale('log')
ax2.set_ylabel('Exposure time (s)', fontsize=35)
ax2.tick_params(axis='both', which='major', labelsize=30, width=1, size=8)
ax2.tick_params(axis='both', which='minor', labelsize=30, width=1, size=4)

for i, tp in enumerate(star_type):
    ax2.text(x[i]-0.2, 10, tp, fontsize=28)

fig.tight_layout()
fig.subplots_adjust(wspace=0)

if args.show:
    plt.show()
else:
    plt.savefig('z_exptime.pdf')
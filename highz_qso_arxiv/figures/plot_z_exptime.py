import numpy as np
import matplotlib.pyplot as plt

import highz_qso_arxiv.plot as hzp

from IPython import embed

# redshift = np.array([6.5, 7, 7.5, 7.7])
# exptime_mag_22 = np.array([15.5, 189, 2763.2, 11586])

snr_thresh = 6
snr = np.array([15.2255819 , 18.4389655 , 20.86029965, 11.36261624, 13.44583139,
       20.16298293, 13.65737015, 13.01396451, 12.36416336, 19.4181692 ,
       21.5769301 ,  8.739841  ,  8.90679418,  9.64735524, 12.57169699,
       10.80387453,  7.02037961,  6.35158582,  4.63115404,  3.63421048,
        2.47087022,  2.33758011])

exptime_mag_22 = 300*(snr_thresh/snr)**2
redshift = np.array([5.5, 5.6, 5.7, 5.8, 5.9, 6. , 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7,
                     6.8, 6.9, 7. , 7.1, 7.2, 7.3, 7.4, 7.5, 7.6])

# redshift = np.array([6. , 6.5, 7. , 7.5, 7.6])
# exptime_mag_22 = np.array([11.80678548, 10.31007746, 41.12275227, 786.2150739, 878.43216636])
exptime_mag_21_5 = exptime_mag_22 * 10**(-0.5*2/2.5)
exptime_mag_22_5 = exptime_mag_22 * 10**(0.5*2/2.5)
exptime_mag_23 = exptime_mag_22 * 10**(1*2/2.5)

fig, ax = plt.subplots(figsize=(10, 8))
ax.plot(redshift, exptime_mag_21_5, color='C0')
ax.scatter(redshift, exptime_mag_21_5, label=r'$m_{J}=21.5$', color='C0')

ax.plot(redshift, exptime_mag_22, color='C1')
ax.scatter(redshift, exptime_mag_22, marker='s', s=50, label=r'$m_{J}=22$', color='C1')

ax.plot(redshift, exptime_mag_22_5, color='C2')
ax.scatter(redshift, exptime_mag_22_5, label=r'$m_{J}=22.5$', color='C2')

ax.plot(redshift, exptime_mag_23, color='C3')
ax.scatter(redshift, exptime_mag_23, label=r'$m_{J}=23$', color='C3')

ax.set_xlabel('Redshift', fontsize=20)
ax.set_ylabel('Exposure time (s)', fontsize=20)
ax.set_yscale('log')
ax.set_xticks([5.5, 6., 6.5, 7., 7.5])
ax.set_title(f'LRIS: S/N={snr_thresh}', fontsize=24)
xmin, xmax = ax.get_xlim()
ax.hlines(300, xmin, xmax, linestyle='--', color='k', label='efficiency=0.5')
ax.fill_between([5.,8], 0, 300, color='grey', alpha=0.3)
ax.set_xlim(xmin, xmax)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.legend(fontsize=20)

plt.show()
# plt.savefig('z_exptime.pdf')
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as pl
from density_contour import density_contour
#from astroML.plotting import scatter_contour

data = pyfits.getdata('decalsdr2_z6_sdss_wise_nir_wise.fits',1)
ind = np.where(data['japermag3']>0)
data = data[ind]
zj,jw1 = data['zj'],data['jw1']

allcand = pyfits.getdata('allcand.fits',1)
zjallcand,jw1allcand = allcand['zj'],allcand['jw1']

#cand = pyfits.getdata('decalsdr2_z6_sdss_wise_nir_wise_point_zjw1_z6_sdssphoto_cand_good.fits',1)
#zjcand,jw1cand = cand['zj'],cand['jw1']
cand = pyfits.getdata('tractor_dr3_fall_grdrop_NIR_ALLWISE_jw1g15l35_izg2.fits',1)
#w1snr = cand['flux_w1']*np.sqrt(cand['flux_ivar_w1'])
#indcand = np.where(w1snr>5)
#cand = cand[indcand]
zjcand,jw1cand = cand['dz_ap3']-cand['japermag3'],cand['japermag3']-cand['w1mpro']

cand65 = pyfits.getdata('decalsdr2_las_z65_candidate_spring.fits',1)
obs = cand65['obs']
zjcand65,jw1cand65 = cand65['zj'],cand65['jw1']
indobs = np.where(obs>0)

simqso = pyfits.getdata('suwz48simqso.fits',0)
z = simqso['z']
szsdss = simqso['synMag'].T[4]
sz,sj,sw1 = simqso['synMag'].T[7],simqso['synMag'].T[9]-0.938,simqso['synMag'].T[12]-2.699

indrs = np.linspace(5.45,7.25,19)
srs = indrs[:-1]+0.05
szjsdss = np.zeros(len(srs))
szj,sjw1 = np.zeros(len(srs)),np.zeros(len(srs))

for i in range(len(srs)):
    indi = np.where((z>indrs[i])&(z<indrs[i+1]))
    szjsdss[i] = np.median(szsdss[indi]-sj[indi])
    szj[i] = np.median(sz[indi]-sj[indi])
    sjw1[i] = np.median(sj[indi]-sw1[indi])
    print srs[i],szj[i],sjw1[i]


### Dwarfs
pars = pyfits.open('MLTQ_DECaLS.fits')
datas = pars[1].data
mzs,mJs,W1s,W2s = 22.5-2.5*np.log10(datas['flux_ap2z']),datas['JAperMag3'],datas['w1mpro'],datas['w2mpro']
type = datas['redshift']
indm = np.where((type=='M')&(mzs-mJs<2.5))
indl = np.where(type=='L')
indt = np.where(type=='T')

parq = pyfits.open('MLTQ_DECaLS_PUBQSO.fits')
dataq = parq[1].data
indz6q = np.where(dataq['redshift']<6.5)
indz65q = np.where((dataq['redshift']>6.5) & (dataq['redshift']<7))
indz7q = np.where(dataq['redshift']>7)
mzq,mJq,W1q,W2q = 22.5-2.5*np.log10(dataq['flux_ap2z']),dataq['JAperMag3'],dataq['w1mpro'],dataq['w2mpro']


pl.figure(figsize=(8,6))
ax=pl.subplot(111)
density_contour(zj,jw1,80,80,ax=ax,color='0.7',alpha=0.5)
#scatter_contour(zj,jw1, threshold=200, log_counts=False, ax=ax,histogram2d_args=dict(bins=80),plot_args=dict(marker=',', linestyle='none', color='black'),contour_args=dict(cmap=pl.cm.bone))

pl.scatter(zjallcand,jw1allcand,marker='.',s=10,color='0.7')
#pl.scatter(zj,jw1,marker='.',s=0.8,color='0.8')
#pl.scatter(zjcand,jw1cand,marker='.',s=20.0,color='k')
#pl.scatter(zjcand65,jw1cand65,marker='.',s=20.0,color='k')
#pl.scatter(zjcand65[indobs],jw1cand65[indobs],marker='.',s=20.0,color='r')

## plot dwarfs
pl.scatter(mzs[indm]-mJs[indm],mJs[indm]-W1s[indm],marker='o',s=50,edgecolor='k',facecolor='none')
pl.scatter(mzs[indl]-mJs[indl],mJs[indl]-W1s[indl],marker='^',s=30,edgecolor='k',facecolor='none')
pl.scatter(mzs[indt]-mJs[indt],mJs[indt]-W1s[indt],marker='s',s=30,edgecolor='k',facecolor='none')

pl.scatter(mzq[indz6q]-mJq[indz6q],mJq[indz6q]-W1q[indz6q],marker='*',s=100,edgecolor='b',facecolor='none')
pl.scatter(mzq[indz65q]-mJq[indz65q],mJq[indz65q]-W1q[indz65q],marker='*',s=100,edgecolor='orange',facecolor='orange')
pl.scatter(mzq[indz7q]-mJq[indz7q],mJq[indz7q]-W1q[indz7q],marker='*',s=100,edgecolor='r',facecolor='r')

pl.text(1.4,0.7,'M',fontsize=18,color='k')
pl.text(2.6,2.85,'L',fontsize=18,color='k')
pl.text(3.0,0.9,'T',fontsize=18,color='k')

pl.text(0.2,3.1,'$z\lesssim6.5$',fontsize=16,color='k')
pl.text(1.6,3.1,'$z\gtrsim6.5$',fontsize=16,color='k')
pl.text(3.1,3.1,'$z>7$',fontsize=16,color='k')

## Plot New Quasars

pl.scatter([0.6836478089788969,0.7317010881496273],[2.6281350000000003,2.72795],marker='*',s=300,edgecolor='b',facecolor='none')
pl.scatter([2.2366994383733036],[2.3708101373291015],marker='*',s=300,edgecolor='orange',facecolor='orange')

pl.plot(szj,sjw1,color='g',lw=1.5)
pl.scatter(szj,sjw1,marker='o',s=10.0,color='g')
pl.scatter([szj[5],szj[10],szj[15]],[sjw1[5],sjw1[10],sjw1[15]],marker='o',s=30.0,color='g')

#pl.scatter(szjsdss,sjw1,marker='.',s=50.0,color='g')
#pl.plot(szjsdss,sjw1,'g-',lw=1.5)

# z~6
pl.plot([-1,1.5,1.5,-1],[1.5,1.5,3.,3.],'--',color='b',lw=2)

# z~6.5
pl.plot([1.5,2.0,2.5,1.5,1.5],[1.5,1.5,3.,3.,1.5],'--',color='orange',lw=2)

# z~7
pl.plot([3.0,4.0,4.0,3.0,3.0],[1.5,1.5,3.,3.,1.5],'--',color='red',lw=2)

# J-W1>1.5 & J-W1<3.0
#pl.plot([-10,2.0],[1.5,1.5],'--',color='gray',lw=2)
#pl.plot([-10,2.5],[3.0,3.0],'--',color='gray',lw=2)

# J-W1>3.0*(z-J-1.5)
#pl.plot([2.0,2.5],[1.5,3.0],'--',color='gray',lw=2)
# z-J<1.5
#pl.plot([1.5,1.5],[1.5,3.0],'--',color='gray',lw=2)

## z>7
#pl.plot([3,3],[1.5,3.0],'--',color='gray',lw=2)
#pl.plot([3,5],[1.5,1.5],'--',color='gray',lw=2)
#pl.plot([3,5],[3.0,3.0],'--',color='gray',lw=2)

#pl.fill([-1,1.5,1.5,-1],[1.5,1.5,3.,3.],color='0.5',alpha=0.3)
#pl.fill([1.5,2.0,2.5,1.5],[1.5,1.5,3.,3.],color='0.5',alpha=0.5)
#pl.fill([3.0,4.0,4.0,3.0],[1.5,1.5,3.,3.],color='0.5',alpha=0.8)


# J-W1>1.3
#pl.plot([-10,1.43],[1.2,1.2],'k--')
# J-w1>3(z-J-1)
#pl.plot([1.4,2],[1.2,3],'k--')

pl.xlim([0.1,3.8])
pl.ylim([-0.1,4.2])
pl.xlabel('$z_{AB}-J_{VEGA}$',fontsize=14)
pl.ylabel('$J_{VEGA}-W1_{VEGA}$',fontsize=14)
pl.savefig('../zjjw1.eps')
pl.show()


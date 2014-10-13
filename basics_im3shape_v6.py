import numpy as np
import pylab as pl
from scipy.stats import binned_statistic

def avg_err(x, values, bins):
    avg = binned_statistic(x, values, bins=bins, statistic='mean')[0]
    std = binned_statistic(x, values, bins=bins, statistic='std')[0]
    N = binned_statistic(x, values, bins=bins, statistic='count')[0]
    return avg, std/np.sqrt(N)

def titles_avg(tiles, ra, dec, e1, e2, ofile):
    unique_tiles = np.unique(tiles)

    tra, tdec, tN, te1, te2, te1e, te2e = [], [], [], [], [], [], []    
    for t in unique_tiles:
        #print t
        c = (tiles == t)
        ra1 = ra[c]
        dec1 = dec[c]
        e11 = e1[c]
        e22 = e2[c]

        tra.append(ra1.mean())
        tdec.append(dec1.mean())
        tN.append(ra1.shape[0])

        te1.append(e11.mean())
        te2.append(e22.mean())
        te1e.append(e11.std()/np.sqrt(e11.shape[0]))
        te2e.append(e22.std()/np.sqrt(e22.shape[0]))

    np.savez(ofile, tiles=unique_tiles, ra=tra, dec=tdec, N=tN, e1=te1, e2=te2, e1e=te1e, e2e=te2e)

def plot_tiles(tile_fname):
    #<e1> and <e2> in each tiles
    #np.unique(im3['tilename']).shape = 672
    ftiles = np.load(tile_fname)
    tiles = ftiles['tiles']
    tra = ftiles['ra']
    tdec = ftiles['dec']
    te1 = ftiles['e1']
    te2 = ftiles['e2']
    te1e = ftiles['e1e']
    te2e = ftiles['e2e']
    tN = ftiles['N']
 
    pl.figure(2, figsize=(7,7))
    ax = pl.subplot(111)
    pl.hist(te1, 500, label='e1', color='k', histtype='step')
    pl.hist(te2, 500, label='e2', color='r', histtype='step')
    pl.legend(loc=0)
    pl.xlim([-0.05, 0.05])
    pl.text(0.05, 0.6, r'$\langle e1 \rangle =%.3e \pm %.3e$'%(np.mean(te1), np.std(te1)), transform=ax.transAxes, size=12)
    pl.text(0.05, 0.5, r'$\langle e2 \rangle =%.3e \pm %.3e$'%(np.mean(te2), np.std(te2)), transform=ax.transAxes, size=12)
    pl.xlabel(r'$\langle e \rangle$ in tiles')
    pl.ylabel('Number')

    pl.figure(3)
    unique_tiles = tiles.shape[0]
    pl.errorbar(range(unique_tiles), te1, te1e, c='k', label='e1', ls='')
    pl.errorbar(range(unique_tiles), te2, te2e, c='r', label='e2', ls='')
    pl.legend(loc=0)
    pl.xlabel('Tile number')
    pl.ylabel(r'$\langle e \rangle$')

    pl.figure(4, figsize=(15, 7))
    pl.subplot(231)
    pl.scatter(tra, tdec, c=tN, edgecolor=None, marker='s', s=30)
    pl.colorbar()
    pl.title('Number of galaxies in tiles')
    pl.axis([60, 93, -62, -42])
    pl.ylabel('DEC')
    pl.xticks(visible=False)

    pl.subplot(232)
    pl.scatter(tra, tdec, c=te1, edgecolors=None, vmin=-0.01, vmax=0.01, marker='s', s=30)
    pl.colorbar()
    pl.title(r'<$e_1$> in tiles')
    pl.axis([60, 93, -62, -42])
    pl.xticks(visible=False)
    pl.yticks(visible=False)

    pl.subplot(233)
    pl.scatter(tra, tdec, c=te2, edgecolor=None, vmin=-0.01, vmax=0.01, marker='s', s=30)
    pl.colorbar()
    pl.title(r'<$e_2$> in tiles')
    pl.axis([60, 93, -62, -42])
    pl.xlabel('RA')
    pl.xticks(visible=False)
    pl.yticks(visible=False)

    pl.subplot(234)
    pl.scatter(tra, tdec, c=te1e, edgecolor=None, vmin=.001, vmax=0.005, marker='s', s=30)
    pl.colorbar()
    pl.title(r'<$e_1^2$> in tiles')
    pl.axis([60, 93, -62, -42])
    pl.xlabel('RA')
    pl.ylabel('DEC')
    pl.subplot(235)
    pl.scatter(tra, tdec, c=te2e, edgecolor=None, vmin=.001, vmax=0.005, marker='s', s=30)
    pl.colorbar()
    pl.title(r'<$e_2^2$> in tiles')
    pl.axis([60, 93, -62, -42])
    pl.xlabel('RA')
    pl.yticks(visible=False)


def plot_maskfrac(maskfrac, e1, e2):
    #<e1> and <e2> as a function of mask_fraction
    mask_bins = np.linspace(0, maskfrac.max(), 10)
    maskb = (mask_bins[1:] + mask_bins[:-1])/2.

    #e1 and e2 as function of radius
    e1b, e1be = avg_err(maskfrac, e1, mask_bins)
    e2b, e2be = avg_err(maskfrac, e2, mask_bins)

    pl.figure(5)
    pl.errorbar(maskb, e1b, e1be, c='k', label='e1', ls='--')
    pl.errorbar(maskb, e2b, e2be, c='r', label='e2', ls='--')
    pl.xlabel('Mask fraction')
    pl.ylabel(r'$\langle e \rangle$')
    pl.legend(loc=0)

def plot_rad(radius, e1, e2, urad, nrad):
    #Radial bins
    rad_bins = np.concatenate(([0], np.logspace(-.09, np.log10(urad), nrad)))
    rad = (rad_bins[1:] + rad_bins[:-1])/2.

    #e1 and e2 as function of radius
    e1b, e1be = avg_err(radius, e1, rad_bins)
    e2b, e2be = avg_err(radius, e2, rad_bins)

    pl.figure(1)
    pl.errorbar(rad, e1b, e1be, c='k', label='e1')
    pl.errorbar(rad, e2b, e2be, c='r', label='e2')
    pl.xlabel('Radius (arcsec)')
    pl.ylabel(r'$\langle e \rangle$')
    pl.legend(loc=0)
    pl.show()


if __name__=='__main__':

    #Reading v6 with MODEST_CLASS=1 (i.e. galaxies only)
    im3 = fits.open('im3shape_v6_r_only_galaxies_Aug3.fits')[1].data
    #Only objects with ERROR_FLAG=0
    e = (im3['error_flag'] == 0) & (im3['snr'] > 0)
    im3 = im3[e]

    titles_avg(im3['TILENAME'], im3['RA'], im3['DEC'], im3['E1'], im3['E2'], 'tiles_all.npz')

    plot_rad(im3['RADIUS'], im3['E1'], im3['E2'])
    plot_maskfrac(im3['MASK_FRACTION'], im3['e1'], im3['e2'])


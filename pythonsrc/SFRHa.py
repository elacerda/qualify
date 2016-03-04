#!/usr/bin/python
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# 
# Calculate the k for SFR(Halpha) = k . L(Halpha)
# 
#     Lacerda@Saco - 9/Jul/2014
#     
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
import numpy as np
from pystarlight.util.base import StarlightBase
from CALIFAUtils.scripts import SFR_parametrize
from CALIFAUtils.scripts import SFR_parametrize_trapz
from CALIFAUtils.plots import add_subplot_axes

#useTrapz = True
useTrapz = False
plot = True
#plot = False
outputImgSuffix = 'pdf'

if plot:
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from matplotlib.ticker import MultipleLocator  #, MaxNLocator
    import matplotlib.gridspec as gridspec
    
    mpl.rcParams['font.size']       = 20
    mpl.rcParams['axes.labelsize']  = 20
    mpl.rcParams['axes.titlesize']  = 22
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16 
    mpl.rcParams['font.family']     = 'serif'
    mpl.rcParams['font.serif']      = 'Times New Roman'
    
    plotConf__Z = [
        dict( c = 'b', lw = 0.5),
        dict( c = 'g', lw = 0.5),
        dict( c = 'r', lw = 0.5),
        dict( c = 'y', lw = 0.5),
        dict( c = 'k', lw = 2.),
        dict( c = 'c', lw = 0.5),
    ] 

bases = [ 'Padova1994.chab', 'Padova1994.salp', 'Padova2000.chab', 'Padova2000.salp' ]
#bases = [ 'Padova2000.salp' ]

baseFile    = '/Users/lacerda/LOCAL/data/Base.bc03.h5'

if __name__ == '__main__':
    for i, b in enumerate(bases):
        base        = StarlightBase(baseFile, b, hdf5 = True)
        
        max_yr      = base.ageBase[-1]
        max_yr      = 1e7
        mask        = base.l_ssp <= 912         # Angstrom
        
        f_ssp   = base.f_ssp[:,:,mask]
        l       = base.l_ssp[mask]
        
        if useTrapz:
            SFRp = SFR_parametrize_trapz
        else:
            SFRp = SFR_parametrize
        
        qh__Zt, Nh__Zt, Nh__Z, k_SFR__Z = SFRp(f_ssp, l, base.ageBase, max_yr)
        
        print b + ':'
        for i, Z in enumerate(base.metBase):
            age98 = base.ageBase[np.where(Nh__Zt[i] / Nh__Zt[i, -1] <= 0.98)][-1]
            age99 = base.ageBase[np.where(Nh__Zt[i] / Nh__Zt[i, -1] <= 0.99)][-1] 
            print '\tZ=%.4f N_H(%.2fMyr)=%e age(N_H(98%%))=%.2f age(N_H(99%%))=%.2f Myr k_SFR=%.2f Msun/yr' % (Z, max_yr / 1e6, Nh__Z[i], age98 / 1e6, age99 / 1e6, k_SFR__Z[i])
            
        if plot is False:
            continue
        else:
            f = plt.figure()
            gs = gridspec.GridSpec(2, 2)
            f.set_size_inches(10,10)
            ax1 = plt.subplot(gs[0, 0])
            ax2 = plt.subplot(gs[0, 1])
            ax3 = plt.subplot(gs[1, :])
    
            subpos = [0.58,  0.20, 0.55, 0.35]
            subax = add_subplot_axes(ax2, subpos)
            
            cmap = plt.get_cmap('spectral_r')
            for iZ, Z in enumerate(base.metBase):
                c = cmap(float(iZ) / base.nMet)
                ax1.plot(np.ma.log10(base.ageBase), Nh__Zt[iZ, :] / 1e60,
                         c = c, lw = plotConf__Z[iZ]['lw'],  
                         label = r'Z $=\ %.2f Z_\odot$' % (Z / base.metBase[4]))
                ax2.plot(np.ma.log10(base.ageBase), Nh__Zt[iZ, :] / Nh__Zt[iZ, -1], 
                         c = c, lw = plotConf__Z[iZ]['lw'], 
                         label = r'Z $=\ %.2f Z_\odot$' % (Z / base.metBase[4]))
                subax.plot(np.ma.log10(base.ageBase), Nh__Zt[iZ, :] / Nh__Zt[iZ, -1], 
                           c = c, lw = plotConf__Z[iZ]['lw'], 
                           label = r'Z $=\ %.2f Z_\odot$' % (Z / base.metBase[4]))
                ax3.plot(np.ma.log10(base.ageBase), np.ma.log10(qh__Zt[iZ, :]),
                         c = c, lw = plotConf__Z[iZ]['lw'], 
                         label = r'Z $=\ %.2f Z_\odot$' % (Z / base.metBase[4]))
                    
            ax3.legend(loc = 1, fontsize=14, frameon=False)
             
            ax2.axhline(y = 0.95, ls = '--', c = 'k')
            subax.axhline(y = 0.95, ls = '--', c = 'k')
            
            ax2.set_ylim([0, 1.1])
            subax.set_xlim([6.4, 7.2])
            subax.set_ylim([0.80, 1.05])
    
            subax.xaxis.set_major_locator(MultipleLocator(0.5))
            subax.xaxis.set_minor_locator(MultipleLocator(0.25))
            subax.yaxis.set_major_locator(MultipleLocator(0.1))
            subax.yaxis.set_minor_locator(MultipleLocator(0.05))
            #subax.xaxis.grid(which='minor')
            #subax.yaxis.grid(which='minor')
            
            ax1.xaxis.set_major_locator(MultipleLocator(1))
            ax1.xaxis.set_minor_locator(MultipleLocator(0.25))
            ax1.yaxis.set_major_locator(MultipleLocator(1))
            ax1.yaxis.set_minor_locator(MultipleLocator(0.25))
            #ax1.xaxis.grid(which='major')
            #ax1.yaxis.grid(which='major')
            
            ax2.xaxis.set_major_locator(MultipleLocator(1))
            ax2.xaxis.set_minor_locator(MultipleLocator(0.25))
            ax2.yaxis.set_major_locator(MultipleLocator(0.25))
            ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
            #ax2.xaxis.grid(which='major')
            #ax2.yaxis.grid(which='major')
            
            ax3.xaxis.set_major_locator(MultipleLocator(1))
            ax3.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax3.yaxis.set_major_locator(MultipleLocator(2))
            ax3.yaxis.set_minor_locator(MultipleLocator(1))
            #ax3.xaxis.grid(which='minor')
            #ax3.yaxis.grid(which='minor')
    
            ax1.set_xlabel(r'$\log\ t\ [yr]$')
            ax1.set_ylabel(r'$\mathcal{N}_H(t)\ [10^{60}\ \gamma\ M_\odot{}^{-1}]$')
            ax2.set_xlabel(r'$\log\ t\ [yr]$')
            ax2.set_ylabel(r'$\mathcal{N}_H(t)/\mathcal{N}_H$')
            ax3.set_xlabel(r'$\log\ t\ [yr]$')
            ax3.set_ylabel(r'$\log\ q_H [s^{-1} M_\odot{}^{-1}]$')
                
            f.tight_layout()
            f.savefig('Nh_logt_metBase_%s.%s' % (b.replace('.', '_'), outputImgSuffix))
            plt.close(f)

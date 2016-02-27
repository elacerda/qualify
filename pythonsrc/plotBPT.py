#!/usr/bin/python
#
# Lacerda@Corrego - 11/Feb/2016
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from pycasso import fitsQ3DataCube
from matplotlib.ticker import MaxNLocator

CALIFASuperFits='/Users/lacerda/CALIFA/gal_fits/v20_q050.d15a/K0277_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.fits'
EmLinesFits='/Users/lacerda/CALIFA/rgb-gas/v20_q050.d15a/K0277_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.fits'
# Carregando arquivos FITS
K = fitsQ3DataCube(CALIFASuperFits)
K.loadEmLinesDataCube(EmLinesFits)
# Agora todos as informacoes sobre as linhas de
# emissao estao instanciadas em K.EL

# Indices dos vetores aonde estao armazenados os
# fluxos de cada linha
Ha_obs__z = K.EL.flux[K.EL.lines.index('6563'), :]
Hb_obs__z = K.EL.flux[K.EL.lines.index('4861'), :]
N2_obs__z = K.EL.flux[K.EL.lines.index('6583'), :]
O3_obs__z = K.EL.flux[K.EL.lines.index('5007'), :]
# Razao entre os fluxos de N2/Ha e O3/Hb
N2Ha__z = np.log10(N2_obs__z) - np.log10(Ha_obs__z)
O3Hb__z = np.log10(O3_obs__z) - np.log10(Hb_obs__z)

# Grafico
f = plt.figure()
f.set_size_inches(5, 4)
ax = f.gca()
plt.axis([-0.6, 0.3, -1.5, 1])
sc = ax.scatter(N2Ha__z, O3Hb__z, c = K.zoneDistance_HLR, 
                cmap = 'viridis_r', marker = 'o', s = 10, 
                alpha = 0.8, edgecolor = 'none')
ax.set_xlabel(r'$\log\ [NII]/H\alpha$', fontsize = 15)
ax.set_ylabel(r'$\log\ [OIII]/H\beta$', fontsize = 15)
cb = plt.colorbar(sc, ticks = [0, .5, 1, 1.5, 2])
cb.set_label('R [HLR]', fontsize = 15)
ax.xaxis.set_major_locator(MaxNLocator(4))
ax.yaxis.set_major_locator(MaxNLocator(4))
f.subplots_adjust(left = 0.2, right = 0.90, 
                  top = 0.95, bottom = 0.15)
f.savefig('%s-BPT.pdf' % K.califaID)
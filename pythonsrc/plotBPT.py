#!/usr/bin/python
#
# Lacerda@Corrego - 11/Feb/2016
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from pycasso import fitsQ3DataCube
from matplotlib.ticker import MaxNLocator
from CALIFAUtils.lines import Lines
from CALIFAUtils.plots import plot_text_ax

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
SNHa = Ha_obs__z / K.EL.eflux[K.EL.lines.index('6563'), :]
SNHb = Hb_obs__z / K.EL.eflux[K.EL.lines.index('4861'), :]
SNN2 = N2_obs__z / K.EL.eflux[K.EL.lines.index('6583'), :]
SNO3 = O3_obs__z / K.EL.eflux[K.EL.lines.index('5007'), :]
# Razao entre os fluxos de N2/Ha e O3/Hb
N2Ha__z = np.log10(N2_obs__z) - np.log10(Ha_obs__z)
O3Hb__z = np.log10(O3_obs__z) - np.log10(Hb_obs__z)

mask = (SNHb > 3) & (SNO3 > 3) & (SNHa > 3) & (SNN2 > 3)  

# Grafico
f = plt.figure()
f.set_size_inches(5, 4)
ax = f.gca()
plt.axis([-1, 0.5, -1, 1])
sc = ax.scatter(N2Ha__z[mask], O3Hb__z[mask], c = K.zoneDistance_HLR[mask], 
                cmap = 'viridis_r', marker = 'o', s = 30, 
                #alpha = 0.6, 
                edgecolor = 'none')
l = np.array((-1,0.5))
trans_angle = plt.gca().transData.transform_angles(np.array((np.arctan(1.01) * 180/np.pi,)),
                                               l.reshape((1, 2)))[0]
ax.set_xlabel(r'$\log\ [NII]/H\alpha$', fontsize = 15)
ax.set_ylabel(r'$\log\ [OIII]/H\beta$', fontsize = 15)
cb = plt.colorbar(sc, ticks = [0, .5, 1, 1.5, 2])
cb.set_label('R [HLR]', fontsize = 15)
L = Lines()
plot_text_ax(ax, 'S06', 0.25, 0.02, '14', 'bottom', 'left', 'k')
plot_text_ax(ax, 'K01', 0.95, 0.02, '14', 'bottom', 'right', 'k')
plot_text_ax(ax, 'CF10', 0.92, 0.98, '14', 'top', 'right', 'k', rotation = 44.62)
ax.plot(L.x['S06'], L.y['S06'], 'k--', label = 'S06')
ax.plot(L.x['K01'], L.y['K01'], 'k--', label = 'K01')
ax.plot(L.x['CF10'], L.y['CF10'], 'k--', label = 'CF10')
ax.xaxis.set_major_locator(MaxNLocator(4))
ax.yaxis.set_major_locator(MaxNLocator(4))
f.subplots_adjust(left = 0.2, right = 0.90, 
                  top = 0.95, bottom = 0.15)
f.savefig('%s-BPT.pdf' % K.califaID)
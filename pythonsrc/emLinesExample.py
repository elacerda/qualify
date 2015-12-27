import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from pycasso import fitsQ3DataCube

# Matplotlib font type and size configurations
mpl.rcParams['font.size']       = 14
mpl.rcParams['axes.labelsize']  = 14
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13 
mpl.rcParams['font.family']     = 'sans serif'
mpl.rcParams['font.serif']      = 'Times'

gal = 'K0277'
CALIFASuperFits='/Users/lacerda/CALIFA/gal_fits/v20_q050.d15a/%s_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.fits' % gal
EmLinesFits='/Users/lacerda/CALIFA/rgb-gas/v20_q050.d15a/%s_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.fits' % gal

# Carregando arquivos FITS
K = fitsQ3DataCube(CALIFASuperFits)
K.loadEmLinesDataCube(EmLinesFits)
# Agora todos as informacoes sobre as linhas de
# emissao estao instanciadas em K.EL

# Indices dos vetores aonde estao armazenados os
# fluxos de cada linha
i_Ha = K.EL.lines.index('6563')
i_Hb = K.EL.lines.index('4861')
i_O3 = K.EL.lines.index('5007')
i_N2 = K.EL.lines.index('6583')
Ha_obs__z = K.EL.flux[i_Ha, :]
Hb_obs__z = K.EL.flux[i_Hb, :]
N2_obs__z = K.EL.flux[i_N2, :]
O3_obs__z = K.EL.flux[i_O3, :]

# Razao entre os fluxos de N2/Ha e O3/Hb
N2Ha__z = np.log10(N2_obs__z) - np.log10(Ha_obs__z)
O3Hb__z = np.log10(O3_obs__z) - np.log10(Hb_obs__z)

# Grafico
f = plt.figure()
ax = f.gca()
sc = ax.scatter(N2Ha__z, O3Hb__z, c = K.zoneDistance_HLR,
           cmap = 'viridis', vmax = 2, vmin = 0,
           marker = 'o', s = 10, alpha = 0.8, edgecolor = 'none')
ax.set_xlabel(r'$\log\ [NII]/H\alpha$')
ax.set_ylabel(r'$\log\ [OIII]/H\beta$')
cb = plt.colorbar(sc)
cb.set_label('Radius [HLR]')
plt.axis([-1, 0.5, -1.5, 1])
f.savefig('%s-BPT.pdf' % K.califaID)

# There are some documentation inside the source code of the EmLinesDataCube class in 
# ../src/pycasso/fitsdatacube.py
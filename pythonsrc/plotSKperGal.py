#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
import sys
import time
import numpy as np
import argparse as ap
import CALIFAUtils as C
import matplotlib as mpl
from matplotlib import pyplot as plt
from CALIFAUtils.objects import runstats
from matplotlib.ticker import MaxNLocator
from CALIFAUtils.plots import plot_text_ax
from matplotlib.pyplot import MultipleLocator
from CALIFAUtils.plots import DrawHLRCircle
from CALIFAUtils.plots import plotOLSbisectorAxis
from CALIFAUtils.scripts import calc_running_stats
from CALIFAUtils.plots import DrawHLRCircleInSDSSImage

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# mpl.rcParams['font.size'] = 16
# mpl.rcParams['axes.labelsize'] = 16
# mpl.rcParams['axes.titlesize'] = 18
# mpl.rcParams['xtick.labelsize'] = 12
# mpl.rcParams['ytick.labelsize'] = 12 
# mpl.rcParams['font.family'] = 'serif'
# mpl.rcParams['font.serif'] = 'Times New Roman'
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'


def verify_files(K, califaID, EL = True, GP = True):
    if K is None:
        print '<<< %s galaxy: miss files' % califaID
        return False
    if EL == True and K.EL is None:
        print '<<< %s galaxy: miss EmLines files' % califaID
        return False
        if K.EL.flux[0, :].sum() == 0.:
            print '<<< %s EmLines FITS problem' % califaID
            return False
    if GP is True and K.GP._hdulist is None:
        print '<<< %s galaxy: miss gasprop file' % califaID
        return False
    # Problem in FITS file
    return True    

def print_args(args):
    for k, v in args.__dict__.iteritems():
        print k, v 

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'hdf5' : None,
        'califaID' : 'K0073',
        'itSF' : 11,
        'vrunstr' : 'v20_q050.d15a',
        'rmin' : None,
        'slice_gals' : None,
        'bamin' : 0,
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--califaID', '-g',
                        metavar = 'FILE',
                        type = str,
                        default = default['califaID'])
    parser.add_argument('--itSF', '-T',
                        help = 'age index',
                        metavar = '',
                        type = int,
                        default = default['itSF'])
    parser.add_argument('--rmin', '-R',
                        help = 'min R (HLR)',
                        metavar = 'HLR',
                        type = float,
                        default = default['rmin'])
    parser.add_argument('--vrunstr',
                        metavar = 'RUNCODE',
                        type = str,
                        default = default['vrunstr'])
    parser.add_argument('--slice_gals', '-S',
                        metavar = 'FILE',
                        type = str,
                        default = default['slice_gals'])
    parser.add_argument('--bamin', '-B',
                        help = 'min b/a',
                        metavar = '',
                        type = float,
                        default = default['bamin'])

    return parser.parse_args()

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


Zsun = 0.019

if __name__ == '__main__':
    t_init_gal = time.clock()

    args = parser_args()
    debug = args.debug
    h5fname = args.hdf5
    galName = args.califaID
    iT = args.itSF
    
    H = C.H5SFRData(args.hdf5)
    tSF__T = H.tSF__T
    ageMyr = tSF__T[iT] / 1e6
    
    if debug:
        print 'califaID: ', galName
        print 'h5fname: ', h5fname
        print 'iTSF: (%.2fMyr)' % (iT, tSF__T[iT] / 1e6)
    
    if (len(np.where(H.califaIDs == galName)[0]) == 0):
        exit('<<< plot: %s: no data.' % galName)
    
    maskradius = args.rmin
    fnamesuffix = '.pdf'
    
    if maskradius is None:
        minR = H.Rbin__r[0]
        namepref = ''
        Rbin__r = H.Rbin__r
        RbinCenter__r = H.RbinCenter__r
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        minR = maskradius
        Rbin__r = H.Rbin__r[H.Rbin__r >= maskradius]
        RbinCenter__r = H.RbinCenter__r[H.RbinCenter__r >= maskradius]
        namepref = '_R%.1f' % args.rmin
        maskRadiusOk__g = (H.zone_dist_HLR__g >= maskradius) & (H.zone_dist_HLR__g <= H.Rbin__r[-1]) 
        maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * (H.RbinCenter__r >= maskradius)).T

    if args.slice_gals is None:
        N_gals = H.N_gals
        gals_slice__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        gals_slice__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
        gals_slice__integr = np.ones(H.califaIDs_all.shape, dtype = np.bool)
        gals_txt = ''
    else:
        l_gals, _ = C.sort_gals(args.slice_gals)
        gals_slice__g, N_gals = H.get_mask_zones_list(l_gals, return_ngals = True)
        gals_slice__rg, N_gals_r = H.get_mask_radius_list(l_gals, return_ngals = True)
        gals_slice__integr, N_gals_integr = H.get_mask_integrated_list(l_gals, return_ngals = True)
        gals_txt = (args.slice_gals).split('/')[-1].split('.')[:-1]
        if len(gals_txt) > 1:
            gals_txt = '.'.join(gals_txt)
        fnamesuffix = '_%s%s' % ('_'.join(gals_txt), fnamesuffix)
    
    ba_max = args.bamin
    mask_GAL__g = np.bitwise_or(np.zeros_like(H.integrated_EW_Ha__g, dtype = np.bool), np.less(H.ba_GAL__g, ba_max))
    mask_GAL__g = np.bitwise_or(mask_GAL__g, ~gals_slice__integr)   
    
    mask__g = np.bitwise_or(np.ma.log10(H.SFRSD__Tg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Tg[iT]).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.SFRSD_Ha__g * 1e6).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.tau_V_neb__g).mask)
    mask__g = np.bitwise_or(mask__g, H.logO3N2_M13__g.mask)
    #mask__g = np.bitwise_or(mask__g, np.less(H.EW_Ha__g, 3.))
    mask__g = np.bitwise_or(mask__g, np.less(H.reply_arr_by_zones(H.ba_GAL__g), ba_max))
    mask__g = np.bitwise_or(mask__g, ~maskRadiusOk__g)
    mask__g = np.bitwise_or(mask__g, ~gals_slice__g)

    mask__rg = mask__rg = np.bitwise_or(~maskRadiusOk__rg, np.less(H.reply_arr_by_radius(H.ba_GAL__g), ba_max))
    mask__rg = np.bitwise_or(mask__rg, ~gals_slice__rg)

    #mask__rg = ~maskRadiusOk__rg
    
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    
    print NgalsOkZones, NgalsOkRbins

    # global
    xOkMin = H.xOkMin
    tauVOkMin = H.tauVOkMin
    tauVNebOkMin = H.tauVNebOkMin
    tauVNebErrMax = H.tauVNebErrMax
    RbinCenter__r = H.RbinCenter__r
    # ALL gal
    ##stellar
    x_Y__g = np.ma.masked_array(H.x_Y__Tg[iT], mask = mask__g)
    x_Y__rg = np.ma.masked_array(H.x_Y__Trg[iT], mask = mask__rg)
    aSFRSD__rg = np.ma.masked_array(H.aSFRSD__Trg[iT], mask = mask__rg)
    tau_V__rg = np.ma.masked_array(H.tau_V__Trg[iT], mask = mask__rg)
    ##nebular
    aSFRSD_Ha__rg = np.ma.masked_array(H.aSFRSD_Ha__Trg[iT], mask = mask__rg)
    tau_V_neb__rg = np.ma.masked_array(H.tau_V_neb__Trg[iT], mask = mask__rg)
    # one gal
    ##stellar
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # tau_V__z = getattr(H, '%s_tau_V__Tg' % galName)[iT]
    # atau_V__r = getattr(H, '%s_tau_V__Trg' % galName)[iT]
    # SFRSD__z = getattr(H, '%s_SFRSD__Tg' % galName)[iT]
    # aSFRSD__r = getattr(H, '%s_aSFRSD__Trg' % galName)[iT]
    # x_Y__z = getattr(H, '%s_x_Y__Tg' % galName)[iT]
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    tau_V__z = H.get_prop_gal(np.ma.masked_array(H.tau_V__Tg[iT], mask = mask__g), galName)
    atau_V__r = H.get_prop_gal(tau_V__rg, galName)
    SFRSD__z = H.get_prop_gal(np.ma.masked_array(H.SFRSD__Tg[iT], mask = mask__g), galName)
    aSFRSD__r = H.get_prop_gal(aSFRSD__rg, galName)
    x_Y__z = H.get_prop_gal(np.ma.masked_array(H.x_Y__Tg[iT], mask = mask__g), galName)
    x_Y__r = H.get_prop_gal(np.ma.masked_array(H.x_Y__Trg[iT], mask = mask__rg), galName)
    zoneDist_HLR__z = H.get_prop_gal(np.ma.masked_array(H.zone_dist_HLR__g, mask = mask__g), galName)
    ##nebular
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # EW_Ha__z = getattr(H, '%s_EW_Ha__g' % galName)
    # EW_Hb__z = getattr(H, '%s_EW_Hb__g' % galName)
    # tau_V_neb__z = getattr(H, '%s_tau_V_neb__g' % galName)
    # atau_V_neb__r = getattr(H, '%s_tau_V_neb__rg' % galName)
    # tau_V_neb_err__z = getattr(H, '%s_tau_V_neb_err__g' % galName)
    # SFRSD_Ha__z = getattr(H, '%s_SFRSD_Ha__g' % galName)
    # aSFRSD_Ha__r = getattr(H, '%s_aSFRSD_Ha__rg' % galName)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    EW_Ha__z = H.get_prop_gal(np.ma.masked_array(H.EW_Ha__g, mask = mask__g), galName)
    EW_Hb__z = H.get_prop_gal(np.ma.masked_array(H.EW_Hb__g, mask = mask__g), galName)
    tau_V_neb__z = H.get_prop_gal(np.ma.masked_array(H.tau_V_neb__g, mask = mask__g), galName)
    atau_V_neb__r = H.get_prop_gal(tau_V_neb__rg, galName)
    tau_V_neb_err__z = H.get_prop_gal(np.ma.masked_array(H.tau_V__Tg[iT], mask = mask__g), galName)
    SFRSD_Ha__z = H.get_prop_gal(np.ma.masked_array(H.tau_V__Tg[iT], mask = mask__g), galName)
    aSFRSD_Ha__r = H.get_prop_gal(aSFRSD_Ha__rg, galName)
    F_obs_Ha__z = H.get_prop_gal(np.ma.masked_array(H.F_obs_Ha__g, mask = mask__g), galName)
    
    imgDir = C.paths.califa_work_dir + 'images/'
    K = C.read_one_cube(galName, EL = True, v_run = args.vrunstr)
    galaxyImgFile = imgDir + galName + '.jpg'
    
    # Setup elliptical-rings geometry
    pa, ba = K.getEllipseParams()
    K.setGeometry(pa, ba)
    
    tipos, tipo, tipo_m, tipo_p = C.get_morfologia(galName)
    
    #stellar
    tau_V__yx = K.zoneToYX(tau_V__z, extensive = False)
    SFRSD__yx = K.zoneToYX(SFRSD__z, extensive = False)
    x_Y__yx = K.zoneToYX(x_Y__z, extensive = False)
    #nebular
    tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
    SFRSD_Ha__yx = K.zoneToYX(SFRSD_Ha__z, extensive = False)
    EW_Ha__yx = K.zoneToYX(EW_Ha__z, extensive = False)
    EW_Hb__yx = K.zoneToYX(EW_Hb__z, extensive = False)
    tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z, extensive = False)
    F_obs_Ha__yx = K.zoneToYX(F_obs_Ha__z, extensive = True)
    #mixed
    deltaTau__z = tau_V_neb__z - tau_V__z 
    deltaTau__yx = K.zoneToYX(deltaTau__z, extensive = False)
    
    t_calc = time.clock()
    print 'calc: elapsed time: %.2f' % (t_calc - t_init_gal)
    
    NRows = 3
    NCols = 4
    
    default_kwargs_imshow = dict(origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'viridis')
    
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     f, axArr = plt.subplots(NRows, NCols)
#     f.set_size_inches((NCols * 5.34, NRows * 5.))
#     
#     for ax in f.axes:
#         ax.set_axis_off()
#     
#     age = tSF__T[iT]
#     
#     tauvmin, tauvmax = -1.5, 0.5
#     SFRSDmin, SFRSDmax = -3.5, 0 
# 
#     ax = axArr[0, 0]
#     ax.set_axis_on()
#     galimg = plt.imread(galaxyImgFile)[::-1, :, :]
#     plt.setp(ax.get_xticklabels(), visible = False)
#     plt.setp(ax.get_yticklabels(), visible = False)
#     ax.imshow(galimg, origin = 'lower')
#     DrawHLRCircleInSDSSImage(ax, K.HLR_pix, pa, ba)
#     
#     ax = axArr[0, 1]
#     ax.set_axis_on()
#     xlabel = r'EW(H$\alpha$) [$\AA$]'
#     im = ax.imshow(EW_Ha__yx, vmax = 20, vmin = 3, **default_kwargs_imshow)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     sigma_dev = 1
#     ax = axArr[0, 2]
#     ax.set_axis_on()
#     xlabel = r'$F_{obs}(H\alpha)$ [erg cm${}^{-2}$ s${}^{-1}$]'
#     vmax = F_obs_Ha__yx.mean() + sigma_dev * F_obs_Ha__yx.std()
#     vmin = 0
#     im = ax.imshow(F_obs_Ha__yx, vmax = vmax, vmin = vmin, **default_kwargs_imshow)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     ax = axArr[0, 3]
#     ax.set_axis_on()
#     xlabel = r'$\epsilon\tau_V^{neb}$'
#     kwargs_imshow = default_kwargs_imshow.copy()
#     im = ax.imshow(tau_V_neb_err__yx, vmax = 1, **default_kwargs_imshow)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     default_rs_kwargs = dict(smooth = True, sigma = 1.2)
#     
#     ax = axArr[1, 0]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V__rg.flatten())
#     y = np.ma.log10(aSFRSD__rg.flatten() * 1e6)
#     xm, ym = C.ma_mask_xyz(x, y)
#     xlabel = r'$\log\ \tau_V^{\star}(R)$'
#     ylabel = r'$\log\ \langle \Sigma_{SFR}^\star(t_\star, R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$'
#     xlim = [tauvmin, tauvmax]
#     ylim = [SFRSDmin, SFRSDmax]
#     sc = ax.scatter(xm, ym, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.3)
#     rs_kwargs = default_rs_kwargs.copy()
#     rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs) 
#     ax.plot(rs.xS, rs.yS, 'k', lw = 2)
#     ax.plot(rs.xPrc[0], rs.yPrc[0], 'k--', lw = 2)
#     ax.plot(rs.xPrc[1], rs.yPrc[1], 'k--', lw = 2)
#     ax.plot(rs.xPrc[2], rs.yPrc[2], 'k--', lw = 2)
#     ax.plot(rs.xPrc[3], rs.yPrc[3], 'k--', lw = 2)
#     a_ols, b_ols, sigma_a_ols, sigma_b_ols = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.98, pos_y = 0.08, fs = 14)
#     rs.OLS_bisector()
#     plotOLSbisectorAxis(ax, rs.OLS_median_slope, rs.OLS_median_intercept, 
#                         pos_x = 0.98, pos_y = 0.02, fs = 14, c = 'b', OLS = True, x_rms = xm, y_rms = ym)
#     ##########################
#     x = np.ma.log10(atau_V__r)
#     y = np.ma.log10(aSFRSD__r * 1e6)
#     xm, ym = C.ma_mask_xyz(x, y)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
#     # vmax = xr.mean() + 2. * xr.std()
#     # vmin = xr.mean() - 2. * xr.std()
#     # ax.scatter(xm, ym, c = xr, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = vmax, vmin = vmin)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     sc = ax.scatter(xm, ym, c = H.RbinCenter__r, cmap = 'viridis', marker = 'o', s = 50, edgecolor = 'black', vmax = 2, vmin = minR)
#     cb = f.colorbar(ax = ax, mappable = sc, use_gridspec = True)
#     cb.set_label(r'R [HLR]')
#     ##########################
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
#     ax.xaxis.set_major_locator(MultipleLocator(0.5))
#     ax.xaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.yaxis.set_major_locator(MultipleLocator(0.5))
#     ax.yaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.grid(which = 'major')
#     
#     sigma_dev = 3.
# 
#     ax = axArr[1, 1]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V__yx)
#     y = np.ma.log10(SFRSD__yx * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
#     ylim = [SFRSDmin, SFRSDmax]
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # vmax = y.mean() + sigma_dev * y.std()
#     # vmin = y.mean() - sigma_dev * y.std()
#     # im = ax.imshow(ym, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin, vmax = vmax)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     im = ax.imshow(y, vmin = SFRSDmin, vmax = SFRSDmax, **default_kwargs_imshow)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     txt = 'not masked: %d' % (~mask).sum()
#     plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
#     
#     ax = axArr[1, 2]
#     ax.set_axis_on()
#     xlabel = r'$\log\ \tau_V^\star$' 
#     im = ax.imshow(xm, vmax = tauvmax, vmin = tauvmin, **default_kwargs_imshow)
#     #im = ax.imshow(x, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     ax = axArr[1, 3]
#     ax.set_axis_on()
#     xlabel = r'$x_Y [\%]$'
#     vmin = xOkMin * 100.
#     im = ax.imshow(100. * x_Y__yx, vmin = vmin, vmax = 50, **default_kwargs_imshow)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
#     
#     ax = axArr[2, 0]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V_neb__rg.flatten())
#     y = np.ma.log10(aSFRSD_Ha__rg.flatten() * 1e6)
#     xm, ym = C.ma_mask_xyz(x, y)
#     xlabel = r'$\log\ \tau_V^{neb}(R)$'
#     ylabel = r'$\log\ \langle \Sigma_{SFR}^{neb}(R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$' 
#     xlim = [tauvmin, tauvmax]
#     ylim = [SFRSDmin, SFRSDmax]
#     sc = ax.scatter(xm, ym, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.3)
#     rs_kwargs = default_rs_kwargs.copy()
#     rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs) 
#     ax.plot(rs.xS, rs.yS, 'k', lw = 2)
#     ax.plot(rs.xPrc[0], rs.yPrc[0], 'k--', lw = 2)
#     ax.plot(rs.xPrc[1], rs.yPrc[1], 'k--', lw = 2)
#     ax.plot(rs.xPrc[2], rs.yPrc[2], 'k--', lw = 2)
#     ax.plot(rs.xPrc[3], rs.yPrc[3], 'k--', lw = 2)
#     a_ols, b_ols, sigma_a_ols, sigma_b_ols = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.98, pos_y = 0.08, fs = 14)
#     rs.OLS_bisector()
#     plotOLSbisectorAxis(ax, rs.OLS_median_slope, rs.OLS_median_intercept, 
#                         pos_x = 0.98, pos_y = 0.02, fs = 14, c = 'b', OLS = True, x_rms = xm, y_rms = ym)
#     ##########################
#     x = np.ma.log10(atau_V_neb__r)
#     y = np.ma.log10(aSFRSD_Ha__r * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
#     # vmax = xr.mean() + 2. * xr.std()
#     # vmin = xr.mean() - 2. * xr.std()
#     # ax.scatter(xm, ym, c = xr, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = vmax, vmin = vmin)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     sc = ax.scatter(xm, ym, c = H.RbinCenter__r, cmap = 'viridis', marker = 'o', s = 50, edgecolor = 'black', vmax = 2, vmin = minR)
#     cb = f.colorbar(ax = ax, mappable = sc, use_gridspec = True)
#     cb.set_label(r'R [HLR]')
#     ##########################
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
#     ax.xaxis.set_major_locator(MultipleLocator(0.5))
#     ax.xaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.yaxis.set_major_locator(MultipleLocator(0.5))
#     ax.yaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.grid(which = 'major')
# 
#     ax = axArr[2, 1]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V_neb__yx)
#     y = np.ma.log10(SFRSD_Ha__yx * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     label = r'$\log\ \Sigma_{SFR}^{neb} [M_\odot yr^{-1} kpc^{-2}]$' 
#     ylim = [SFRSDmin, SFRSDmax]
#     im = ax.imshow(y, vmin = SFRSDmin, vmax = SFRSDmax, **default_kwargs_imshow)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     txt = 'not masked: %d' % (~mask).sum()
#     plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#     ax.set_title(label, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
#     
#     ax = axArr[2, 2]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V_neb__yx)
#     y = np.ma.log10(SFRSD_Ha__yx * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     label = r'$\log\ \tau_V^{neb}$' 
#     im = ax.imshow(xm, vmin = tauvmin, vmax = tauvmax, **default_kwargs_imshow)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(label, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     ax = axArr[2, 3]
#     ax.set_axis_on()
#     label = r'$\delta\ \tau_V$'
#     im = ax.imshow(deltaTau__yx, **default_kwargs_imshow)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     mask = deltaTau__yx.mask
#     txt = 'not masked: %d' % (~mask).sum()
#     plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#     ax.set_title(label, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # x = np.ma.log10(tau_V__Tz[iT])
#     # y = np.ma.log10(SFRSD__Tz[iT] * 1e6)
#     # mask = x.mask | y.mask
#     # xm = np.ma.masked_array(x, mask = mask)
#     # ym = np.ma.masked_array(y, mask = mask)
#     # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
#     # xr__yx = K.zoneToYX(xr, extensive = False)
#     # xlabel = r'xr'
#     # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # # vmax = xr__yx.mean() + 2. * xr__yx.std()
#     # # vmin = xr__yx.mean() - 2. * xr__yx.std()
#     # # im = ax.imshow(xr__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = vmax, vmin = vmin)
#     # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # im = ax.imshow(xr__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
#     # DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     # txt = 'Nz: %d' % K.N_zone 
#     # plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#     # txt = 'masked: %d' % mask.sum() 
#     # plot_text_ax(ax, txt, 0.02, 0.92, 14, 'top', 'left')
#     # ax.set_title(xlabel, y=-0.15)
#     # ax.grid()
#     # f.colorbar(ax = ax, mappable = im, use_gridspec = True)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# 
#     f.suptitle(r'%s - morph:%s  b/a:%.2f  age:%.2fMyr  $x_Y$(min):%.0f%%' % (galName, tipos, ba, ageMyr, xOkMin * 100.), fontsize = 24)
#     f.subplots_adjust(left = 0.07, bottom = 0.1, right = 0.99, wspace = 0.1, top = 0.9)
#     f.savefig('%s_mosaic%s' % (galName, fnamesuffix))
#     plt.close(f)
#     t_plot = time.clock()
#     print 'plot: elapsed time: %.2f' % (t_plot - t_calc)
#     print 'total: elapsed time: %.2f' % (t_plot - t_init_gal)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    f = plt.figure()
    NCols = 3
    NRows = 1
    page_size_inches = [15, 5]
    f.set_size_inches(page_size_inches)
    f.set_dpi(100)
    grid_shape = (NRows, NCols)

    ax = plt.subplot2grid(grid_shape, (0, 0))
    ax.set_axis_on()
    galimg = plt.imread(galaxyImgFile)[::-1, :, :]
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)
    ax.set_xlabel(galName, labelpad = 20)
    ax.imshow(galimg, origin = 'lower')
    DrawHLRCircleInSDSSImage(ax, K.HLR_pix, pa, ba)
    
    ax = plt.subplot2grid(grid_shape, (0, 1))
    ax.set_axis_on()
    ax.set_title(r'$\mathrm{frac\c\~ao\ de\ luz\ provenientes\ de\ populac\c\~oes\ jovens}$ ($x_Y$)', y = 1.08)
    vmin = xOkMin * 100.
    im = ax.imshow(100. * x_Y__yx, vmin = vmin, vmax = 50, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'viridis')
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    ax.set_xlim(0, 65)
    ax.set_ylim(-5, 70)
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.grid()
    cb = f.colorbar(ax = ax, mappable = im)
    cb.ax.yaxis.set_major_locator(MaxNLocator(5))
    
    ax = plt.subplot2grid(grid_shape, (0, 2))
    ax.set_axis_on()
    ylabel = r'$x_Y[\%]$'
    ax.set_xlabel('R [HLR]')
    ax.set_ylabel(ylabel, labelpad = 15)
    ax.set_xlim(minR, 2)
    ax.set_ylim(0, 50)
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    xm, ym = C.ma_mask_xyz(zoneDist_HLR__z, 100. * x_Y__z)
    rs = C.runstats(xm.compressed(), ym.compressed(), smooth = True, frac = 0.05, sigma = 1.2)
    ax.scatter(zoneDist_HLR__z, 100. * x_Y__z, marker = 'o', s = 10, edgecolor = 'none', c = '0.8', label = 'zones')
    ax.plot(RbinCenter__r, 100. * x_Y__r, 'k-', lw = 2, label = 'radial bins')
    ax.plot(rs.xS, rs.yS, 'k--', label = 'smoothed median')

    f.subplots_adjust(left = 0.03, bottom = 0.15, right = 0.95, wspace = 0.25, top = 0.85)
    f.savefig('%s_xY_radialProfile%s' % (galName, fnamesuffix))
    plt.close(f)
    

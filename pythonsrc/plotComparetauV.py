#!/opt/local/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import sys
import numpy as np
import CALIFAUtils as C
import matplotlib as mpl
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
from os.path import basename
import argparse as ap
from CALIFAUtils.lines import Lines
from matplotlib import pyplot as plt
#from matplotlib.pyplot import MultipleLocator
from matplotlib.pyplot import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
from CALIFAUtils.plots import density_contour
from CALIFAUtils.plots import plotOLSbisectorAxis, plot_text_ax

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
    default = {
        'debug' : False,
        'hdf5' : None,
        'output' : 'output.pdf',
        'maskradius' : None,
        'slice_gals' : None,
        'dryrun' : False,
        'bamin' : 0,
        'itSF' : 11,
    }
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--dryrun',
                        action = 'store_true',
                        default = default['dryrun'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--output', '-o',
                        help = 'Name of the output PDF file.',
                        metavar = 'FILENAME',
                        type = str,
                        default = default['output'])
    parser.add_argument('--slice_gals', '-S',
                        metavar = 'FILE',
                        type = str,
                        default = default['slice_gals'])
    parser.add_argument('--maskradius', '-R',
                        help = 'initial RDisc value in HLR',
                        metavar = 'NUM',
                        type = float,
                        default = default['maskradius'])
    parser.add_argument('--bamin', '-B',
                        help = 'min b/a',
                        metavar = '',
                        type = float,
                        default = default['bamin'])
    parser.add_argument('--itSF', '-T',
                        help = 'SF age index.',
                        metavar = '',
                        type = int,
                        default = default['itSF'])
    return parser.parse_args()

if __name__ == '__main__':
    args = parser_args()
    
    C.debug_var(args.debug, args = args)
    
    h5file = args.hdf5
    H = C.H5SFRData(h5file)
    
    minR = args.maskradius
    #zlim = [RNuc, 3]
    A4Size_inches = [ 8.3, 8.3 ]
    #A4Size_inches = [ 8.267, 11.692 ]
    LetterSize_inches = [ 8.5, 11 ]
    
    tSF__T = H.tSF__T
    iT = args.itSF
    tSF = tSF__T[iT]
    zlim = [minR, H.Rbin__r[-1]]
    # fix integrated_tau_V_neb mask
    aux = np.less(H.integrated_tau_V_neb__g, 0)
    H.integrated_tau_V_neb__g[aux] = np.ma.masked
    fnamesuffix = '.pdf'
    
    if minR is None:
        ticks_r = [0, 1.5, 3]
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        ticks_r = [0.7, 1.85, 3]
        maxR = H.Rbin__r[-1]
        maskRadiusOk__g = (H.zone_dist_HLR__g >= minR) & (H.zone_dist_HLR__g <= maxR) 
        maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * ((H.RbinCenter__r >= minR) & (H.RbinCenter__r <= maxR))).T
        fnamesuffix = '_maskradius%s' % fnamesuffix
        
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
    
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    
    print NgalsOkZones, NgalsOkRbins, H.N_zones__g.sum() - mask__g.astype(int).sum(), (H.N_gals * H.NRbins) - mask__rg.astype(int).sum() 
    
    axial_ratio = 0.13
    axial_ratio_sq = axial_ratio ** 2.0
    
    # There is 2 galaxies with ba < 0.14. Those galaxies will assume mu = 0
    mu_GAL = np.ma.where(H.ba_GAL__g < 0.14, 0.05, ((H.ba_GAL__g ** 2.0 - axial_ratio_sq) / (1 - axial_ratio_sq)) ** 0.5)
    mu_GAL = np.ma.where(H.ba_GAL__g < 0.13, 0.0, ((H.ba_GAL__g ** 2.0 - axial_ratio_sq) / (1 - axial_ratio_sq)) ** 0.5)
    
    alpha = 1.
    
    # zones
    mu_GAL__g = H.reply_arr_by_zones(mu_GAL)
    D_tau__g = np.ma.masked_array(H.tau_V_neb__g - H.tau_V__Tg[iT], mask = mask__g)
    R_tau__g = np.ma.masked_array(H.tau_V_neb__g / H.tau_V__Tg[iT], mask = mask__g)
    aux = 1. - (alpha * H.x_Y__Tg[iT])
    tau_BC__g = np.ma.masked_array(D_tau__g / aux, mask = mask__g)
    tau_ISM__g = np.ma.masked_array(mu_GAL__g * (H.tau_V__Tg[iT] - (alpha * H.x_Y__Tg[iT] * H.tau_V_neb__g)) / aux, mask = mask__g)
    rho__g = np.ma.masked_array(tau_BC__g / tau_ISM__g, mask = mask__g)
    D_tau_xY__g = np.ma.masked_array(aux * tau_BC__g, mask = mask__g)
    tau_O__g = np.ma.masked_array(tau_ISM__g / mu_GAL__g, mask = mask__g)
    tau_Y__g = np.ma.masked_array(tau_O__g + alpha * tau_BC__g, mask = mask__g)
    aux = mu_GAL__g * alpha * rho__g
    R_tau_xY__g = np.ma.masked_array((1 + aux) / (1 + H.x_Y__Tg[iT] * aux), mask = mask__g)

    # radial bins (HLR)
    mu_GAL__rg = H.reply_arr_by_radius(mu_GAL)
    D_tau__rg = np.ma.masked_array(H.tau_V_neb__Trg[iT] - H.tau_V__Trg[iT], mask = mask__rg)
    R_tau__rg = np.ma.masked_array(H.tau_V_neb__Trg[iT] / H.tau_V__Trg[iT], mask = mask__rg)
    aux = 1. - (alpha * H.x_Y__Trg[iT])
    tau_BC__rg = np.ma.masked_array(D_tau__rg / aux, mask = mask__rg)
    tau_ISM__rg = np.ma.masked_array(mu_GAL__rg * (H.tau_V__Trg[iT] - (alpha * H.x_Y__Trg[iT] * H.tau_V_neb__Trg[iT])) / aux, mask = mask__rg)
    rho__rg = np.ma.masked_array(tau_BC__rg / tau_ISM__rg, mask = mask__rg)
    D_tau_xY__rg = np.ma.masked_array(aux * tau_BC__rg, mask = mask__rg)
    tau_O__rg = np.ma.masked_array(tau_ISM__rg / mu_GAL__rg, mask = mask__rg)
    tau_Y__rg = np.ma.masked_array(tau_O__rg + alpha * tau_BC__rg, mask = mask__rg)
    aux = mu_GAL__rg * alpha * rho__rg
    R_tau_xY__rg = np.ma.masked_array((1 + aux) / (1 + H.x_Y__Trg[iT] * aux), mask = mask__rg)
 
    # integrated
    integrated_D_tau = np.ma.masked_array(H.integrated_tau_V_neb__g - H.integrated_tau_V__g, mask = mask_GAL__g)
    integrated_R_tau = np.ma.masked_array(H.integrated_tau_V_neb__g / H.integrated_tau_V__g, mask = mask_GAL__g)
    aux = 1. - (alpha * H.integrated_x_Y__Tg[iT])
    aux2 = (alpha * H.integrated_x_Y__Tg[iT] * H.integrated_tau_V_neb__g)
    integrated_tau_BC = np.ma.masked_array(integrated_D_tau / aux, mask = mask_GAL__g)
    integrated_tau_ISM = np.ma.masked_array(mu_GAL * (H.integrated_tau_V__g - aux2) / aux, mask = mask_GAL__g)
    integrated_rho = np.ma.masked_array(integrated_tau_BC / integrated_tau_ISM, mask = mask_GAL__g)
    integrated_D_tau_xY = np.ma.masked_array(aux * integrated_tau_BC, mask = mask_GAL__g)
    integrated_tau_O = np.ma.masked_array(integrated_tau_ISM / mu_GAL, mask = mask_GAL__g)
    integrated_tau_Y = np.ma.masked_array(integrated_tau_O + alpha * integrated_tau_BC, mask = mask_GAL__g)
    aux = mu_GAL * alpha * integrated_rho
    integrated_R_tau_xY = np.ma.masked_array((1 + aux) / (1 + H.integrated_x_Y__Tg[iT] * aux), mask = mask_GAL__g)
    
    default_rs_kwargs = dict(smooth = True, sigma = 4, debug = True, frac = 0.02, gs_prc = True)
    default_sc_kwargs = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = '')
    default_plot_rs_kwargs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.9)
    default_figure_kwargs = dict(figsize = (10, 8), dpi = 100)
    default_scatter_kwargs = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.45, label = '')
    default_ols_plot_kwargs = dict(c = 'r', ls = '--', lw = 2, label = 'OLS')
    default_ols_kwargs = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 6, rms = True, text = True, kwargs_plot = default_ols_plot_kwargs)
    default_zbins_kwargs = dict(
        ols = True,
        zmask = None,
        debug = True,
        write_N = True,
        spearmanr = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        kwargs_rs = default_rs_kwargs,
        kwargs_ols = default_ols_kwargs,
        kwargs_scatter = default_sc_kwargs,
        kwargs_figure = default_figure_kwargs,
        kwargs_plot_rs = default_plot_rs_kwargs,
        kwargs_ols_plot = default_ols_plot_kwargs,
    )    
    if args.dryrun: sys.exit(1)

    #with PdfPages('CompareTauV_%.2fMyr%s%s' % ((tSF/1e6), basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix)) as pdf:
    NCols = 2
    NRows = 2
    f = plt.figure()
    page_size_inches = [10,10]
    f.set_size_inches(page_size_inches)
    f.set_dpi(100)
    grid_shape = (NRows, NCols)
    
    bins = (30, 30)
    cmap = 'Blues'
    
    ols_kwargs = default_ols_kwargs.copy()
    ols_kwargs.update(dict(
        va = 'top',
        ha = 'right', 
        pos_x = 0.98, 
        fs = 15, 
        rms = True, 
        text = True, 
        pos_y = 0.98, 
        kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
    )
    ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    xran = [0, 1.5]
    yran = [0, 3]
    xm, ym = C.ma_mask_xyz(x = H.tau_V__Tg[iT], y = H.tau_V_neb__g, mask = mask__g) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
    # X, Y = np.meshgrid(xedges, yedges)
    # im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range = [xran, yran], colors = [ 'b', 'y', 'r' ])
    kw_text = dict(pos_x = 0.01, pos_y = 0.99, fs = 8, va = 'top', ha = 'left', c = 'k')
    plot_text_ax(ax, '%s' % xm.count(), **kw_text)
    ax.scatter(xm, ym, marker = 'o', s = 10, edgecolor = 'none', c = '0.8', alpha = 0.45)
    ax.set_xlim(xran)
    ax.set_ylim(yran)
    ax.set_title('zonas')
    rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
    #for i in xrange(len(rs.xPrcS)):
    #   ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
    ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
    ax.set_xlabel(r'$\tau_V^\star$') 
    ax.set_ylabel(r'$\tau_V^{neb}$')
    rs.poly1d()
    p = [ rs.poly1d_slope, rs.poly1d_intercep ]
    ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
    plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
    ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.9, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
    plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    
    ols_kwargs.update(dict(
        va = 'top',
        ha = 'right', 
        pos_x = 0.98, 
        fs = 15, 
        rms = True, 
        text = True, 
        pos_y = 0.98, 
        kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
    )
    ax = plt.subplot2grid(grid_shape, loc = (0, 1))
    xran = [0, 1.5]
    yran = [0, 3]
    xm, ym = C.ma_mask_xyz(x = H.tau_V__Trg[iT], y = H.tau_V_neb__Trg[iT], mask = mask__rg) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
    # X, Y = np.meshgrid(xedges, yedges)
    # im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range = [xran, yran], colors = [ 'b', 'y', 'r' ])
    kw_text = dict(pos_x = 0.01, pos_y = 0.99, fs = 8, va = 'top', ha = 'left', c = 'k')
    plot_text_ax(ax, '%s' % xm.count(), **kw_text)
    ax.scatter(xm, ym, marker = 'o', s = 10, edgecolor = 'none', c = '0.8', alpha = 0.45)
    ax.set_xlim(xran)
    ax.set_ylim(yran)
    ax.set_title('bins radiais')
    rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
    #for i in xrange(len(rs.xPrcS)):
    #    ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
    ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
    ax.set_xlabel(r'$\tau_V^\star$') 
    ax.set_ylabel(r'$\tau_V^{neb}$')
    rs.poly1d()
    p = [ rs.poly1d_slope, rs.poly1d_intercep ]
    ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
    plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
    ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.9, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
    plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))

    ols_kwargs.update(dict(
        va = 'top',
        ha = 'right', 
        pos_x = 0.98, 
        fs = 15, 
        rms = True, 
        text = True, 
        pos_y = 0.98, 
        kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
    )
    ax = plt.subplot2grid(grid_shape, loc = (1, 0))
    xran = [0, 1.5]
    yran = [0, 3]
    xm, ym = C.ma_mask_xyz(x = H.integrated_tau_V__g, y = H.integrated_tau_V_neb__g, mask = mask_GAL__g) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
    # X, Y = np.meshgrid(xedges, yedges)
    # im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range = [xran, yran], colors = [ 'b', 'y', 'r' ])
    kw_text = dict(pos_x = 0.01, pos_y = 0.99, fs = 8, va = 'top', ha = 'left', c = 'k')
    plot_text_ax(ax, '%s' % xm.count(), **kw_text)
    ax.scatter(xm, ym, marker = 'o', s = 10, edgecolor = 'none', c = '0.8')
    ax.set_xlim(xran)
    ax.set_ylim(yran)
    ax.set_title('GAL')
    rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
    #for i in xrange(len(rs.xPrcS)):
    #    ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
    ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
    ax.set_xlabel(r'$\tau_V^\star$') 
    ax.set_ylabel(r'$\tau_V^{neb}$')
    rs.poly1d()
    p = [ rs.poly1d_slope, rs.poly1d_intercep ]
    ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
    plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
    ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.9, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
    plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))

    sum_Ha__g = H.sum_prop_gal(H.F_obs_Ha__g, mask_zones = mask__g)
    sum_Hb__g = H.sum_prop_gal(H.F_obs_Hb__g, mask_zones = mask__g)
    tau_V_neb_resolved__g = 2.886 * (np.ma.log(sum_Ha__g/sum_Hb__g) - np.ma.log(2.86))

    ols_kwargs.update(dict(
        va = 'top',
        ha = 'right', 
        pos_x = 0.98, 
        fs = 15, 
        rms = True, 
        text = True, 
        pos_y = 0.98, 
        kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
    )
    ax = plt.subplot2grid(grid_shape, loc = (1, 1))
    xran = [0, 1.5]
    yran = [0, 3]
    xm, ym = C.ma_mask_xyz(x = H.integrated_tau_V__g, y = tau_V_neb_resolved__g, mask = mask_GAL__g) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
    # X, Y = np.meshgrid(xedges, yedges)
    # im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range = [xran, yran], colors = [ 'b', 'y', 'r' ])
    kw_text = dict(pos_x = 0.01, pos_y = 0.99, fs = 8, va = 'top', ha = 'left', c = 'k')
    plot_text_ax(ax, '%s' % xm.count(), **kw_text)
    ax.scatter(xm, ym, marker = 'o', s = 10, edgecolor = 'none', c = '0.8')
    ax.set_xlim(xran)
    ax.set_ylim(yran)
    ax.set_title('GAL RES.')
    rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
    #for i in xrange(len(rs.xPrcS)):
    #    ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
    ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
    ax.set_xlabel(r'$\tau_V^\star$') 
    ax.set_ylabel(r'$\tau_V^{neb}$')
    rs.poly1d()
    p = [ rs.poly1d_slope, rs.poly1d_intercep ]
    ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
    plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
    ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.9, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
    plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    
    f.subplots_adjust(bottom = 0.1, top = 0.95, hspace = 0.3, wspace = 0.25, right = 0.95, left = 0.1)
    #pdf.savefig(f)
    f.savefig('CompareTauV_%.2fMyr%s%s' % ((tSF/1e6), basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix))
    plt.close(f)

        
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         NCols = 3
#         NRows = 2
#         f = plt.figure(figsize = (15, 8), dpi = 100)
#         grid_shape = (NRows, NCols)
# 
#         bins = (30, 30)
#         cmap = 'Blues'
#         
#         ols_kwargs = default_ols_kwargs.copy()
#         ols_kwargs.update(dict(
#             va = 'top',
#             ha = 'right', 
#             pos_x = 0.98, 
#             fs = 14, 
#             rms = True, 
#             text = True, 
#             pos_y = 0.98, 
#             kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
#         )
#         ax = plt.subplot2grid(grid_shape, loc = (0, 0))
#         xran = [0, 1.5]
#         yran = [0, 3]
#         xm, ym = C.ma_mask_xyz(x = H.tau_V__Tg[iT], y = H.tau_V_neb__g, mask = mask__g) 
#         h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
#         X, Y = np.meshgrid(xedges, yedges)
#         im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
#         ax.set_xlim(xran)
#         ax.set_ylim(yran)
#         rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
#         for i in xrange(len(rs.xPrcS)):
#             ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
#         ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
#         ax.set_xlabel(r'$\tau_V^\star$') 
#         ax.set_ylabel(r'$\tau_V^{neb}$')
#         p = np.polyfit(xm.compressed(), ym.compressed(), 1) 
#         ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
#         plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
#         ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.9, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
#         plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
#         ax.xaxis.set_major_locator(MaxNLocator(3))
#         ax.yaxis.set_major_locator(MaxNLocator(3))
# 
#         ols_kwargs.update(dict(
#             va = 'top',
#             ha = 'right', 
#             pos_x = 0.98, 
#             fs = 10, 
#             rms = True, 
#             text = True, 
#             pos_y = 0.98, 
#             kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
#         )
#         ax = plt.subplot2grid(grid_shape, loc = (0, 1))
#         ax.set_title('zonas')
#         xran = [0, 1]
#         yran = [0, 10]
#         xm, ym = C.ma_mask_xyz(x = H.x_Y__Tg[iT], y = R_tau__g, mask = mask__g) 
#         h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
#         X, Y = np.meshgrid(xedges, yedges)
#         im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
#         rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
#         for i in xrange(len(rs.xPrcS)):
#             ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
#         ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
#         ax.set_xlabel(r'$x_Y$')
#         ax.set_ylabel(r'$\mathcal{R}_\tau$') 
#         ax.set_xlim(xran)
#         ax.set_ylim(yran)
#         p = np.polyfit(xm.compressed(), ym.compressed(), 1) 
#         ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
#         plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
#         ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.93, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
#         plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
#         ax.xaxis.set_major_locator(MaxNLocator(4))
#         ax.yaxis.set_major_locator(MaxNLocator(3))
#         ax.axhline(y = 1, ls = '-.', c = 'k')
# 
#         ols_kwargs.update(dict(
#             va = 'top',
#             ha = 'right', 
#             pos_x = 0.98, 
#             fs = 10, 
#             rms = True, 
#             text = True, 
#             pos_y = 0.98, 
#             kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
#         )
#         ax = plt.subplot2grid(grid_shape, loc = (0, 2))
#         xran = [0, 1]
#         yran = [-1, 2]
#         xm, ym = C.ma_mask_xyz(x = H.x_Y__Tg[iT], y = D_tau__g, mask = mask__g) 
#         h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
#         X, Y = np.meshgrid(xedges, yedges)
#         im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
#         ax.set_xlim(xran)
#         ax.set_ylim(yran)
#         rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
#         for i in xrange(len(rs.xPrcS)):
#             ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
#         ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
#         ax.set_xlabel(r'$x_Y$')
#         ax.set_ylabel(r'$\mathcal{D}_\tau$') 
#         p = np.polyfit(xm.compressed(), ym.compressed(), 1) 
#         ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
#         plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
#         ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.93, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
#         plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
#         ax.xaxis.set_major_locator(MaxNLocator(4))
#         ax.yaxis.set_major_locator(MaxNLocator(3))
#         ax.axhline(y = 0, ls = '-.', c = 'k')
#         
#         ols_kwargs.update(dict(
#             va = 'top',
#             ha = 'right', 
#             pos_x = 0.98, 
#             fs = 10, 
#             rms = True, 
#             text = True, 
#             pos_y = 0.98, 
#             kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
#         )
#         ax = plt.subplot2grid(grid_shape, loc = (1, 0))
#         xran = [0, 1.5]
#         yran = [0, 3]
#         xm, ym = C.ma_mask_xyz(x = H.tau_V__Trg[iT], y = H.tau_V_neb__Trg[iT], mask = mask__rg) 
#         h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
#         X, Y = np.meshgrid(xedges, yedges)
#         im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
#         ax.set_xlim(xran)
#         ax.set_ylim(yran)
#         rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
#         for i in xrange(len(rs.xPrcS)):
#             ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
#         ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
#         ax.set_xlabel(r'$\tau_V^\star$') 
#         ax.set_ylabel(r'$\tau_V^{neb}$')
#         p = np.polyfit(xm.compressed(), ym.compressed(), 1) 
#         ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
#         plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
#         ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.93, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
#         plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
#         ax.xaxis.set_major_locator(MaxNLocator(3))
#         ax.yaxis.set_major_locator(MaxNLocator(3))
# 
#         ols_kwargs.update(dict(
#             va = 'top',
#             ha = 'right', 
#             pos_x = 0.98, 
#             fs = 10, 
#             rms = True, 
#             text = True, 
#             pos_y = 0.98, 
#             kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
#         )
#         ax = plt.subplot2grid(grid_shape, loc = (1, 1))
#         ax.set_title('bins radiais')
#         xran = [0, 1]
#         yran = [0, 10]
#         xm, ym = C.ma_mask_xyz(x = H.x_Y__Trg[iT], y = R_tau__rg, mask = mask__rg) 
#         h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
#         X, Y = np.meshgrid(xedges, yedges)
#         im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
#         ax.set_xlim(xran)
#         ax.set_ylim(yran)
#         rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
#         for i in xrange(len(rs.xPrcS)):
#             ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
#         ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
#         ax.set_xlabel(r'$x_Y$')
#         ax.set_ylabel(r'$\mathcal{R}_\tau$') 
#         p = np.polyfit(xm.compressed(), ym.compressed(), 1) 
#         ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
#         plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
#         ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.93, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
#         plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
#         ax.xaxis.set_major_locator(MaxNLocator(4))
#         ax.yaxis.set_major_locator(MaxNLocator(3))
#         ax.axhline(y = 1, ls = '-.', c = 'k')
# 
#         ols_kwargs.update(dict(
#             va = 'top',
#             ha = 'right', 
#             pos_x = 0.98, 
#             fs = 10, 
#             rms = True, 
#             text = True, 
#             pos_y = 0.98, 
#             kwargs_plot = dict(c = 'k', ls = '--', lw = 2)),
#         )
#         ax = plt.subplot2grid(grid_shape, loc = (1, 2))
#         xran = [0, 1]
#         yran = [-1, 2]
#         xm, ym = C.ma_mask_xyz(x = H.x_Y__Trg[iT], y = D_tau__rg, mask = mask__rg) 
#         h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
#         X, Y = np.meshgrid(xedges, yedges)
#         im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
#         ax.set_xlim(xran)
#         ax.set_ylim(yran)
#         rs = C.runstats(xm.compressed(), ym.compressed(), **default_rs_kwargs)
#         for i in xrange(len(rs.xPrcS)):
#             ax.plot(rs.xPrc[i], rs.yPrc[i], 'k--', lw = 1)
#         ax.plot(rs.xS, rs.yS, 'k-', lw = 2)
#         ax.set_xlabel(r'$x_Y$')
#         ax.set_ylabel(r'$\mathcal{D}_\tau$') 
#         p = np.polyfit(xm.compressed(), ym.compressed(), 1) 
#         ols_kwargs.update(dict(c = 'k', label = 'poly1d', OLS = True, x_rms = xm.compressed()))
#         plotOLSbisectorAxis(ax, p[0], p[1], **ols_kwargs)
#         ols_kwargs.update(OLS = None, c = 'r', label = 'OLS', pos_y = 0.93, x_rms = xm.compressed(), kwargs_plot = dict(c = 'r', ls = '--', lw = 2))    
#         plotOLSbisectorAxis(ax, xm.compressed(), ym.compressed(), **ols_kwargs)
#         ax.xaxis.set_major_locator(MaxNLocator(4))
#         ax.yaxis.set_major_locator(MaxNLocator(3))
#         ax.axhline(y = 0, ls = '-.', c = 'k')
# 
#         f.tight_layout()
#         pdf.savefig(f)
#         plt.close(f)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
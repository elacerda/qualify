#!/usr/bin/python
#
# Lacerda@UFSC_232 - 21/Jul/2015
#
import sys
import numpy as np
import argparse as ap
import CALIFAUtils as C
import matplotlib as mpl
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
from scipy import stats as st
from matplotlib import pyplot as plt
from CALIFAUtils.objects import runstats
from matplotlib.ticker import MaxNLocator
from CALIFAUtils.scripts import OLS_bisector
from matplotlib.backends.backend_pdf import PdfPages
from CALIFAUtils.plots import plot_text_ax, plotOLSbisectorAxis
from CALIFAUtils.plots import density_contour

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'itSF' : 11,
        'bamin' : 0,
        'hdf5' : None,
        'debug' : False,
        'dryrun' : False,
        'maskradius' : None,
        'slice_gals' : None,
        'radialbins' : False,
        'integrated' : False,
        'output' : 'out.pdf',
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--dryrun', '-C',
                        action = 'store_true',
                        default = default['dryrun'])
    parser.add_argument('--radialbins',
                        action = 'store_true',
                        default = default['radialbins'])
    parser.add_argument('--integrated',
                        action = 'store_true',
                        default = default['integrated'])
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
    paths = C.CALIFAPaths()
    
    H = C.H5SFRData(args.hdf5)
    iT = args.itSF
    
    minR = args.maskradius

    if minR is None:
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
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

    default_rs_kwargs = dict(smooth = True, sigma = 1.2, debug = True, frac = 0.1, gs_prc = True)
    default_sc_kwargs = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = '')
    default_plot_rs_kwargs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.9)
    default_figure_kwargs = dict(figsize = (10, 8), dpi = 100)
    default_scatter_kwargs = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.45, label = '')
    default_ols_plot_kwargs = dict(c = 'k', ls = '--', lw = 2, label = 'OLS')
    default_ols_kwargs = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 10, rms = True, text = True, kwargs_plot = default_ols_plot_kwargs)

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
    
    
    axial_ratio = 0.13
    axial_ratio_sq = axial_ratio ** 2.0
    
    # There is 2 galaxies with ba < 0.14. Those galaxies will assume mu = 0
    mu_GAL = np.ma.where(H.ba_GAL__g < 0.14, 0.05, ((H.ba_GAL__g ** 2.0 - axial_ratio_sq) / (1 - axial_ratio_sq)) ** 0.5)
    
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

    if args.dryrun or args.debug: sys.exit()     
    
    with PdfPages('zones_%s' % args.output) as pdf:
        ##########################
        ######### ZONES ##########
        ##########################
        ######### PAGE 1 #########
        ##########################
        NRows = 3
        NCols = 4
        f, axArr = plt.subplots(NRows, NCols)
        page_size_inches = (NCols * 3, NRows * 2)
        f.set_size_inches(page_size_inches)
        f.suptitle(r'$N_{gals}$: %d  $N_{zones}$: %d (inside plot: mean, median, std, max, min)' % (H.N_gals, H.N_zones__g.sum()))
        for ax in f.axes:
            ax.set_axis_off()
        props = [ 
            H.tau_V__Tg[iT],
            H.tau_V_neb__g,
            tau_ISM__g,
            tau_BC__g,
            rho__g,
            mu_GAL__g,
            tau_Y__g,
            tau_O__g,
            D_tau__g,
            R_tau__g,
            D_tau_xY__g,
            R_tau_xY__g,
        ]
        props_label = [
            r'$\tau_V^\star$',
            r'$\tau_V^{neb}$',
            r'$\tau_{ISM}$',
            r'$\tau_{BC}$',
            r'$\rho$',
            r'$\mu$',
            r'$\tau_Y$',
            r'$\tau_O$',
            r'$\mathcal{D}_\tau$',
            r'$\mathcal{R}_\tau$',
            r'$\mathcal{D}_\tau (x_Y)$',
            r'$\mathcal{R}_\tau (x_Y)$',
        ]
        props_limits = [
            [0, 2],
            [0, 3.5],
            [-0.5, 1.],
            [-1, 2],
            [0, 10],
            [0, 1],
            [0, 10],
            [0, 10],
            [-1.5, 2.0],
            [0, 10],
            [-1.5, 2.0],
            [0, 10],
        ]
        props_bins = [
            30,
            30,
            30,
            30,
            30,
            30,
            30,
            30,
            30,
            30,
            30,
            30,
        ]
        for i, p in enumerate(props):
            if isinstance(p, np.ma.core.MaskedArray):
                aux = np.ma.masked_array(p, mask = (p.mask | mask__g))
                x = aux.compressed()
            else:
                x = p[~mask__g]
            ax = f.axes[i]
            ax.set_axis_on()
            ax.hist(x, bins = props_bins[i], range = props_limits[i])
            ax.set_xlabel(props_label[i])
            ax.set_xlim(props_limits[i])
            txt = r'%.2f' % np.mean(x)
            kw_text = dict(pos_x = 0.96, pos_y = 0.96, fs = 8, va = 'top', ha = 'right', c = 'k')
            plot_text_ax(ax, txt, **kw_text)
            txt = r'%.2f' % np.median(x)
            kw_text = dict(pos_x = 0.96, pos_y = 0.88, fs = 8, va = 'top', ha = 'right', c = 'k')
            plot_text_ax(ax, txt, **kw_text)
            txt = r'%.2f' % np.std(x)
            kw_text = dict(pos_x = 0.96, pos_y = 0.80, fs = 8, va = 'top', ha = 'right', c = 'k')
            plot_text_ax(ax, txt, **kw_text)
            txt = r'%.2f' % np.max(x)
            kw_text = dict(pos_x = 0.96, pos_y = 0.72, fs = 8, va = 'top', ha = 'right', c = 'k')
            plot_text_ax(ax, txt, **kw_text)
            txt = r'%.2f' % np.min(x)
            kw_text = dict(pos_x = 0.96, pos_y = 0.64, fs = 8, va = 'top', ha = 'right', c = 'k')
            plot_text_ax(ax, txt, **kw_text)
 
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
 
        ##########################
        ######### PAGE 2 #########
        ##########################
        default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
        default_rs_kwargs = dict(smooth = True, sigma = 1.2, frac = 0.05)
        default_im_kwargs = dict(interpolation = 'nearest', origin = 'lower', aspect = 'auto', cmap = mpl.cm.spectral)
          
        f = plt.figure()
        NRows = 2
        NCols = 4
        page_size_inches = (NCols * 3, NRows * 2)
        f.set_size_inches(page_size_inches)
        f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
        for ax in f.axes:
            ax.set_axis_off()
        grid_shape = (NRows, NCols)
  
        # tau_V vs tau_V_neb
        ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
        xm, ym = C.ma_mask_xyz(H.tau_V__Tg[iT], H.tau_V_neb__g, mask = mask__g)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        density_contour(xm.compressed(), ym.compressed(), 30, 30, ax)
        ax.set_xlim([0, 2])
        ax.set_ylim([0, 3.5])
        plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
        plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
        Nprc = len(rs.xPrc)
        for i in range(Nprc):
            ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
        ax.plot(rs.xMedian, rs.yMedian, '-k')
        ax.set_xlabel(r'$\tau_V^{\star}$')
        ax.set_ylabel(r'$\tau_V^{neb}$')
        ax.grid()
  
        # tau_V_neb vs ln(EWHa/EWHb)
        ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
        y1 = np.ma.log(H.EW_Ha__g)
        y2 = np.ma.log(H.EW_Hb__g)
        y1m, y2m = C.ma_mask_xyz(y1, y2, mask = mask__g)
        ym = y1m - y2m
        xm, ym = C.ma_mask_xyz(H.tau_V_neb__g, ym, mask = mask__g)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        density_contour(xm.compressed(), ym.compressed(), 30, 30, ax)
        ax.set_xlim([0, 3.5])
        ax.set_ylim([0, 2.8])
        plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
        plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
        Nprc = len(rs.xPrc)
        for i in range(Nprc):
            ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
        ax.plot(rs.xS, rs.yS, '-k')
        xlabel = r'$\tau_V^{neb}$'
        ylabel = r'$\ln\ W_{H\alpha} / W_{H\beta}$'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid()
  
        ### End page
        f.subplots_adjust(hspace = 0.2, wspace = 0.4, bottom = 0.2)
        pdf.savefig(f)
        plt.close(f)
  
        ##########################
        ######### PAGE 3 #########
        ##########################
        f = plt.figure()
        NRows = 4
        NCols = 4
        page_size_inches = (NCols * 3, NRows * 2)
        f.set_size_inches(page_size_inches)
        f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
        for ax in f.axes:
            ax.set_axis_off()
        grid_shape = (NRows, NCols)
  
        ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
        xm, ym = C.ma_mask_xyz(H.x_Y__Tg[iT], D_tau__g, mask = mask__g)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        density_contour(xm.compressed(), ym.compressed(), 30, 30, ax)
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 2.5)
        plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
        plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
        Nprc = len(rs.xPrc)
        for i in range(Nprc):
            ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
        ax.plot(rs.xS, rs.yS, '-k')
        xlabel = r'$x_Y$ [%]'
        ylabel = r'$\mathcal{D}_\tau$'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid()
  
        ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        xm, ym, zm = C.ma_mask_xyz(H.x_Y__Tg[iT], D_tau__g, z = mu_GAL__g, mask = mask__g)
        mask = []
        mask.append(np.less_equal(zm, 0.25))
        mask.append(np.bitwise_and(np.greater(zm, 0.25), np.less_equal(zm, 0.5)))
        mask.append(np.bitwise_and(np.greater(zm, 0.5), np.less_equal(zm, 0.75)))
        mask.append(np.bitwise_and(np.greater(zm, 0.75), np.less_equal(zm, 1.)))
        zticks = [0, 0.25, 0.5, 0.75, 1.0]
        #ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 2.5)
        zcmap = mpl.cm.ScalarMappable()
        zcmap.set_cmap(mpl.cm.spectral_r)
        norm = mpl.colors.Normalize(vmin = 0., vmax = 1.)
        zcmap.set_norm(norm)
        for i, msk in enumerate(mask):
            j = i + 1
            XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
            rs = runstats(XM.compressed(), YM.compressed(), nBox = 20, **rs_kwargs)
            c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
            ax.plot(rs.xS, rs.yS, label = r'$%.2f < \mu \leq %.2f$' % (zticks[i], zticks[j]), c = c)
        ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
        xlabel = r'$x_Y$ [%]'
        ax.set_xlabel(xlabel)
        ax.grid()
  
        ax = plt.subplot2grid(grid_shape, loc = (2, 0), colspan = 2, rowspan = 2)
        xm, ym = C.ma_mask_xyz(mu_GAL__g, D_tau__g, mask = mask__g)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 2.5)
        rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        density_contour(xm.compressed(), ym.compressed(), 30, 30, ax)
        plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
        plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
        Nprc = len(rs.xPrc)
        for i in range(Nprc):
            ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
        ax.plot(rs.xS, rs.yS, '-k')
        xlabel = r'$\mu$'
        ylabel = r'$\mathcal{D}_\tau$'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid()
  
        ax = plt.subplot2grid(grid_shape, loc = (2, 2), colspan = 2, rowspan = 2)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        xm, ym, zm = C.ma_mask_xyz(z = H.x_Y__Tg[iT], y = D_tau__g, x = mu_GAL__g, mask = mask__g)
        mask = []
        mask.append(np.less_equal(zm, 0.15))
        mask.append(np.bitwise_and(np.greater(zm, 0.15), np.less_equal(zm, 0.3)))
        mask.append(np.bitwise_and(np.greater(zm, 0.3), np.less_equal(zm, 0.45)))
        mask.append(np.bitwise_and(np.greater(zm, 0.45), np.less_equal(zm, 0.6)))
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 2.5)
        zticks = [0, 0.15, 0.3, 0.45, 0.6]
        #ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        zcmap = mpl.cm.ScalarMappable()
        zcmap.set_cmap(mpl.cm.spectral_r)
        norm = mpl.colors.Normalize(vmin = 0., vmax = 0.6)
        zcmap.set_norm(norm)
        for i, msk in enumerate(mask):
            j = i + 1
            XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
            rs = runstats(XM.compressed(), YM.compressed(), nBox = 20, **rs_kwargs)
            c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
            ax.plot(rs.xS, rs.yS, label = r'$%.2f < x_Y \leq %.2f$' % (zticks[i], zticks[j]), c = c)
        ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
        xlabel = r'$\mu$'
        ax.set_xlabel(xlabel)
        ax.grid()
          
        ### End page
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
  
        ##########################
        ######### PAGE 4 #########
        ##########################
        f = plt.figure()
        NRows = 4
        NCols = 4
        page_size_inches = (NCols * 3, NRows * 2)
        f.set_size_inches(page_size_inches)
        f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
        for ax in f.axes:
            ax.set_axis_off()
        grid_shape = (NRows, NCols)
  
        ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
        xm, ym = C.ma_mask_xyz(H.x_Y__Tg[iT], R_tau__g, mask = mask__g)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        density_contour(xm.compressed(), ym.compressed(), 30, 30, ax)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 7)
        plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
        plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
        Nprc = len(rs.xPrc)
        for i in range(Nprc):
            ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
        ax.plot(rs.xS, rs.yS, '-k')
        xlabel = r'$x_Y$'
        ylabel = r'$\mathcal{R}_\tau$'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid()
  
        ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        xm, ym, zm = C.ma_mask_xyz(H.x_Y__Tg[iT], R_tau__g, z = mu_GAL__g, mask = mask__g)
        mask = []
        mask.append(np.less_equal(zm, 0.25))
        mask.append(np.bitwise_and(np.greater(zm, 0.25), np.less_equal(zm, 0.5)))
        mask.append(np.bitwise_and(np.greater(zm, 0.5), np.less_equal(zm, 0.75)))
        mask.append(np.bitwise_and(np.greater(zm, 0.75), np.less_equal(zm, 1.)))
        zticks = [0, 0.25, 0.5, 0.75, 1.0]
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 7)
        #ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        zcmap = mpl.cm.ScalarMappable()
        zcmap.set_cmap(mpl.cm.spectral_r)
        norm = mpl.colors.Normalize(vmin = 0., vmax = 1.)
        zcmap.set_norm(norm)
        for i, msk in enumerate(mask):
            j = i + 1
            XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
            rs = runstats(XM.compressed(), YM.compressed(), nBox = 20, **rs_kwargs)
            c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
            ax.plot(rs.xS, rs.yS, label = r'$%.2f < \mu \leq %.2f$' % (zticks[i], zticks[j]), c = c)
        ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
        xlabel = r'$x_Y$'
        ax.set_xlabel(xlabel)
        ax.grid()
  
        ax = plt.subplot2grid(grid_shape, loc = (2, 0), colspan = 2, rowspan = 2)
        xm, ym = C.ma_mask_xyz(mu_GAL__g, R_tau__g, mask = mask__g)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        density_contour(xm.compressed(), ym.compressed(), 30, 30, ax)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 7)
        plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
        plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
        Nprc = len(rs.xPrc)
        for i in range(Nprc):
            ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
        ax.plot(rs.xS, rs.yS, '-k')
        xlabel = r'$\mu$'
        ylabel = r'$\mathcal{R}_\tau$'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid()
  
        ax = plt.subplot2grid(grid_shape, loc = (2, 2), colspan = 2, rowspan = 2)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        sc_kwargs['s'] = 5
        sc_kwargs['alpha'] = 0.3
        xm, ym, zm = C.ma_mask_xyz(z = H.x_Y__Tg[iT], y = R_tau__g, x = mu_GAL__g, mask = mask__g)
        mask = []
        mask.append(np.less_equal(zm, 0.15))
        mask.append(np.bitwise_and(np.greater(zm, 0.15), np.less_equal(zm, 0.3)))
        mask.append(np.bitwise_and(np.greater(zm, 0.3), np.less_equal(zm, 0.45)))
        mask.append(np.bitwise_and(np.greater(zm, 0.45), np.less_equal(zm, 0.6)))
        zticks = [0, 0.15, 0.3, 0.45, 0.6]
        #ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 7)
        zcmap = mpl.cm.ScalarMappable()
        zcmap.set_cmap(mpl.cm.spectral_r)
        norm = mpl.colors.Normalize(vmin = 0., vmax = 0.6)
        zcmap.set_norm(norm)
        for i, msk in enumerate(mask):
            j = i + 1
            XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
            rs = runstats(XM.compressed(), YM.compressed(), nBox = 20, **rs_kwargs)
            c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
            ax.plot(rs.xS, rs.yS, label = r'$%.2f < x_Y \leq %.2f$' % (zticks[i], zticks[j]), c = c)
        ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
        xlabel = r'$\mu$'
        ax.set_xlabel(xlabel)
        ax.grid()
          
        ### End page
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
  
        ##########################
        ######### PAGE 5 #########
        ##########################
        f = plt.figure()
        NRows = 2
        NCols = 4
        page_size_inches = (NCols * 3, NRows * 2)
        f.set_size_inches(page_size_inches)
        f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
        for ax in f.axes:
            ax.set_axis_off()
        grid_shape = (NRows, NCols)
  
        ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
        xm, ym = C.ma_mask_xyz(x = np.ma.log10(H.tau_V_neb__g), y = np.ma.log10(H.SFRSD_Ha__g) + 6, mask = mask__g)
        rs_kwargs = default_rs_kwargs.copy()
        rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        counts, xe, ye, im = ax.hist2d(xm.compressed(), ym.compressed(), bins = 30, cmap = mpl.cm.Blues)
        density_contour(xm.compressed(), ym.compressed(), 30, 30, ax)
        ax.set_xlim([-1.5, .5])
        ax.set_ylim([-3, 0])
        plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
        plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
        ax.plot(rs.xS, rs.yS, '--*')
        ax.legend(loc = 'upper left', fontsize = 8, frameon = False)
        xlabel = r'$\tau_V^{neb}$'
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        f.colorbar(im)
        ax.grid()
      
        ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
        xm, ym = C.ma_mask_xyz(x = np.ma.log10(tau_BC__g), y = np.ma.log10(H.SFRSD_Ha__g) + 6, mask = mask__g)
        rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        counts, xe, ye, im = ax.hist2d(xm.compressed(), ym.compressed(), bins = 30, cmap = mpl.cm.Blues)
        density_contour(xm.compressed(), ym.compressed(), 30, 30, ax)
        #counts, xe, ye, im = ax.hist2d(zm.compressed(), np.ma.log10(ym).compressed(), bins = 30, range = np.asarray([[-1, 2], [3, 6.1]]), cmap = mpl.cm.Blues)
        ax.set_xlim([-1.5, .5])
        ax.set_ylim([-3, 0])
        plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
        plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
        ax.plot(rs.xS, rs.yS, '--*')
        ax.legend(loc = 'upper left', fontsize = 8, frameon = False)
        xlabel = r'$\log \tau_{BC}$' 
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        f.colorbar(im)
        ax.grid()
  
        ### End page
        f.subplots_adjust(hspace = 0.2, wspace = 0.4, bottom = 0.2)
        pdf.savefig(f)
        plt.close(f)
        
    ##########################
    ## Physical properties ###
    ##########################
    with PdfPages('zones_physprop_%s' % args.output) as pdf:
        default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
        default_rs_kwargs = dict(smooth = True, sigma = 1.2, frac = 0.05)
        default_im_kwargs = dict(interpolation = 'nearest', origin = 'lower', aspect = 'auto', cmap = mpl.cm.spectral)

        xaxis = dict(
            logMcorSD = dict(v = np.ma.log10(H.McorSD__Tg[iT]), label = r'$\log \Sigma_\star$'),
            logLintHa = dict(v = np.ma.log10(H.L_int_Ha__g), label = r'$\log L_{H\alpha}^{int}$'),
            logSFRSD = dict(v = np.ma.log10(H.SFRSD_Ha__g) + 6, label = r'$\log \Sigma_{SFR}^{H\alpha}$'),
            logSSFR = dict(v = np.ma.log10(H.SFR_Ha__g / H.Mcor__g), label = r'$\log$ SSFR'),
            #logArea = dict(v = np.ma.log10(H.zone_area_pc2__g), label = r'$\log A_{zone}\ [pc^2]$'),
            logZ = dict(v = H.logO3N2_M13__g, label = r'$12\ +\ \log (O/H)$'),
        )
        yaxis = dict(
            tauVSyn = dict(v = H.tau_V__Tg[iT], label = r'$\tau_V^\star$', lim = [0, 2]),
            tauVNeb = dict(v = H.tau_V_neb__g, label = r'$\tau_V^{neb}$', lim = [0, 3.5]),
            tauO = dict(v = tau_O__g, label = r'$\tau_O$', lim = [0, 4]),
            tauY = dict(v = tau_Y__g, label = r'$\tau_Y$', lim = [0, 4]),
            tauBC = dict(v = tau_BC__g, label = r'$\tau_{BC}$', lim = [-1, 2]),
            tauISM = dict(v = tau_ISM__g, label = r'$\tau_{ISM}$', lim = [-0.5, 1]),
        )
        yord = ['tauO', 'tauY', 'tauISM', 'tauBC', 'tauVNeb', 'tauVSyn']
        #xord = ['logMcorSD', 'logSFRSD', 'logSSFR', 'logArea', 'logLintHa', 'logZ']
        xord = ['logMcorSD', 'logSFRSD', 'logSSFR', 'logLintHa', 'logZ']
        NRows = len(yaxis)
        NCols = len(xaxis)
        grid_shape = (NRows, NCols)
        ax_size = (2, 2)
        f = plt.figure()
        page_size_inches = (grid_shape[1] * ax_size[0], grid_shape[0] * ax_size[1])
        f.set_size_inches(page_size_inches)
        row, col = 0, 0
        for yk in yord:
            for xk in xord:
                ax = plt.subplot2grid(grid_shape, loc = (row, col))
                x = xaxis[xk]
                y = yaxis[yk]
                xm, ym = C.ma_mask_xyz(x['v'], y['v'], mask = mask__g)
                rs_kwargs = default_rs_kwargs.copy()
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                rhoSpearman, pvalSpearman = st.spearmanr(rs.x, rs.y)
                p = np.polyfit(rs.xS, rs.yS, 1)
                a, b = p
                ax.hist2d(rs.x, rs.y, bins = 30, cmap = mpl.cm.Blues, range = [[rs.x.min(), rs.x.max()], y['lim']])
                kw_text = dict(pos_x = 0.98, pos_y = 0.98, fs = 10, va = 'top', ha = 'right', c = 'k')
                plot_text_ax(ax, 'R:%.3f' % rhoSpearman, **kw_text)
                kw_text['pos_y'] = 0.88
                plot_text_ax(ax, 'a:%.2f' % a, **kw_text)
                kw_text['pos_y'] = 0.78
                plot_text_ax(ax, '$y_{rms}:%.2f$' % (rs.y - (a * rs.x + b)).std(), **kw_text)
                ax.plot(rs.xS, rs.yS, '-*k')
                Nprc = len(rs.xPrc)
                ax.yaxis.set_major_locator(MaxNLocator(6, prune = 'both'))
                ax.xaxis.set_major_locator(MaxNLocator(6, prune = 'both'))
                for i in range(Nprc):
                    ax.plot(rs.xbinCenter, rs.yPrc[i], '--k')
                ax.set_ylim(y['lim'])
                if col == 0: 
                    ax.set_ylabel(y['label'])
                else:
                    plt.setp(ax.get_yticklabels(), visible = False)
                if row == (NRows - 1): 
                    ax.set_xlabel(x['label'])
                else:
                    plt.setp(ax.get_xticklabels(), visible = False)
                col += 1
            row += 1
            col = 0
        f.subplots_adjust(hspace = 0, wspace = 0)
        pdf.savefig(f)
        plt.close(f)

        if args.radialbins is True:
            with PdfPages('radialbins_%s' % args.output) as pdf:
                ##########################
                ##### #RADIAL BINS #######
                ##########################
                ######### PAGE 1 #########
                ##########################
                NRows = 3
                NCols = 4
                f, axArr = plt.subplots(NRows, NCols)
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d  $N_{zones}$: %d (inside plot: mean, median, std, max, min)' % (H.N_gals, H.N_zones__g.sum()))
                for ax in f.axes:
                    ax.set_axis_off()
                props = [ 
                    H.tau_V__Trg[iT],
                    H.tau_V_neb__Trg[iT],
                    tau_ISM__rg,
                    tau_BC__rg,
                    rho__rg,
                    mu_GAL__rg,
                    tau_Y__rg,
                    tau_O__rg,
                    D_tau__rg,
                    R_tau__rg,
                    D_tau_xY__rg,
                    R_tau_xY__rg,
                ]
                props_label = [
                    r'$\tau_V^\star$',
                    r'$\tau_V^{neb}$',
                    r'$\tau_{ISM}$',
                    r'$\tau_{BC}$',
                    r'$\rho$',
                    r'$\mu$',
                    r'$\tau_Y$',
                    r'$\tau_O$',
                    r'$\mathcal{D}_\tau$',
                    r'$\mathcal{R}_\tau$',
                    r'$\mathcal{D}_\tau (x_Y)$',
                    r'$\mathcal{R}_\tau (x_Y)$',
                ]
                props_limits = [
                    [0, 2],
                    [0, 3.5],
                    [-0.5, 3],
                    [-1, 2],
                    [0, 10],
                    [0, 1],
                    [0, 10],
                    [0, 10],
                    [-1.5, 2.0],
                    [0, 10],
                    [-1.5, 2.0],
                    [0, 10],
                ]
                props_bins = [
                    30,
                    30,
                    30,
                    30,
                    30,
                    30,
                    30,
                    30,
                    30,
                    30,
                    30,
                    30,
                ]
                for i, p in enumerate(props):
                    if isinstance(p, np.ma.core.MaskedArray):
                        aux = np.ma.masked_array(p, mask = (p.mask | mask__rg))
                        x = aux.compressed()
                    else:
                        x = p[~mask__rg]
                    ax = f.axes[i]
                    ax.set_axis_on()
                    ax.hist(x, bins = props_bins[i], range = props_limits[i])
                    ax.set_xlabel(props_label[i])
                    ax.set_xlim(props_limits[i])
                    txt = r'%.2f' % np.mean(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.96, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
                    txt = r'%.2f' % np.median(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.88, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
                    txt = r'%.2f' % np.std(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.80, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
                    txt = r'%.2f' % np.max(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.72, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
                    txt = r'%.2f' % np.min(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.64, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
        
                f.subplots_adjust(hspace = 0.4, wspace = 0.4)
                pdf.savefig(f)
                plt.close(f)
        
                default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                default_rs_kwargs = dict(smooth = True, sigma = 1.2, overlap = 0.4)
                default_im_kwargs = dict(interpolation = 'nearest', origin = 'lower', aspect = 'auto', cmap = mpl.cm.spectral)
        
                ##########################
                ######### PAGE 2 #########
                ##########################
                f = plt.figure()
                NRows = 2
                NCols = 4
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
                for ax in f.axes:
                    ax.set_axis_off()
                grid_shape = (NRows, NCols)
        
                # tau_V vs tau_V_neb
                ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(H.tau_V__Trg[iT], H.tau_V_neb__Trg[iT], mask = mask__rg)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                ax.set_xlim([0, 2])
                ax.set_ylim([0, 3.5])
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.yPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xMedian, rs.yMedian, '-k')
                ax.set_xlabel(r'$\tau_V^{\star}$')
                ax.set_ylabel(r'$\tau_V^{neb}$')
                ax.grid()
        
                # tau_V_neb vs ln(EWHa/EWHb)
                ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
                y1 = np.ma.log(H.EW_Ha__Trg[iT])
                y2 = np.ma.log(H.EW_Hb__Trg[iT])
                y1m, y2m = C.ma_mask_xyz(y1, y2, mask = mask__rg)
                ym = y1m - y2m
                xm, ym = C.ma_mask_xyz(H.tau_V_neb__Trg[iT], ym, mask = mask__rg)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.set_xlim([0, 3.5])
                ax.set_ylim([0, 2.8])
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$\tau_V^{neb}$'
                ylabel = r'$\ln\ W_{H\alpha} / W_{H\beta}$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ### End page
                f.subplots_adjust(hspace = 0.2, wspace = 0.4, bottom = 0.2)
                pdf.savefig(f)
                plt.close(f)
        
                ##########################
                ######### PAGE 3 #########
                ##########################
                f = plt.figure()
                NRows = 4
                NCols = 4
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
                for ax in f.axes:
                    ax.set_axis_off()
                grid_shape = (NRows, NCols)
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(H.x_Y__Trg[iT], D_tau__rg, mask = mask__rg)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                ax.set_xlim(0, 0.6)
                ax.set_ylim(-0.5, 2.5)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$x_Y$ [%]'
                ylabel = r'$\mathcal{D}_\tau$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                xm, ym, zm = C.ma_mask_xyz(H.x_Y__Trg[iT], D_tau__rg, z = mu_GAL__rg, mask = mask__rg)
                mask = []
                mask.append(np.less_equal(zm, 0.25))
                mask.append(np.bitwise_and(np.greater(zm, 0.25), np.less_equal(zm, 0.5)))
                mask.append(np.bitwise_and(np.greater(zm, 0.5), np.less_equal(zm, 0.75)))
                mask.append(np.bitwise_and(np.greater(zm, 0.75), np.less_equal(zm, 1.)))
                zticks = [0, 0.25, 0.5, 0.75, 1.0]
                ax.set_xlim(0, 0.6)
                ax.set_ylim(-0.5, 2.5)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                zcmap = mpl.cm.ScalarMappable()
                zcmap.set_cmap(mpl.cm.spectral_r)
                norm = mpl.colors.Normalize(vmin = 0., vmax = 1.)
                zcmap.set_norm(norm)
                for i, msk in enumerate(mask):
                    j = i + 1
                    XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
                    if len(XM.compressed()) > 3 and len(YM.compressed()) > 3:
                        rs = runstats(XM.compressed(), YM.compressed(), **rs_kwargs)
                        c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
                        ax.plot(rs.xS, rs.yS, label = r'$%.2f < \mu \leq %.2f$' % (zticks[i], zticks[j]), c = c)
                ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
                xlabel = r'$x_Y$ [%]'
                ax.set_xlabel(xlabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (2, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(mu_GAL__rg, D_tau__rg, mask = mask__rg)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                ax.set_xlim(0, 1)
                ax.set_ylim(-0.5, 2.5)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$\mu$'
                ylabel = r'$\mathcal{D}_\tau$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (2, 2), colspan = 2, rowspan = 2)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                xm, ym, zm = C.ma_mask_xyz(z = H.x_Y__Trg[iT], y = D_tau__rg, x = mu_GAL__rg, mask = mask__rg)
                mask = []
                mask.append(np.less_equal(zm, 0.15))
                mask.append(np.bitwise_and(np.greater(zm, 0.15), np.less_equal(zm, 0.3)))
                mask.append(np.bitwise_and(np.greater(zm, 0.3), np.less_equal(zm, 0.45)))
                mask.append(np.bitwise_and(np.greater(zm, 0.45), np.less_equal(zm, 0.6)))
                zticks = [0, 0.15, 0.3, 0.45, 0.6]
                ax.set_xlim(0, 1)
                ax.set_ylim(-0.5, 2.5)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                zcmap = mpl.cm.ScalarMappable()
                zcmap.set_cmap(mpl.cm.spectral_r)
                norm = mpl.colors.Normalize(vmin = 0., vmax = 0.6)
                zcmap.set_norm(norm)
                for i, msk in enumerate(mask):
                    j = i + 1
                    XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
                    if len(XM.compressed()) > 3 and len(YM.compressed()) > 3:
                        rs = runstats(XM.compressed(), YM.compressed(), **rs_kwargs)
                        c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
                        ax.plot(rs.xS, rs.yS, label = r'$%.2f < x_Y \leq %.2f$' % (zticks[i], zticks[j]), c = c)
                ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
                xlabel = r'$\mu$'
                ax.set_xlabel(xlabel)
                ax.grid()
                
                ### End page
                f.subplots_adjust(hspace = 0.4, wspace = 0.4)
                pdf.savefig(f)
                plt.close(f)
        
                ##########################
                ######### PAGE 4 #########
                ##########################
                f = plt.figure()
                NRows = 4
                NCols = 4
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
                for ax in f.axes:
                    ax.set_axis_off()
                grid_shape = (NRows, NCols)
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(H.x_Y__Trg[iT], R_tau__rg, mask = mask__rg)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.set_xlim(0, 0.6)
                ax.set_ylim(0, 7)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$x_Y$'
                ylabel = r'$\mathcal{R}_\tau$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                xm, ym, zm = C.ma_mask_xyz(H.x_Y__Trg[iT], R_tau__rg, z = mu_GAL__rg, mask = mask__rg)
                mask = []
                mask.append(np.less_equal(zm, 0.25))
                mask.append(np.bitwise_and(np.greater(zm, 0.25), np.less_equal(zm, 0.5)))
                mask.append(np.bitwise_and(np.greater(zm, 0.5), np.less_equal(zm, 0.75)))
                mask.append(np.bitwise_and(np.greater(zm, 0.75), np.less_equal(zm, 1.)))
                zticks = [0, 0.25, 0.5, 0.75, 1.0]
                ax.set_xlim(0, 0.6)
                ax.set_ylim(0, 7)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                zcmap = mpl.cm.ScalarMappable()
                zcmap.set_cmap(mpl.cm.spectral_r)
                norm = mpl.colors.Normalize(vmin = 0., vmax = 1.)
                zcmap.set_norm(norm)
                for i, msk in enumerate(mask):
                    j = i + 1
                    XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
                    if len(XM.compressed()) > 3 and len(YM.compressed()) > 3:
                        rs = runstats(XM.compressed(), YM.compressed(), **rs_kwargs)
                        c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
                        ax.plot(rs.xS, rs.yS, label = r'$%.2f < \mu \leq %.2f$' % (zticks[i], zticks[j]), c = c)
                ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
                xlabel = r'$x_Y$'
                ax.set_xlabel(xlabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (2, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(mu_GAL__rg, R_tau__rg, mask = mask__rg)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 7)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$\mu$'
                ylabel = r'$\mathcal{R}_\tau$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (2, 2), colspan = 2, rowspan = 2)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                xm, ym, zm = C.ma_mask_xyz(z = H.x_Y__Trg[iT], y = R_tau__rg, x = mu_GAL__rg, mask = mask__rg)
                mask = []
                mask.append(np.less_equal(zm, 0.15))
                mask.append(np.bitwise_and(np.greater(zm, 0.15), np.less_equal(zm, 0.3)))
                mask.append(np.bitwise_and(np.greater(zm, 0.3), np.less_equal(zm, 0.45)))
                mask.append(np.bitwise_and(np.greater(zm, 0.45), np.less_equal(zm, 0.6)))
                zticks = [0, 0.15, 0.3, 0.45, 0.6]
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 7)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                zcmap = mpl.cm.ScalarMappable()
                zcmap.set_cmap(mpl.cm.spectral_r)
                norm = mpl.colors.Normalize(vmin = 0., vmax = 0.6)
                zcmap.set_norm(norm)
                for i, msk in enumerate(mask):
                    j = i + 1
                    XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
                    if len(XM.compressed()) > 3 and len(YM.compressed()) > 3:
                        rs = runstats(XM.compressed(), YM.compressed(), **rs_kwargs)
                        c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
                        ax.plot(rs.xS, rs.yS, label = r'$%.2f < x_Y \leq %.2f$' % (zticks[i], zticks[j]), c = c)
                ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
                xlabel = r'$\mu$'
                ax.set_xlabel(xlabel)
                ax.grid()
                
                ### End page
                f.subplots_adjust(hspace = 0.4, wspace = 0.4)
                pdf.savefig(f)
                plt.close(f)
        
                ##########################
                ######### PAGE 5 #########
                ##########################
                f = plt.figure()
                NRows = 2
                NCols = 4
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
                for ax in f.axes:
                    ax.set_axis_off()
                grid_shape = (NRows, NCols)
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(x = np.ma.log10(H.tau_V_neb__Trg[iT]), y = np.ma.log10(H.aSFRSD_Ha__Trg[iT]) + 6., mask = mask__rg)
                rs_kwargs = default_rs_kwargs.copy()
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                counts, xe, ye, im = ax.hist2d(xm.compressed(), ym.compressed(), bins = 30, cmap = mpl.cm.Blues)
                ax.set_xlim([-1.5, .5])
                ax.set_ylim([-3, 0])
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
                ax.plot(rs.xS, rs.yS, '--*')
                ax.legend(loc = 'upper left', fontsize = 8, frameon = False)
                xlabel = r'$\log\ \tau_V^{neb}$'
                ylabel = r'$\log\ \Sigma_{SFR}^{neb} (R)\ [M_\odot yr^{-1} kpc^{-2}]$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                f.colorbar(im)
                ax.grid()
            
                ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(x = np.ma.log10(tau_BC__rg), y = np.ma.log10(H.aSFRSD_Ha__Trg[iT]) + 6.0, mask = mask__rg)
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                counts, xe, ye, im = ax.hist2d(xm.compressed(), ym.compressed(), bins = 30, cmap = mpl.cm.Blues)
                ax.set_xlim([-1.5, .5])
                ax.set_ylim([-3, 0])
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
                ax.plot(rs.xS, rs.yS, '--*')
                ax.legend(loc = 'upper left', fontsize = 8, frameon = False)
                xlabel = r'$\log\ \tau_{BC}$' 
                ylabel = r'$\log\ \Sigma_{SFR}^{neb} (R)\ [M_\odot yr^{-1} kpc^{-2}]$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                f.colorbar(im)
                ax.grid()
        
                ### End page
                f.subplots_adjust(hspace = 0.2, wspace = 0.4, bottom = 0.2)
                pdf.savefig(f)
                plt.close(f)

        if args.integrated is True:
            with PdfPages('integrated_%s' % args.output) as pdf:
                ##########################
                ####### INTEGRATED #######
                ##########################
                ######### PAGE 1 #########
                ##########################
                NRows = 3
                NCols = 4
                f, axArr = plt.subplots(NRows, NCols)
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d (inside plot: mean, median, std, max, min)' % H.N_gals)
                for ax in f.axes:
                    ax.set_axis_off()
                props = [ 
                    H.integrated_tau_V__g,
                    H.integrated_tau_V_neb__g,
                    integrated_tau_ISM,
                    integrated_tau_BC,
                    integrated_rho,
                    mu_GAL,
                    integrated_tau_Y,
                    integrated_tau_O,
                    integrated_D_tau,
                    integrated_R_tau,
                    integrated_D_tau_xY,
                    integrated_R_tau_xY,
                ]
                props_label = [
                    r'$\tau_V^\star$',
                    r'$\tau_V^{neb}$',
                    r'$\tau_{ISM}$',
                    r'$\tau_{BC}$',
                    r'$\rho$',
                    r'$\mu$',
                    r'$\tau_Y$',
                    r'$\tau_O$',
                    r'$\mathcal{D}_\tau$',
                    r'$\mathcal{R}_\tau$',
                    r'$\mathcal{D}_\tau (x_Y)$',
                    r'$\mathcal{R}_\tau (x_Y)$',
                ]
                props_limits = [
                    [0, 2],
                    [0, 3.5],
                    [-0.5, 3],
                    [-1, 2],
                    [0, 10],
                    [0, 1],
                    [0, 10],
                    [0, 10],
                    [-1.5, 2.0],
                    [0, 10],
                    [-1.5, 2.0],
                    [0, 10],
                ]
                props_bins = [
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                ]
                for i, p in enumerate(props):
                    if isinstance(p, np.ma.core.MaskedArray):
                        aux = np.ma.masked_array(p, mask = (p.mask|mask_GAL__g))
                        x = aux.compressed()
                    else:
                        x = p[~mask_GAL__g]
                    ax = f.axes[i]
                    ax.set_axis_on()
                    ax.hist(x, bins = props_bins[i], range = props_limits[i])
                    ax.set_xlabel(props_label[i])
                    ax.set_xlim(props_limits[i])
                    txt = r'%.2f' % np.mean(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.96, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
                    txt = r'%.2f' % np.median(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.88, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
                    txt = r'%.2f' % np.std(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.80, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
                    txt = r'%.2f' % np.max(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.72, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
                    txt = r'%.2f' % np.min(x)
                    kw_text = dict(pos_x = 0.96, pos_y = 0.64, fs = 8, va = 'top', ha = 'right', c = 'k')
                    plot_text_ax(ax, txt, **kw_text)
        
                f.subplots_adjust(hspace = 0.4, wspace = 0.4)
                pdf.savefig(f)
                plt.close(f)
        
                default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                default_rs_kwargs = dict(smooth = True, sigma = 1.2, frac = 0.05)
                default_im_kwargs = dict(interpolation = 'nearest', origin = 'lower', aspect = 'auto', cmap = mpl.cm.spectral)
        
                ##########################
                ######### PAGE 2 #########
                ##########################
                f = plt.figure()
                NRows = 2
                NCols = 4
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
                for ax in f.axes:
                    ax.set_axis_off()
                grid_shape = (NRows, NCols)
        
                # tau_V vs tau_V_neb
                ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(H.integrated_tau_V__g, H.integrated_tau_V_neb__g, mask = mask_GAL__g)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.set_xlim([0, 2])
                ax.set_ylim([0, 3.5])
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.04, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.12, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xMedian, rs.yMedian, '-k')
                ax.set_xlabel(r'$\tau_V^{\star}$')
                ax.set_ylabel(r'$\tau_V^{neb}$')
                ax.grid()
        
                # tau_V_neb vs ln(EWHa/EWHb)
                ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
                y1 = np.ma.log(H.integrated_EW_Ha__g)
                y2 = np.ma.log(H.integrated_EW_Hb__g)
                y1m, y2m = C.ma_mask_xyz(y1, y2, mask = mask_GAL__g)
                ym = y1m - y2m
                xm, ym = C.ma_mask_xyz(H.integrated_tau_V_neb__g, ym, mask = mask_GAL__g)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.set_xlim([0, 3.5])
                ax.set_ylim([0, 2.8])
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.04, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.12, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$\tau_V^{neb}$'
                ylabel = r'$\ln\ W_{H\alpha} / W_{H\beta}$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ### End page
                f.subplots_adjust(hspace = 0.2, wspace = 0.4, bottom = 0.2)
                pdf.savefig(f)
                plt.close(f)
        
                ##########################
                ######### PAGE 3 #########
                ##########################
                f = plt.figure()
                NRows = 4
                NCols = 4
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
                for ax in f.axes:
                    ax.set_axis_off()
                grid_shape = (NRows, NCols)
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(H.integrated_x_Y__Tg[iT], integrated_D_tau, mask = mask_GAL__g)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.set_xlim(0, 0.4)
                ax.set_ylim(-0.5, 2.5)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$x_Y$ [%]'
                ylabel = r'$\mathcal{D}_\tau$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                xm, ym, zm = C.ma_mask_xyz(H.integrated_x_Y__Tg[iT], integrated_D_tau, z = mu_GAL, mask = mask_GAL__g)
                mask = []
                mask.append(np.less_equal(zm, 0.25))
                mask.append(np.bitwise_and(np.greater(zm, 0.25), np.less_equal(zm, 0.5)))
                mask.append(np.bitwise_and(np.greater(zm, 0.5), np.less_equal(zm, 0.75)))
                mask.append(np.bitwise_and(np.greater(zm, 0.75), np.less_equal(zm, 1.)))
                #mask.append(np.bitwise_and(np.greater(zm, 0.8), np.less_equal(zm, 1.)))
                ax.set_xlim(0, 0.4)
                ax.set_ylim(-0.5, 2.5)
                zticks = [0, 0.25, 0.5, 0.75, 1.0]
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                zcmap = mpl.cm.ScalarMappable()
                zcmap.set_cmap(mpl.cm.spectral_r)
                norm = mpl.colors.Normalize(vmin = 0., vmax = 1.)
                zcmap.set_norm(norm)
                for i, msk in enumerate(mask):
                    j = i + 1
                    XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
                    if len(XM.compressed()) > 3 and len(YM.compressed()) > 3:
                        rs = runstats(XM.compressed(), YM.compressed(), **rs_kwargs)
                        c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
                        ax.plot(rs.xS, rs.yS, label = r'$%.2f < \mu \leq %.2f$' % (zticks[i], zticks[j]), c = c)
                ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
                xlabel = r'$x_Y$ [%]'
                ax.set_xlabel(xlabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (2, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(mu_GAL, integrated_D_tau, mask = mask_GAL__g)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                ax.set_xlim(0, 1)
                ax.set_ylim(-0.5, 2.5)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$\mu$'
                ylabel = r'$\mathcal{D}_\tau$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (2, 2), colspan = 2, rowspan = 2)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                xm, ym, zm = C.ma_mask_xyz(z = H.integrated_x_Y__Tg[iT], y = integrated_D_tau, x = mu_GAL, mask = mask_GAL__g)
                mask = []
                mask.append(np.less_equal(zm, 0.15))
                mask.append(np.bitwise_and(np.greater(zm, 0.15), np.less_equal(zm, 0.3)))
                mask.append(np.bitwise_and(np.greater(zm, 0.3), np.less_equal(zm, 0.45)))
                mask.append(np.bitwise_and(np.greater(zm, 0.45), np.less_equal(zm, 0.6)))
                ax.set_xlim(0, 1)
                ax.set_ylim(-0.5, 2.5)
                zticks = [0, 0.15, 0.3, 0.45, 0.6]
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                zcmap = mpl.cm.ScalarMappable()
                zcmap.set_cmap(mpl.cm.spectral_r)
                norm = mpl.colors.Normalize(vmin = 0., vmax = 0.6)
                zcmap.set_norm(norm)
                for i, msk in enumerate(mask):
                    j = i + 1
                    XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
                    if len(XM.compressed()) > 3 and len(YM.compressed()) > 3:
                        rs = runstats(XM.compressed(), YM.compressed(), **rs_kwargs)
                        c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
                        ax.plot(rs.xS, rs.yS, label = r'$%.2f < x_Y \leq %.2f$' % (zticks[i], zticks[j]), c = c)
                ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
                xlabel = r'$\mu$'
                ax.set_xlabel(xlabel)
                ax.grid()
                
                ### End page
                f.subplots_adjust(hspace = 0.4, wspace = 0.4)
                pdf.savefig(f)
                plt.close(f)
        
                ##########################
                ######### PAGE 4 #########
                ##########################
                f = plt.figure()
                NRows = 4
                NCols = 4
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
                for ax in f.axes:
                    ax.set_axis_off()
                grid_shape = (NRows, NCols)
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(H.integrated_x_Y__Tg[iT], integrated_R_tau, mask = mask_GAL__g)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.set_xlim(0, 0.4)
                ax.set_ylim(0, 7)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$x_Y$'
                ylabel = r'$\mathcal{R}_\tau$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                xm, ym, zm = C.ma_mask_xyz(H.integrated_x_Y__Tg[iT], integrated_R_tau, z = mu_GAL, mask = mask_GAL__g)
                mask = []
                mask.append(np.less_equal(zm, 0.25))
                mask.append(np.bitwise_and(np.greater(zm, 0.25), np.less_equal(zm, 0.5)))
                mask.append(np.bitwise_and(np.greater(zm, 0.5), np.less_equal(zm, 0.75)))
                mask.append(np.bitwise_and(np.greater(zm, 0.75), np.less_equal(zm, 1.)))
                zticks = [0, 0.25, 0.5, 0.75, 1.0]
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                ax.set_xlim(0, 0.4)
                ax.set_ylim(0, 7)
                zcmap = mpl.cm.ScalarMappable()
                zcmap.set_cmap(mpl.cm.spectral_r)
                norm = mpl.colors.Normalize(vmin = 0., vmax = 1.)
                zcmap.set_norm(norm)
                for i, msk in enumerate(mask):
                    j = i + 1
                    XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
                    if len(XM.compressed()) > 3 and len(YM.compressed()) > 3:
                        rs = runstats(XM.compressed(), YM.compressed(), **rs_kwargs)
                        c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
                        ax.plot(rs.xS, rs.yS, label = r'$%.2f < \mu \leq %.2f$' % (zticks[i], zticks[j]), c = c)
                ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
                xlabel = r'$x_Y$'
                ax.set_xlabel(xlabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (2, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(mu_GAL, integrated_R_tau, mask = mask_GAL__g)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 7)
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.96, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.88, fs = 8, c = 'r', rms = True, label = 'OLS', text = True))
                Nprc = len(rs.xPrc)
                for i in range(Nprc):
                    ax.plot(rs.xPrc[i], rs.yPrc[i], '--k')
                ax.plot(rs.xS, rs.yS, '-k')
                xlabel = r'$\mu$'
                ylabel = r'$\mathcal{R}_\tau$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.grid()
        
                ax = plt.subplot2grid(grid_shape, loc = (2, 2), colspan = 2, rowspan = 2)
                rs_kwargs = default_rs_kwargs.copy()
                sc_kwargs = default_sc_kwargs.copy()
                sc_kwargs['s'] = 10
                sc_kwargs['alpha'] = 0.8
                xm, ym, zm = C.ma_mask_xyz(z = H.integrated_x_Y__Tg[iT], y = integrated_R_tau, x = mu_GAL, mask = mask_GAL__g)
                mask = []
                mask.append(np.less_equal(zm, 0.15))
                mask.append(np.bitwise_and(np.greater(zm, 0.15), np.less_equal(zm, 0.3)))
                mask.append(np.bitwise_and(np.greater(zm, 0.3), np.less_equal(zm, 0.45)))
                mask.append(np.bitwise_and(np.greater(zm, 0.45), np.less_equal(zm, 0.6)))
                zticks = [0, 0.15, 0.3, 0.45, 0.6]
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 7)
                ax.scatter(xm, ym, c = '0.5', **sc_kwargs)
                zcmap = mpl.cm.ScalarMappable()
                zcmap.set_cmap(mpl.cm.spectral_r)
                norm = mpl.colors.Normalize(vmin = 0., vmax = 0.6)
                zcmap.set_norm(norm)
                for i, msk in enumerate(mask):
                    j = i + 1
                    XM, YM = C.ma_mask_xyz(xm, ym, mask = ~msk)
                    if len(XM.compressed()) > 3 and len(YM.compressed()) > 3:
                        rs = runstats(XM.compressed(), YM.compressed(), **rs_kwargs)
                        c = zcmap.to_rgba(0.5 * (zticks[i] + zticks[j]))
                        ax.plot(rs.xS, rs.yS, label = r'$%.2f < x_Y \leq %.2f$' % (zticks[i], zticks[j]), c = c)
                ax.legend(loc = 'upper right', fontsize = 8, frameon = False)
                xlabel = r'$\mu$'
                ax.set_xlabel(xlabel)
                ax.grid()
                
                ### End page
                f.subplots_adjust(hspace = 0.4, wspace = 0.4)
                pdf.savefig(f)
                plt.close(f)
                
                ##########################
                ######### PAGE 5 #########
                ##########################
                f = plt.figure()
                NRows = 2
                NCols = 4
                page_size_inches = (NCols * 3, NRows * 2)
                f.set_size_inches(page_size_inches)
                f.suptitle(r'$N_{gals}$: %d  $N_{zones}$:  %d' % (H.N_gals, H.N_zones__g.sum()))
                for ax in f.axes:
                    ax.set_axis_off()
                grid_shape = (NRows, NCols)
        
                ax = plt.subplot2grid(grid_shape, loc = (0, 0), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(x = np.ma.log10(H.integrated_tau_V_neb__g), y = np.ma.log10(H.integrated_SFRSD_Ha__g) + 6., mask = mask_GAL__g)
                rs_kwargs = default_rs_kwargs.copy()
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                counts, xe, ye, im = ax.hist2d(xm.compressed(), ym.compressed(), bins = 10, cmap = mpl.cm.Blues)
                ax.set_xlim([-1.5, .5])
                ax.set_ylim([-3, 0])
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
                ax.plot(rs.xS, rs.yS, '--*')
                ax.legend(loc = 'upper left', fontsize = 8, frameon = False)
                xlabel = r'$\tau_V^{neb}$'
                ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                f.colorbar(im)
                ax.grid()
            
                ax = plt.subplot2grid(grid_shape, loc = (0, 2), colspan = 2, rowspan = 2)
                xm, ym = C.ma_mask_xyz(x = np.ma.log10(integrated_tau_BC), y = np.ma.log10(H.integrated_SFRSD_Ha__g) + 6., mask = mask_GAL__g)
                rs = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
                counts, xe, ye, im = ax.hist2d(xm.compressed(), ym.compressed(), bins = 10, cmap = mpl.cm.Blues)
                ax.set_xlim([-1.5, .5])
                ax.set_ylim([-3, 0])
                plotOLSbisectorAxis(ax, rs.xS, rs.yS, **dict(pos_x = 0.96, pos_y = 0.08, fs = 8, c = 'r', rms = True, label = 'OLS(tend)', text = True))
                plotOLSbisectorAxis(ax, rs.x, rs.y, **dict(pos_x = 0.96, pos_y = 0.16, fs = 8, c = 'k', rms = True, label = 'OLS', text = True))
                ax.plot(rs.xS, rs.yS, '--*')
                ax.legend(loc = 'upper left', fontsize = 8, frameon = False)
                xlabel = r'$\tau_{BC}$' 
                ylabel = r'$\log\ \Sigma_{SFR}^{neb} (R)\ [M_\odot yr^{-1} kpc^{-2}]$'
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                f.colorbar(im)
                ax.grid()
        
                ### End page
                f.subplots_adjust(hspace = 0.2, wspace = 0.4, bottom = 0.2)
                pdf.savefig(f)
                plt.close(f)

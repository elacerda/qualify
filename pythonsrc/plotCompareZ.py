#!/opt/local/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import sys
import numpy as np
import CALIFAUtils as C
import matplotlib as mpl
from matplotlib.cm import cmap_d
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
from CALIFAUtils.plots import plotOLSbisectorAxis, plot_text_ax, next_row_col
from scipy.stats import spearmanr

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
    
    default_rs_kwargs = dict(smooth = True, sigma = 1.2, frac = 0.07)
    default_sc_kwargs = dict(marker = 'o', edgecolor = 'none')
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
    f = plt.figure()
    NCols = 3
    NRows = 2
    page_size_inches = [15, 8]
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
    
    y = H.integrated_logO3N2_M13__g - 8.69
    #yran = [0, 3]
    row, col = 0, 0
    for iU, tZ in enumerate(H.tZ__U):
        ax = plt.subplot2grid(grid_shape, loc = (row, col))
        #xran = [0, 1.5]
        xm, ym = C.ma_mask_xyz(H.alogZ_mass_GAL__Ug[iU], y = y, mask = mask_GAL__g)
        ax.scatter(xm, ym, c = '0.8', **default_sc_kwargs)
        #density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range = [xran, yran], colors = [ 'b', 'y', 'r' ])
        kw_text = dict(pos_x = 0.01, pos_y = 0.99, fs = 8, va = 'top', ha = 'left', c = 'k')
        plot_text_ax(ax, '%s' % xm.count(), **kw_text)
        rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
        ax.plot(rs.xS, rs.yS, 'k-')
        ax.set_title('%.2f Ganos' % (tZ / 1e9))
        ax.set_xlim(-2., 0.4)
        ax.set_ylim(-0.7, 0.4)
        ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--')
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        if col == 0:
            ax.set_ylabel(r'$\log\ \left(\frac{(O/H)}{(O/H)_\odot}\right)$')
        if row == 1:
            ax.set_xlabel(r'$\langle \log\ Z_\star \rangle_M$ [$Z_\odot$]')
        row, col = next_row_col(row, col, NRows, NCols)
    f.subplots_adjust(bottom = 0.15, top = 0.95, hspace = 0.3, wspace = 0.3, right = 0.95, left = 0.1)
    #pdf.savefig(f)
    fname = 'CompareZ_%s%s' % (basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix)
    f.savefig(fname)
    plt.close(f)

    ######################
    ######################
    ######################
    f = plt.figure()

    NCols = 3
    NRows = 2
    
    lw = 2
    bins = 30
    cmap_radii = 'jet_r'
    
    page_size_inches = [15, 8]
    f.set_size_inches(page_size_inches)
    f.set_dpi(100)
    grid_shape = (NRows, NCols)
    # zones
    ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    x = np.ma.log10(H.McorSD__Tg[iT])
    xm, ym, zm = C.ma_mask_xyz(y = H.alogZ_mass__Ug[-1], x = x, z = H.zone_dist_HLR__g, mask = mask__g)    
    #density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range = [xran, yran], colors = [ 'b', 'y', 'r' ])
    #counts, xedges, yedges, im = ax.hist2d(xm.compressed(), ym.compressed(),  cmin = 0.001 * xm.count(), bins = bins, cmap = 'Blues')
    sc = ax.scatter(xm.compressed(), ym.compressed(), c = zm.compressed(), s = 2, vmin = minR, vmax = H.RbinFin, cmap = cmap_radii, alpha = 0.4, **default_sc_kwargs)
    cax = plt.axes([0.92, 0.1, 0.025, 0.8])
    cb = f.colorbar(sc, cax = cax, ticks = ticks_r)
    cb.ax.set_yticklabels(ticks_r)
    #cb.ax.yaxis.set_major_locator(MaxNLocator(4))
    cmap = plt.get_cmap('viridis')
    for iU, tZ in enumerate(H.tZ__U):
        c = cmap(float(iU) / H.N_U)
        xm, ym = C.ma_mask_xyz(x = x, y = H.alogZ_mass__Ug[iU], mask = mask__g)
        rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
        rs.OLS_bisector()
        #print tZ/1e9, rs.OLS_median_intercept, rs.OLS_median_slope, rs.OLS_intercept, rs.OLS_slope
        ax.plot(rs.xS, rs.yS, c = c, lw = lw)
        kw_text = dict(pos_x = 0.01, pos_y = 0.99 - (iU * 0.08), fs = 15, va = 'top', ha = 'left', c = c)
        plot_text_ax(ax, '%.2f Ga' % (tZ / 1e9), **kw_text)
        row, col = next_row_col(row, col, NRows, NCols)
    xm, ym = C.ma_mask_xyz(y = H.logO3N2_M13__g - 8.69, x = x, mask = mask__g)    
    rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
    kw_text = dict(transform = False, pos_x = rs.xS[0] - 0.02, pos_y = rs.yS[0], fs = 12, va = 'center', ha = 'right', c = 'k')
    plot_text_ax(ax, r'$Z_{neb}$' % (tZ / 1e9), **kw_text)
    ax.plot(rs.xS, rs.yS, 'k--')
    ax.set_ylabel(r'metallicity [$Z_\odot$]')
    ax.set_xlabel(r'$\log\ \mu_\star$ [$M_\odot pc^{-2}$]')
    ax.set_xlim(0.,4.5)
    ax.set_ylim(-2,0.4)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.set_title('zonas')

    # radii
    ax = plt.subplot2grid(grid_shape, loc = (0, 1))
    x = np.ma.log10(H.McorSD__Trg[iT])
    xm, ym = C.ma_mask_xyz(y = H.alogZ_mass__Urg[-1], x = x, mask = mask__rg)    
    #ax.hist2d(xm.compressed(), ym.compressed(), bins = bins, cmap = 'Blues')
    ax.scatter(xm, ym, c = H.Rtoplot(), cmap = cmap_radii, vmin = minR, alpha = 0.4, s = 4, **default_sc_kwargs)
    #cb = plt.colorbar(sc)
    #cb.ax.yaxis.set_label('R [HLR]')
    cmap = plt.get_cmap('viridis')
    kw_text = dict(pos_x = 0.01, pos_y = 0.99, fs = 12, va = 'top', ha = 'left', c = 'k')
    for iU, tZ in enumerate(H.tZ__U):
        c = cmap(float(iU) / H.N_U)
        xm, ym = C.ma_mask_xyz(x = x, y = H.alogZ_mass__Urg[iU], mask = mask__rg)
        rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
        rs.OLS_bisector()
        #print tZ/1e9, rs.OLS_median_intercept, rs.OLS_median_slope, rs.OLS_intercept, rs.OLS_slope
        ax.plot(rs.xS, rs.yS, c = c, lw = lw)
        row, col = next_row_col(row, col, NRows, NCols)
    xm, ym = C.ma_mask_xyz(y = H.logO3N2_M13__Trg[iT] - 8.69, x = x, mask = mask__rg)    
    rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
    kw_text = dict(transform = False, pos_x = rs.xS[0] - 0.02, pos_y = rs.yS[0], fs = 12, va = 'center', ha = 'right', c = 'k')
    plot_text_ax(ax, r'$Z_{neb}$' % (tZ / 1e9), **kw_text)
    ax.plot(rs.xS, rs.yS, 'k--')
    ax.set_xlabel(r'$\langle\log\ \mu_\star\rangle_R$ [$M_\odot pc^{-2}$]')
    ax.set_xlim(0., 4.5)
    ax.set_ylim(-2,0.4)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.set_title('bins radiais')

    # integrate
    ax = plt.subplot2grid(grid_shape, loc = (0, 2))
    x = np.ma.log10(H.McorSD_GAL__g)
    xm, ym = C.ma_mask_xyz(x = x, y = H.alogZ_mass_GAL__Ug[-1], mask = mask_GAL__g)
    ax.scatter(xm, ym, c = '0.8', **default_sc_kwargs)    
    cmap = plt.get_cmap('viridis')
    for iU, tZ in enumerate(H.tZ__U):
        c = cmap(float(iU) / H.N_U)
        xm, ym = C.ma_mask_xyz(x = x, y = H.alogZ_mass_GAL__Ug[iU], mask = mask_GAL__g)    
        rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
        rs.OLS_bisector()
        #print tZ/1e9, rs.OLS_median_intercept, rs.OLS_median_slope, rs.OLS_intercept, rs.OLS_slope
        ax.plot(rs.xS, rs.yS, c = c, lw = lw)
        row, col = next_row_col(row, col, NRows, NCols)
    xm, ym = C.ma_mask_xyz(y = H.integrated_logO3N2_M13__g - 8.69, x = x, mask = mask_GAL__g)
    rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
    kw_text = dict(transform = False, pos_x = rs.xS[0] - 0.02, pos_y = rs.yS[0], fs = 12, va = 'center', ha = 'right', c = 'k')
    plot_text_ax(ax, '$Z_{neb}$' % (tZ / 1e9), **kw_text)
    ax.plot(rs.xS, rs.yS, 'k--')
    ax.set_xlabel(r'$\langle\log\ \mu_\star\rangle_{GAL}$ [$M_\odot$]')
    ax.set_xlim(0, 4.5)
    ax.set_ylim(-2,0.4)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.set_title('GAL')
    print H.RbinFin
    # zones
    ax = plt.subplot2grid(grid_shape, loc = (1, 0))
    y = H.logO3N2_M13__g - 8.69
    xm, ym = C.ma_mask_xyz(x = H.alogZ_mass__Ug[-1], y = y, mask = mask__g)
    #density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range = [xran, yran], colors = [ 'b', 'y', 'r' ])
    #ax.hist2d(xm.compressed(), ym.compressed(), bins = bins, cmap = 'Blues')
    ax.scatter(xm, ym, c = H.zone_dist_HLR__g, s = 2, vmin = minR, vmax = H.RbinFin, cmap = cmap_radii, alpha = 0.4, **default_sc_kwargs)
    #cb = plt.colorbar(sc)
    #cb.ax.yaxis.set_major_locator(MaxNLocator(4))
    cmap = plt.get_cmap('viridis')
    for iU, tZ in enumerate(H.tZ__U):
        c = cmap(float(iU) / H.N_U)
        xm, ym = C.ma_mask_xyz(x = H.alogZ_mass__Ug[iU], y = y, mask = mask__g)
        rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
        rs.OLS_bisector()
        print tZ/1e9, rs.OLS_median_intercept, rs.OLS_median_slope, rs.OLS_intercept, rs.OLS_slope
        ax.plot(rs.xS, rs.yS, c = c, lw = lw)
        kw_text = dict(pos_x = 0.01 + (iU * 0.17), pos_y = 0.01, fs = 15, va = 'bottom', ha = 'left', c = c)
        plot_text_ax(ax, '%.3f' % spearmanr(xm.compressed(), ym.compressed())[0], **kw_text)
        row, col = next_row_col(row, col, NRows, NCols)
    ax.set_xlabel(r'$\langle \log\ Z_\star \rangle_M$ [$Z_\odot$]')
    ax.set_ylabel(r'$\log\ \left(\frac{(O/H)}{(O/H)_\odot}\right)$')
    ax.set_xlim(-2, 0.4)
    ax.set_ylim(-0.7, 0)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    #ax.set_title('zonas')

    # radii
    ax = plt.subplot2grid(grid_shape, loc = (1, 1))
    y = H.logO3N2_M13__Trg[iT] - 8.69    
    xm, ym = C.ma_mask_xyz(x = H.alogZ_mass__Urg[-1], y = y, mask = mask__rg)    
    #ax.hist2d(xm.compressed(), ym.compressed(), bins = bins, cmap = 'Blues')
    ax.scatter(xm, ym, c = H.Rtoplot(), cmap = cmap_radii, vmin = minR, alpha = 0.4, s = 4, **default_sc_kwargs)
    #cb = plt.colorbar(sc)
    #cb.ax.yaxis.set_major_locator(MaxNLocator(4))
    cmap = plt.get_cmap('viridis')
    kw_text = dict(pos_x = 0.01, pos_y = 0.99, fs = 12, va = 'top', ha = 'left', c = 'k')
    for iU, tZ in enumerate(H.tZ__U):
        c = cmap(float(iU) / H.N_U)
        xm, ym = C.ma_mask_xyz(x = H.alogZ_mass__Urg[iU], y = y, mask = mask__rg)
        rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
        rs.OLS_bisector()
        print tZ/1e9, rs.OLS_median_intercept, rs.OLS_median_slope, rs.OLS_intercept, rs.OLS_slope
        ax.plot(rs.xS, rs.yS, c = c, lw = lw)
        kw_text = dict(pos_x = 0.01 + (iU * 0.17), pos_y = 0.01, fs = 15, va = 'bottom', ha = 'left', c = c)
        plot_text_ax(ax, '%.3f' % spearmanr(xm.compressed(), ym.compressed())[0], **kw_text)
        row, col = next_row_col(row, col, NRows, NCols)
    ax.set_xlabel(r'$\langle \log\ Z_\star \rangle_M(R)$ [$Z_\odot$]')
    ax.set_xlim(-2, 0.4)
    ax.set_ylim(-0.7, 0)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    #ax.set_title('bins radiais')

    # integrate
    ax = plt.subplot2grid(grid_shape, loc = (1, 2))
    y = H.integrated_logO3N2_M13__g - 8.69
    xm, ym = C.ma_mask_xyz(x = H.alogZ_mass_GAL__Ug[-1], y = y, mask = mask_GAL__g)
    ax.scatter(xm, ym, c = '0.8', **default_sc_kwargs)    
    cmap = plt.get_cmap('viridis')
    for iU, tZ in enumerate(H.tZ__U):
        c = cmap(float(iU) / H.N_U)
        xm, ym = C.ma_mask_xyz(x = H.alogZ_mass_GAL__Ug[iU], y = y, mask = mask_GAL__g)    
        rs = C.runstats(xm.compressed(), ym.compressed(), debug = args.debug, **default_rs_kwargs)
        rs.OLS_bisector()
        print tZ/1e9, rs.OLS_median_intercept, rs.OLS_median_slope, rs.OLS_intercept, rs.OLS_slope
        ax.plot(rs.xS, rs.yS, c = c, lw = lw)
        kw_text = dict(pos_x = 0.01 + (iU * 0.17), pos_y = 0.01, fs = 15, va = 'bottom', ha = 'left', c = c)
        plot_text_ax(ax, '%.3f' % spearmanr(xm.compressed(), ym.compressed())[0], **kw_text)
        row, col = next_row_col(row, col, NRows, NCols)
    ax.set_xlabel(r'$\langle \log\ Z_\star \rangle_M^{GAL}$ [$Z_\odot$]')
    ax.set_xlim(-2, 0.4)
    ax.set_ylim(-0.7, 0)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    #ax.set_title('GAL')
    
    f.subplots_adjust(bottom = 0.15, top = 0.9, hspace = 0.4, wspace = 0.25, right = 0.9, left = 0.1)
    fname = 'stellarmuZR_%s%s' % (basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix)
    f.savefig(fname)
    plt.close(f)
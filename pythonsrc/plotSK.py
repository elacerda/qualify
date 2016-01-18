#!/opt/local/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import sys
import numpy as np
import CALIFAUtils as C
import matplotlib as mpl
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['axes.titlesize'] = 10
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
from os.path import basename
import argparse as ap
from CALIFAUtils.lines import Lines
from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib.backends.backend_pdf import PdfPages
from CALIFAUtils.plots import plotOLSbisectorAxis, plot_text_ax

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'hdf5' : None,
        'output' : '',
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
    parser.add_argument('--slice_gals', '-S',
                        metavar = 'FILE',
                        type = str,
                        default = default['slice_gals'])
    parser.add_argument('--maskradius', '-R',
                        help = 'initial RDisc value in HLR',
                        metavar = 'NUM',
                        type = float,
                        default = default['maskradius'])
    parser.add_argument('--output', '-o',
                        help = 'Name of the output PDF file.',
                        metavar = 'FILENAME',
                        type = str,
                        default = default['output'])
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
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        maxR = H.Rbin__r[-1]
        maxR = 3
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

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # mask__rg = np.bitwise_or(np.ma.log10(H.aSFRSD__Trg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Trg[iT]).mask)
    # mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.aSFRSD_Ha__rg * 1e6).mask)
    # mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.tau_V_neb__rg).mask)
    # mask__rg = np.bitwise_or(mask__rg, H.O_O3N2_M13__rg.mask)
    # #mask__rg = np.bitwise_or(mask__rg, np.less(H.EW_Ha__rg, 3.))
    # mask__rg = np.bitwise_or(mask__rg, np.less(H.reply_arr_by_radius(H.ba_GAL__g), ba_max))
    # mask__rg = np.bitwise_or(mask__rg, ~maskRadiusOk__rg)
    # mask__rg = np.bitwise_or(mask__rg, ~gals_slice__rg)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))

    default_kwargs_rs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 50)
    default_sc_kwargs = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = '')
    default_kwargs_zbins = dict(
        debug = True,
        zlabel = r'R (HLR)',
        zmask = None,
        #xlim = [-1.5, 0.5],
        #ylim = [-3.5, 1],
        zlim = zlim,
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        kwargs_rs = default_kwargs_rs,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.9),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.45, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 6, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        write_N = True,
    )
    kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS')
    kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 6, rms = True, text = True, kwargs_plot = kwargs_ols_plot)
            
    with PdfPages('SK_%.2fMyr%s%s' % ((tSF/1e6), basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix)) as pdf:
        ############## BPT ##############
        ############## BPT ##############
        ############## BPT ##############
        
        NRows = 2
        NCols = 2
        f = plt.figure()
        page_size_inches = np.asarray(A4Size_inches)
        #page_size_inches[1] /= 2.
        f.set_size_inches(page_size_inches.tolist())
        grid_shape = (NRows, NCols)

        if minR is None:
            suptitle = r'NGals:%d  Nzones:%d  NRbins:%d  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (~mask__g).sum(), (~mask__rg).sum(), (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        else:
            suptitle = r'NGals:%d  Nzones:%d  NRbins:%d  R > %.1fHLR  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (~mask__g).sum(), (~mask__rg).sum(), minR, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    
        f.suptitle(suptitle, fontsize = 9)
        
        EW_Ha__g = H.EW_Ha__g
        EW_Ha__rg = H.EW_Ha__Trg[iT]
        Ha__g = H.F_int_Ha__g
        Ha__rg = H.F_int_Ha__Trg[iT]
        Hb__g = H.F_int_Hb__g
        Hb__rg = H.F_int_Hb__Trg[iT]
        O3__g = H.F_int_O3__g
        O3__rg = H.F_int_O3__Trg[iT]
        N2__g = H.F_int_N2__g
        N2__rg = H.F_int_N2__Trg[iT]
        
        logWHa__g = np.ma.masked_array(np.ma.log10(EW_Ha__g), mask = mask__g)
        logWHa__rg = np.ma.masked_array(np.ma.log10(EW_Ha__rg), mask = mask__rg)
        logO3Hb__g = np.ma.masked_array(np.ma.log10(O3__g) - np.ma.log10(Hb__g), mask = mask__g)
        logO3Hb__rg = np.ma.masked_array(np.ma.log10(O3__rg) - np.ma.log10(Hb__rg), mask = mask__rg)
        logN2Ha__g = np.ma.masked_array(np.ma.log10(N2__g) - np.ma.log10(Ha__g), mask = mask__g)
        logN2Ha__rg = np.ma.masked_array(np.ma.log10(N2__rg) - np.ma.log10(Ha__rg), mask = mask__rg)

        l = Lines()
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        for line in l.linesbpt:
            ax.plot(l.x[line], l.y[line], label = line)
        ax.text(-0.7, -0.5, 'SF')
        ax.text(-0.3, 0.6, 'Seyfert')
        ax.text(0.15, -0.5, 'LINER')
        C.plot_zbins(
            debug = True,
            cb = False,
            write_N = True,
            zmask = None,
            f = f,
            ax = ax,
            x = logN2Ha__g, 
            y = logO3Hb__g, 
            z = H.zone_dist_HLR__g, 
            zlim = zlim,
            xlim = (-1.5, 0.5),
            ylim = (-1.5, 1.0),
            x_major_locator = 1.,
            x_minor_locator = 0.2,
            y_major_locator = 1.,
            y_minor_locator = 0.2,
            kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.45, label = ''),
            xlabel = r'$\log [NII]/H\alpha$',
            ylabel = r'$\log [OIII]/H\beta$',
            zlabel = r'R (HLR)',
        )

        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        ax.set_axis_on()
        for line in l.linesbpt:
            ax.plot(l.x[line], l.y[line], label = line)
        ax.text(-0.7, -0.5, 'SF')
        ax.text(-0.3, 0.6, 'Seyfert')
        ax.text(0.15, -0.5, 'LINER')
        C.plot_zbins(
            debug = True,
            cb = True,
            write_N = True,
            zmask = None,
            f = f,
            ax = ax,
            x = logN2Ha__rg, 
            y = logO3Hb__rg, 
            z = H.Rtoplot(), 
            zlim = zlim,
            xlim = (-1.5, 0.5),
            ylim = (-1.5, 1.0),
            x_major_locator = 1.,
            x_minor_locator = 0.2,
            y_major_locator = 1.,
            y_minor_locator = 0.2,
            kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.45, label = ''),
            xlabel = r'$\log [NII]/H\alpha(R)$',
            ylabel = r'$\log [OIII]/H\beta(R)$',
            zlabel = r'R (HLR)',
        )

        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        ax.plot((-0.4, -0.4), (np.log10(3), 3), 'k-')
        ax.plot((-0.4, 0.5), np.ma.log10([6, 6]), 'k-')
        ax.axhline(y = np.log10(3), c = 'k')
        p = [np.log10(0.5/5.0), np.log10(0.5)]
        xini = (np.log10(3.) - p[1]) / p[0]
        ax.plot((xini, 0.), np.polyval(p, [xini, 0.]), 'k:')
        ax.plot((0, 0.5), np.log10([0.5, 0.5]), 'k:')
        ax.text(-1.4, 0.75, 'SF')
        ax.text(0.17, 0.9, 'sAGN')
        ax.text(0.15, 0.55, 'wAGN')
        ax.text(0.3, 0.0, 'RG')
        ax.text(-0.8, 0, 'PG')        
        C.plot_zbins(
            debug = True,
            cb = False,
            write_N = True,
            zmask = None,
            f = f,
            ax = ax,
            x = logN2Ha__g, 
            y = logWHa__g, 
            z = H.zone_dist_HLR__g, 
            xlim = (-1.5, 0.5),
            ylim = (-0.5, 3),
            zlim = zlim,
            x_major_locator = 1.,
            x_minor_locator = 0.2,
            kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.45, label = ''),
            xlabel = r'$\log [NII]/H\alpha$',
            ylabel = r'$\log WH\alpha$',
            zlabel = r'R (HLR)',
        )

        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        ax.set_axis_on()
        ax.plot((-0.4, -0.4), (np.log10(3.), 3.), 'k-')
        ax.plot((-0.4, 0.5),  np.ma.log10([6., 6.]), 'k-')
        ax.axhline(np.log10(3), c = 'k')
        p = [np.log10(0.5/5.0), np.log10(0.5)]
        xini = (np.log10(3.) - p[1]) / p[0]
        ax.plot((xini, 0.), np.polyval(p, [xini, 0.]), 'k:')
        ax.plot((0, 0.5), np.log10([0.5, 0.5]), 'k:')
        ax.text(-1.4, 0.75, 'SF')
        ax.text(0.13, 0.9, 'sAGN')
        ax.text(0.11, 0.55, 'wAGN')
        ax.text(0.25, 0.0, 'RG')
        ax.text(-0.8, 0, 'PG')        
        kw = C.plot_zbins(
            debug = True,
            cb = True,
            return_kwargs = True, 
            write_N = True,
            zmask = None,
            f = f,
            ax = ax,
            x = logN2Ha__rg, 
            y = logWHa__rg, 
            z = H.Rtoplot(), 
            xlim = (-1.5, 0.5),
            ylim = (-0.5, 3),
            zlim = zlim,
            x_major_locator = 1.,
            x_minor_locator = 0.2,
            kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.45, label = ''),
            xlabel = r'$\log [NII]/H\alpha(R)$',
            ylabel = r'$\log WH\alpha(R)$',
            zlabel = r'R (HLR)',
        )

        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)

        ############## Syn ##############
        ############## Syn ##############
        ############## Syn ##############

        NRows = 3
        NCols = 4
        f = plt.figure()
        kwargs_zbins = default_kwargs_zbins.copy()
        kwargs_zbins.update(f = f)
        #f, axArr = plt.subplots(NRows, NCols)
        #page_size_inches = (NCols * 3, NRows * 1.5)
        page_size_inches = A4Size_inches
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)

        #f, axArr = plt.subplots(1, 3)
        #f.set_dpi(100)
        #f.set_size_inches(15, 8)
        f.suptitle(suptitle, fontsize = 9)
         
        tau_V_norm__rg = H.tau_V__Trg[iT] / H.tau_V_oneHLR__Tg[iT]
        aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / H.aSFRSD_oneHLR__Tg[iT]

        x1, y1 = C.ma_mask_xyz(np.ma.log10(H.tau_V__Tg[iT]), np.ma.log10(H.SFRSD__Tg[iT] * 1e6), mask = mask__g)
        x2, y2 = C.ma_mask_xyz(np.ma.log10(H.tau_V__Trg[iT]), np.ma.log10(H.aSFRSD__Trg[iT] * 1e6), mask = mask__rg) 
        x3, y3 = C.ma_mask_xyz(np.ma.log10(tau_V_norm__rg), np.ma.log10(aSFRSD_norm__rg), mask = mask__rg) 
        x4, y4 = C.ma_mask_xyz(np.ma.log10(H.tau_V_oneHLR__Tg[iT]), np.ma.log10(H.aSFRSD_oneHLR__Tg[iT] * 1e6), mask = mask_GAL__g)
        x5, y5 = C.ma_mask_xyz(np.ma.log10(H.integrated_tau_V__g), np.ma.log10(H.integrated_SFRSD__Tg[iT] * 1e6), mask = mask_GAL__g)
        
        x1label = r'$\log\ \tau_V^{\star}$'
        y1label = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x2label = r'$\log\ \tau_V^{\star}(R)$'
        y2label = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x3label = r'$\log\ \frac{\tau_V^\star(R)}{\tau_V^\star(@1HLR)}$'
        y3label = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'
        x4label = r'$\log\ \tau_V^\star(@1HLR)$'
        y4label = r'$\log\ \Sigma_{SFR}^\star(@1HLR)$'
        x5label = r'$\log\ \tau_V^\star (int)$'
        y5label = r'$\log\ \Sigma_{SFR}^\star (int)$'
    
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 0.5])
        kwargs_zbins.update(ylim = [-3.5, 1])
        #kwargs_zbins.update(ylim = [-6.05, 2.05])
        kwargs_zbins.update(cb = False)
        C.plot_zbins(x = x1, y = y1, xlabel = x1label, ylabel = y1label, z = H.zone_dist_HLR__g, **kwargs_zbins)
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        C.plot_zbins(x = x2.flatten(), y = y2.flatten(),  z = H.Rtoplot(x2.shape).flatten(), xlabel = x2label, ylabel = y2label, **kwargs_zbins)    

        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 1])
        kwargs_zbins.update(ylim = [-2.5, 2])
        kwargs_zbins.update(cb = True)
        C.plot_zbins(x = x3.flatten(), y = y3.flatten(), z = H.Rtoplot(x3.shape).flatten(), xlabel = x3label, ylabel = y3label, **kwargs_zbins)
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 3))
        ax.set_axis_on()
        ax.scatter(x4, y4, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x4.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax2 = ax.twinx()
        ax.set_xlabel(x4label)
        ax2.set_ylabel(y4label)
        plt.setp(ax2.get_yticklabels(), visible = False)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x4, y4, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (1, 3))
        ax.set_axis_on()
        ax.scatter(x5, y5, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x5.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax.set_xlabel(x5label)
        ax.set_ylabel(y5label)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x5, y5, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (2, 3))
        ax.set_axis_on()
        ax.plot(H.Rtoplot()[:, 0], x2.mean(axis = 1), label = x2label)
        ax.plot(H.Rtoplot()[:, 0], x3.mean(axis = 1), label = x3label)
        ax.plot(H.Rtoplot()[:, 0], y2.mean(axis = 1), label = y2label.replace('\ [M_\odot yr^{-1} kpc^{-2}]', ''))
        ax.plot(H.Rtoplot()[:, 0], y3.mean(axis = 1), label = y3label)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylim(-3.5, 0.5)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.legend(bbox_to_anchor = (1.0, 1.15), fontsize = 6, frameon = False, ncol = 4)
        ax.grid()
        
        Nbins = 20
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        ax.hist(x1.compressed(), bins = Nbins)
        ax.set_xlabel(x1label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        ax.set_axis_on()
        ax.hist(x2.compressed(), bins = Nbins)
        ax.set_xlabel(x2label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        ax.set_axis_on()
        ax.hist(x3.compressed(), bins = Nbins)
        ax.set_xlabel(x3label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        ax.hist(y1.compressed(), bins = Nbins)
        ax.set_xlabel(y1label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 1))
        ax.set_axis_on()
        ax.hist(y2.compressed(), bins = Nbins)
        ax.set_xlabel(y2label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 2))
        ax.set_axis_on()
        ax.hist(y3.compressed(), bins = Nbins)
        ax.set_xlabel(y3label)

        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
        
        ############## Neb ##############
        ############## Neb ##############
        ############## Neb ##############
        NRows = 3
        NCols = 4
        f = plt.figure()
        kwargs_zbins = default_kwargs_zbins.copy()
        kwargs_zbins.update(f = f)
        #f, axArr = plt.subplots(NRows, NCols)
        #page_size_inches = (NCols * 3, NRows * 1.5)
        page_size_inches = A4Size_inches
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
    
        f.suptitle(suptitle, fontsize = 9)

        tau_V_neb_norm__rg = H.tau_V_neb__Trg[iT] / H.tau_V_neb_oneHLR__Tg[iT]
        aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__Trg[iT] / H.aSFRSD_Ha_oneHLR__Tg[iT]
     
        x1, y1 = C.ma_mask_xyz(np.ma.log10(H.tau_V_neb__g), np.ma.log10(H.SFRSD_Ha__g * 1e6), mask = mask__g) 
        x2, y2 = C.ma_mask_xyz(np.ma.log10(H.tau_V_neb__Trg[iT]), np.ma.log10(H.aSFRSD_Ha__Trg[iT] * 1e6), mask = mask__rg) 
        x3, y3 = C.ma_mask_xyz(np.ma.log10(tau_V_neb_norm__rg), np.ma.log10(aSFRSD_Ha_norm__rg), mask = mask__rg) 
        x4, y4 = C.ma_mask_xyz(np.ma.log10(H.tau_V_neb_oneHLR__Tg[iT]), np.ma.log10(H.aSFRSD_Ha_oneHLR__Tg[iT] * 1e6), mask = mask_GAL__g)
        x5, y5 = C.ma_mask_xyz(np.ma.log10(H.integrated_tau_V_neb__g), np.ma.log10(H.integrated_SFRSD_Ha__g * 1e6), mask = mask_GAL__g)
      
        x1label = r'$\log\ \tau_V^{neb}$'
        y1label = r'$\log\ \Sigma_{SFR}^{H\alpha}\ [M_\odot yr^{-1} kpc^{-2}]$'
        x2label = r'$\log\ \tau_V^{neb}(R)$'
        y2label = r'$\log\ \Sigma_{SFR}^{H\alpha}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x3label = r'$\log\ \frac{\tau_V^{neb}(R)}{\tau_V^{neb}(@1HLR)}$'
        y3label = r'$\log\ \frac{\Sigma_{SFR}^{H\alpha}(R)}{\Sigma_{SFR}^{H\alpha}(@1HLR)}$'
        x4label = r'$\log\ \tau_V^{neb}(@1HLR)$'
        y4label = r'$\log\ \Sigma_{SFR}^{H\alpha}(@1HLR)$'
        x5label = r'$\log\ \tau_V^{neb} (int)$'
        y5label = r'$\log\ \Sigma_{SFR}^{H\alpha} (int)$'
      
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 0.5])
        kwargs_zbins.update(ylim = [-3.5, 1])
        kwargs_zbins.update(cb = False)
        C.plot_zbins(x = x1, y = y1, z = H.zone_dist_HLR__g, xlabel = x1label, ylabel = y1label, **kwargs_zbins)    
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        C.plot_zbins(x = x2.flatten(), y = y2.flatten(), xlabel = x2label, ylabel = y2label, z = H.Rtoplot(x2.shape).flatten(), **kwargs_zbins)    
         
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 1])
        kwargs_zbins.update(ylim = [-2.5, 2])
        kwargs_zbins.update(cb = True)
        C.plot_zbins(x = x3.flatten(), y = y3.flatten(), xlabel = x3label, ylabel = y3label, z = H.Rtoplot(x3.shape).flatten(), **kwargs_zbins)
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 3))
        ax.set_axis_on()
        ax.scatter(x4, y4, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x4.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax2 = ax.twinx()
        ax.set_xlabel(x4label)
        ax2.set_ylabel(y4label)
        plt.setp(ax2.get_yticklabels(), visible = False)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x4, y4, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (1, 3))
        ax.set_axis_on()
        ax.scatter(x5, y5, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x5.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax.set_xlabel(x5label)
        ax.set_ylabel(y5label)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x5, y5, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (2, 3))
        ax.set_axis_on()
        ax.plot(H.Rtoplot()[:, 0], x2.mean(axis = 1), label = x2label)
        ax.plot(H.Rtoplot()[:, 0], x3.mean(axis = 1), label = x3label)
        ax.plot(H.Rtoplot()[:, 0], y2.mean(axis = 1), label = y2label.replace('\ [M_\odot yr^{-1} kpc^{-2}]', ''))
        ax.plot(H.Rtoplot()[:, 0], y3.mean(axis = 1), label = y3label)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylim(-3.5, 0.5)
        ax.set_ylim(-3.5, 0.5)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.legend(bbox_to_anchor = (1.0, 1.15), fontsize = 6, frameon = False, ncol = 4)
        #ax.legend(fontsize = 8, frameon = False)
        ax.grid()

        Nbins = 20
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        ax.hist(x1.compressed(), bins = Nbins)
        ax.set_xlabel(x1label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        ax.set_axis_on()
        ax.hist(x2.compressed(), bins = Nbins)
        ax.set_xlabel(x2label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        ax.set_axis_on()
        ax.hist(x3.compressed(), bins = Nbins)
        ax.set_xlabel(x3label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        ax.hist(y1.compressed(), bins = Nbins)
        ax.set_xlabel(y1label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 1))
        ax.set_axis_on()
        ax.hist(y2.compressed(), bins = Nbins)
        ax.set_xlabel(y2label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 2))
        ax.set_axis_on()
        ax.hist(y3.compressed(), bins = Nbins)
        ax.set_xlabel(y3label)

        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
      
        ############## Mixed ##############
        ############## Mixed ##############
        ############## Mixed ############## 
        NRows = 3
        NCols = 4
        f = plt.figure()
        kwargs_zbins = default_kwargs_zbins.copy()
        kwargs_zbins.update(f = f)
        #f, axArr = plt.subplots(NRows, NCols)
        #page_size_inches = (NCols * 3, NRows * 1.5)
        page_size_inches = A4Size_inches
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)

        f.suptitle(suptitle, fontsize = 9)
           
        tau_V_neb_norm__rg = H.tau_V_neb__Trg[iT] / H.tau_V_neb_oneHLR__Tg[iT]
        aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / H.aSFRSD_oneHLR__Tg[iT]
       
        x1, y1 = C.ma_mask_xyz(np.ma.log10(H.tau_V_neb__g), np.ma.log10(H.SFRSD__Tg[iT] * 1e6), mask = mask__g) 
        x2, y2 = C.ma_mask_xyz(np.ma.log10(H.tau_V_neb__Trg[iT]), np.ma.log10(H.aSFRSD__Trg[iT] * 1e6), mask = mask__rg) 
        x3, y3 = C.ma_mask_xyz(np.ma.log10(tau_V_neb_norm__rg), np.ma.log10(aSFRSD_norm__rg), mask = mask__rg) 
        x4, y4 = C.ma_mask_xyz(np.ma.log10(H.tau_V_neb_oneHLR__Tg[iT]), np.ma.log10(H.aSFRSD_oneHLR__Tg[iT] * 1e6), mask = mask_GAL__g)
        x5, y5 = C.ma_mask_xyz(np.ma.log10(H.integrated_tau_V_neb__g), np.ma.log10(H.integrated_SFRSD__Tg[iT] * 1e6), mask = mask_GAL__g)
 
        x1label = r'$\log\ \tau_V^{neb}$'
        y1label = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x2label = r'$\log\ \tau_V^{neb}(R)$'
        y2label = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x3label = r'$\log\ \frac{\tau_V^{neb}(R)}{\tau_V^{neb}(@1HLR)}$'
        y3label = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'
        x4label = r'$\log\ \tau_V^{neb}(@1HLR)$'
        y4label = r'$\log\ \Sigma_{SFR}^\star(@1HLR)$'
        x5label = r'$\log\ \tau_V^{neb} (int)$'
        y5label = r'$\log\ \Sigma_{SFR}^\star (int)$'

        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 0.5])
        kwargs_zbins.update(ylim = [-3.5, 1])
        kwargs_zbins.update(cb = False)
        C.plot_zbins(x = x1, y = y1, xlabel = x1label, ylabel = y1label, z = H.zone_dist_HLR__g, **kwargs_zbins)
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        C.plot_zbins(x = x2.flatten(), y = y2.flatten(), z = H.Rtoplot(x2.shape).flatten(), xlabel = x2label, ylabel = y2label, **kwargs_zbins)
         
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 1])
        kwargs_zbins.update(ylim = [-2.5, 2])
        kwargs_zbins.update(cb = True)
        C.plot_zbins(x = x3.flatten(), y = y3.flatten(), z = H.Rtoplot(x3.shape).flatten(), xlabel = x3label, ylabel = y3label, **kwargs_zbins)    
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 3))
        ax.set_axis_on()
        ax.scatter(x4, y4, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x4.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax2 = ax.twinx()
        ax.set_xlabel(x4label)
        ax2.set_ylabel(y4label)
        plt.setp(ax2.get_yticklabels(), visible = False)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x4, y4, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (1, 3))
        ax.set_axis_on()
        ax.scatter(x5, y5, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x5.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax.set_xlabel(x5label)
        ax.set_ylabel(y5label)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x5, y5, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (2, 3))
        ax.set_axis_on()
        ax.plot(H.Rtoplot()[:, 0], x2.mean(axis = 1), label = x2label)
        ax.plot(H.Rtoplot()[:, 0], x3.mean(axis = 1), label = x3label)
        ax.plot(H.Rtoplot()[:, 0], y2.mean(axis = 1), label = y2label.replace('\ [M_\odot yr^{-1} kpc^{-2}]', ''))
        ax.plot(H.Rtoplot()[:, 0], y3.mean(axis = 1), label = y3label)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylim(-3.5, 0.5)
        ax.set_ylim(-3.5, 0.5)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.legend(bbox_to_anchor = (1.0, 1.15), fontsize = 6, frameon = False, ncol = 4)
        #ax.legend(fontsize = 8, frameon = False)
        ax.grid()
      
        Nbins = 20
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        ax.hist(x1.compressed(), bins = Nbins)
        ax.set_xlabel(x1label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        ax.set_axis_on()
        ax.hist(x2.compressed(), bins = Nbins)
        ax.set_xlabel(x2label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        ax.set_axis_on()
        ax.hist(x3.compressed(), bins = Nbins)
        ax.set_xlabel(x3label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        ax.hist(y1.compressed(), bins = Nbins)
        ax.set_xlabel(y1label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 1))
        ax.set_axis_on()
        ax.hist(y2.compressed(), bins = Nbins)
        ax.set_xlabel(y2label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 2))
        ax.set_axis_on()
        ax.hist(y3.compressed(), bins = Nbins)
        ax.set_xlabel(y3label)

        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
      
        ############## Mixed ##############
        ############## Mixed ##############
        ############## Mixed ############## 
        NRows = 3
        NCols = 4
        f = plt.figure()
        kwargs_zbins = default_kwargs_zbins.copy()
        kwargs_zbins.update(f = f)
        #f, axArr = plt.subplots(NRows, NCols)
        #page_size_inches = (NCols * 3, NRows * 1.5)
        page_size_inches = A4Size_inches
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        
        f.suptitle(suptitle, fontsize = 9)

        tau_V_norm__rg = H.tau_V__Trg[iT] / H.tau_V_oneHLR__Tg[iT]
        aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__Trg[iT] / H.aSFRSD_Ha_oneHLR__Tg[iT]
       
        x1, y1 = C.ma_mask_xyz(np.ma.log10(H.tau_V__Tg[iT]), np.ma.log10(H.SFRSD_Ha__g * 1e6), mask = mask__g) 
        x2, y2 = C.ma_mask_xyz(np.ma.log10(H.tau_V__Trg[iT]), np.ma.log10(H.aSFRSD_Ha__Trg[iT] * 1e6), mask = mask__rg) 
        x3, y3 = C.ma_mask_xyz(np.ma.log10(tau_V_norm__rg), np.ma.log10(aSFRSD_Ha_norm__rg), mask = mask__rg) 
        x4, y4 = C.ma_mask_xyz(np.ma.log10(H.tau_V_oneHLR__Tg[iT]), np.ma.log10(H.aSFRSD_Ha_oneHLR__Tg[iT] * 1e6), mask = mask_GAL__g)
        x5, y5 = C.ma_mask_xyz(np.ma.log10(H.integrated_tau_V__g), np.ma.log10(H.integrated_SFRSD_Ha__g * 1e6), mask = mask_GAL__g)
      
        x1label = r'$\log\ \tau_V^{\star}$'
        y1label = r'$\log\ \Sigma_{SFR}^{H\alpha}\ [M_\odot yr^{-1} kpc^{-2}]$'
        x2label = r'$\log\ \tau_V^{\star}(R)$'
        y2label = r'$\log\ \Sigma_{SFR}^{H\alpha}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x3label = r'$\log\ \frac{\tau_V^\star(R)}{\tau_V^\star(@1HLR)}$'
        y3label = r'$\log\ \frac{\Sigma_{SFR}^{H\alpha}(R)}{\Sigma_{SFR}^{H\alpha}(@1HLR)}$'
        x4label = r'$\log\ \tau_V^\star(@1HLR)$'
        y4label = r'$\log\ \Sigma_{SFR}^{H\alpha}(@1HLR)$'
        x5label = r'$\log\ \tau_V^\star (int)$'
        y5label = r'$\log\ \Sigma_{SFR}^{H\alpha} (int)$'
       
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 0.5])
        kwargs_zbins.update(ylim = [-3.5, 1])
        kwargs_zbins.update(cb = False)
        C.plot_zbins(x = x1, y = y1, xlabel = x1label, ylabel = y1label,z = H.zone_dist_HLR__g, **kwargs_zbins)    
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        C.plot_zbins(x = x2.flatten(), y = y2.flatten(), z = H.Rtoplot(x2.shape).flatten(), xlabel = x2label, ylabel = y2label, **kwargs_zbins)
         
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 1])
        kwargs_zbins.update(ylim = [-2.5, 2])
        kwargs_zbins.update(cb = True)
        C.plot_zbins(x = x3.flatten(), y = y3.flatten(), z = H.Rtoplot(x3.shape).flatten(), xlabel = x3label, ylabel = y3label,**kwargs_zbins)
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 3))
        ax.set_axis_on()
        ax.scatter(x4, y4, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x4.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax2 = ax.twinx()
        ax.set_xlabel(x4label)
        ax2.set_ylabel(y4label)
        plt.setp(ax2.get_yticklabels(), visible = False)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x4, y4, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (1, 3))
        ax.set_axis_on()
        ax.scatter(x5, y5, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x5.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax.set_xlabel(x5label)
        ax.set_ylabel(y5label)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x5, y5, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (2, 3))
        ax.set_axis_on()
        ax.plot(H.Rtoplot()[:, 0], x2.mean(axis = 1), label = x2label)
        ax.plot(H.Rtoplot()[:, 0], x3.mean(axis = 1), label = x3label)
        ax.plot(H.Rtoplot()[:, 0], y2.mean(axis = 1), label = y2label.replace('\ [M_\odot yr^{-1} kpc^{-2}]', ''))
        ax.plot(H.Rtoplot()[:, 0], y3.mean(axis = 1), label = y3label)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylim(-3.5, 0.5)
        ax.set_ylim(-3.5, 0.5)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.legend(bbox_to_anchor = (1.0, 1.15), fontsize = 6, frameon = False, ncol = 4)
        #ax.legend(fontsize = 8, frameon = False)
        ax.grid()

        Nbins = 20
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        ax.hist(x1.compressed(), bins = Nbins)
        ax.set_xlabel(x1label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        ax.set_axis_on()
        ax.hist(x2.compressed(), bins = Nbins)
        ax.set_xlabel(x2label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        ax.set_axis_on()
        ax.hist(x3.compressed(), bins = Nbins)
        ax.set_xlabel(x3label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        ax.hist(y1.compressed(), bins = Nbins)
        ax.set_xlabel(y1label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 1))
        ax.set_axis_on()
        ax.hist(y2.compressed(), bins = Nbins)
        ax.set_xlabel(y2label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 2))
        ax.set_axis_on()
        ax.hist(y3.compressed(), bins = Nbins)
        ax.set_xlabel(y3label)

        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
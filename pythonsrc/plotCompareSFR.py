#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014 
#                 - 10/Sept/2015
#
from matplotlib.ticker import MultipleLocator
import sys
import CALIFAUtils as C
from scipy import stats as st
#from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.ticker import MultipleLocator
#from matplotlib.ticker import MaxNLocator
from CALIFAUtils.plots import plotOLSbisectorAxis
from CALIFAUtils.plots import plot_linreg_params
from CALIFAUtils.plots import plot_text_ax
from CALIFAUtils.scripts import OLS_bisector
#from CALIFAUtils.plots import plot_zbins
#from CALIFAUtils.objects import runstats
from os.path import basename
from matplotlib import pyplot as plt
import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys
from numpy.ma.extras import mask_cols

mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

def mask_zones_iT(iT, H, args, maskRadiusOk, gals_slice):
    mask__g = np.bitwise_or(np.ma.log10(H.SFRSD__Tg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Tg[iT]).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.SFRSD_Ha__g * 1e6).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.tau_V_neb__g).mask)
    mask__g = np.bitwise_or(mask__g, H.logO3N2_M13__g.mask)
    mask__g = np.bitwise_or(mask__g, np.less(H.reply_arr_by_zones(H.ba_GAL__g), args.bamin))
    mask__g = np.bitwise_or(mask__g, ~maskRadiusOk__g)
    mask__g = np.bitwise_or(mask__g, ~gals_slice__g)
    #mask__g = np.bitwise_or(mask__g, np.less(H.EW_Ha__g, 3.))
    return mask__g

def mask_radius_iT(iT, H, args, maskRadiusOk, gals_slice):
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # mask__rg = np.bitwise_or(np.ma.log10(H.aSFRSD__Trg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Trg[iT]).mask)
    # mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.aSFRSD_Ha__rg * 1e6).mask)
    # mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.tau_V_neb__rg).mask)
    # mask__rg = np.bitwise_or(mask__rg, H.O_O3N2_M13__rg.mask)
    # mask__rg = np.bitwise_or(mask__rg, np.less(H.reply_arr_by_radius(H.ba_GAL__g), args.bamin))
    # mask__rg = np.bitwise_or(mask__rg, ~maskRadiusOk)
    # mask__rg = np.bitwise_or(mask__rg, ~gals_slice)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    mask__rg = mask__rg = np.bitwise_or(~maskRadiusOk, np.less(H.reply_arr_by_radius(H.ba_GAL__g), args.bamin))
    mask__rg = np.bitwise_or(mask__rg, ~gals_slice)
    #mask__rg = np.bitwise_or(mask__rg, np.less(H.EW_Ha__rg, 3.))
    return mask__rg

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

    return parser.parse_args()

if __name__ == '__main__':
    args = parser_args()
    
    C.debug_var(args.debug, args = args)
    
    H = C.H5SFRData(args.hdf5)
    
    minR = 0
    fnamesuffix = '.png'
    
    if args.maskradius is None:
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        minR = args.maskradius
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

    ##########################
    ######### MASKS ##########
    ##########################
    ba_max = args.bamin
    mask_GAL__g = np.bitwise_or(np.zeros_like(H.integrated_EW_Ha__g, dtype = np.bool), np.less(H.ba_GAL__g, ba_max))
    mask_GAL__g = np.bitwise_or(mask_GAL__g, ~gals_slice__integr)   
    
    txt_suptitle = r'$\Longrightarrow$ %s  NGals:%d  $x_Y$(min):%.0f%%  $\tau_V^\star $(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f ' % (gals_txt, N_gals, H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)    
    #tSF__T              = H.tSF__T[0:20]
    tSF__T = H.tSF__T[0:20]
    RRange = H.RRange
    RColor = H.RColor
    SFR__Tg = H.get_data_h5('SFR__Tg')
    SFR_Ha__g = H.get_data_h5('SFR_Ha__g')
    aSFRSD__Trg = H.get_data_h5('aSFRSD__Trg')
    aSFRSD_Ha__Trg = H.get_data_h5('aSFRSD_Ha__Trg')
    SFRSD__Tg = H.get_data_h5('SFRSD__Tg')
    SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
    dist_zone__g = H.get_data_h5('dist_zone__g')
    
    ols_kwargs = dict(pos_x = 0.96, y_pos = 0.07, fs = 8)
     
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(300)
    f.set_size_inches(11.69, 8.27) 
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ SFR_\star(t_\star)\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    f.suptitle(txt_suptitle, fontsize = 11)
    filename = 'SFR_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(SFR__Tg[iT])
            y = np.ma.log10(SFR_Ha__g)
            mask__g = mask_zones_iT(iT, H, args, maskRadiusOk__g, gals_slice__g)
            xm, ym = C.ma_mask_xyz(x, y, mask = mask__g)
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())        
            xran = [-6, 0]
            yran = [-6, 0]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1
            
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # if iT < (H.N_T - 1):
    #     for i in range(iT, H.N_T):
    #         x = np.ma.log10(SFR__Tg[i])
    #         y = np.ma.log10(SFR_Ha__g)
    #         mask__g = mask_zones_iT(iT, H, args, maskRadiusOk__g, gals_slice__g)
    #         xm, ym = C.ma_mask_xyz(x, y, mask = mask__g)
    #         a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
    #         b2[i] = (ym - xm).mean()
    #         Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
    #         Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    xlabel = r'$\log\ t_\star$ [yr]'  
    x = np.log10(tSF__T)
    plot_linreg_params(a, x, xlabel, 'slope',
                       'SFR_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel, r'intercept [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'SFR_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel, r'intercept2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'SFR_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel, 'Rs', 
                       'SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(Rp, x, xlabel, 'Rp', 
                       'SFR_Rp_age%s' % fnamesuffix, 1., 16) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # plot_linreg_params(sigma, x, xlabel, 
    #                    r'$\sigma$', 'SFR_linregress_sigma_age%s' % fnamesuffix)
    # plot_linreg_params(r**2., x, xlabel, 
    #                    r'$r^2$', 'SFR_linregress_sqrcorr_age%s' % fnamesuffix, 1., 16) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    Rs_SFR = np.copy(Rs)
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27) 
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical') 
    filename = 'SFRSD_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)  
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(SFRSD__Tg[iT] * 1e6)
            y = np.ma.log10(SFRSD_Ha__g * 1e6)
            mask__g = mask_zones_iT(iT, H, args, maskRadiusOk__g, gals_slice__g)
            xm, ym = C.ma_mask_xyz(x, y, mask = mask__g)
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
            xran = [-3.5, 1]
            yran = [-3.5, 1]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # if iT < (H.N_T - 1):
    #     for i in range(iT, H.N_T):
    #         x = np.ma.log10(SFRSD__Tg[i] * 1e6)
    #         y = np.ma.log10(SFRSD_Ha__g * 1e6)
    #         mask__g = mask_zones_iT(iT, H, args, maskRadiusOk__g, gals_slice__g)
    #         xm, ym = C.ma_mask_xyz(x, y, mask = mask__g)
    #         a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
    #         b2[i] = (ym - xm).mean()
    #         Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
    #         Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel, 'slope', 
                       'SFRSD_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel, r'intercept [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'SFRSD_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel, r'intercept2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'SFRSD_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel, 'Rs', 
                       'SFRSD_Rs_age%s' % fnamesuffix, 1., 16)
    plot_linreg_params(Rp, x, xlabel, 'Rp', 
                       'SFRSD_Rp_age%s' % fnamesuffix, 1., 16)
    Rs_SFRSD = np.copy(Rs) 
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    pos_y_ini = 0.38
    pos_step = 0.09
    Rfontsize = 12
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27)
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')   
    filename = 'aSFRSD_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)
    NAxes = len(f.axes)
    iT = 0
    jump = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)          
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(aSFRSD__Trg[iT] * 1e6)
            y = np.ma.log10(aSFRSD_Ha__Trg[iT] * 1e6)
            mask__rg = mask_radius_iT(iT, H, args, maskRadiusOk__rg, gals_slice__rg)
            xm, ym = C.ma_mask_xyz(x, y, mask = mask__rg)
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
            xran = [-3.5, 1]
            yran = [-3.5, 1]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # if iT < (H.N_T - 1):
    #     for i in range(iT, H.N_T):
    #         x = np.ma.log10(aSFRSD__Trg[i] * 1e6)
    #         y = np.ma.log10(aSFRSD_Ha__rg * 1e6)
    #         mask__rg = mask_radius_iT(iT, H, args, maskRadiusOk__rg, gals_slice__rg)
    #         xm, ym = C.ma_mask_xyz(x, y, mask = mask__rg)
    #         a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
    #         b2[i] = (ym - xm).mean()
    #         Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
    #         Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
             
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel, 'slope', 
                       'aSFRSD_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel, r'intercept [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'aSFRSD_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel, r'intercept2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'aSFRSD_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel, 'Rs', 
                       'aSFRSD_Rs_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(Rp, x, xlabel, 'Rp', 
                       'aSFRSD_Rp_age%s' % fnamesuffix, 1., 16)
    Rs_aSFRSD = np.copy(Rs) 
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    pos_y_ini = 0.38
    pos_step = 0.09
    Rfontsize = 12
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27)
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'
    ylabel = r'$\log\ \frac{\Sigma_{SFR}^{neb}(R)}{\Sigma_{SFR}^{neb}(@1HLR)}$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    fnamesuftmp = '_norm%s' % fnamesuffix   
    filename = 'aSFRSD_linregress_report%s' % fnamesuftmp
    C.debug_var(args.debug, filename = filename)
    NAxes = len(f.axes)
    iT = 0
    jump = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)          
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j]
            age = tSF__T[iT]
            aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / H.aSFRSD_oneHLR__Tg[iT]
            aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__Trg[iT] / H.aSFRSD_Ha_oneHLR__Tg[iT]
            xran = [-1.5, 1.5]
            yran = [-1.5, 1.5]
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # binsx = np.linspace(-4.5, 1., 51)
            # binsy = np.linspace(min(ym),max(ym), 51)
            # density_contour(xm, ym, binsx, binsy, ax=ax)
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            x = np.ma.log10(aSFRSD_norm__rg)
            y = np.ma.log10(aSFRSD_Ha_norm__rg)
            mask__rg = mask_radius_iT(iT, H, args, maskRadiusOk__rg, gals_slice__rg)
            xm, ym = C.ma_mask_xyz(x, y, mask = mask__rg)
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.6)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # if iT < (H.N_T - 1):
    #     for i in range(iT, H.N_T):
    #         aSFRSD_norm__rg = H.aSFRSD__Trg[i] / H.aSFRSD_oneHLR__Tg[i]
    #         aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__rg / H.aSFRSD_Ha_oneHLR__g
    #         x = np.ma.log10(aSFRSD_norm__rg)
    #         y = np.ma.log10(aSFRSD_Ha_norm__rg)
    #         mask__rg = mask_radius_iT(iT, H, args, maskRadiusOk__rg, gals_slice__rg)
    #         xm, ym = C.ma_mask_xyz(x, y, mask = mask__rg)
    #         a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
    #         b2[i] = (ym - xm).mean()
    #         Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
    #         Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel, 'slope', 
                       'aSFRSD_linregress_slope_age%s' % fnamesuftmp, 1., 16) 
    plot_linreg_params(b, x, xlabel, r'intercept [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'aSFRSD_linregress_intercep_age%s' % fnamesuftmp, 0., 16)
    plot_linreg_params(b2, x, xlabel, r'intercept2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'aSFRSD_linregress_intercep2_age%s' % fnamesuftmp, 0., 16)
    plot_linreg_params(Rs, x, xlabel, 'Rs', 
                       'aSFRSD_Rs_age%s' % fnamesuftmp, 1., 16) 
    plot_linreg_params(Rp, x, xlabel, 'Rp', 
                       'aSFRSD_Rp_age%s' % fnamesuftmp, 1., 16)
    Rs_aSFRSD_norm = np.copy(Rs) 
    ############# Integrated #############
    ############# Integrated #############
    ############# Integrated #############
    integrated_SFR__Tg = H.integrated_SFR__Tg
    integrated_SFR_Ha__g = H.integrated_SFR_Ha__g
    integrated_SFRSD__Tg = H.integrated_SFRSD__Tg
    integrated_SFRSD_Ha__g = H.integrated_SFRSD_Ha__g
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(300)
    f.set_size_inches(11.69, 8.27) 
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ SFR_\star^{int}(t_\star)\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}^{int}\ [M_\odot yr^{-1}]$'
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    filename = 'integrated_SFR_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)   
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(integrated_SFR__Tg[iT])
            y = np.ma.log10(integrated_SFR_Ha__g)
            xm, ym = C.ma_mask_xyz(x, y, mask = mask_GAL__g)
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'integrated SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())        
            xran = [-5, 2]
            yran = [-5, 2]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # if iT < (H.N_T - 1):
    #     for i in range(iT, H.N_T):
    #         x = np.ma.log10(integrated_SFR__Tg[i])
    #         y = np.ma.log10(integrated_SFR_Ha__g)
    #         xm, ym = C.ma_mask_xyz(x, y, mask = mask_GAL__g)
    #         a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
    #         b2[i] = (ym - xm).mean()
    #         Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
    #         Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel, 'slope', 
                       'integrated_SFR_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel, r'intercept [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'integrated_SFR_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel, r'intercept2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'integrated_SFR_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel, 'Rs', 
                       'integrated_SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(Rp, x, xlabel, 'Rp', 
                       'integrated_SFR_Rp_age%s' % fnamesuffix, 1., 16)
    Rs_SFR_int = np.copy(Rs) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # plot_linreg_params(sigma, x, xlabel, 
    #                    r'$\sigma$', 'SFR_linregress_sigma_age%s' % fnamesuffix)
    # plot_linreg_params(r**2., x, xlabel, 
    #                    r'$r^2$', 'SFR_linregress_sqrcorr_age%s' % fnamesuffix, 1., 16) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27) 
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(int, t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(int)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    filename = 'integrated_SFRSD_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)   
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)  
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(integrated_SFRSD__Tg[iT] * 1e6)
            y = np.ma.log10(integrated_SFRSD_Ha__g * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = mask_GAL__g)
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'integrated SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
            xran = [-5, 0]
            yran = [-5, 0]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # if iT < (H.N_T - 1):
    #     for i in range(iT, H.N_T):
    #         x = np.ma.log10(integrated_SFRSD__Tg[i] * 1e6)
    #         y = np.ma.log10(integrated_SFRSD_Ha__g * 1e6)
    #         xm, ym = C.ma_mask_xyz(x, y, mask = mask_GAL__g)
    #         a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
    #         b2[i] = (ym - xm).mean()
    #         Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
    #         Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel, 'slope', 
                       'integrated_SFRSD_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel, r'intercept [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'integrated_SFRSD_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel, r'intercept2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 
                       'integrated_SFRSD_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel, 'Rs', 
                       'integrated_SFRSD_Rs_age%s' % fnamesuffix, 1., 16)
    plot_linreg_params(Rp, x, xlabel,'Rp', 
                       'integrated_SFRSD_Rp_age%s' % fnamesuffix, 1., 16)
    Rs_SFRSD_int = np.copy(Rs)

    kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS')
    ols_kwargs = dict(c = 'r', pos_x = 0.98, pos_y = 0.01, fs = 10, rms = True, text = True, kwargs_plot = kwargs_ols_plot)

    from matplotlib.pyplot import MaxNLocator
    f = plt.figure()
    f.set_size_inches(7, 10)
    iT = 11
    tSF = tSF__T[iT]
    txt_suptitle = r'$\mathrm{Correlac\c\~ao\ da\ taxa\ de\ formac\c\~ao\ estelar}$'
    xlabel = r'$\log\ t_\star$ [yr]'
    ylabel = r'$R_s$'
    grid_shape = (4,2)
    ax = plt.subplot2grid(grid_shape, loc = (0,0), rowspan = 2, colspan = 2)
    ax.set_title(txt_suptitle)
    ax.plot(np.log10(tSF__T), Rs_SFR, 'b--', label = 'SFR')
    ax.plot(np.log10(tSF__T), Rs_SFRSD, 'r--', label = r'$\Sigma_{SFR}$')
    ax.plot(np.log10(tSF__T), Rs_aSFRSD, 'k-', label = r'$\Sigma_{SFR}(R)$')
    ax.plot(np.log10(tSF__T), Rs_aSFRSD_norm, 'g--', label = r'$\frac{\Sigma_{SFR}}{\Sigma_{SFR}(@1HLR)}$')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(6,8)
    ax.set_ylim(0,1)
    ax.legend(bbox_to_anchor = (0.3, 1), fontsize = 12, frameon = False)
    xlim_inf, xlim_sup = ax.get_xlim()
    ylim_inf, ylim_sup = ax.get_ylim()
    x_pos = np.log10(tSF__T[iT])
    y_pos = Rs_aSFRSD[iT]
    arrow_size_x = 0.2
    arrow_size_y = 0.15
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.) 
    ax.annotate('%.2f Myr' % (tSF / 1e6),
        xy = (x_pos, y_pos), xycoords = 'data',
        xytext = (x_pos + arrow_size_x, y_pos + arrow_size_y),
        textcoords = 'data',
        verticalalignment = 'top', horizontalalignment = 'left',
        bbox = textbox,
        arrowprops = dict(
                        arrowstyle = '->',
                        color = 'k',
                        connectionstyle = 'angle3,angleA=0,angleB=90',
                    ),
    )
    
    bins = (30, 30)
    cmap = 'Blues'
    ax = plt.subplot2grid(grid_shape, loc = (2,0))
    xlabel = r'$\log\ SFR_\star(t_\star)\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
    x = np.ma.log10(SFR__Tg[iT])
    y = np.ma.log10(SFR_Ha__g)
    mask__g = mask_zones_iT(iT, H, args, maskRadiusOk__g, gals_slice__g)
    xm, ym = C.ma_mask_xyz(x, y, mask = mask__g)
    age = tSF__T[iT]
    C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
    #print 'SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())        
    xran = [-6, 0]
    yran = [-6, 0]
    h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
    X, Y = np.meshgrid(xedges, yedges)
    im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
    ax.set_xlim(xran)
    ax.set_ylim(yran)
    #scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
    a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
    b2 = (ym - xm).mean()
    Y2 = xm + b2
    Yrms = (ym - Y2).std()
    ax.plot(ax.get_xlim(), np.asarray(ax.get_xlim()) + b2, c = 'b', ls = '--', lw = 2)
    #ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
    txt = r'(1.000, %.3f, %.3f)' % (b2, Yrms)
    plot_text_ax(ax, txt, ols_kwargs['pos_x'], ols_kwargs['pos_y'] + ols_kwargs['fs'] / 100., ols_kwargs['fs'], 'bottom', 'right', color = 'b')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--', lw = 2)
        
    ax = plt.subplot2grid(grid_shape, loc = (2,1))
    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$' 
    x = np.ma.log10(SFRSD__Tg[iT] * 1e6)
    y = np.ma.log10(SFRSD_Ha__g * 1e6)
    mask__g = mask_zones_iT(iT, H, args, maskRadiusOk__g, gals_slice__g)
    xm, ym = C.ma_mask_xyz(x, y, mask = mask__g)
    age = tSF__T[iT]
    C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
    #print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
    xran = [-3.5, 1]
    yran = [-3.5, 1]
    ax.set_xlim(xran)
    ax.set_ylim(yran)
    h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
    X, Y = np.meshgrid(xedges, yedges)
    im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
    #scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
    a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
    b2 = (ym - xm).mean()
    Y2 = xm + b2
    Yrms = (ym - Y2).std()
    ax.plot(ax.get_xlim(), np.asarray(ax.get_xlim()) + b2, c = 'b', ls = '--', lw = 2)
    #ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
    txt = r'(1.000, %.3f, %.3f)' % (b2, Yrms)
    plot_text_ax(ax, txt, ols_kwargs['pos_x'], ols_kwargs['pos_y'] + ols_kwargs['fs'] / 100., ols_kwargs['fs'], 'bottom', 'right', color = 'b')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--', lw = 2)
        
    ax = plt.subplot2grid(grid_shape, loc = (3,0))
    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    x = np.ma.log10(aSFRSD__Trg[iT] * 1e6)
    y = np.ma.log10(aSFRSD_Ha__Trg[iT] * 1e6)
    mask__rg = mask_radius_iT(iT, H, args, maskRadiusOk__rg, gals_slice__rg)
    xm, ym = C.ma_mask_xyz(x, y, mask = mask__rg)
    age = tSF__T[iT]
    C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
    xran = [-3.5, 1]
    yran = [-3.5, 1]
    h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
    X, Y = np.meshgrid(xedges, yedges)
    im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
    ax.set_xlim(xran)
    ax.set_ylim(yran)
    #scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
    a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
    b2 = (ym - xm).mean()
    Y2 = xm + b2
    Yrms = (ym - Y2).std()
    ax.plot(ax.get_xlim(), np.asarray(ax.get_xlim()) + b2, c = 'b', ls = '--', lw = 2)
    #ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
    txt = r'(1.000, %.3f, %.3f)' % (b2, Yrms)
    plot_text_ax(ax, txt, ols_kwargs['pos_x'], ols_kwargs['pos_y'] + ols_kwargs['fs'] / 100., ols_kwargs['fs'], 'bottom', 'right', color = 'b')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--', lw = 2)

    ax = plt.subplot2grid(grid_shape, loc = (3,1))
    xlabel = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'
    ylabel = r'$\log\ \frac{\Sigma_{SFR}^{neb}(R)}{\Sigma_{SFR}^{neb}(@1HLR)}$' 
    aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / H.aSFRSD_oneHLR__Tg[iT]
    aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__Trg[iT] / H.aSFRSD_Ha_oneHLR__Tg[iT]
    xran = [-1.5, 1.5]
    yran = [-1.5, 1.5]
    ax.set_xlim(xran)
    ax.set_ylim(yran)
    x = np.ma.log10(aSFRSD_norm__rg)
    y = np.ma.log10(aSFRSD_Ha_norm__rg)
    mask__rg = mask_radius_iT(iT, H, args, maskRadiusOk__rg, gals_slice__rg)
    xm, ym = C.ma_mask_xyz(x, y, mask = mask__rg)
    h, xedges, yedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = bins, range = [xran, yran])
    X, Y = np.meshgrid(xedges, yedges)
    im = ax.pcolormesh(X, Y, h.T, cmap = cmap)
    #scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.6)
    a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
    b2 = (ym - xm).mean()
    Y2 = xm + b2
    Yrms = (ym - Y2).std()
    ax.plot(ax.get_xlim(), np.asarray(ax.get_xlim()) + b2, c = 'b', ls = '--', lw = 2)
    #ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
    txt = r'(1.000, %.3f, %.3f)' % (b2, Yrms)
    plot_text_ax(ax, txt, ols_kwargs['pos_x'], ols_kwargs['pos_y'] + ols_kwargs['fs'] / 100., ols_kwargs['fs'], 'bottom', 'right', color = 'b')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--', lw = 2)

    f.subplots_adjust(hspace = 0.4, wspace = 0.4)
    f.savefig('Rs_allSFR.pdf')
    plt.close(f)
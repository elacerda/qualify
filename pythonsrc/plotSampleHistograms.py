#!/usr/bin/python
#
# Lacerda@Granada - 02/Dec/2015
#
#from matplotlib.ticker import MultipleLocator
from matplotlib.backends.backend_pdf import PdfPages
from CALIFAUtils.plots import plot_text_ax, next_row_col
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Polygon
from matplotlib import pyplot as plt
from os.path import basename
import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

default_sc_kwargs = dict(marker = 'o', s = 10, edgecolor = 'none', cmap = 'Spectral')
default_rs_kwargs = dict(smooth = True, sigma = 1.2, frac = 0.02)
default_ols_plot_kwargs = dict(c = 'r', ls = '--', lw = 2, label = 'OLS')
default_ols_kwargs = dict(c = 'r', pos_x = 0.98, pos_y = 0.01, fs = 15, rms = True, text = True, kwargs_plot = default_ols_plot_kwargs)

A4Size_inches = [ 8.267, 11.692 ]
LetterSize_inches = [ 8.5, 11 ]

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'hdf5' : None,
        'output' : None,
        'itSF' : 11,
        'itZ' : -1,
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
                        default = default['hdf5'])
    parser.add_argument('--itSF', '-T',
                        help = 'SF age index.',
                        metavar = '',
                        type = int,
                        default = default['itSF'])
    parser.add_argument('--itZ', '-U',
                        help = 'Stellar Metallicity age index.',
                        metavar = '',
                        type = int,
                        default = default['itZ'])
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
    #Hallsample = C.H5SFRData('/Users/lacerda/dev/astro/EmissionLines/SFR_all416gal_rgbcuts_3HLR.h5')
    iT = args.itSF
    iU = args.itZ

    fnamesuffix = '.pdf'
    if args.output is not None:
        fnamesuffix = '%s%s' % (args.output, fnamesuffix)
    minR = args.maskradius

    if minR is None:
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        maxR = H.Rbin__r[-1]
        maskRadiusOk__g = (H.zone_dist_HLR__g >= minR)
        maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * (H.RbinCenter__r >= minR)).T
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
    
    print ((~mask__g).sum()),((~mask__rg).sum()), NgalsOkZones, NgalsOkRbins
 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # count = 0
    # fraclim = 5
    # gals_ok = []
    #  
    # for gal in np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]):
    #     igal_all = H.califaIDs_all.tolist().index(gal)
    #     igal = H.califaIDs.tolist().index(gal)
    #     mask__z = H.get_prop_gal(mask__g, gal)
    #     mask__r = H.get_prop_gal(mask__rg, gal)
    #     NzOk = (~(mask__z)).astype(int).sum()
    #     Nz = len(mask__z)
    #     NzOkfrac = 100.*NzOk / Nz 
    #       
    #     if (NzOk < fraclim):# or (NrOk < fraclim):
    #         count += 1
    #         pref = '>>>'
    #     else:
    #         gals_ok.append(gal)
    #         pref = ''
    #     print pref, gal, Nz, mask__z.sum(), NzOk, '%.2f%%'%NzOkfrac
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    D_gals = {
        gal : {
            'morfType' : H.morfType_GAL__g[H.califaIDs_all.tolist().index(gal)], 
            'Mcor' : H.Mcor_GAL__g[H.califaIDs_all.tolist().index(gal)], 
            'McorSD' : H.McorSD_GAL__g[H.califaIDs_all.tolist().index(gal)], 
            'atflux' : H.at_flux_GAL__g[H.califaIDs_all.tolist().index(gal)]
        } for gal in np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g])
    }
    
    props = { 
        'Mcor' : dict(v = np.ma.log10(H.Mcor__g), 
                      vm = np.ma.masked_array(np.ma.log10(H.Mcor__g), mask = mask__g), 
                      v_int = np.ma.log10(H.Mcor_GAL__g), 
                      vm_int = np.ma.masked_array(np.ma.log10(H.Mcor_GAL__g), mask = mask_GAL__g), 
                      label = r'$\log\ M_\star$ [ $M_\odot$ ]'),
        'McorSD' : dict(v = np.ma.log10(H.McorSD__Tg[iT]),
                        vm = np.ma.masked_array(np.ma.log10(H.McorSD__Tg[iT]), mask = mask__g), 
                        v_int = np.ma.log10(H.McorSD_GAL__g),
                        vm_int = np.ma.masked_array(np.ma.log10(H.McorSD_GAL__g), mask = mask_GAL__g), 
                        label = r'$\log\ \mu_\star$ [ $M_\odot\ pc^{-2}$ ]'), 
        'at_flux' : dict(v = H.at_flux__g, 
                         vm = np.ma.masked_array(H.at_flux__g, mask = mask__g), 
                         v_int = H.at_flux_GAL__g, 
                         vm_int = np.ma.masked_array(H.at_flux_GAL__g, mask = mask_GAL__g), 
                         label = r'$\langle\log\ t_\star\rangle$ [ anos ]'),
        'ba' : dict(v = H.reply_arr_by_zones(H.ba_GAL__g), 
                    vm = np.ma.masked_array(H.reply_arr_by_zones(H.ba_GAL__g), mask = mask__g), 
                    v_int = H.ba_GAL__g, 
                    vm_int = np.ma.masked_array(H.ba_GAL__g, mask = mask_GAL__g), 
                    label = r'$b/a$'),
        'Zneb' : dict(v = H.logO3N2_M13__g - 8.69, 
                      vm = np.ma.masked_array(H.logO3N2_M13__g - 8.69, mask = mask__g), 
                      v_int = H.integrated_logO3N2_M13__g - 8.69, 
                      vm_int = np.ma.masked_array(H.integrated_logO3N2_M13__g - 8.69, mask = mask_GAL__g), 
                      label = r'$12\ +\ \log (O/H)$ [ $Z_\odot$ ]'),
        'tauV' : dict(v = H.tau_V__Tg[iT], 
                      vm = np.ma.masked_array(H.tau_V__Tg[iT], mask = mask__g), 
                      v_int = H.integrated_tau_V__g, 
                      vm_int = np.ma.masked_array(H.integrated_tau_V__g, mask = mask_GAL__g), 
                      label = r'$\tau_V^\star$'),
        'tauVNeb' : dict(v = H.tau_V_neb__g, 
                         vm = np.ma.masked_array(H.tau_V_neb__g, mask = mask__g), 
                         v_int = H.integrated_tau_V_neb__g, 
                         vm_int = np.ma.masked_array(H.integrated_tau_V_neb__g, mask = mask_GAL__g), 
                         label = r'$\tau_V^{neb}$'),
        'alogZmass' : dict(v = H.alogZ_mass__Ug[iU], 
                           vm = np.ma.masked_array(H.alogZ_mass__Ug[iU], mask = mask__g), 
                           v_int = H.alogZ_mass_GAL__Ug[iU], 
                           vm_int = np.ma.masked_array(H.alogZ_mass_GAL__Ug[iU], mask = mask_GAL__g), 
                           label =  r'$\langle \log\ Z_\star \rangle_M (R)$ (t < %.2f Gyr) [ $Z_\odot$ ]' % (H.tZ__U[iU] / 1e9)),
        'xY' : dict(v = H.x_Y__Tg[iT], 
                    vm = np.ma.masked_array(H.x_Y__Tg[iT], mask = mask__g), 
                    v_int = H.integrated_x_Y__Tg[iT], 
                    vm_int = np.ma.masked_array(H.integrated_x_Y__Tg[iT], mask = mask_GAL__g), 
                    label =  r'$x_Y$'),
    }

    ######################### Morf #########################
    ######################### Morf #########################
    ######################### Morf #########################
    
    maskSaSab = (H.morfType_GAL__g < 2) & (H.morfType_GAL__g >= 0)
    maskSb = (H.morfType_GAL__g == 2.)
    maskSbc = (H.morfType_GAL__g == 3.)
    maskScScd = (H.morfType_GAL__g < 6.) & (H.morfType_GAL__g >= 4.)
    maskSdSmIrr = (H.morfType_GAL__g >= 6) 
    mask_morf = [ maskSaSab, maskSb, maskSbc, maskScScd, maskSdSmIrr  ]
    types_morf = [ 0., 1., 2., 3., 4. ]
    Ntypes = len(types_morf)
    color_morf = [ 'orange', 'green', '#00D0C9', '#0076C9', 'blue' ]
    label_morf = [ '', 'Sa', 'Sb', 'Sbc', 'Sc', 'Sd', '' ]
    tagname_morf = [ l.replace(' + ', '') for l in label_morf ]
    tickpos = np.arange(types_morf[0] - 1,types_morf[-1] + 2, 1)
    
    if args.dryrun: sys.exit()
    
    f = plt.figure()
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    suptitle_R = r'Nzones:%d(%d gals) $t_{SF}$:%.2fMyr ' % ((~mask__g).sum(), NgalsOkZones, (H.tSF__T[iT] / 1e6))
    NRows = 2
    NCols = 2
    page_size_inches = [10, 8]
    f.set_size_inches(page_size_inches)
    f.set_dpi(100)
    grid_shape = (NRows, NCols)
    #f.suptitle(suptitle_R, fontsize = 15)
    row, col = 0, 0
    props_order = [ 'Mcor', 'at_flux', 'xY', 'ba', ] # 'alogZmass', 'tauV', 'tauVNeb', 'Zneb',  ]
    mode = 'integrated_mask'
    for k in props_order:
        zones_tot = 0
        ax = plt.subplot2grid(grid_shape, loc = (row, col))
        p = props[k]
        Nboxes = Ntypes
        if mode.find('integrated') >= 0:
            if mode.find('mask') >= 0:
                v = p['vm_int']
            else:
                v = p['v_int']
            boxes_data = []
            for i_mt, msk in enumerate(mask_morf):
                vm = v[msk].compressed()
                boxes_data.append(vm)
                msk_aux = msk
                if isinstance(v, np.ma.MaskedArray):
                    msk_aux = np.bitwise_and(msk_aux, ~(v.mask))
                Nmsk = msk_aux.astype(int).sum()
                print mode, k, types_morf[i_mt], np.percentile(vm, [ 50, 68, 95, 99 ]), vm.mean(), vm.std()    
                if (row == 0) and (col == 0):
                    plot_text_ax(ax, '%d' % Nmsk, (i_mt + 1) / 6., 1.02, 20, 'bottom', 'center', color_morf[i_mt])
        else:
            if mode.find('mask') >= 0:
                if k == 'ba':
                    v = p['vm_int']
                    boxes_data = [ (v[m]).compressed() if isinstance(v, np.ma.MaskedArray) else v[m] for m in mask_morf ]
                else:
                    v = p['vm']
                    boxes_data = [ (v[H.reply_arr_by_zones(m)]).compressed() if isinstance(v, np.ma.MaskedArray) 
                        else v[H.reply_arr_by_zones(m)] for m in mask_morf 
                    ]
            else:
                if k == 'ba':
                    v = p['v_int']
                    boxes_data = [ (v[m]).compressed() if isinstance(v, np.ma.MaskedArray) else v[m] for m in mask_morf ]
                else:
                    v = p['v']
                    boxes_data = [ (v[H.reply_arr_by_zones(m)]).compressed() if isinstance(v, np.ma.MaskedArray) 
                        else v[H.reply_arr_by_zones(m)] for m in mask_morf 
                    ]
            for i_mt, msk in enumerate(mask_morf):
                if k == 'ba':
                    if mode.find('mask') >= 0:
                        v = p['vm']
                    else:
                        v = p['v']
                msk_aux = np.bitwise_and(H.reply_arr_by_zones(msk), ~(v.mask))    
                vm = v[msk_aux].compressed()
                Nmsk = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[msk_aux]))
                N_zones_bin = msk_aux.astype(int).sum() 
                zones_tot += N_zones_bin 
                print mode, k, types_morf[i_mt], np.percentile(vm, [ 50, 68, 95, 99 ]), vm.mean(), vm.std()
                if (row == 0) and (col == 0):
                    plot_text_ax(ax, '%d' % Nmsk, (i_mt + 1) / 6., 0.02, 20, 'bottom', 'center', color_morf[i_mt])
                    plot_text_ax(ax, '%d' % N_zones_bin, (i_mt + 1) / 6., 1.02, 20, 'bottom', 'center', color_morf[i_mt])
        bp = plt.boxplot(boxes_data, sym='+', whis = [0.3,99.7], positions = types_morf)
        plt.setp(bp['boxes'], color='black')
        plt.setp(bp['whiskers'], color='black')
        plt.setp(bp['fliers'], color='red', marker='+')
        medians = list(range(Nboxes))
        for i, m in enumerate(mask_morf):
            box = bp['boxes'][i]
            boxX = []
            boxY = []
            for j in range(5):
                xb = box.get_xdata()[j]
                yb = box.get_ydata()[j]
                boxX.append(xb)
                boxY.append(yb)
            boxCoords = list(zip(boxX, boxY))
            boxPolygon = Polygon(boxCoords, facecolor=color_morf[i])
            ax.add_patch(boxPolygon)            
            med = bp['medians'][i]
            medianX = []
            medianY = []
            for j in range(2):
                medianX.append(med.get_xdata()[j])
                medianY.append(med.get_ydata()[j])
                ax.plot(medianX, medianY, 'k')
                medians[i] = medianY[0]
            ax.plot([np.average(med.get_xdata())], [np.average(boxes_data[i])],
                     color='w', marker='*', markersize = 15, markeredgecolor='k')
        m_int_find = mode.find('integrated')
        if k == 'McorSD':
            ax.set_ylim(0,5)
        if k == 'Mcor':
            ax.set_ylim(8.5,12)
        if k == 'ba':
            ax.set_ylim(0, 1)    
        elif k == 'xY':
            if m_int_find >= 0:
                ax.set_ylim(0, .5)
            else:
                ax.set_ylim(0, 1.)
            #ax.set_ylim(0, .5)
        elif k == 'at_flux':
            if m_int_find >= 0:
                ax.set_ylim(8, 10)
            else:
                ax.set_ylim(6, 10)
        ax.yaxis.set_major_locator(MaxNLocator(4))
        #ax.yaxis.label.set_size(15)
        ax.set_ylabel(p['label'])
        ax.set_xticks(tickpos)
        ax.set_xticklabels(label_morf)
        ax.xaxis.set_ticks_position('none')
        #ax.tick_params(axis='both', which='major', labelsize=12)
        for j, l in enumerate(ax.get_xticklabels()):
            if j == 0 or j == 6:
                continue
            l.set_color(color_morf[j - 1])
        #ax.label_params(axis = 'both', size = 10)
        if col == (NCols - 1): 
            col = 0 
            row += 1
            if (row == NRows):
                row = 0
        else: 
            col += 1
        print k, zones_tot
    f.subplots_adjust(bottom = 0.05, top = 0.95, hspace = 0.2, wspace = 0.28, right = 0.98, left = 0.1)
    plt.close(f)
    fname = 'sample%s' % fnamesuffix
    f.savefig(fname)
    
    ##################################
    Hall = H = C.H5SFRData('/Users/lacerda/dev/astro/EmissionLines/SFR_all305spiral_qualificacao.h5')
    mask_morf_all = (Hall.morfType_GAL__g >= 0)
    mask_all__g = np.less(Hall.reply_arr_by_zones(H.ba_GAL__g), ba_max)
    mask_all__rg = np.less(H.reply_arr_by_radius(Hall.ba_GAL__g), ba_max)
    mask_GAL_all__g = np.bitwise_or(np.zeros_like(Hall.integrated_EW_Ha__g, dtype = np.bool), np.less(Hall.ba_GAL__g, ba_max))

    NgalsOkZones = len(np.unique(Hall.reply_arr_by_zones(Hall.califaIDs)[~mask_all__g]))  
    NgalsOkRbins = len(np.unique(Hall.reply_arr_by_radius(Hall.califaIDs_all)[~mask_all__rg]))
    print ((~mask_all__g).sum()),((~mask_all__rg).sum()), NgalsOkZones, NgalsOkRbins

    props_all = { 
        'Mcor' : dict(v = np.ma.log10(Hall.Mcor__g), 
                      vm = np.ma.masked_array(np.ma.log10(Hall.Mcor__g), mask = mask_all__g), 
                      v_int = np.ma.log10(Hall.Mcor_GAL__g), 
                      vm_int = np.ma.masked_array(np.ma.log10(Hall.Mcor_GAL__g), mask = mask_GAL_all__g), 
                      label = r'$\log\ M_\star$ [ $M_\odot$ ]'),
        'McorSD' : dict(v = np.ma.log10(Hall.McorSD__Tg[iT]),
                        vm = np.ma.masked_array(np.ma.log10(Hall.McorSD__Tg[iT]), mask = mask_all__g), 
                        v_int = np.ma.log10(Hall.McorSD_GAL__g),
                        vm_int = np.ma.masked_array(np.ma.log10(Hall.McorSD_GAL__g), mask = mask_GAL_all__g), 
                        label = r'$\log\ \mu_\star$ [ $M_\odot\ pc^{-2}$ ]'), 
        'at_flux' : dict(v = Hall.at_flux__g, 
                         vm = np.ma.masked_array(Hall.at_flux__g, mask = mask_all__g), 
                         v_int = Hall.at_flux_GAL__g, 
                         vm_int = np.ma.masked_array(Hall.at_flux_GAL__g, mask = mask_GAL_all__g), 
                         label = r'$\langle\log\ t_\star\rangle$ [ anos ]'),
        'ba' : dict(v = Hall.reply_arr_by_zones(Hall.ba_GAL__g), 
                    vm = np.ma.masked_array(Hall.reply_arr_by_zones(Hall.ba_GAL__g), mask = mask_all__g), 
                    v_int = Hall.ba_GAL__g, 
                    vm_int = np.ma.masked_array(Hall.ba_GAL__g, mask = mask_GAL_all__g), 
                    label = r'$b/a$'),
        'xY' : dict(v = Hall.x_Y__Tg[iT], 
                    vm = np.ma.masked_array(Hall.x_Y__Tg[iT], mask = mask_all__g), 
                    v_int = Hall.integrated_x_Y__Tg[iT], 
                    vm_int = np.ma.masked_array(Hall.integrated_x_Y__Tg[iT], mask = mask_GAL_all__g), 
                    label =  r'$x_Y$'),
        'Zneb' : dict(v = Hall.logO3N2_M13__g - 8.69, 
                      vm = np.ma.masked_array(Hall.logO3N2_M13__g - 8.69, mask = mask_all__g), 
                      v_int = Hall.integrated_logO3N2_M13__g - 8.69, 
                      vm_int = np.ma.masked_array(Hall.integrated_logO3N2_M13__g - 8.69, mask = mask_GAL_all__g), 
                      label = r'$12\ +\ \log (O/H)$ [ $Z_\odot$ ]'),
        'tauV' : dict(v = Hall.tau_V__Tg[iT], 
                      vm = np.ma.masked_array(Hall.tau_V__Tg[iT], mask = mask_all__g), 
                      v_int = Hall.integrated_tau_V__g, 
                      vm_int = np.ma.masked_array(Hall.integrated_tau_V__g, mask = mask_GAL_all__g), 
                      label = r'$\tau_V^\star$'),
        'tauVNeb' : dict(v = Hall.tau_V_neb__g, 
                         vm = np.ma.masked_array(Hall.tau_V_neb__g, mask = mask_all__g), 
                         v_int = Hall.integrated_tau_V_neb__g, 
                         vm_int = np.ma.masked_array(Hall.integrated_tau_V_neb__g, mask = mask_GAL_all__g), 
                         label = r'$\tau_V^{neb}$'),
        'alogZmass' : dict(v = Hall.alogZ_mass__Ug[iU], 
                           vm = np.ma.masked_array(Hall.alogZ_mass__Ug[iU], mask = mask_all__g), 
                           v_int = Hall.alogZ_mass_GAL__Ug[iU], 
                           vm_int = np.ma.masked_array(Hall.alogZ_mass_GAL__Ug[iU], mask = mask_GAL_all__g), 
                           label =  r'$\langle \log\ Z_\star \rangle_M (R)$ (t < %.2f Gyr) [ $Z_\odot$ ]' % (Hall.tZ__U[iU] / 1e9)),
    }

    f = plt.figure()
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    NRows = 2
    NCols = 2
    page_size_inches = [10, 8]
    f.set_size_inches(page_size_inches)
    f.set_dpi(100)
    grid_shape = (NRows, NCols)
    prop_histo = [ 'McorSD', 'at_flux', 'xY', 'ba' ]
    bins = 30
    row, col = 0, 0
    for k in prop_histo:
        p = props[k]
        pall = props_all[k]
        x = p['vm'].compressed()
        xall = pall['vm'][Hall.reply_arr_by_zones(mask_morf_all)].compressed()
        ax = plt.subplot2grid(grid_shape, loc = (row, col))
        ax.set_axis_on()
        ax.set_xlabel(p['label'])
        print k, len(x), len(xall)
        c = 'r'
        pos_x = 0.86
        va = 'top'
        ha = 'right'
        if k == 'at_flux' or k == 'ba':
            pos_x = 0.02
            ha = 'left'
        ax.hist(xall, bins = bins, color = c, align = 'mid', alpha = 0.6, histtype='stepfilled', normed=True, )#, range = props_limits[i])
        txt = r'%.2f' % np.mean(xall)
        kw_text = dict(pos_x = pos_x, pos_y = 0.96, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        txt = r'%.2f' % np.median(xall)
        kw_text = dict(pos_x = pos_x, pos_y = 0.88, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        txt = r'%.2f' % np.std(xall)
        kw_text = dict(pos_x = pos_x, pos_y = 0.80, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        txt = r'%.2f' % np.max(xall)
        kw_text = dict(pos_x = pos_x, pos_y = 0.72, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        txt = r'%.2f' % np.min(xall)
        kw_text = dict(pos_x = pos_x, pos_y = 0.64, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)

        c = 'b'
        pos_x = 0.98
        va = 'top'
        if k == 'at_flux' or k == 'ba':
            pos_x = 0.14
        ax.hist(x, bins = bins, color = c, align = 'mid', alpha = 0.6, histtype='stepfilled', normed=True, )#, range = props_limits[i])
        txt = r'%.2f' % np.mean(x)
        kw_text = dict(pos_x = pos_x, pos_y = 0.96, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        txt = r'%.2f' % np.median(x)
        kw_text = dict(pos_x = pos_x, pos_y = 0.88, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        txt = r'%.2f' % np.std(x)
        kw_text = dict(pos_x = pos_x, pos_y = 0.80, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        txt = r'%.2f' % np.max(x)
        kw_text = dict(pos_x = pos_x, pos_y = 0.72, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        txt = r'%.2f' % np.min(x)
        kw_text = dict(pos_x = pos_x, pos_y = 0.64, fs = 10, va = va, ha = ha, c = c)
        plot_text_ax(ax, txt, **kw_text)
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)

        if k == 'McorSD':
            ax.set_xlim(0,5)
        elif k == 'Mcor':
            ax.set_xlim(8.5,12)
        elif k == 'ba':
            ax.set_xlim(0, 1)    
        elif k == 'xY':
            ax.set_xlim(0, 1)
        elif k == 'at_flux':
            ax.set_xlim(6, 10)
            
        if col == 0:
            ax.set_ylabel(r'$\mathrm{frac\c\~ao}$')
            
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))            
        row, col = next_row_col(row, col, NRows, NCols)
    f.subplots_adjust(bottom = 0.1, top = 0.95, hspace = 0.28, wspace = 0.28, right = 0.95, left = 0.12)
    fname = 'sample_histo_%s' % fnamesuffix        
    f.savefig(fname)
    plt.close(f)
         
    
    

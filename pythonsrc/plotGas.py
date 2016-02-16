#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
from CALIFAUtils.scripts import get_CALIFAID_by_NEDName
from matplotlib.backends.backend_pdf import PdfPages
from CALIFAUtils.scripts import mask_radius_iT
from CALIFAUtils.scripts import mask_zones_iT
from CALIFAUtils.plots import plot_histo_axis
from CALIFAUtils.plots import density_contour
from matplotlib.ticker import MaxNLocator
from CALIFAUtils.plots import plot_text_ax
from CALIFAUtils.objects import runstats
from matplotlib import pyplot as plt
from os.path import basename
import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

A4Size_inches = [ 8.267, 11.692 ]
LetterSize_inches = [ 8.5, 11 ]

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'hdf5' : None,
        'output' : None,
        'itSF' : 11,
        'itZ' :-1,
        'maskradius' : None,
        'slice_gals' : None,
        'dryrun' : False,
        'bamin' : 0,
        'scatter': False,
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--scatter',
                        action = 'store_true',
                        default = default['scatter'])
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
    
    h5file = args.hdf5
    H = C.H5SFRData(h5file)
    iT = args.itSF
    tSF = H.tSF__T[iT]
    iU = args.itZ

    fnamesuffix = '.pdf'
    if args.output is not None:
        fnamesuffix = '%s%s' % (args.output, fnamesuffix)
    minR = args.maskradius

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
    
    mask__g = mask_zones_iT(iT, H, args, maskRadiusOk__g, gals_slice__g)
    mask__rg = mask_radius_iT(iT, H, args, maskRadiusOk__rg, gals_slice__rg)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # mask__g = np.bitwise_or(mask__g, np.ma.log10(H.SFRSD_Ha__g * 1e6).mask)
    # mask__g = np.bitwise_or(mask__g, np.ma.log10(H.tau_V_neb__g).mask)
    # mask__g = np.bitwise_or(mask__g, H.logO3N2_M13__g.mask)
    # #mask__g = np.bitwise_or(mask__g, np.less(H.EW_Ha__g, 3.))
    # mask__g = np.bitwise_or(mask__g, np.less(H.reply_arr_by_zones(H.ba_GAL__g), ba_max))
    # mask__g = np.bitwise_or(mask__g, ~maskRadiusOk__g)
    # mask__g = np.bitwise_or(mask__g, ~gals_slice__g)
    # 
    # mask__rg = np.bitwise_or(~maskRadiusOk__rg, np.less(H.reply_arr_by_radius(H.ba_GAL__g), ba_max))
    # mask__rg = np.bitwise_or(mask__rg, ~gals_slice__rg)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    
    print ((~mask__g).sum()), ((~mask__rg).sum()), NgalsOkZones, NgalsOkRbins
    
    Area_GAL__g = (H.Mcor_GAL__g / H.McorSD_GAL__g)
    dustdim = 0.2 # md / rhod

    gasnames = [ 'BRTauVStar', 'BRTauVNeb', 'SKSFRSDStar', 'SKSFRSDHa' ]
    gas_dict = {}
    
    ######################
    # SK Law 
    ######################
    SK_zero = 1.6e-4
    SK_slope = 1.4
    aux = 1e6 ** (1. / SK_slope)
    SK_SigmaGas__g = aux * (H.SFRSD__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__g = aux * (H.SFRSD_Ha__g / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas__rg = aux * (H.aSFRSD__Trg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__rg = aux * (H.aSFRSD_Ha__Trg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_oneHLR__g = aux * (H.aSFRSD_oneHLR__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha_oneHLR__g = aux * (H.aSFRSD_Ha_oneHLR__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_integrated_SigmaGas = aux * (H.integrated_SFRSD__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_integrated_SigmaGas_Ha = aux * (H.integrated_SFRSD_Ha__g / SK_zero) ** (1. / SK_slope) 
    SK_DGR__g = dustdim * H.tau_V__Tg[iT] / SK_SigmaGas__g
    SK_DGR_Ha__g = dustdim * H.tau_V_neb__g / SK_SigmaGas_Ha__g
    SK_DGR__rg = dustdim * H.tau_V__Trg[iT] / SK_SigmaGas__rg
    SK_DGR_Ha__rg = dustdim * H.tau_V_neb__Trg[iT] / SK_SigmaGas_Ha__rg
    SK_DGR_oneHLR__g = dustdim * H.tau_V_oneHLR__Tg[iT] / SK_SigmaGas_oneHLR__g
    SK_DGR_Ha_oneHLR__g = dustdim * H.tau_V_neb_oneHLR__Tg[iT] / SK_SigmaGas_Ha_oneHLR__g
    SK_integrated_DGR = dustdim * H.integrated_tau_V__g / SK_integrated_SigmaGas
    SK_integrated_DGR_Ha = dustdim * H.integrated_tau_V_neb__g / SK_integrated_SigmaGas_Ha
    SK_GSR__g = SK_SigmaGas__g / H.McorSD__Tg[iT]
    SK_GSR_Ha__g = SK_SigmaGas_Ha__g / H.McorSD__Tg[iT]
    SK_GSR__rg = SK_SigmaGas__rg / H.McorSD__Trg[iT]
    SK_GSR_Ha__rg = SK_SigmaGas_Ha__rg / H.McorSD__Trg[iT]
    SK_GSR_oneHLR__g = SK_SigmaGas_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    SK_GSR_Ha_oneHLR__g = SK_SigmaGas_Ha_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    SK_integrated_GSR = SK_integrated_SigmaGas / H.McorSD_GAL__g
    SK_integrated_GSR_Ha = SK_integrated_SigmaGas_Ha / H.McorSD_GAL__g
    SK_f_gas__g = 1. / (1. + 1. / SK_GSR__g)
    SK_f_gas_Ha__g = 1. / (1. + 1. / SK_GSR_Ha__g)
    SK_f_gas__rg = 1. / (1. + 1. / SK_GSR__rg)
    SK_f_gas_Ha__rg = 1. / (1. + 1. / SK_GSR_Ha__rg)
    SK_f_gas_oneHLR__g = 1. / (1. + 1. / SK_GSR_oneHLR__g)
    SK_f_gas_Ha_oneHLR__g = 1. / (1. + 1. / SK_GSR_Ha_oneHLR__g)
    SK_integrated_f_gas = 1. / (1. + 1. / SK_integrated_GSR)
    SK_integrated_f_gas_Ha = 1. / (1. + 1. / SK_integrated_GSR_Ha)
    gas_dict['SKSFRSDStar'] = dict(
        c = 'g',
        SigmaGas__g = SK_SigmaGas__g,
        SigmaGas__rg = SK_SigmaGas__rg,
        integrated_SigmaGas = SK_integrated_SigmaGas,
        SigmaGas_oneHLR = SK_SigmaGas_oneHLR__g,
        DGR__g = SK_DGR__g,
        DGR__rg = SK_DGR__rg,
        integrated_DGR = SK_integrated_DGR,
        DGR_oneHLR = SK_DGR_oneHLR__g,
        GSR__g = SK_GSR__g,
        GSR__rg = SK_GSR__rg,
        integrated_GSR = SK_integrated_GSR,
        GSR_oneHLR = SK_GSR_oneHLR__g,
        f_gas__g = SK_f_gas__g,
        f_gas__rg = SK_f_gas__rg,
        integrated_f_gas = SK_integrated_f_gas,
        f_gas_oneHLR = SK_f_gas_oneHLR__g,
    )
    gas_dict['SKSFRSDHa'] = dict(
        c = 'k',
        SigmaGas__g = SK_SigmaGas_Ha__g,
        SigmaGas__rg = SK_SigmaGas_Ha__rg,
        integrated_SigmaGas = SK_integrated_SigmaGas_Ha,
        SigmaGas_oneHLR = SK_SigmaGas_Ha_oneHLR__g,
        DGR__g = SK_DGR_Ha__g,
        DGR__rg = SK_DGR_Ha__rg,
        integrated_DGR = SK_integrated_DGR_Ha,
        DGR_oneHLR = SK_DGR_Ha_oneHLR__g,
        GSR__g = SK_GSR_Ha__g,
        GSR__rg = SK_GSR_Ha__rg,
        integrated_GSR = SK_integrated_GSR_Ha,
        GSR_oneHLR = SK_GSR_Ha_oneHLR__g,
        f_gas__g = SK_f_gas_Ha__g,
        f_gas__rg = SK_f_gas_Ha__rg,
        integrated_f_gas = SK_integrated_f_gas_Ha,
        f_gas_oneHLR = SK_f_gas_Ha_oneHLR__g,
    ) 
    ######################
    
    ######################
    #Brinchmann
    ######################
    DGR_conv_lim_sup = 1.1e-2
    DGR_conv_lim_inf = 5.3e-3
    DGR_interval = np.array([DGR_conv_lim_inf, DGR_conv_lim_sup])
    DGR_cte = DGR_interval.mean()
    OHSunBrinch_inv = 1 / (10.**(8.82 - 12))
    ######################
    # from OLS Bisector
    import h5py
    try:
        h5 = h5py.File('/Users/lacerda/dev/latex/qualify/pythonsrc/OHcalibMPAM13.h5', 'r')
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # x = 10**(h5['logOH_M13'].value - 8.69)
        # y = 10**(h5['logOH_MPA'].value - 8.89)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        x = (h5['logOH_M13'].value - 8.69)
        y = (h5['logOH_MPA'].value - 8.89)
        rs = C.runstats(x, y, **dict(smooth = True, sigma = 1.2, frac = 0.005, gs_prc = True, OLS = True))
        #p = np.polyfit(rs.xS, rs.yS, 3)
        rs.OLS_bisector()
        p = [ rs.OLS_median_slope, rs.OLS_median_intercept ]
        C.debug_var(args.debug, third_order_coeff = p)
    except IOError:
        #p = [-3.8954267, 8.32488463, -2.57174644, 0.35947907]
        p = [ 2.01634637207, 0.425233054925 ]
    aux = (H.logO3N2_M13__g - 8.69).compressed()
    aux_mask = H.logO3N2_M13__g.mask
    BR_DGR_up__g = np.ma.masked_all(aux_mask.shape)
    BR_DGR_down__g = np.ma.copy(BR_DGR_up__g)
    BR_DGR__g = np.ma.copy(BR_DGR_up__g)
    BR_DGR_up__g[~aux_mask] = DGR_conv_lim_sup * 10**np.polyval(p, aux)
    BR_DGR_down__g[~aux_mask] = DGR_conv_lim_inf * 10**np.polyval(p, aux)
    BR_DGR__g[~aux_mask] = DGR_cte * 10**np.polyval(p, aux)
    aux = (H.logO3N2_M13__Trg[iT] - 8.69).compressed()
    aux_mask = H.logO3N2_M13__Trg[iT].mask
    BR_DGR_up__rg = np.ma.masked_all(aux_mask.shape)
    BR_DGR_down__rg = np.ma.copy(BR_DGR_up__rg)
    BR_DGR__rg = np.ma.copy(BR_DGR_up__rg)
    BR_DGR_up__rg[~aux_mask] = DGR_conv_lim_sup * 10**np.polyval(p, aux)
    BR_DGR_down__rg[~aux_mask] = DGR_conv_lim_inf * 10**np.polyval(p, aux)
    BR_DGR__rg[~aux_mask] = DGR_cte * 10**np.polyval(p, aux)
    aux = (H.logO3N2_M13_oneHLR__Tg[iT] - 8.69).compressed()
    aux_mask = H.logO3N2_M13_oneHLR__Tg[iT].mask
    BR_DGR_up_oneHLR__g = np.ma.masked_all(aux_mask.shape)
    BR_DGR_down_oneHLR__g = np.ma.copy(BR_DGR_up_oneHLR__g)
    BR_DGR_oneHLR__g = np.ma.copy(BR_DGR_up_oneHLR__g)
    BR_DGR_up_oneHLR__g[~aux_mask] = DGR_conv_lim_sup * 10**np.polyval(p, aux)
    BR_DGR_down_oneHLR__g[~aux_mask] = DGR_conv_lim_inf * 10**np.polyval(p, aux)
    BR_DGR_oneHLR__g[~aux_mask] = DGR_cte * 10**np.polyval(p, aux)
    aux = (H.integrated_logO3N2_M13__g[iT] - 8.69)
    aux_mask = np.ones_like(aux, dtype = np.bool_)
    BR_integrated_DGR_up = np.ma.masked_all(aux_mask.shape)
    BR_integrated_DGR_down = np.ma.copy(BR_integrated_DGR_up)
    BR_integrated_DGR = np.ma.copy(BR_integrated_DGR_up)
    BR_integrated_DGR_up[~aux_mask] = DGR_conv_lim_sup * 10**np.polyval(p, aux)
    BR_integrated_DGR_down[~aux_mask] = DGR_conv_lim_inf * 10**np.polyval(p, aux)
    BR_integrated_DGR[~aux_mask] = DGR_cte * 10**np.polyval(p, aux)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # p = np.array([1.87, -6.98])
    # BR_OHBrinch__g = np.ma.masked_all((H.logO3N2_M13__g.shape))
    # BR_OHBrinch__g[~H.logO3N2_M13__g.mask] = np.polyval(p, H.logO3N2_M13__g.compressed()) 
    # BR_OHBrinch__rg = np.ma.masked_all((H.logO3N2_M13__Trg[iT].shape))
    # BR_OHBrinch__rg[~H.logO3N2_M13__Trg[iT].mask] = np.polyval(p, H.logO3N2_M13__Trg[iT].compressed())
    # BR_OHBrinch_oneHLR__g = np.ma.masked_all((H.logO3N2_M13_oneHLR__Tg[iT].shape))
    # BR_OHBrinch_oneHLR__g[~H.logO3N2_M13_oneHLR__Tg[iT].mask] = np.polyval(p, H.logO3N2_M13_oneHLR__Tg[iT].compressed())
    # BR_integrated_OHBrinch = np.ma.masked_all((H.integrated_logO3N2_M13__g.shape))
    # BR_integrated_OHBrinch[~H.integrated_logO3N2_M13__g.mask] = np.polyval(p, H.integrated_logO3N2_M13__g.compressed())
    # BR_DGR_up__g = DGR_conv_lim_sup * (10 ** (BR_OHBrinch__g - 12) * OHSunBrinch_inv)
    # BR_DGR_up__rg = DGR_conv_lim_sup * (10 ** (BR_OHBrinch__rg - 12) * OHSunBrinch_inv)
    # BR_DGR_down__g = DGR_conv_lim_inf * (10 ** (BR_OHBrinch__g - 12) * OHSunBrinch_inv)
    # BR_DGR_down__rg = DGR_conv_lim_inf * (10 ** (BR_OHBrinch__rg - 12) * OHSunBrinch_inv)
    # BR_DGR__g = DGR_cte * (10 ** (BR_OHBrinch__g - 12) * OHSunBrinch_inv)
    # BR_DGR__rg = DGR_cte * (10 ** (BR_OHBrinch__rg - 12) * OHSunBrinch_inv)
    # BR_DGR_oneHLR__g = DGR_cte * (10 ** (BR_OHBrinch_oneHLR__g - 12) * OHSunBrinch_inv)
    # BR_DGR_down_oneHLR__g = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_oneHLR__g - 12) * OHSunBrinch_inv)
    # BR_DGR_up_oneHLR__g = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_oneHLR__g - 12) * OHSunBrinch_inv)
    # BR_integrated_DGR = DGR_cte * (10 ** (BR_integrated_OHBrinch - 12) * OHSunBrinch_inv)
    # BR_integrated_DGR_up = DGR_conv_lim_sup * (10 ** (BR_integrated_OHBrinch - 12) * OHSunBrinch_inv)
    # BR_integrated_DGR_down = DGR_conv_lim_inf * (10 ** (BR_integrated_OHBrinch - 12) * OHSunBrinch_inv)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    BR_SigmaGas_up__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_up__g
    BR_SigmaGas_up__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_up__rg
    BR_SigmaGas_Ha_up__g = dustdim * H.tau_V_neb__g / BR_DGR_up__g
    BR_SigmaGas_Ha_up__rg = dustdim * H.tau_V_neb__Trg[iT] / BR_DGR_up__rg
    BR_SigmaGas_down__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_down__g
    BR_SigmaGas_down__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_down__rg
    BR_SigmaGas_Ha_down__g = dustdim * H.tau_V_neb__g / BR_DGR_down__g
    BR_SigmaGas_Ha_down__rg = dustdim * H.tau_V_neb__Trg[iT] / BR_DGR_down__rg
    BR_SigmaGas__g = dustdim * H.tau_V__Tg[iT] / BR_DGR__g
    BR_SigmaGas__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR__rg
    BR_SigmaGas_Ha__g = dustdim * H.tau_V_neb__g / BR_DGR__g
    BR_SigmaGas_Ha__rg = dustdim * H.tau_V_neb__Trg[iT] / BR_DGR__rg
    BR_SigmaGas_oneHLR__g = dustdim * H.tau_V_oneHLR__Tg[iT] / BR_DGR_oneHLR__g
    BR_SigmaGas_up_oneHLR__g = dustdim * H.tau_V_oneHLR__Tg[iT] / BR_DGR_up_oneHLR__g
    BR_SigmaGas_down_oneHLR__g = dustdim * H.tau_V_oneHLR__Tg[iT] / BR_DGR_down_oneHLR__g
    BR_SigmaGas_Ha_oneHLR__g = dustdim * H.tau_V_neb_oneHLR__Tg[iT] / BR_DGR_oneHLR__g
    BR_SigmaGas_Ha_up_oneHLR__g = dustdim * H.tau_V_neb_oneHLR__Tg[iT] / BR_DGR_up_oneHLR__g
    BR_SigmaGas_Ha_down_oneHLR__g = dustdim * H.tau_V_neb_oneHLR__Tg[iT] / BR_DGR_down_oneHLR__g
    BR_integrated_SigmaGas = dustdim * H.integrated_tau_V__g / BR_integrated_DGR
    BR_integrated_SigmaGas_up = dustdim * H.integrated_tau_V__g / BR_integrated_DGR_up
    BR_integrated_SigmaGas_down = dustdim * H.integrated_tau_V__g / BR_integrated_DGR_down
    BR_integrated_SigmaGas_Ha = dustdim * H.integrated_tau_V_neb__g / BR_integrated_DGR
    BR_integrated_SigmaGas_Ha_up = dustdim * H.integrated_tau_V_neb__g / BR_integrated_DGR_up
    BR_integrated_SigmaGas_Ha_down = dustdim * H.integrated_tau_V_neb__g / BR_integrated_DGR_down
    BR_GSR_up__g = BR_SigmaGas_up__g / H.McorSD__Tg[iT]
    BR_GSR_up__rg = BR_SigmaGas_up__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_up__g = BR_SigmaGas_Ha_up__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_up__rg = BR_SigmaGas_Ha_up__rg / H.McorSD__Trg[iT]
    BR_GSR_down__g = BR_SigmaGas_down__g / H.McorSD__Tg[iT]
    BR_GSR_down__rg = BR_SigmaGas_down__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_down__g = BR_SigmaGas_Ha_down__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_down__rg = BR_SigmaGas_Ha_down__rg / H.McorSD__Trg[iT]
    BR_GSR__g = BR_SigmaGas__g / H.McorSD__Tg[iT]
    BR_GSR__rg = BR_SigmaGas__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha__g = BR_SigmaGas_Ha__g / H.McorSD__Tg[iT]
    BR_GSR_Ha__rg = BR_SigmaGas_Ha__rg / H.McorSD__Trg[iT]
    BR_GSR_oneHLR__g = BR_SigmaGas_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_up_oneHLR__g = BR_SigmaGas_up_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_down_oneHLR__g = BR_SigmaGas_down_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_Ha_oneHLR__g = BR_SigmaGas_Ha_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_Ha_up_oneHLR__g = BR_SigmaGas_Ha_up_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_Ha_down_oneHLR__g = BR_SigmaGas_Ha_down_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_integrated_GSR = BR_integrated_SigmaGas / H.McorSD_GAL__g
    BR_integrated_GSR_up = BR_integrated_SigmaGas_up / H.McorSD_GAL__g
    BR_integrated_GSR_down = BR_integrated_SigmaGas_down / H.McorSD_GAL__g
    BR_integrated_GSR_Ha = BR_integrated_SigmaGas_Ha / H.McorSD_GAL__g
    BR_integrated_GSR_Ha_up = BR_integrated_SigmaGas_Ha_up / H.McorSD_GAL__g
    BR_integrated_GSR_Ha_down = BR_integrated_SigmaGas_Ha_down / H.McorSD_GAL__g
    BR_f_gas_up__g = 1. / (1. + 1. / BR_GSR_up__g)
    BR_f_gas_up__rg = 1. / (1. + 1. / BR_GSR_up__rg)
    BR_f_gas_Ha_up__g = 1. / (1. + 1. / BR_GSR_Ha_up__g)
    BR_f_gas_Ha_up__rg = 1. / (1. + 1. / BR_GSR_Ha_up__rg)
    BR_f_gas_down__g = 1. / (1. + 1. / BR_GSR_down__g)
    BR_f_gas_down__rg = 1. / (1. + 1. / BR_GSR_down__rg)
    BR_f_gas_Ha_down__g = 1. / (1. + 1. / BR_GSR_Ha_down__g)
    BR_f_gas_Ha_down__rg = 1. / (1. + 1. / BR_GSR_Ha_down__rg)
    BR_f_gas__g = 1. / (1. + 1. / BR_GSR__g)
    BR_f_gas__rg = 1. / (1. + 1. / BR_GSR__rg)
    BR_f_gas_Ha__g = 1. / (1. + 1. / BR_GSR_Ha__g)
    BR_f_gas_Ha__rg = 1. / (1. + 1. / BR_GSR_Ha__rg)
    BR_f_gas_oneHLR__g = 1. / (1. + 1. / BR_GSR_oneHLR__g)
    BR_f_gas_up_oneHLR__g = 1. / (1. + 1. / BR_GSR_up_oneHLR__g)
    BR_f_gas_down_oneHLR__g = 1. / (1. + 1. / BR_GSR_down_oneHLR__g)
    BR_f_gas_Ha_oneHLR__g = 1. / (1. + 1. / BR_GSR_Ha_oneHLR__g)
    BR_f_gas_Ha_up_oneHLR__g = 1. / (1. + 1. / BR_GSR_Ha_up_oneHLR__g)
    BR_f_gas_Ha_down_oneHLR__g = 1. / (1. + 1. / BR_GSR_Ha_down_oneHLR__g)
    BR_integrated_f_gas = 1. / (1. + 1. / BR_integrated_GSR)
    BR_integrated_f_gas_up = 1. / (1. + 1. / BR_integrated_GSR_up)
    BR_integrated_f_gas_down = 1. / (1. + 1. / BR_integrated_GSR_down)
    BR_integrated_f_gas_Ha = 1. / (1. + 1. / BR_integrated_GSR_Ha)
    BR_integrated_f_gas_Ha_up = 1. / (1. + 1. / BR_integrated_GSR_Ha_up)
    BR_integrated_f_gas_Ha_down = 1. / (1. + 1. / BR_integrated_GSR_Ha_down)
    gas_dict['BRTauVStar'] = dict(
        c = 'b',
        SigmaGas__g = BR_SigmaGas__g,
        SigmaGas__rg = BR_SigmaGas__rg,
        integrated_SigmaGas = BR_integrated_SigmaGas,
        SigmaGas_oneHLR = BR_SigmaGas_oneHLR__g,
        DGR__g = BR_DGR__g,
        DGR__rg = BR_DGR__rg,
        integrated_DGR = BR_integrated_DGR,
        DGR_oneHLR = BR_DGR_oneHLR__g,
        GSR__g = BR_GSR__g,
        GSR__rg = BR_GSR__rg,
        integrated_GSR = BR_integrated_GSR,
        GSR_oneHLR = BR_GSR_oneHLR__g,
        f_gas__g = BR_f_gas__g,
        f_gas__rg = BR_f_gas__rg,
        integrated_f_gas = BR_integrated_f_gas,
        f_gas_oneHLR = BR_f_gas_oneHLR__g,
        up = dict(
            SigmaGas__g = BR_SigmaGas_up__g,
            SigmaGas__rg = BR_SigmaGas_up__rg,
            integrated_SigmaGas = BR_integrated_SigmaGas_up,
            SigmaGas_oneHLR = BR_SigmaGas_up_oneHLR__g,
            DGR__g = BR_DGR_up__g,
            DGR__rg = BR_DGR_up__rg,
            integrated_DGR = BR_integrated_DGR_up,
            DGR_oneHLR = BR_DGR_up_oneHLR__g,
            GSR__g = BR_GSR_up__g,
            GSR__rg = BR_GSR_up__rg,
            integrated_GSR = BR_integrated_GSR_up,
            GSR_oneHLR = BR_GSR_up_oneHLR__g,
            f_gas__g = BR_f_gas_up__g,
            f_gas__rg = BR_f_gas_up__rg,
            integrated_f_gas = BR_integrated_f_gas_up,
            f_gas_oneHLR = BR_f_gas_up_oneHLR__g,
        ),
        down = dict(
            SigmaGas__g = BR_SigmaGas_down__g,
            SigmaGas__rg = BR_SigmaGas_down__rg,
            integrated_SigmaGas = BR_integrated_SigmaGas_down,
            SigmaGas_oneHLR = BR_SigmaGas_down_oneHLR__g,
            DGR__g = BR_DGR_down__g,
            DGR__rg = BR_DGR_down__rg,
            integrated_DGR = BR_integrated_DGR_down,
            DGR_oneHLR = BR_DGR_down_oneHLR__g,
            GSR__g = BR_GSR_down__g,
            GSR__rg = BR_GSR_down__rg,
            integrated_GSR = BR_integrated_GSR_down,
            GSR_oneHLR = BR_GSR_down_oneHLR__g,
            f_gas__g = BR_f_gas_down__g,
            f_gas__rg = BR_f_gas_down__rg,
            integrated_f_gas = BR_integrated_f_gas_down,
            f_gas_oneHLR = BR_f_gas_down_oneHLR__g,
        ),
    )
    gas_dict['BRTauVNeb'] = dict(
        c = 'r',
        SigmaGas__g = BR_SigmaGas_Ha__g,
        SigmaGas__rg = BR_SigmaGas_Ha__rg,
        integrated_SigmaGas = BR_integrated_SigmaGas_Ha,
        SigmaGas_oneHLR = BR_SigmaGas_Ha_oneHLR__g,
        DGR__g = BR_DGR__g,
        DGR__rg = BR_DGR__rg,
        integrated_DGR = BR_integrated_DGR,
        DGR_oneHLR = BR_DGR_oneHLR__g,
        GSR__g = BR_GSR_Ha__g,
        GSR__rg = BR_GSR_Ha__rg,
        integrated_GSR = BR_integrated_GSR_Ha,
        GSR_oneHLR = BR_GSR_Ha_oneHLR__g,
        f_gas__g = BR_f_gas_Ha__g,
        f_gas__rg = BR_f_gas_Ha__rg,
        integrated_f_gas = BR_integrated_f_gas_Ha,
        f_gas_oneHLR = BR_f_gas_Ha_oneHLR__g,
        up = dict(
            SigmaGas__g = BR_SigmaGas_Ha_up__g,
            SigmaGas__rg = BR_SigmaGas_Ha_up__rg,
            integrated_SigmaGas = BR_integrated_SigmaGas_Ha_up,
            SigmaGas_oneHLR = BR_SigmaGas_Ha_up_oneHLR__g,
            DGR__g = BR_DGR_up__g,
            DGR__rg = BR_DGR_up__rg,
            integrated_DGR = BR_integrated_DGR_up,
            DGR_oneHLR = BR_DGR_up_oneHLR__g,
            GSR__g = BR_GSR_Ha_up__g,
            GSR__rg = BR_GSR_Ha_up__rg,
            integrated_GSR = BR_integrated_GSR_Ha_up,
            GSR_oneHLR = BR_GSR_Ha_up_oneHLR__g,
            f_gas__g = BR_f_gas_Ha_up__g,
            f_gas__rg = BR_f_gas_Ha_up__rg,
            integrated_f_gas = BR_integrated_f_gas_Ha_up,
            f_gas_oneHLR = BR_f_gas_Ha_up_oneHLR__g,
        ),
        down = dict(
            SigmaGas__g = BR_SigmaGas_Ha_down__g,
            SigmaGas__rg = BR_SigmaGas_Ha_down__rg,
            integrated_SigmaGas = BR_integrated_SigmaGas_Ha_down,
            SigmaGas_oneHLR = BR_SigmaGas_Ha_down_oneHLR__g,
            DGR__g = BR_DGR_down__g,
            DGR__rg = BR_DGR_down__rg,
            integrated_DGR = BR_integrated_DGR_down,
            DGR_oneHLR = BR_DGR_down_oneHLR__g,
            GSR__g = BR_GSR_Ha_down__g,
            GSR__rg = BR_GSR_Ha_down__rg,
            integrated_GSR = BR_integrated_GSR_Ha_down,
            GSR_oneHLR = BR_GSR_Ha_down_oneHLR__g,
            f_gas__g = BR_f_gas_Ha_down__g,
            f_gas__rg = BR_f_gas_Ha_down__rg,
            integrated_f_gas = BR_integrated_f_gas_Ha_down,
            f_gas_oneHLR = BR_f_gas_Ha_down_oneHLR__g,
        ),
    )

    ######################
    # from cubic polynomial fit
    #p_cubic = np.array([-4.91783872, 122.48149162, -1014.51941088, 2803.24285985])
    ######################
    
    ####################################
    #### Mass from CALIFA Datafiles ####
    ####################################
    ####################################
    '''
    FILE: M_H1_CALIFA.csv - Miguel A. Perez
    delimiter = ','
    comment = '#'
    columns:
        1 - CALIFA No
        2 - NED Name
        3 - Distance
        4 - RA(J2000.0)
        5 - DEC(J2000.0)
        6 - Sobs
        7 - Scor
        8 - Sabs
        9 - e_Sabs,
        10 - M(HI)abs
        11 - e_M(HI)abs
        
    FILE: CALIFA_HI_angel.dat - Angel R. Lopez-Sanchez
    delimiter = ','
    comments = ';'
    columns:
        1 - NED Name
        2 - redshift
        3 - log(Ms)
        4 - e_log(Ms)
        5 - 12+log(O/H)
        6 - e_12+log(O/H)
        7 - log(SFR)
        8 - e_log(SFR)
        9 - log(Mg = 1.32 MHI)
        10 - e_log(Mg = 1.32 MHI)
    '''
    dirs = C.CALIFAPaths()
    file_miguel = '%sM_HI_CALIFA.csv' % dirs.califa_work_dir
    dtype_miguel = np.dtype([('califaID', '|S5'), ('M_HI', np.float)])
    read_miguel = np.loadtxt(file_miguel, delimiter = ',', usecols = [0, 9], dtype = dtype_miguel)
    map_miguel = {}
    for i, g in enumerate(read_miguel['califaID']):
        map_miguel[g] = i
    aux = set(H.califaIDs.tolist())
    gals_miguel_intersect = sorted([g for g in map_miguel.keys() if g in aux])
    gals_miguel_slice = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
    integrated_M_HI_miguel__g = np.ma.masked_all(H.califaIDs_all.shape, dtype = np.float_)
    for g in gals_miguel_intersect:
        i = H.califaIDs_all.tolist().index(g)
        i_r = map_miguel[g]
        #print g, i, i_r
        gals_miguel_slice[i] = True
        integrated_M_HI_miguel__g[i] = read_miguel['M_HI'][i_r]
    integrated_f_gas_miguel = 1. / (1. + 1 / ((integrated_M_HI_miguel__g) / H.Mcor_GAL__g))
    #integrated_f_gas_miguel = 1./(1. + 1/((1.32 * integrated_M_HI_miguel__g) / H.Mcor_GAL__g))
    #integrated_f_gas_angel = 1./(1. + 1/((10. ** integrated_log_M_gas_angel__g) / H.Mcor_GAL__g))

    default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
    default_rs_kwargs = dict(smooth = True, sigma = 1.2, frac = 0.05, gs_prc = True)
    
    if args.dryrun is True:
        sys.exit('dryrun')
    
    if args.output is None:
        if args.maskradius is not None:
            output = 'gas_maskRadius%.1f.pdf' % args.maskradius
        else:
            output = 'gas.pdf'
    else:
        output = args.output
          
#    with PdfPages(output) as pdf:
    ##########################
    ######### PAGE 1 #########
    ##########################
    NRows = 1
    NCols = 3
    f = plt.figure()
    #page_size_inches = A4Size_inches[::-1]
    page_size_inches = [NCols * 4, NRows * 4]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
     
    ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    rs_kwargs = default_rs_kwargs.copy()
    sc_kwargs = default_sc_kwargs.copy()
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_DGR__rg, mask = mask__rg)
    rs_BR_DGR = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    m_aux = mask__rg | BR_DGR_up__rg.mask | BR_DGR_up__rg.mask
    xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_DGR_up__rg, mask = m_aux)
    rs_BR_DGR_up = runstats(xm.compressed(), yup.compressed(), **rs_kwargs)
    xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_DGR_down__rg, mask = m_aux)
    rs_BR_DGR_down = runstats(xm.compressed(), ydown.compressed(), **rs_kwargs)
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = SK_DGR__rg, mask = mask__rg) 
    rs_SK_DGR = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = SK_DGR_Ha__rg, mask = mask__rg) 
    rs_SK_DGR_Ha = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    if args.scatter is True:
        ax.scatter(rs_BR_DGR.x, rs_BR_DGR.y, c = 'b', **sc_kwargs)
        ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, c = 'g', **sc_kwargs)
        ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, c = 'g', **sc_kwargs)
    ax.plot(rs_BR_DGR.xS, rs_BR_DGR.yS, '.-', c = 'b')
    ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'g')
    ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'black')
    ax.fill_between(rs_BR_DGR_up.xS, rs_BR_DGR_up.yS, rs_BR_DGR_down.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
    ax.set_xlim(minR, H.RbinFin)
    #ax.set_ylim(0, .5)
    ax.set_xlabel(r'R [HLR]')
    ax.set_ylabel(r'$\delta_{DGR}$')
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    #ax.grid()
    #plt.setp(ax.get_xticklabels(), visible = False)
      
    #ax = plt.subplot2grid(grid_shape, loc = (1, 0))
    ax = plt.subplot2grid(grid_shape, loc = (0, 1))
    ax.set_axis_on()
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas__rg), mask = mask__rg) 
    rs_BR_SigmaGas = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_up__rg), mask = mask__rg) 
    rs_BR_SigmaGas_up = runstats(xm.compressed(), yup.compressed(), **rs_kwargs)
    xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_down__rg), mask = yup.mask) 
    rs_BR_SigmaGas_down = runstats(xm.compressed(), ydown.compressed(), **rs_kwargs)
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_Ha__rg), mask = mask__rg)
    rs_BR_SigmaGas_Ha = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_Ha_up__rg), mask = mask__rg) 
    rs_BR_SigmaGas_Ha_up = runstats(xm.compressed(), yup.compressed(), **rs_kwargs)
    xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_Ha_down__rg), mask = yup.mask) 
    rs_BR_SigmaGas_Ha_down = runstats(xm.compressed(), ydown.compressed(), **rs_kwargs)
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(SK_SigmaGas__rg), mask = mask__rg)
    rs_SK_SigmaGas = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(SK_SigmaGas_Ha__rg), mask = mask__rg)
    rs_SK_SigmaGas_Ha = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    if args.scatter is True:
        ax.scatter(rs_BR_SigmaGas.x, rs_BR_SigmaGas.y, c = 'b', **sc_kwargs)
        ax.scatter(rs_BR_SigmaGas_Ha.x, rs_BR_SigmaGas_Ha.y, c = 'r', **sc_kwargs)
        ax.scatter(rs_SK_SigmaGas.x, rs_SK_SigmaGas.y, c = 'g', **sc_kwargs)
        ax.scatter(rs_SK_SigmaGas_Ha.x, rs_SK_SigmaGas_Ha.y, c = 'black', **sc_kwargs)
    ax.plot(rs_BR_SigmaGas.xS, rs_BR_SigmaGas.yS, '.-', c = 'b')
    ax.plot(rs_BR_SigmaGas_Ha.xS, rs_BR_SigmaGas_Ha.yS, '.-', c = 'r')
    ax.plot(rs_SK_SigmaGas.xS, rs_SK_SigmaGas.yS, '.-', c = 'g')
    ax.plot(rs_SK_SigmaGas_Ha.xS, rs_SK_SigmaGas_Ha.yS, '.-', c = 'black')
    ax.fill_between(rs_BR_SigmaGas_up.xS, rs_BR_SigmaGas_up.yS, rs_BR_SigmaGas_down.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
    ax.fill_between(rs_BR_SigmaGas_Ha_up.xS, rs_BR_SigmaGas_Ha_up.yS, rs_BR_SigmaGas_Ha_down.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
    ax.set_xlim(minR, H.RbinFin)
    ax.set_ylim(0.4, 2)
    #ax.legend(loc = 'upper right')
    ax.yaxis.set_major_locator(MaxNLocator(4))
    #ax.grid()
    ax.set_xlabel(r'R [HLR]')
    #plt.setp(ax.get_xticklabels(), visible = False)
    ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [M${}_\odot$ pc${}^{-2}$]')
     
    #ax = plt.subplot2grid(grid_shape, loc = (2, 0))
    ax = plt.subplot2grid(grid_shape, loc = (0, 2))
    ax.set_axis_on()
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas__rg, mask = mask__rg) 
    rs_BR_f_gas = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_up__rg, mask = mask__rg) 
    rs_BR_f_gas_up = runstats(xm.compressed(), yup.compressed(), **rs_kwargs)
    xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_down__rg, mask = yup.mask) 
    rs_BR_f_gas_down = runstats(xm.compressed(), ydown.compressed(), **rs_kwargs)
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_Ha__rg, mask = mask__rg)
    rs_BR_f_gas_Ha = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_Ha_up__rg, mask = mask__rg) 
    rs_BR_f_gas_Ha_up = runstats(xm.compressed(), yup.compressed(), **rs_kwargs)
    xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_Ha_down__rg, mask = yup.mask) 
    rs_BR_f_gas_Ha_down = runstats(xm.compressed(), ydown.compressed(), **rs_kwargs)
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = SK_f_gas__rg, mask = mask__rg)
    rs_SK_f_gas = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = SK_f_gas_Ha__rg, mask = mask__rg)
    rs_SK_f_gas_Ha = runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    if args.scatter is True:
        ax.scatter(rs_BR_f_gas.x, rs_BR_f_gas.y, c = 'b', **sc_kwargs)
        ax.scatter(rs_BR_f_gas_Ha.x, rs_BR_f_gas_Ha.y, c = 'r', **sc_kwargs)
        ax.scatter(rs_SK_f_gas.x, rs_SK_f_gas.y, c = 'g', **sc_kwargs)
        ax.scatter(rs_SK_f_gas_Ha.x, rs_SK_f_gas_Ha.y, c = 'black', **sc_kwargs)
    ax.plot(rs_BR_f_gas.xS, rs_BR_f_gas.yS, '.-', c = 'b', label = r'BR13 $\tau_V^\star$')
    ax.plot(rs_BR_f_gas_Ha.xS, rs_BR_f_gas_Ha.yS, '.-', c = 'r', label = r'BR13 $\tau_V^{neb}$')
    ax.plot(rs_SK_f_gas.xS, rs_SK_f_gas.yS, '.-', c = 'g', label = r'KS syn.')
    ax.plot(rs_SK_f_gas_Ha.xS, rs_SK_f_gas_Ha.yS, '.-', c = 'black', label = r'KS neb.')
    ax.fill_between(rs_BR_f_gas_up.xS, rs_BR_f_gas_up.yS, rs_BR_f_gas_down.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
    ax.fill_between(rs_BR_f_gas_Ha_up.xS, rs_BR_f_gas_Ha_up.yS, rs_BR_f_gas_Ha_down.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
    ax.legend(bbox_to_anchor = (0.99, 0.99), fontsize = 12, frameon = False, ncol = 2)
    ax.set_xlim(minR, H.RbinFin)
    ax.set_ylim(0, .4)
    ax.set_xlabel(r'R [HLR]')
    ax.set_ylabel(r'$f_{\mathrm{gas}}$')
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    #ax.grid()
                 
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     ax = plt.subplot2grid(grid_shape, loc = (0, NCols - 2), rowspan = NRows)
#     ax.set_axis_off()
#     txt = r'NGals:%d  tSF:%.2f Myr  $(\frac{m_d}{\sigma_d})$ = %.1f' % (N_gals, (H.tSF__T[iT] / 1e6), dustdim)
#     kw_text = dict(pos_x = 0, pos_y = 0.65, fs = 12, va = 'bottom', ha = 'left', c = 'k')
#     plot_text_ax(ax, txt, **kw_text)
#     txt = r'$x_Y$(min):%.0f%%  $\tau_V^\star $(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
#     kw_text = dict(pos_x = 0, pos_y = 0.60, fs = 12, va = 'bottom', ha = 'left', c = 'k')
#     plot_text_ax(ax, txt, **kw_text)
# 
#     txt = r'$\delta_{DGR} = (%.2f\times 10^{-3} - %.2f\times 10^{-2}) (\frac{O/H}{(O/H)_\odot})$' % (DGR_conv_lim_inf / 1e-3, DGR_conv_lim_sup / 1e-2)
#     kw_text = dict(pos_x = 0, pos_y = 0.55, fs = 12, va = 'bottom', ha = 'left', c = 'b')
#     plot_text_ax(ax, txt, **kw_text)
#     txt = r'$\delta_{DGR} = (\frac{m_d}{\sigma_d}) \times (\frac{\tau_V^\star}{\Sigma_{gas}})$'
#     kw_text = dict(pos_x = 0, pos_y = 0.50, fs = 12, va = 'bottom', ha = 'left', c = 'g')
#     plot_text_ax(ax, txt, **kw_text)
#     txt = r'$\delta_{DGR} = (\frac{m_d}{\sigma_d}) \times (\frac{\tau_V^{neb}}{\Sigma_{gas}})$'
#     kw_text = dict(pos_x = 0, pos_y = 0.44, fs = 12, va = 'bottom', ha = 'left', c = 'k')
#     plot_text_ax(ax, txt, **kw_text)
# 
#     txt = r'$\Sigma_{gas} = (\frac{m_d}{\sigma_d})\times (\frac{\tau_V^\star}{\delta_{DGR}})$'
#     kw_text = dict(pos_x = 0, pos_y = 0.34, fs = 12, va = 'bottom', ha = 'left', c = 'b')
#     plot_text_ax(ax, txt, **kw_text)
#     txt = r'$\Sigma_{gas} = (\frac{m_d}{\sigma_d})\times (\frac{\tau_V^{neb}}{\delta_{DGR}})$'
#     kw_text = dict(pos_x = 0, pos_y = 0.28, fs = 12, va = 'bottom', ha = 'left', c = 'r')
#     plot_text_ax(ax, txt, **kw_text)
#     txt = r'$\Sigma_{gas} = (\frac{\Sigma_{SFR}^\star}{1.6\times 10^{-4}})^{\frac{1}{1.4}}$' 
#     kw_text = dict(pos_x = 0, pos_y = 0.22, fs = 12, va = 'bottom', ha = 'left', c = 'k')
#     plot_text_ax(ax, txt, **kw_text)
#     txt = r'$\Sigma_{gas} = (\frac{\Sigma_{SFR}^{H\alpha}}{1.6\times 10^{-4}})^{\frac{1}{1.4}}$' 
#     kw_text = dict(pos_x = 0, pos_y = 0.16, fs = 12, va = 'bottom', ha = 'left', c = 'g')
#     plot_text_ax(ax, txt, **kw_text)
# 
#     txt = r'$f_{gas}\ =\ \frac{\Sigma_{gas}}{\mu_\star + \Sigma_{gas}}$'
#     kw_text = dict(pos_x = 0, pos_y = 0.08, fs = 15, va = 'bottom', ha = 'left', c = 'k')
#     plot_text_ax(ax, txt, **kw_text)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    f.subplots_adjust(hspace = 0.3, wspace = 0.3, left = 0.1, right = 0.95, top = 0.95, bottom = 0.2)
    f.savefig('gas_R_%.2fMyr_%s%s' % ((tSF/1e6), basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix))
    plt.close(f)

    ##########################
    ######### PAGE 2 #########
    ##########################
    NRows = 3
    NCols = 3
    f = plt.figure()
    #page_size_inches = A4Size_inches[::-1]
    page_size_inches = [NCols * 4, NRows * 4]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    yaxis = [ 'alogSFRSDHakpcR', 'alogSFRSDkpcR', 'atfluxR', 'logO3N2M13R', 'xYR', 'logMcorSDR', 'alogZmassR', 'logtauVR', 'logtauVNebR' ]
    row, col = 0, 0
    bins = (20,20)
    for yk in yaxis:
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        ax = plt.subplot2grid(grid_shape, loc = (row, col))
        _, ydict = H.get_plot_dict(iT, -1, yk)
        x = H.Rtoplot()
        y = ydict['v']
        xm, ym = C.ma_mask_xyz(x = x, y = y, mask = mask__rg)
        #rs_kwargs['sigma'] = 1.2
        rs = runstats(xm.compressed(), ym.compressed(), xbin = H.Rbin__r[H.Rbin__r >= 0.7], **rs_kwargs)
        for i, p in enumerate(rs.xPrc):
            ax.plot(rs.xPrcS[i], rs.yPrcS[i], 'k--', lw = 2.)
        #if args.scatter is True:
        density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range = [[minR, H.RbinFin], ydict['lim']], colors = [ 'b', 'y', 'r' ])
        ax.scatter(rs.x, rs.y, marker = 'o', c = 'b', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
        ax.plot(rs.xS, rs.yS, '.-', c = 'k')
        ax.set_xlim(minR, H.RbinFin)
        ax.set_ylim(ydict['lim'])
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylabel(ydict['label'])
        ax.yaxis.set_major_locator(MaxNLocator(4))
        #ax.grid()
        #plt.setp(ax.get_yticklabels(), visible = False)
        if col == (NCols - 1):
            row += 1
            col = 0
        else:
            col += 1
    f.subplots_adjust(hspace = 0.3, wspace = 0.3, left = 0.1, right = 0.95, top = 0.95)
    f.savefig('propGas_R_%.2fMyr_%s%s' % ((tSF/1e6), basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix))
    plt.close(f)

    NRows = 3
    NCols = 3
    f = plt.figure()
    #page_size_inches = A4Size_inches[::-1]
    page_size_inches = [NCols * 4, NRows * 4]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    row, col = 0, 0
    bins = (20,20)
    for yk in yaxis:
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        ax = plt.subplot2grid(grid_shape, loc = (row, col))
        _, ydict = H.get_plot_dict(iT, -1, yk)
        y = ydict['v']
        mask_aux = (mask__rg | y.mask)
        plot_histo_axis(ax, np.ma.masked_array(y, mask = mask_aux).compressed())
        ax.set_xlim(ydict['lim'])
        ax.set_xlabel(ydict['label'])
        ax.xaxis.set_major_locator(MaxNLocator(4))
        plt.setp(ax.get_yticklabels(), visible = False)
        if col == (NCols - 1):
            row += 1
            col = 0
        else:
            col += 1
    f.subplots_adjust(hspace = 0.3, wspace = 0.3, left = 0.05, right = 0.95, top = 0.95)
    f.savefig('Histo_propGas_R_%.2fMyr_%s%s' % ((tSF/1e6), basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix))
    plt.close(f)

    ##########################
    ######### PAGE 3 #########
    ##########################
    NRows = 2
    NCols = 3
    f = plt.figure()
    #f, axArr = plt.subplots(NRows, NCols)
    #f.set_size_inches((NCols * 3, NRows * 1.5))
    #page_size_inches = (NCols * 3, NRows * 2.5)
    #page_size_inches = A4Size_inches[::-1]
    page_size_inches = [NCols * 4, NRows * 4]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    xaxis = [ 'atfluxR', 'xYR', 'logMcorSDR', 'baR', 'alogZmassR', 'logO3N2M13R' ]
    #xaxis = [ 'atfluxR', 'xYR', 'baR', 'alogZmassR', 'logMcorSDR' ]
    row, col = 0, 0
    for xk in xaxis:
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        ax = plt.subplot2grid(grid_shape, loc = (row, col))
        _, xdict = H.get_plot_dict(iT, -1, xk)
        x = xdict['v']
        xm, ym_BR_DGR = C.ma_mask_xyz(x = x, y = BR_DGR__rg, mask = mask__rg)
        rs_BR_DGR = runstats(xm.compressed(), ym_BR_DGR.compressed(), **rs_kwargs)
        xm, ym_BR_DGR_up = C.ma_mask_xyz(x = x, y = BR_DGR_up__rg, mask = mask__rg)
        rs_BR_DGR_up = runstats(xm.compressed(), ym_BR_DGR_up.compressed(), **rs_kwargs)
        xm, ym_BR_DGR_down = C.ma_mask_xyz(x = x, y = BR_DGR_down__rg, mask = mask__rg)
        rs_BR_DGR_down = runstats(xm.compressed(), ym_BR_DGR_down.compressed(), **rs_kwargs)
        xm, ym_SK_DGR = C.ma_mask_xyz(x = x, y = SK_DGR__rg, mask = mask__rg)
        rs_SK_DGR = runstats(xm.compressed(), ym_SK_DGR.compressed(), **rs_kwargs)
        xm, ym_SK_DGR_Ha = C.ma_mask_xyz(x = x, y = SK_DGR_Ha__rg, mask = mask__rg)
        rs_SK_DGR_Ha = runstats(xm.compressed(), ym_SK_DGR_Ha.compressed(), **rs_kwargs)
        if args.scatter is True:
            ax.scatter(rs_BR_DGR.x, rs_BR_DGR.y, marker = 'o', c = 'b', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, marker = 'o', c = 'g', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, marker = 'o', c = 'black', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
        ax.plot(rs_BR_DGR.xS, rs_BR_DGR.yS, '.-', c = 'b')
        ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'g')
        ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'black')
        ax.fill_between(rs_BR_DGR_up.xS, rs_BR_DGR_up.yS, rs_BR_DGR_down.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
        ax.yaxis.set_major_locator(MaxNLocator(4))
        #ax.grid()
        ax.set_xlabel(xdict['label'])
        ax.set_ylim(0, 16e-3)
        ax.set_xlim(xdict['lim'])
        if col != 0:
            plt.setp(ax.get_yticklabels(), visible = False)
        if col == 0:
            ax.set_ylabel(r'$\delta_{DGR}$')
        if col == (NCols - 1):
            row += 1
            col = 0
            prune = None
        else:
            col += 1
            prune = 'upper'
        ax.xaxis.set_major_locator(MaxNLocator(5, prune = prune))
    f.subplots_adjust(hspace = 0.3, wspace = 0, left = 0.1, right = 0.95, top = 0.95)
    f.savefig('props_DGR_%.2fMyr_%s%s' % ((tSF/1e6), basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix))
    plt.close(f)

    ##########################
    ######### PAGE 4 #########
    ##########################
    NRows = 2
    NCols = 3
    f = plt.figure()
    #f, axArr = plt.subplots(NRows, NCols)
    #f.set_size_inches((NCols * 3, NRows * 1.5))
    #page_size_inches = (NCols * 3, NRows * 2.5)
    #page_size_inches = A4Size_inches[::-1]
    page_size_inches = [NCols * 4, NRows * 4]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    xaxis = [ 'atfluxR', 'xYR', 'logMcorSDR', 'baR', 'alogZmassR', 'logO3N2M13R' ]
    #xaxis = [ 'atfluxR', 'xYR', 'baR', 'alogZmassR', 'logMcorSDR', ]
    row, col = 0, 0
    for xk in xaxis:
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        ax = plt.subplot2grid(grid_shape, loc = (row, col))
        _, xdict = H.get_plot_dict(iT, -1, xk)
        x = xdict['v']
        xm, ym_BR_SigmaGas = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas__rg), mask = mask__rg)
        rs_BR_SigmaGas = runstats(xm.compressed(), ym_BR_SigmaGas.compressed(), **rs_kwargs)
        xm, ym_BR_SigmaGas_up = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_up__rg), mask = mask__rg)
        rs_BR_SigmaGas_up = runstats(xm.compressed(), ym_BR_SigmaGas_up.compressed(), **rs_kwargs)
        xm, ym_BR_SigmaGas_down = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_down__rg), mask = mask__rg)
        rs_BR_SigmaGas_down = runstats(xm.compressed(), ym_BR_SigmaGas_down.compressed(), **rs_kwargs)
        xm, ym_BR_SigmaGas_Ha = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha__rg), mask = mask__rg)
        rs_BR_SigmaGas_Ha = runstats(xm.compressed(), ym_BR_SigmaGas_Ha.compressed(), **rs_kwargs)
        xm, ym_BR_SigmaGas_Ha_up = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_up__rg), mask = mask__rg)
        rs_BR_SigmaGas_Ha_up = runstats(xm.compressed(), ym_BR_SigmaGas_Ha_up.compressed(), **rs_kwargs)
        xm, ym_BR_SigmaGas_Ha_down = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_down__rg), mask = mask__rg)
        rs_BR_SigmaGas_Ha_down = runstats(xm.compressed(), ym_BR_SigmaGas_Ha_down.compressed(), **rs_kwargs)
        xm, ym_SK_SigmaGas = C.ma_mask_xyz(x = x, y = np.log10(SK_SigmaGas__rg), mask = mask__rg)
        rs_SK_SigmaGas = runstats(xm.compressed(), ym_SK_SigmaGas.compressed(), **rs_kwargs)
        xm, ym_SK_SigmaGas_Ha = C.ma_mask_xyz(x = x, y = np.log10(SK_SigmaGas_Ha__rg), mask = mask__rg)
        rs_SK_SigmaGas_Ha = runstats(xm.compressed(), ym_SK_SigmaGas_Ha.compressed(), **rs_kwargs)
        if args.scatter is True:
            ax.scatter(rs_BR_SigmaGas.x, rs_BR_SigmaGas.y, marker = 'o', c = 'b', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.scatter(rs_BR_SigmaGas_Ha.x, rs_BR_SigmaGas_Ha.y, marker = 'o', c = 'r', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.scatter(rs_SK_SigmaGas.x, rs_SK_SigmaGas.y, marker = 'o', c = 'g', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.scatter(rs_SK_SigmaGas_Ha.x, rs_SK_SigmaGas_Ha.y, marker = 'o', c = 'black', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
        ax.plot(rs_BR_SigmaGas.xS, rs_BR_SigmaGas.yS, '.-', c = 'b')
        ax.plot(rs_BR_SigmaGas_Ha.xS, rs_BR_SigmaGas_Ha.yS, '.-', c = 'r')
        ax.plot(rs_SK_SigmaGas.xS, rs_SK_SigmaGas.yS, '.-', c = 'g')
        ax.plot(rs_SK_SigmaGas_Ha.xS, rs_SK_SigmaGas_Ha.yS, '.-', c = 'black')
        ax.fill_between(rs_BR_SigmaGas_up.xS, rs_BR_SigmaGas_up.yS, rs_BR_SigmaGas_down.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
        ax.fill_between(rs_BR_SigmaGas_Ha_up.xS, rs_BR_SigmaGas_Ha_up.yS, rs_BR_SigmaGas_Ha_down.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)            
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        #ax.grid()
        ax.set_xlabel(xdict['label'])
        ax.set_ylim(0.4, 2)
        ax.set_xlim(xdict['lim'])
        if col != 0:
            plt.setp(ax.get_yticklabels(), visible = False)
        if col == 0:
            ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [M${}_\odot$ pc${}^{-2}$]')
        if col == (NCols - 1):
            row += 1
            col = 0
            prune = None
        else:
            col += 1
            prune = 'upper'
        ax.xaxis.set_major_locator(MaxNLocator(5, prune = prune))
    f.subplots_adjust(hspace = 0.3, wspace = 0, left = 0.1, right = 0.95, top = 0.95)
    f.savefig('props_SigmaGas_%.2fMyr_%s%s' % ((tSF/1e6), basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix))
    plt.close(f)

    ##########################
    ######### PAGE 5 #########
    ##########################
    NRows = 2
    NCols = 3
    f = plt.figure()
    #f, axArr = plt.subplots(NRows, NCols)
    #f.set_size_inches((NCols * 3, NRows * 1.5))
    #page_size_inches = (NCols * 3, NRows * 2.5)
    #page_size_inches = A4Size_inches[::-1]
    page_size_inches = [NCols * 4, NRows * 4]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    xaxis = [ 'atfluxR', 'xYR', 'logMcorSDR', 'baR', 'alogZmassR', 'logO3N2M13R' ]
    #xaxis = [ 'atfluxR', 'xYR', 'baR', 'alogZmassR', 'logMcorSDR' ]
    row, col = 0, 0
    for xk in xaxis:
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        ax = plt.subplot2grid(grid_shape, loc = (row, col))
        _, xdict = H.get_plot_dict(iT, -1, xk)
        x = xdict['v']
        xm, ym_BR_f_gas = C.ma_mask_xyz(x = x, y = BR_f_gas__rg, mask = mask__rg)
        rs_BR_f_gas = runstats(xm.compressed(), ym_BR_f_gas.compressed(), **rs_kwargs)
        xm, ym_BR_f_gas_up = C.ma_mask_xyz(x = x, y = BR_f_gas_up__rg, mask = mask__rg)
        rs_BR_f_gas_up = runstats(xm.compressed(), ym_BR_f_gas_up.compressed(), **rs_kwargs)
        xm, ym_BR_f_gas_down = C.ma_mask_xyz(x = x, y = BR_f_gas_down__rg, mask = mask__rg)
        rs_BR_f_gas_down = runstats(xm.compressed(), ym_BR_f_gas_down.compressed(), **rs_kwargs)
        xm, ym_BR_f_gas_Ha = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha__rg, mask = mask__rg)
        rs_BR_f_gas_Ha = runstats(xm.compressed(), ym_BR_f_gas_Ha.compressed(), **rs_kwargs)
        xm, ym_BR_f_gas_Ha_up = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_up__rg, mask = mask__rg)
        rs_BR_f_gas_Ha_up = runstats(xm.compressed(), ym_BR_f_gas_Ha_up.compressed(), **rs_kwargs)
        xm, ym_BR_f_gas_Ha_down = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_down__rg, mask = mask__rg)
        rs_BR_f_gas_Ha_down = runstats(xm.compressed(), ym_BR_f_gas_Ha_down.compressed(), **rs_kwargs)
        xm, ym_SK_f_gas = C.ma_mask_xyz(x = x, y = SK_f_gas__rg, mask = mask__rg)
        rs_SK_f_gas = runstats(xm.compressed(), ym_SK_f_gas.compressed(), **rs_kwargs)
        xm, ym_SK_f_gas_Ha = C.ma_mask_xyz(x = x, y = SK_f_gas_Ha__rg, mask = mask__rg)
        rs_SK_f_gas_Ha = runstats(xm.compressed(), ym_SK_f_gas_Ha.compressed(), **rs_kwargs)
        if args.scatter is True:
            ax.scatter(rs_BR_f_gas.x, rs_BR_f_gas.y, marker = 'o', c = 'b', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.scatter(rs_BR_f_gas_Ha.x, rs_BR_f_gas_Ha.y, marker = 'o', c = 'r', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.scatter(rs_SK_f_gas.x, rs_SK_f_gas.y, marker = 'o', c = 'g', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.scatter(rs_SK_f_gas_Ha.x, rs_SK_f_gas_Ha.y, marker = 'o', c = 'black', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
        ax.plot(rs_BR_f_gas.xS, rs_BR_f_gas.yS, '.-', c = 'b')
        ax.plot(rs_BR_f_gas_Ha.xS, rs_BR_f_gas_Ha.yS, '.-', c = 'r')
        ax.plot(rs_SK_f_gas.xS, rs_SK_f_gas.yS, '.-', c = 'g')
        ax.plot(rs_SK_f_gas_Ha.xS, rs_SK_f_gas_Ha.yS, '.-', c = 'black')
        ax.fill_between(rs_BR_f_gas_up.xS, rs_BR_f_gas_up.yS, rs_BR_f_gas_down.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
        ax.fill_between(rs_BR_f_gas_Ha_up.xS, rs_BR_f_gas_Ha_up.yS, rs_BR_f_gas_Ha_down.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)            
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        #ax.grid()
        ax.set_xlabel(xdict['label'])
        ax.set_ylim(0, 0.4)
        ax.set_xlim(xdict['lim'])
        if col != 0:
            plt.setp(ax.get_yticklabels(), visible = False)
        if col == 0:
            ax.set_ylabel(r'$f_{\mathrm{gas}}$')
        if col == (NCols - 1):
            row += 1
            col = 0
            prune = None
        else:
            col += 1
            prune = 'upper'
        ax.xaxis.set_major_locator(MaxNLocator(5, prune = prune))
    f.subplots_adjust(hspace = 0.3, wspace = 0, left = 0.1, right = 0.95, top = 0.95)
    f.savefig('props_fGas_%.2fMyr_%s%s' % ((tSF/1e6), basename(h5file).replace('SFR_', '').replace('.h5', ''), fnamesuffix))
    plt.close(f)


#!/usr/bin/env python
#



import os, sys, re, json, numpy, astropy, matplotlib, subprocess
#matplotlib.use('Qt5Agg')

from astropy.table import Table
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
import numpy as np
from pprint import pprint

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)

sys.path.append('/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code')

from setup_matplotlib import setup_matplotlib; setup_matplotlib()

from calc_galaxy_luminosity_function import (calc_IR_LF_Gruppioni2013)

ln = np.log






# 
# Set output name
# 
output_name = 'Plot_IR_LF_bin_by_z'



# 
# User Setting
# 
#obs_area = 7200*u.arcmin*u.arcmin # 1.4*1.4*u.deg*u.deg
#obs_area = 1.4*1.4*u.deg*u.deg
obs_area = 1.5546582999901375*u.deg*u.deg
print('obs_area = %s [%s]'%(obs_area.to(u.arcmin*u.arcmin).value, obs_area.to(u.arcmin*u.arcmin).unit))
print('obs_area = %s [%s]'%(obs_area.to(u.steradian).value, obs_area.to(u.steradian).unit))
#print('obs_area = %s [%s]'%(7200 * 3600 / 4.25451703e10, 'steradian')) # checked consistent



# 
# Read data points
# 
tb = Table.read('datatable_generated_galaxies_with_coordinates.fits')
#print(tb.colnames)
#print(tb['MSTAR'].data.shape)
#print(tb['MSTAR'][0][0], tb['SFR'][0][0])
data_lgMstar = tb['lgMstar'].data.flatten()
data_Mstar = 10**data_lgMstar
data_lgSFR = tb['lgSFR'].data.flatten()
data_SFR = 10**data_lgSFR
data_lgLIR = data_lgSFR + 10.0 # Chabrier IMF, LIR = SFR_IR * 1e10
data_redshift = tb['z'].data.flatten()
#sys.exit()



# 
# apply IRX=IR/UV
# 
do_apply_IRX_Schreiber2017 = False
do_apply_IRX_Whitaker2017 = True
# read from user command line input
for i in range(len(sys.argv)):
    if sys.argv[i] == '-IRX-Schreiber2017':
        do_apply_IRX_Schreiber2017 = True
        do_apply_IRX_Whitaker2017 = False
    if sys.argv[i] == '-IRX-Whitaker2017':
        do_apply_IRX_Schreiber2017 = False
        do_apply_IRX_Whitaker2017 = True
# 
if np.count_nonzero([do_apply_IRX_Schreiber2017, do_apply_IRX_Whitaker2017]) > 1:
    print('Error! There should have only one True value for [do_apply_IRX_Schreiber2017, do_apply_IRX_Whitaker2017]!')
    sys.exit()
# 
if do_apply_IRX_Schreiber2017:
    # Schreiber2017 Eq. 13 on page 8
    IRX = (0.45 * data_redshift + 0.35) * (data_lgMstar-10.5) + 1.2 # log10(LIR/LUV)
    IRX[(data_redshift>3.0)] = (0.45 * 3.0 + 0.35) * (data_lgMstar[(data_redshift>3.0)]-10.5) + 1.2 # log10(LIR/LUV)
    SFR_total = data_SFR
    SFR_IR = SFR_total * 1.0/(1.0+(1.0/(10**IRX)))
    SFR_UV = SFR_total * 1.0/(1.0+(10**IRX))
    if os.path.isfile('dump_SFR_IR_SFR_UV_IRX_Schreiber2017.fits'):
        os.remove('dump_SFR_IR_SFR_UV_IRX_Schreiber2017.fits')
    if not os.path.isfile('dump_SFR_IR_SFR_UV_IRX_Schreiber2017.fits'):
        dump_data_table = Table( {'z':data_redshift, 'lgMstar':data_lgMstar, 'SFR':data_SFR, 'SFR_IR':SFR_IR, 'SFR_UV':SFR_UV, 'IRX':IRX} )
        dump_data_table.write('dump_SFR_IR_SFR_UV_IRX_Schreiber2017.fits', overwrite=True)
        print('Written "%s"'%('dump_SFR_IR_SFR_UV_IRX_Schreiber2017.fits'))
    data_SFR = SFR_IR
    output_name = output_name + '_applied_IRX_Schreiber2017'
# 
if do_apply_IRX_Whitaker2017:
    # Whitaker2017 Eq. 1 on page 3
    a = 1.96e9
    b = -2.277
    SFR_total = data_SFR # SFR_total is SFR(UV+IR)
    f_obscured = 1.0 / (1.0 + a * np.exp(b * data_lgMstar)) # f_obscured = SFR(IR)/SFR(UV+IR)
    # Whitaker2017 Eq. 2 on page 5
    #zlo = np.array([0.02, 0.5, 1.0, 1.5, 2.0])
    #zhi = np.array([0.05, 1.0, 1.5, 2.0, 2.5])
    #a = np.array([1.234, 2.092, 3.917, 6.806, 10.701])
    #b = np.array([-2.858, -2.384, -2.596, -2.673, -2.516])
    #SFR_total = data_SFR # SFR_total is SFR(UV+IR)
    #f_obscured = 1.0 / (1.0 + a * np.exp(b * np.log10(SFR_total))) # f_obscured = SFR(IR)/SFR(UV+IR)
    # 
    SFR_IR = SFR_total * f_obscured
    SFR_UV = SFR_total * (1.0-f_obscured)
    IRX = np.log10(SFR_IR/SFR_UV)
    if os.path.isfile('dump_SFR_IR_SFR_UV_IRX_Whitaker2017.fits'):
        os.remove('dump_SFR_IR_SFR_UV_IRX_Whitaker2017.fits')
    if not os.path.isfile('dump_SFR_IR_SFR_UV_IRX_Whitaker2017.fits'):
        dump_data_table = Table( {'z':data_redshift, 'lgMstar':data_lgMstar, 'SFR':data_SFR, 'SFR_IR':SFR_IR, 'SFR_UV':SFR_UV, 'f_obscured':f_obscured, 'IRX':IRX} )
        dump_data_table.write('dump_SFR_IR_SFR_UV_IRX_Whitaker2017.fits', overwrite=True)
        print('Written "%s"'%('dump_SFR_IR_SFR_UV_IRX_Whitaker2017.fits'))
    data_SFR = SFR_IR
    output_name = output_name + '_applied_IRX_Whitaker2017'



# 
# set L_IR to SFR_IR only
# 
L_IR = ((data_SFR*1e10)) # Lsun
lgL_IR_min = 8.0
lgL_IR_max = 14.0



# 
# set label_MS from current folder name
# 
label_MS = ''
if re.match(r'.*_using_([a-zA-Z]+)([0-9]+)_MS(_.*|)', os.path.basename(os.getcwd())):
    label_MS = 'with %s MS'%(re.sub(r'.*_using_([a-zA-Z]+)([0-9]+)_MS(_.*|)', r'\1+\2', os.path.basename(os.getcwd())))










z_edges = [0.02, 0.25, 0.50, 0.75, 1.00, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
ncol = 4
nrow = int(np.ceil((len(z_edges)-1) / float(ncol)))

fig = plt.figure(figsize=(15.0, 9.0))
fig.subplots_adjust(left=0.09, right=0.97, bottom=0.10, top=0.98, wspace=0.20, hspace=0.30)

# 
# loop z bin
for i in range(len(z_edges)-1):
    # 
    print('z %s - %s'%(z_edges[i], z_edges[i+1]))
    # 
    z = (z_edges[i]+z_edges[i+1])/2.0
    # 
    # 
    # calc comving colume
    comoving_volume = ((cosmo.comoving_volume(z_edges[i+1]) - cosmo.comoving_volume(z_edges[i])) / (4.0*np.pi*u.steradian) * obs_area.to(u.steradian))
    print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    # 
    differntial_z_list = np.linspace(z_edges[i], z_edges[i+1], num=10, endpoint=True)
    comoving_volume = np.sum((cosmo.differential_comoving_volume(differntial_z_list[1:]) * np.diff(differntial_z_list) * obs_area.to(u.steradian)))
    print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    # 
    # 
    # 
    # ax
    ax = fig.add_subplot(nrow, ncol, i+1)
    # 
    # xylim
    #ax.set_xlim([-1.0, 4.0])
    #ax.set_ylim([-7.5, 0.5])
    ax.set_yscale('log')
    # 
    # 
    # prepare legend
    legend_handles = []
    legend_labels = []
    # 
    # 
    # 
    # select and plot data histogram
    data_selection = np.logical_and.reduce((data_redshift >= z_edges[i], data_redshift < z_edges[i+1], data_SFR > 0))
    bin_hists, bin_edges = np.histogram(np.log10(L_IR[data_selection]), bins=np.arange(lgL_IR_min, lgL_IR_max, 0.2) )
    bin_centers = (bin_edges[0:-1]+bin_edges[1:])/2.0
    bin_widths = (bin_edges[1:]-bin_edges[0:-1])
    bin_hists = bin_hists / (bin_edges[1]-bin_edges[0]) / comoving_volume.value # Mpc-3 dex-1
    bin_mask = (bin_hists>0)
    plot_data_1 = ax.bar(bin_centers[bin_mask], bin_hists[bin_mask], width=bin_widths[bin_mask], alpha=0.6, align='edge' )
    label_data_1 = 'Mock galaxies\' LF' + '\n  ' + r'from SMF$\times$SFMS'
    if label_MS != '':
        label_data_1 += '\n  ' + label_MS
    legend_handles.append(plot_data_1)
    legend_labels.append(label_data_1)
    # 
    # 
    # 
    # plot data points from Sanders2003 AJ Table 6 (z=0)
    #<TODO>#
    # 
    # 
    # 
    # fix ylim
    ax.set_ylim(ax.get_ylim())
    # 
    # 
    # 
    # plot model line
    try:
        lgL = np.arange(lgL_IR_min, lgL_IR_max, 0.2)
        lgPhi = calc_IR_LF_Gruppioni2013(z, lgL)
        #print(lgL, lgPhi)
        Phi = 10**lgPhi
        plot_data_2 = ax.plot(lgL, Phi, ls='dashed', lw=2.5, color='C1')
        label_data_2 = 'Gruppioni+2013 LF' # Gruppioni+2013 used Herschel PACS 160um in four deep fields and 
        legend_handles.append(plot_data_2[0])
        legend_labels.append(label_data_2)
    except:
        legend_handles.append(plot_data_2[0])
        legend_labels.append(label_data_2)
    # 
    # 
    # 
    # label current bin
    ax.text(0.97, 0.97, r'$z = %0.1f - %0.1f$'%(z_edges[i], z_edges[i+1]), transform=ax.transAxes, fontsize=16, va='top', ha='right')
    # 
    # ax tick params
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.yscale('log')
    #plt.xlabel(r'$\log_{10} (\mathrm{IR})$', fontsize=20, labelpad=2)
    plt.xlabel(r'$\log_{10} (L_{\mathrm{IR}} / \mathrm{[L_{\odot}]})$', fontsize=17, labelpad=2)
    # 
    # ylabel
    if i%ncol == 0 and int(i/ncol) == int((nrow-1)/2):
        plt.ylabel(r'$\Phi \ [\mathrm{Mpc^{-3}\,dex^{-1}}]$', fontsize=20, labelpad=20)
    # 
    # legend
    if i+1 == len(z_edges)-1:
        plot_legend1 = plt.legend(\
                        legend_handles, 
                        legend_labels, 
                        fontsize=16, loc='center left', bbox_to_anchor=(1.1, 0.0, 0.2, 1.0), 
                        #borderpad=0.6, borderaxespad=0.6, handlelength=2.8, 
                       )
        #leg.get_frame().set_edgecolor('#cccccc')
        #leg.get_frame().set_linewidth(2.0)
        #plt.setp(plot_legend1.get_texts(), fontsize='18')
        ax.add_artist(plot_legend1)
    # 
    # show y minor ticks
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=0.1))
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    #ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=99))
    #ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=99))
    #ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    # 
    # xylim
    #ax.set_ylim([ax.get_ylim()[0], ax.get_ylim()[1]*1.5])
    # 
    # grid
    ax.grid(True, ls='dotted', c='#cccccc')


# 
# savefig
fig.savefig('%s.pdf'%(output_name), transparent=True)
print('Output to "%s"!' % ('%s.pdf'%(output_name)) )
os.system('open "%s"' % ('%s.pdf'%(output_name)) )


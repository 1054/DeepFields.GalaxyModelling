#!/usr/bin/env python
#



import os, sys, json, numpy, astropy, matplotlib, subprocess
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

from calc_galaxy_luminosity_function import (calc_IR_250um_LF_Koprowski2017)

ln = np.log






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


L_IR_250um = ((data_SFR*1e10)*3.839e26)/(2.99792458e8/250.0e-6) # W Hz-1
lgL_IR_250um_min = 24.0
lgL_IR_250um_max = 27.0
raise NotImplementedError('To assign SED for each galaxy')










z_edges = [0.02, 0.25, 0.50, 0.75, 1.00, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
ncol = 4
nrow = int(np.ceil((len(z_edges)-1) / float(ncol)))

fig = plt.figure(figsize=(15.0, 10.0))
fig.subplots_adjust(left=0.10, right=0.99, bottom=0.10, top=0.99, wspace=0.20, hspace=0.30)

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
    # select and plot data histogram
    data_selection = np.logical_and.reduce((data_redshift >= z_edges[i], data_redshift < z_edges[i+1], data_SFR > 0))
    bin_hists, bin_edges = np.histogram(np.log10(L_IR_250um[data_selection]), bins=np.arange(lgL_IR_250um_min, lgL_IR_250um_max, 0.2) )
    bin_centers = (bin_edges[0:-1]+bin_edges[1:])/2.0
    bin_widths = (bin_edges[1:]-bin_edges[0:-1])
    bin_mask = (bin_hists>0)
    bin_hists = bin_hists / (bin_edges[1]-bin_edges[0]) / comoving_volume.value # Mpc-3 dex-1
    plot_data_1 = ax.bar(bin_centers[bin_mask], bin_hists[bin_mask], width=bin_widths[bin_mask], alpha=0.6, align='edge' )
    label_data_1 = 'mock galaxies\' LF' + '\n  ' + r'from SMF$\times$SFMS'
    legend_handles.append(plot_data_1)
    legend_labels.append(label_data_1)
    # 
    # 
    # fix ylim
    ax.set_ylim(ax.get_ylim())
    # 
    # 
    # 
    # plot model line
    if z>=0.5 and z<=4.5:
        lgL = np.arange(lgL_IR_250um_min, lgL_IR_250um_max, 0.2)
        lgPhi = calc_IR_250um_LF_Koprowski2017(z, lgL)
        Phi = 10**lgPhi
        plot_data_2 = ax.plot(lgL, Phi, ls='dashed', lw=2.5, color='C1')
        label_data_2 = 'Koprowski+2017 LF'
        legend_handles.append(plot_data_2[0])
        legend_labels.append(label_data_2)
    else:
        pass
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
    plt.xlabel(r'$\log_{10} (L_{\mathrm{250\,{\mu}m,\,rest}} / \mathrm{[W \, Hz^{-1}]})$', fontsize=17, labelpad=2)
    # 
    # ylabel
    if i%ncol == 0 and int(i/ncol) == int((nrow-1)/2):
        plt.ylabel(r'$\Phi \ [\mathrm{Mpc^{-3}\,dex^{-1}}]$', fontsize=18)
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
fig.savefig('Plot_IR_250um_LF_bin_by_z.pdf', transparent=True)
print('Output to "%s"!' % ('Plot_IR_250um_LF_bin_by_z.pdf') )
os.system('open "%s"' % ('Plot_IR_250um_LF_bin_by_z.pdf') )


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
#matplotlib.rcParams['text.usetex'] = True 
#matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
#matplotlib.rcParams['font.sans-serif'] = 'cm'
#matplotlib.rcParams['font.family'] = 'sans-serif'

from calc_galaxy_main_sequence import (
                                        calc_SFR_MS_Speagle2014, 
                                        calc_SFR_MS_Sargent2014, 
                                        calc_SFR_MS_Whitaker2014, 
                                        calc_SFR_MS_Bethermin2015, 
                                        calc_SFR_MS_Schreiber2015, 
                                        calc_SFR_MS_Lee2015, 
                                        calc_SFR_MS_Tomczak2016, 
                                        calc_SFR_MS_Pearson2018, 
                                        calc_SFR_MS_Leslie20180901, 
                                        calc_SFR_MS_Leslie20190111, 
                                        calc_SFR_MS_Leslie20190515, 
                                        calc_SFR_MS_Leslie20190710, 
                                        calc_SFR_MS_Scoville2017, 
                                        )

from calc_galaxy_stellar_mass_function import (calc_SMF_Davidzon2017, calc_SMF_Ilbert2013, calc_SMF_Peng2010, calc_SMF_Wright2018_single_component, calc_SMF_Wright2018_double_component, calc_SMF_dzliu2018)





# 
# Read data points
# 
tb = Table.read('datatable_generated_galaxies_with_coordinates.fits')
data_lgMstar = tb['lgMstar'].data.flatten()
data_Mstar = 10**data_lgMstar
data_lgSFR = tb['lgSFR'].data.flatten()
data_SFR = 10**data_lgSFR
data_redshift = tb['z'].data.flatten()
#sys.exit()




# 
# Read data2 points
# 
tb2 = Table.read('/Users/dzliu/Work/AlmaCosmos/Catalogs/Galaxy_Modelling_SIDES/Mock_cat_Bethermin2017.fits')
data2_Mstar = tb2['MSTAR'].data.flatten()
data2_lgMstar = np.log10(data2_Mstar)
data2_SFR = tb2['SFR'].data.flatten()
data2_lgSFR = np.log10(data2_SFR)
data2_redshift = tb2['REDSHIFT'].data.flatten()

#sys.exit()



# 
# Check MS
# 
plot_label_MS = ''
current_dir = os.path.basename(os.getcwd())
if re.match(r'^.*_using_([a-zA-Z0-9]+)_MS([_].*|)$', current_dir):
    plot_label_MS = '(%s MS)'%(re.sub(r'^.*_using_([a-zA-Z0-9]+)_MS([_].*|)$', r'\1', current_dir))











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
    data_selection = np.logical_and.reduce((data_redshift >= z_edges[i], data_redshift < z_edges[i+1], data_SFR > 0, data_lgMstar > 10.25, data_lgMstar < 10.75))
    bin_hists, bin_edges = np.histogram(data_lgSFR[data_selection], bins=np.arange(-1.0, 4.5, 0.1) )
    bin_centers = (bin_edges[0:-1]+bin_edges[1:])/2.0
    bin_widths = (bin_edges[1:]-bin_edges[0:-1])
    bin_mask = (bin_hists>0)
    plot_data_1 = ax.bar(bin_centers[bin_mask], bin_hists[bin_mask], width=bin_widths[bin_mask], alpha=0.6, align='edge' )
    label_data_1 = 'dzliu model'+'\n'+plot_label_MS+'\n'+'lgMstar=10.25-10.75'
    legend_handles.append(plot_data_1)
    legend_labels.append(label_data_1)
    # 
    # 
    # select and plot data histogram
    data2_selection = np.logical_and.reduce((data2_redshift >= z_edges[i], data2_redshift < z_edges[i+1], data2_SFR > 0, data2_lgMstar > 10.25, data2_lgMstar < 10.75))
    bin_hists, bin_edges = np.histogram(data2_lgSFR[data2_selection], bins=np.arange(-1.0, 4.5, 0.1) )
    bin_centers = (bin_edges[0:-1]+bin_edges[1:])/2.0
    bin_widths = (bin_edges[1:]-bin_edges[0:-1])
    bin_mask = (bin_hists>0)
    plot_data_2 = ax.bar(bin_centers[bin_mask], bin_hists[bin_mask], width=bin_widths[bin_mask], alpha=0.6, align='edge' )
    label_data_2 = 'Bethermin+2017\nlgMstar=10.25-10.75'
    legend_handles.append(plot_data_2)
    legend_labels.append(label_data_2)
    # 
    # label current bin
    ax.text(0.97, 0.97, r'$z = %0.1f - %0.1f$'%(z_edges[i], z_edges[i+1]), transform=ax.transAxes, fontsize=16, va='top', ha='right')
    #ax.text(0.95, 0.15, r'$log_{10} M_{\star} = %0.1f - %0.1f$'%(10.25, 10.75), transform=ax.transAxes, fontsize=16, va='bottom', ha='right')
    # 
    # ax tick params
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.yscale('log')
    plt.xlabel(r'$\log_{10} (\mathrm{SFR})$', fontsize=20, labelpad=2)
    # 
    # ylabel
    if i%ncol == 0 and int(i/ncol) == int((nrow-1)/2):
        plt.ylabel(r'$N$', fontsize=18)
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
fig.savefig('Plot_SFR_Histogram_bin_by_z_comparing_with_Bethermin2017.pdf', transparent=True)
print('Output to "%s"!' % ('Plot_SFR_Histogram_bin_by_z_comparing_with_Bethermin2017.pdf') )
os.system('open "%s"' % ('Plot_SFR_Histogram_bin_by_z_comparing_with_Bethermin2017.pdf') )


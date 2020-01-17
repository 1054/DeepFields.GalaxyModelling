#!/usr/bin/env python
#



import os, sys, json, numpy, astropy, matplotlib
#matplotlib.use('Qt5Agg')

from astropy.table import Table
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






# 
# User Setting
# 





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
data_redshift = tb['z'].data.flatten()
#sys.exit()










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
    lgMstar = np.linspace(8.0, 13.0, num=100, endpoint=True)
    # 
    # SF MS curve
    SFR_MS_Leslie2019 = calc_SFR_MS_Leslie20190710(z, lgMstar)
    lgSFR_MS_Leslie2019 = np.log10(SFR_MS_Leslie2019)
    # 
    # ax
    ax = fig.add_subplot(nrow, ncol, i+1)
    # 
    # xylim
    ax.set_xlim([8.0, 12.0])
    ax.set_ylim([-2.5, 3.5])
    # 
    # plot SF MS curve
    plot_curve_1a = ax.plot(lgMstar, lgSFR_MS_Leslie2019, c='red', ls='solid', solid_capstyle='butt', alpha=0.35, lw=20)
    plot_curve_1b = ax.plot(lgMstar, lgSFR_MS_Leslie2019, c='red', ls='solid', solid_capstyle='butt', alpha=0.95, lw=2)
    plot_curve_1 = (plot_curve_1a[0], plot_curve_1b[0])
    label_curve_1 = r'Leslie+2019'
    # 
    # select and plot data points
    data_selection = np.logical_and(np.logical_and(data_redshift >= z_edges[i], data_redshift < z_edges[i+1]), data_SFR > 0)
    x = np.log10(data_Mstar[data_selection])
    y = np.log10(data_SFR[data_selection])
    bins = [100, 100] # number of bins
    hists, yedges, xedges = np.histogram2d(y, x, bins=bins, range=[ax.get_ylim(), ax.get_xlim()], density=True) # histogram the data
    #c = np.array([ hists[np.argmax(a<=xedges[1:]), np.argmax(b<=yedges[1:])] for a,b in zip(x,y) ]) # Sort the points by density, so that the densest points are plotted last
    #idx = c.argsort()
    #x2, y2, c2 = x[idx], y[idx], c[idx]
    #plot_data_1a = ax.scatter(x2, y2, c=c2, s=0.5, cmap='jet', marker='.')
    #plot_data_1b = ax.contour(hists, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    plot_data_1b = ax.contour((xedges[0:-1]+xedges[1:])/2.0, (yedges[0:-1]+yedges[1:])/2.0, hists, levels=np.logspace(-1.5, 0.0, num=12, endpoint=True), linewidths=0.5, cmap='jet')
    plot_data_1 = plot_data_1b.collections[0]
    #label_data_1 = 'dzliu model with\nLeslie+2019 MS'
    #label_data_1 = 'dzliu model with\nSpeagle+2014 MS'
    label_data_1 = 'dzliu model'
    # 
    # label current bin
    ax.text(0.95, 0.05, r'$z = %0.1f - %0.1f$'%(z_edges[i], z_edges[i+1]), transform=ax.transAxes, fontsize=16, va='bottom', ha='right')
    # 
    # ax tick params
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.yscale('log')
    plt.xlabel(r'$\log_{10} (M_{\star}/\mathrm{M_{\odot}})$', fontsize=20, labelpad=2)
    # 
    # ylabel
    if i%ncol == 0:
        plt.ylabel(r'$\log_{10} \, \mathrm{SFR}_{\mathrm{MS}}$', fontsize=20)
    # 
    # legend
    if i+1 == len(z_edges)-1:
        plot_legend1 = plt.legend(\
                        [ plot_curve_1, plot_data_1 ], 
                        [ label_curve_1, label_data_1 ], 
                        fontsize=18, loc='center left', bbox_to_anchor=(1.2, 0.0, 0.2, 1.0), 
                        #borderpad=0.6, borderaxespad=0.6, handlelength=2.8, 
                       )
        #leg.get_frame().set_edgecolor('#cccccc')
        #leg.get_frame().set_linewidth(2.0)
        #plt.setp(plot_legend1.get_texts(), fontsize='18')
        ax.add_artist(plot_legend1)
    # 
    # show y minor ticks
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
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
fig.savefig('Plot_Mstar_SFR_bin_by_z.pdf', transparent=True)
print('Output to "%s"!' % ('Plot_Mstar_SFR_bin_by_z.pdf') )
os.system('open "%s"' % ('Plot_Mstar_SFR_bin_by_z.pdf') )


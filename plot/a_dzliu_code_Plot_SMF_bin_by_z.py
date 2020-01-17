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

from calc_galaxy_stellar_mass_function import (calc_SMF_Davidzon2017, calc_SMF_Ilbert2013, calc_SMF_Peng2010, calc_SMF_Wright2018_single_component, calc_SMF_Wright2018_double_component, calc_SMF_dzliu2018)





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
data_redshift = tb['z'].data.flatten()
#sys.exit()



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
    lgMstar_edges = np.linspace(8.0, 13.0, num=100, endpoint=True)
    lgMstar_centers = (lgMstar_edges[:-1]+lgMstar_edges[1:])/2.0
    # 
    comoving_volume = ((cosmo.comoving_volume(z_edges[i+1]) - cosmo.comoving_volume(z_edges[i])) / (4.0*np.pi*u.steradian) * obs_area.to(u.steradian))
    print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    # 
    #comoving_z_list = np.linspace(z_edges[i], z_edges[i+1], num=100, endpoint=True)
    #comoving_volume = np.sum((cosmo.comoving_volume(comoving_z_list[1:]) - cosmo.comoving_volume(comoving_z_list[0:-1])) / (4.0*np.pi*u.steradian) * obs_area.to(u.steradian))
    #print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    # 
    differntial_z_list = np.linspace(z_edges[i], z_edges[i+1], num=10, endpoint=True)
    comoving_volume = np.sum((cosmo.differential_comoving_volume(differntial_z_list[1:]) * np.diff(differntial_z_list) * obs_area.to(u.steradian)))
    print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    ##print(cosmo.de_density_scale(z)) # should be 1
    ##print(cosmo._Ogamma0, cosmo._Onu0)
    ##print(cosmo.efunc(z), np.sqrt(0.27*(1.+z)**3 + 0*(1.+z)**2 + 0.73) ) # checked consistent
    ##print(cosmo._hubble_distance, 2.997902458e5/70 ) # checked consistent
    #sys_lumdist_output = subprocess.getoutput("/Users/dzliu/Cloud/Github/Crab.Toolkit.PdBI/bin/lumdist -h0 70 -verbose %s | grep 'lumdist d_L=' | sed -e 's/=/ /g'"%(z))
    ##print(cosmo.angular_diameter_distance(z), sys_lumdist_output.split()[8], '(z = %s)'%(z), sys_lumdist_output ) # 
    #dH_astropy = cosmo._hubble_distance
    #Ez_astropy = cosmo.efunc(z)
    #dA_astropy = cosmo.angular_diameter_distance(z)
    #dH_dzliu = 2.997902458e5/70
    #Ez_dzliu = np.sqrt(0.27*(1.+z)**3 + 0*(1.+z)**2 + 0.73)
    #dA_dzliu = float(sys_lumdist_output.split()[8])
    #print(dH_astropy/Ez_astropy*dA_astropy, dH_dzliu/Ez_dzliu*dA_dzliu) # chekced consistent
    #zp1 = 1.0 + z
    #print(cosmo.differential_comoving_volume(z), dH_astropy*((zp1*dA_astropy)**2)/Ez_astropy, dH_dzliu/Ez_dzliu*dA_dzliu**2*zp1**2) # chekced consistent
    # z = 2.0 - 2.5, area = 0.0006092348395183178 steradian, dVc = 43627725623.05944, 
    # 43627725623.05944 * 0.25 / 10
    # 
    # 
    # ax
    ax = fig.add_subplot(nrow, ncol, i+1)
    # 
    # xylim
    ax.set_xlim([8.0, 12.0])
    ax.set_ylim([-7.5, 0.5])
    # 
    # calc and plot SMF curve
    legend_handles = []
    legend_labels = []
    # 
    # 
    # 
    # plot curve from Davidzon+2017
    try:
        #label_curve_1 = r'SMF dzliu'
        #data_curve_1 = calc_SMF_dzliu2018(z, lgMstar_centers, verbose=False) # lgPhiMstar, per dex
        label_curve_1 = r'SMF Davidzon+2017'
        data_curve_1 = calc_SMF_Davidzon2017(z, lgMstar_centers, galaxy_type='SFG') # lgPhiMstar, per dex
        plot_curve_1 = ax.plot(lgMstar_centers, data_curve_1, c='red', ls='solid', solid_capstyle='butt', alpha=0.8, lw=2)
        legend_handles.append(plot_curve_1[0])
        legend_labels.append(label_curve_1)
    except ValueError as error:
        print(error)
        pass
    # 
    # 
    # 
    # plot curve from dzliu model
    try:
        label_curve_2 = r'SMF renorm to $\int\rho_\mathrm{SFR}$'
        data_curve_2 = calc_SMF_dzliu2018(z, lgMstar_centers, verbose=False, galaxy_type='SFG') # lgPhiMstar, per dex
        plot_curve_2 = ax.plot(lgMstar_centers, data_curve_2, c='green', ls='dotted', solid_capstyle='butt', alpha=0.8, lw=2)
        legend_handles.append(plot_curve_2[0])
        legend_labels.append(label_curve_2)
    except ValueError as error:
        print(error)
        pass
    # 
    # 
    # 
    # select and plot data histogram
    data_selection = np.logical_and(np.logical_and(data_redshift >= z_edges[i], data_redshift < z_edges[i+1]), data_SFR > 0)
    data_count_lgMstar = []
    for j in range(len(lgMstar_edges)-1):
        data_selection_by_lgMstar = np.logical_and(np.logical_and(data_lgMstar>=lgMstar_edges[j], data_lgMstar<lgMstar_edges[j+1]), data_selection)
        data_count_lgMstar.append(np.count_nonzero(data_selection_by_lgMstar) / (lgMstar_edges[j+1]-lgMstar_edges[j])) # Msun per dex
    data_count_lgMstar = np.array(data_count_lgMstar)
    data_Phi_lgMstar = data_count_lgMstar / comoving_volume.value # we are already counting in lgMstar, so the unit will already be dex^{-1}
    data_Phi_lgMstar = np.log10(data_Phi_lgMstar)
    #plot_data_1 = ax.step(lgMstar_centers, data_Phi_lgMstar, where='mid', alpha=0.6 )
    plot_data_1 = ax.fill_between(lgMstar_centers, data_Phi_lgMstar*0.0+ax.get_ylim()[0], data_Phi_lgMstar, step='mid', alpha=0.6 )
    #plot_data_1 = ax.bar(lgMstar_edges[0:-1], data_Phi_lgMstar, width=lgMstar_edges[1:]-lgMstar_edges[0:-1], align='edge', bottom=ax.get_ylim()[0], alpha=0.6 )
    label_data_1 = 'Mock galaxies\' SMF'+'\n  '+r'renorm to $\int\rho_\mathrm{SFR}$'
    if os.path.basename(os.getcwd()).find('SMF_no_renorm') > 0:
        label_data_1 = 'Mock galaxies\' SMF'
    if label_MS != '':
        label_data_1 += '\n  ' + label_MS
    #legend_handles.append(plot_data_1[0])
    legend_handles.append(plot_data_1)
    legend_labels.append(label_data_1)
    # 
    # 
    # 
    # plot data from Davidzon+2017
    #<TODO>
    # 
    # 
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
    if i%ncol == 0 and int(i/ncol) == int((nrow-1)/2):
        plt.ylabel(r'$\log_{10} \, \Phi \,/\, [\mathrm{Mpc^{-3} \, dex^{-1}}]$', fontsize=20, labelpad=20)
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
fig.savefig('Plot_SMF_bin_by_z.pdf', transparent=True)
print('Output to "%s"!' % ('Plot_SMF_bin_by_z.pdf') )
#os.system('open "%s"' % ('Plot_SMF_bin_by_z.pdf') )


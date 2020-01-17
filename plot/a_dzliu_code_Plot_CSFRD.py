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
                                        calc_SFR_MS_Leslie20191212, 
                                        calc_SFR_MS_Scoville2017, 
                                        )

from calc_cosmic_star_formation_rate_density import (calc_CSFRD_Madau2014, convert_age_to_z)





# 
# User Setting
# 
#obs_area = 7200*u.arcmin*u.arcmin # 1.4*1.4*u.deg*u.deg
#obs_area = 1.4*1.4*u.deg*u.deg
obs_area = 1.5546582999901375*u.deg*u.deg
#obs_area = 2.0*u.deg*u.deg
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
# def
# 
def tick_function(X):
    V = cosmo.age(X).value
    return ['%0.1f' % t for t in V]


# 
# fig
# 
fig = plt.figure(figsize=(6.8,4.8))
fig.subplots_adjust(left=0.15, right=0.95, bottom=0.105, top=0.885)
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel('Redshift', fontsize=16, labelpad=1)
ax1.set_ylabel(r'$\log_{10} \, \rho_{\mathrm{SFR}}$ [$\mathrm{M_{\odot}\,yr^{-1}\,Mpc^{-3}}$]', fontsize=17, labelpad=15)
ax1.tick_params(axis='both', labelsize=14)
ax1.tick_params(direction='in', axis='both', which='both')
ax1.tick_params(top=False, right=True, which='both')

my_tick_locations = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
ax1.set_xticks(my_tick_locations)
#ax1.set_xlim([-0.3, np.max(my_tick_locations)])
ax1.set_xlim([-0.2, 5.5])
ax1.set_ylim([-2.0, -0.2])
ax1.grid(True, ls='--', lw=0.25)

#new_tick_locations = convert_age_to_z([13.7, 2.0, 1.0, 0.5, 0.3]) #<20190915><BUGGY># 
new_tick_locations = convert_age_to_z([cosmo.age(0).value, 5.0, 3.0, 2.0, 1.0, 0.7, ])
new_tick_locations = new_tick_locations[np.argwhere(np.logical_and(new_tick_locations >= ax1.get_xlim()[0], new_tick_locations <= ax1.get_xlim()[1])).flatten()]
print('new_tick_locations', new_tick_locations)
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r"Cosmic Age [$\mathrm{Gyr}$]", fontsize=16, labelpad=6)
ax2.minorticks_off()
ax2.grid(None)

# show y minor ticks
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
#ax1.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=99))
#ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=99))

# 
# z
z_edges = np.linspace(0.0, 6.0, num=30, endpoint=True) # np.array([0.02, 0.25, 0.50, 0.75, 1.00, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0])
z_centers = (z_edges[0:-1] + z_edges[1:]) / 2.0
lgPhi_SFR = z_centers * 0.0 - 99

# 
# loop z bin
for i in range(len(z_edges)-1):
    # 
    print('z %s - %s'%(z_edges[i], z_edges[i+1]))
    # 
    z = z_centers[i]
    # 
    #comoving_volume = ((cosmo.comoving_volume(z_edges[i+1]) - cosmo.comoving_volume(z_edges[i])) / (4.0*np.pi*u.steradian) * obs_area.to(u.steradian))
    #print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    # 
    differntial_z_list = np.linspace(z_edges[i], z_edges[i+1], num=10, endpoint=True)
    comoving_volume = np.sum((cosmo.differential_comoving_volume(differntial_z_list[1:]) * np.diff(differntial_z_list) * obs_area.to(u.steradian)))
    print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    # 
    # select and count CSFRD
    data_selection = np.logical_and.reduce((data_redshift >= z_edges[i], data_redshift < z_edges[i+1], data_SFR > 0, data_lgMstar >= 9.0))
    data_Phi_SFR = np.sum(data_SFR[data_selection]) / comoving_volume.value
    lgPhi_SFR[i] = np.log10(data_Phi_SFR)



# 
# 
Phi_SFR_MD14 = calc_CSFRD_Madau2014(z_centers)
lgPhi_SFR_MD14 = np.log10(Phi_SFR_MD14)
ax1.plot(z_centers, lgPhi_SFR_MD14, c='red', ls='solid', solid_capstyle='butt', alpha=0.8, lw=2, label=r'SMF MD14')
# 
# 
plot_label = r'dzliu model (lgMstar$\gtrsim$9.0)'
current_dir = os.path.basename(os.getcwd())
if re.match(r'^.*_using_([a-zA-Z0-9]+)_MS([_].*|)$', current_dir):
    plot_label = plot_label + '\n(%s MS)'%(re.sub(r'^.*_using_([a-zA-Z0-9]+)_MS([_].*|)$', r'\1', current_dir))
ax1.step(z_centers, lgPhi_SFR, where='mid', alpha=0.6, label=plot_label)
# 
# 
#plot_legend1 = plt.legend(\
#                legend_handles, 
#                legend_labels, 
#                fontsize=16, loc='upper right', 
#                #borderpad=0.6, borderaxespad=0.6, handlelength=2.8, 
#               )
#ax1.add_artist(plot_legend1)
ax1.legend(loc='upper left', ncol=2, framealpha=0.5)


# 
# savefig
fig.savefig('Plot_CSFRD.pdf', transparent=True)
print('Output to "%s"!' % ('Plot_CSFRD.pdf') )
os.system('open "%s"' % ('Plot_CSFRD.pdf') )


#!/usr/bin/env python
# 
"""
This code analyzes galaxy evolution prescriptions and generates some plots.
"""
from __future__ import print_function
import os, sys, re, json, shutil
import numpy as np
from astropy.table import Table
script_dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(os.path.join(script_dir, 'python', 'lib', 'crab', 'crabgalaxy'))
#sys.path.append(os.path.join(script_dir, 'python', 'lib', 'crab', 'crabtable'))
sys.path.append(os.getenv('HOME')+'/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code') #<TODO># replace by a3g
#from CrabGalaxy import CrabGalaxy
#from CrabTable import CrabTableReadInfo
from calc_galaxy_stellar_mass_function import (calc_SMF_Davidzon2017, calc_SMF_Ilbert2013, calc_SMF_Peng2010, calc_SMF_Wright2018_single_component, calc_SMF_Wright2018_double_component, calc_SMF_dzliu2018)
from calc_galaxy_main_sequence import (calc_SFR_MS_Speagle2014, calc_SFR_MS_Sargent2014, calc_SFR_MS_Schreiber2015, calc_SFR_MS_Leslie20190710)
from calc_cosmic_star_formation_rate_density import (calc_CSFRD_Madau2014, calc_CSFRD_Liu2018, convert_age_to_z)
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from setup_matplotlib import setup_matplotlib; setup_matplotlib()
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=69.8, Om0=0.3, Tcmb0=2.725) # Freedman 2019ApJ...882...34F

#InfoDict = CrabTableReadInfo(InfoFile, verbose=0)


global calc_SFR_MS
#calc_SFR_MS = calc_SFR_MS_Leslie20190710 ; label_SF_MS = 'Leslie+2019'
#calc_SFR_MS = calc_SFR_MS_Speagle2014 ; label_SF_MS = 'Speagle+2014'
#calc_SFR_MS = calc_SFR_MS_Sargent2014 ; label_SF_MS = 'Sargent+2014'
calc_SFR_MS = calc_SFR_MS_Schreiber2015 ; label_SF_MS = 'Schreiber+2015'


#global SF_MS_scatter
#SF_MS_scatter = 0.3


global calc_CSFRD
calc_CSFRD = calc_CSFRD_Madau2014 ; label_CSFRD = 'MD14'


#global calc_SMF
#calc_SMF = calc_SMF_dzliu2018 ; label_SMF = 'Liu2019'




"""
Function for SF_MS_scatter
"""
def SF_MS_scatter(lgMstar):
    # See https://arxiv.org/pdf/1812.07057.pdf Figure 5
    # It seems their SED fitting derived main sequence scatter is most stable at around 0.29 dex.
    #return np.interp(lgMstar, [6.0, 10.0, 11.0, 14.0], [0.3, 0.3, 0.4, 0.4])
    return np.interp(lgMstar, [6.0, 10.0, 11.0, 14.0], [0.29, 0.29, 0.29, 0.29])



"""
Functions for axis tick display
"""
def lower_ticks_for_CSFRD_vs_opz(opz):
    # input 1+z, return z
    V = opz-1.0
    return ['%0g' % t for t in V]

def upper_ticks_for_CSFRD_vs_opz(opz):
    # input 1+z, return age
    V = cosmo.age(opz-1.0).value
    return ['%0.1f' % t for t in V]

def lower_ticks_for_CSFRD_vs_z(z):
    # input 1+z, return z
    V = z
    return ['%0g' % t for t in V]

def upper_ticks_for_CSFRD_vs_z(z):
    # input 1+z, return age
    V = cosmo.age(z).value
    return ['%0.1f' % t for t in V]






"""
Plot cosmic star formation rate (SFR) density versus redshift diagram
"""
def plot_CSFRD():
    # 
    global calc_CSFRD
    global label_CSFRD
    global calc_SFR_MS
    global label_SF_MS
    #calc_CSFRD = calc_CSFRD_Madau2014
    #calc_CSFRD = calc_CSFRD_Liu2018
    # 
    output_name = 'Plot_CSFRD'
    # 
    set_x_opz = False
    # 
    fig = plt.figure(figsize=(6.8,4.8)) # before 20190722: figsize=(6.8,5.2) # before 20190512: figsize=(7.2,5.8)
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.105, top=0.885)
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_xlabel('Redshift', fontsize=16, labelpad=1)
    ax1.set_ylabel(r'$\log_{10} \, \rho_{\mathrm{SFR}}$ [$\mathrm{M_{\odot}\,yr^{-1}\,Mpc^{-3}}$]', fontsize=17, labelpad=15)
    ax1.tick_params(axis='both', labelsize=14)
    ax1.tick_params(direction='in', axis='both', which='both')
    ax1.tick_params(top=False, right=True, which='both') # top='on', right='on' -- deprecated -- use True/False instead
    ax1.grid(True, ls='--', lw=0.25)
    if set_x_opz:
        # y range
        ax1.set_ylim([-3.0, 0.5])
        # y minor ticks
        #ax1.set_yscale('log')
        #ax1.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=99))
        #ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=99))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(base=0.5))
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.1))
        # bottom redshift axis
        ax1.set_xlim([1.0-0.2, 1.0+10.0])
        ax1.set_xscale('log')
        lower_tick_locations = 1.0 + np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]) # 1+z
        lower_tick_locations = lower_tick_locations[np.argwhere(np.logical_and(lower_tick_locations >= ax1.get_xlim()[0], lower_tick_locations <= ax1.get_xlim()[1])).flatten()]
        ax1.set_xticks(lower_tick_locations)
        ax1.set_xticklabels(lower_ticks_for_CSFRD_vs_opz(1.0+lower_tick_locations))
        # top redshift axis
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        upper_tick_locations = 1.0 + convert_age_to_z([cosmo.age(0).value, 2.0, 1.0, 0.5, 0.3]) # age -> 1+z
        upper_tick_locations = upper_tick_locations[np.argwhere(np.logical_and(upper_tick_locations >= ax2.get_xlim()[0], upper_tick_locations <= ax2.get_xlim()[1])).flatten()]
        print('upper_tick_locations', upper_tick_locations, '(opz)')
        ax2.set_xticks(upper_tick_locations)
        ax2.set_xticklabels(upper_ticks_for_CSFRD_vs_opz(1.0+upper_tick_locations))
    else:
        # y range
        ax1.set_ylim([-2.4, -0.4])
        # y minor ticks
        #ax1.set_yscale('log')
        #ax1.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=99))
        #ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=99))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(base=0.4))
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.05))
        # bottom redshift axis
        ax1.set_xlim([-0.2, 5.5])
        lower_tick_locations = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]) # z
        lower_tick_locations = lower_tick_locations[np.argwhere(np.logical_and(lower_tick_locations >= ax1.get_xlim()[0], lower_tick_locations <= ax1.get_xlim()[1])).flatten()]
        ax1.set_xticks(lower_tick_locations)
        ax1.set_xticklabels(lower_ticks_for_CSFRD_vs_z(lower_tick_locations))
        # top redshift axis
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        upper_tick_locations = convert_age_to_z([cosmo.age(0).value, 5.0, 3.0, 2.0, 1.0, 0.5, ]) # age -> z
        upper_tick_locations = upper_tick_locations[np.argwhere(np.logical_and(upper_tick_locations >= ax2.get_xlim()[0], upper_tick_locations <= ax2.get_xlim()[1])).flatten()]
        #upper_tick_locations = convert_age_to_z([2.0])
        print('upper_tick_locations', upper_tick_locations, '(z)')
        ax2.set_xticks(upper_tick_locations)
        ax2.set_xticklabels(upper_ticks_for_CSFRD_vs_z(upper_tick_locations))
    ax2.set_xlabel(r"Cosmic Age [$\mathrm{Gyr}$]", fontsize=16, labelpad=6)
    ax2.tick_params(axis='x', which='minor', top=False)
    ax2.grid(None)
    # 
    # 
    # plot curve
    if set_x_opz:
        opz_list = np.logspace(np.log10(1.0+0.0), np.log10(1.0+10.0), num=100, endpoint=True)
        z_list = opz_list - 1.0
        CSFRD_list = calc_CSFRD(z_list)
        ax1.plot(opz_list, np.log10(CSFRD_list), label=label_CSFRD + ' (Chabrier IMF)')
    else:
        z_list = np.linspace(0.0, 10.0, num=100, endpoint=True)
        CSFRD_list = calc_CSFRD(z_list)
        ax1.plot(z_list, np.log10(CSFRD_list), label=label_CSFRD + ' (Chabrier IMF)')
    # 
    # plot data point if it exists
    datatable_SMF_SFRD_file = 'datatable_integrating_SMF_SFRD_for_CSFRD_evolution' + '_with_' + re.sub(r'[^a-zA-Z0-9]', r'', label_SF_MS) + '_MS.txt'
    if os.path.isfile(datatable_SMF_SFRD_file):
        tb = Table.read(datatable_SMF_SFRD_file, format='ascii')
        if set_x_opz:
            ax1.plot(1.0+tb['z'], np.log10(tb['PhiSFR']), ls='dashed', label=r'$\int \, \mathrm{SMF_{\sigma}} \times \mathrm{SFR_{MS,%s}}$'%(label_SF_MS.replace('+','')))
        else:
            ax1.plot(tb['z'], np.log10(tb['PhiSFR']), ls='dashed', label=r'$\int \, \mathrm{SMF_{\sigma}} \times \mathrm{SFR_{MS,%s}}$'%(label_SF_MS.replace('+','')))
    # 
    ax1.legend(loc='upper left', framealpha=0.5, borderaxespad=0.6)
    # 
    # savefig
    fig.savefig(output_name+'.pdf', transparent=True)
    print('Output to "%s"' % (output_name+'.pdf'))
    os.system('open "%s"' % (output_name+'.pdf'))
    # 
    shutil.copy(output_name+'.pdf', output_name+'_with_'+re.sub(r'[^a-zA-Z0-9]', r'', label_SF_MS)+'_MS.pdf')






"""
Integrate cosmic SFR density (CSFRD) since redshift +inf to the redshift under examination, 
This should give a consistent cosmic stellar mass density with the integration of stellar mass function over all stellar mass bins at that redshift. 
"""
def integrate_CSFRD_for_rho_Mstar():
    # 
    global calc_CSFRD
    global label_CSFRD
    #global calc_SMF
    #global label_SMF
    # 
    lgMstar_min = 9.0
    lgMstar_max = 13.0
    z_min = 0.0
    z_max = 10.0
    # 
    output_name = 'datatable_integrating_CSFRD_for_rho_Mstar'
    # 
    opz_list = np.logspace(np.log10(1.0+z_min), np.log10(1.0+z_max), num=50, endpoint=True)
    # 
    z_list = opz_list - 1.0
    t_cosmic_age_list = cosmo.age(z_list).value # Gyr
    t_cosmic_time_interval_list = np.abs(t_cosmic_age_list[1:] - t_cosmic_age_list[0:-1]) # Gyr
    CSFRD_list = calc_CSFRD(z_list)
    CSFRD_moving_average_list = (CSFRD_list[1:] + CSFRD_list[0:-1]) / 2.0 # Msun/yr
    #z_moving_average_list = (z_list[1:] + z_list[0:-1]) / 2.0
    CSFRD_rho_Mstar_list = z_list * 0.0 # integrating CSFRD since z_max to the redshift under examination.
    CSFRD_delta_Mstar_list = z_list * 0.0 # integrating CSFRD in each redshift bin at the redshift under examination.
    SMF_rho_Mstar_list = z_list * 0.0 # integrating SMF at the redshift under examination.
    # 
    # loop redshift bins
    for i in range(len(z_list)-2, -1, -1):
        z = z_list[i]
        t_cosmic_age = t_cosmic_age_list[i]
        t_cosmic_time_interval = t_cosmic_time_interval_list[i]
        CSFRD_at_t_cosmic_time = CSFRD_moving_average_list[i]
        print('z = %s, CSFRD = %s'%(z, CSFRD_at_t_cosmic_time))
        # 
        # integrating CSFRD (see "calc_galaxy_stellar_mass_function.py")
        mass_loss_time_scale = 0.3 # Myr, Conroy & Wechsler (2009, bibcode 2009ApJ...696..620C) arxiv PDF page 5 Eq (11). 
        Mstar_loss_frac = 0.05 * np.log(1.0+(t_cosmic_age)/(mass_loss_time_scale*1e-3)) 
                                # see Ilbert et al. (2013) PDF page 11 left middle; Conroy & Wechsler (2009) arxiv PDF page 5 Eq (11). 
                                # see https://arxiv.org/pdf/1404.5299.pdf PDF page 3 Eq (6); Conroy & Wechsler (2009) arxiv PDF page 5 Eq (11). 
        #Mstar_loss_frac = 0.0 # no mass loss at all <TODO>
        CSFRD_delta_Mstar_list[i] = CSFRD_at_t_cosmic_time * t_cosmic_time_interval*1e9 * (1.0 - Mstar_loss_frac) # newly formed star in stellar mass
        # 
        # integrating SMF
        lgMstar = np.linspace(lgMstar_min, lgMstar_max, num=100, endpoint=True)
        Mstar = 10**lgMstar
        lgPhiMstar = calc_SMF_dzliu2018(z, lgMstar, galaxy_type='ALL', verbose=False) # Phi is per dex and already included ln(10), see "calc_galaxy_stellar_mass_function.py"
        SMF_rho_Mstar_list[i] = np.sum(10**lgPhiMstar * (lgMstar[1]-lgMstar[0]) * 10**(lgMstar) )
    # 
    # 
    CSFRD_rho_Mstar_list = np.cumsum(CSFRD_delta_Mstar_list[::-1])[::-1]
    # 
    # output
    output_dict = {}
    output_dict['z'] = z_list
    output_dict['CSFRD_delta_Mstar'] = CSFRD_delta_Mstar_list
    output_dict['CSFRD_rho_Mstar'] = CSFRD_rho_Mstar_list
    output_dict['SMF_rho_Mstar'] = SMF_rho_Mstar_list
    tbout = Table(output_dict)
    tbout['z'].format = '%.6f'
    tbout['CSFRD_delta_Mstar'].format = '%e'
    tbout['CSFRD_rho_Mstar'].format = '%e'
    tbout['SMF_rho_Mstar'].format = '%e'
    tbout.write(output_name+'.txt', format='ascii.fixed_width', delimiter=' ', bookend=True, overwrite=True)
    with open(output_name+'.txt', 'r+') as fp:
        fp.seek(0); fp.write('#')
    print('Output to "%s"!' % (output_name+'.txt') )
    #os.system('open "%s"' % (output_name+'.txt') )






"""
Plot cosmic stellar mass density evolution.
"""
def plot_rho_Mstar():
    # 
    global calc_SFR_MS
    global label_SF_MS
    # 
    tb1 = Table.read('datatable_integrating_CSFRD_for_rho_Mstar.txt', format='ascii')
    tb2 = Table.read('datatable_MD14_rho_Mstar_Salpeter_IMF.txt', format='ascii')
    tb2['lg_rho_Mstar'] = tb2['lg_rho_Mstar'] - np.log10(1.64) # see MD14 ARAA Sect. 5.3, PDF page 50.
    # 
    output_name = 'Plot_CSMD'
    # 
    set_x_opz = False
    # 
    fig = plt.figure(figsize=(6.8,4.8)) # before 20190722: figsize=(6.8,5.2) # before 20190512: figsize=(7.2,5.8)
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.105, top=0.885)
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_xlabel('Redshift', fontsize=16, labelpad=1)
    ax1.set_ylabel(r'$\log_{10} \, \rho_{M_{\star}}$ [$\mathrm{M_{\odot}\,Mpc^{-3}}$]', fontsize=17, labelpad=15)
    ax1.tick_params(axis='both', labelsize=14)
    ax1.tick_params(direction='in', axis='both', which='both')
    ax1.tick_params(top=False, right=True, which='both') # top='on', right='on' -- deprecated -- use True/False instead
    ax1.grid(True, ls='--', lw=0.25)
    if set_x_opz:
        # y range
        ax1.set_ylim([6.0, 9.0])
        # y minor ticks
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(base=0.5))
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.1))
        # bottom redshift axis
        ax1.set_xlim([1.0-0.2, 1.0+10.0])
        ax1.set_xscale('log')
        lower_tick_locations = 1.0 + np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]) # 1+z
        lower_tick_locations = lower_tick_locations[np.argwhere(np.logical_and(lower_tick_locations >= ax1.get_xlim()[0], lower_tick_locations <= ax1.get_xlim()[1])).flatten()]
        ax1.set_xticks(lower_tick_locations)
        ax1.set_xticklabels(lower_ticks_for_CSFRD_vs_opz(1.0+lower_tick_locations))
        # top redshift axis
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        upper_tick_locations = 1.0 + convert_age_to_z([cosmo.age(0).value, 2.0, 1.0, 0.5, 0.3]) # age -> 1+z
        upper_tick_locations = upper_tick_locations[np.argwhere(np.logical_and(upper_tick_locations >= ax2.get_xlim()[0], upper_tick_locations <= ax2.get_xlim()[1])).flatten()]
        print('upper_tick_locations', upper_tick_locations, '(opz)')
        ax2.set_xticks(upper_tick_locations)
        ax2.set_xticklabels(upper_ticks_for_CSFRD_vs_opz(1.0+upper_tick_locations))
    else:
        # y range
        ax1.set_ylim([6.0, 9.0])
        # y minor ticks
        #ax1.yaxis.set_major_locator(ticker.MultipleLocator(base=0.2))
        #ax1.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.01))
        # bottom redshift axis
        ax1.set_xlim([-0.2, 5.5])
        lower_tick_locations = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]) # z
        lower_tick_locations = lower_tick_locations[np.argwhere(np.logical_and(lower_tick_locations >= ax1.get_xlim()[0], lower_tick_locations <= ax1.get_xlim()[1])).flatten()]
        ax1.set_xticks(lower_tick_locations)
        ax1.set_xticklabels(lower_ticks_for_CSFRD_vs_z(lower_tick_locations))
        # top redshift axis
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        upper_tick_locations = convert_age_to_z([cosmo.age(0).value, 5.0, 3.0, 2.0, 1.0, 0.5, ]) # age -> z
        upper_tick_locations = upper_tick_locations[np.argwhere(np.logical_and(upper_tick_locations >= ax2.get_xlim()[0], upper_tick_locations <= ax2.get_xlim()[1])).flatten()]
        print('upper_tick_locations', upper_tick_locations, '(z)')
        ax2.set_xticks(upper_tick_locations)
        ax2.set_xticklabels(upper_ticks_for_CSFRD_vs_z(upper_tick_locations))
    ax2.set_xlabel(r"Cosmic Age [$\mathrm{Gyr}$]", fontsize=16, labelpad=6)
    ax2.tick_params(axis='x', which='minor', top=False)
    ax2.grid(None)
    # 
    # 
    # plot data
    if set_x_opz:
        ax1.plot(1.0+tb1['z'], np.log10(tb1['CSFRD_rho_Mstar']), label=r'$\int \, \mathrm{CSFRD} \, \mathrm{d} t$')
        ax1.plot(1.0+tb1['z'], np.log10(tb1['SMF_rho_Mstar']), ls='dashed', label=r'$\int \, \mathrm{SMF} \, \mathrm{d} M_{\star}$')
    else:
        ax1.plot(tb1['z'], np.log10(tb1['CSFRD_rho_Mstar']), label=r'$\int \, \mathrm{CSFRD} \, \mathrm{d} t$')
        ax1.plot(tb1['z'], np.log10(tb1['SMF_rho_Mstar']), ls='dashed', label=r'$\int \, \mathrm{SMF} \, \mathrm{d} M_{\star}$')
        ax1.errorbar((tb2['zLo']+tb2['zHi'])/2.0, tb2['lg_rho_Mstar'], yerr=(-tb2['eLo'], tb2['eHi']), ls='none', marker='o', markersize=5, alpha=0.5, lw=0.75, label='MD14 collection (Chabrier IMF)')
    # 
    ax1.legend(loc='lower left', fontsize=11, framealpha=0.3)
    # 
    # savefig
    fig.savefig(output_name+'.pdf', transparent=True)
    print('Output to "%s"' % (output_name+'.pdf'))
    os.system('open "%s"' % (output_name+'.pdf'))






"""
Integrate stellar mass function weighted main sequence SFR over all stellar mass bins at each redshift under examination. 
This should give us a consistent cosmic SFR density (CSFRD) as LF-integrated CSFRD at each redshift. 
"""
def integrate_SMF_SFRD_for_CSFRD_evolution():
    # 
    global calc_SFR_MS
    global label_SF_MS
    #global calc_SMF
    #global label_SMF
    # 
    #MS_SFR_center = 0.87 # 1.0 # 0.87 #<TODO># see Bethermin+2017 Sect. 2.6, The distribution of main-sequence galaxies is centered on 0.87 SFR_MS and 5.3 SFR_MS for the starbursts (Schreiber et al. 2015).
    #MS_SFR_meanlg_to_lgmean_correction = 10**(-0.5*(SF_MS_scatter*np.log(10))**2) # see "check_lognormal_mean_scatter.py", see also Padoan & Nordlund 2002; Leroy et al. 2017 (ApJ 835:217)
    #print('MS_SFR_meanlg_to_lgmean_correction = %s'%(MS_SFR_meanlg_to_lgmean_correction))
    # 
    lgMstar_min = 9.0
    lgMstar_max = 13.0
    # 
    output_name = 'datatable_integrating_SMF_SFRD_for_CSFRD_evolution' + '_with_' + re.sub(r'[^a-zA-Z0-9]', r'', label_SF_MS) + '_MS'
    # 
    opz_list = np.logspace(np.log10(1.0+0.0), np.log10(1.0+10.0), num=25, endpoint=True) # z = 1 - 10
    z_list = opz_list - 1.0
    PhiSFR_list = z_list * 0.0
    PhiSFR_list2 = z_list * 0.0
    # 
    # loop redshift bins
    for i in range(len(z_list)):
        z = z_list[i]
        print('z = %s'%(z))
        lgMstar = np.linspace(lgMstar_min, lgMstar_max, num=100, endpoint=True)
        Mstar = 10**lgMstar
        lgPhiMstar = calc_SMF_dzliu2018(z, lgMstar, verbose=False)
        SFR_MS = calc_SFR_MS(z, lgMstar)
        SSFR_MS = SFR_MS / Mstar
        MS_SFR_meanlg_to_lgmean_correction = np.exp(-0.5*(SF_MS_scatter(lgMstar)*np.log(10))**2) # see "check_lognormal_mean_scatter.py", see also Padoan & Nordlund 2002; Leroy et al. 2017 (ApJ 835:217)
        print('SF_MS_scatter(lgMstar)', SF_MS_scatter(lgMstar[0]), SF_MS_scatter(lgMstar[-1]))
        print('MS_SFR_meanlg_to_lgmean_correction', MS_SFR_meanlg_to_lgmean_correction[0], MS_SFR_meanlg_to_lgmean_correction[-1])
        PhiSFR_list[i] = np.sum(10**lgPhiMstar * SFR_MS * MS_SFR_meanlg_to_lgmean_correction * (lgMstar[1]-lgMstar[0]) ) # Phi is per dex and already included ln(10), see "calc_galaxy_stellar_mass_function.py"
    # 
    # output
    tbout = Table( {'z': z_list, 'PhiSFR': PhiSFR_list} )
    #tbout = Table( {'z': z_list, 'PhiSFR': PhiSFR_list, 'PhiSFR2': PhiSFR_list2} )
    tbout.write(output_name+'.txt', format='ascii.fixed_width', delimiter=' ', bookend=True, overwrite=True)
    with open(output_name+'.txt', 'r+') as fp:
        fp.seek(0); fp.write('#')
    print('Output to "%s"!' % (output_name+'.txt') )
    #os.system('open "%s"' % (output_name+'.txt') )






"""
Plot stellar mass function weighted cosmic SFR per stellar mass bin distribution
"""
def plot_SMF_SFRD():
    # 
    global calc_SFR_MS
    global label_SF_MS
    # 
    output_name = 'Plot_SMF_SFRD'
    z_list = [0.05, 0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 9.5]
    panel_nx = 5
    panel_ny = 2
    fig = plt.figure(figsize=(15.0,5.6))
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.12, top=0.97, wspace=0.02, hspace=0.12)
    for i in range(len(z_list)):
        z = z_list[i]
        lgMstar = np.linspace(7.5, 13.0, num=100, endpoint=True)
        # 
        # Use SMF-1 and multiply SF MS SFR
        try:
            lgPhiMstar1 = calc_SMF_Davidzon2017(z, lgMstar)
            lgPhiSFR1 = np.log10(10**lgPhiMstar1 * calc_SFR_MS(z, lgMstar))
        except ValueError as error:
            print(error)
            lgPhiMstar1 = lgMstar * 0 - 99
            lgPhiSFR1 = lgMstar * 0 - 99
        # 
        # Use SMF-2 and multiply SF MS SFR
        try:
            lgPhiMstar2 = calc_SMF_Ilbert2013(z, lgMstar)
            lgPhiSFR2 = np.log10(10**lgPhiMstar2 * calc_SFR_MS(z, lgMstar))
        except ValueError as error:
            print(error)
            lgPhiMstar2 = lgMstar * 0 - 99
            lgPhiSFR2 = lgMstar * 0 - 99
        # 
        # Use SMF-9 and multiply SF MS SFR
        lgPhiMstar9 = calc_SMF_dzliu2018(z, lgMstar, verbose=False)
        lgPhiSFR9 = np.log10(10**lgPhiMstar9 * calc_SFR_MS(z, lgMstar))
        # 
        # ax
        ax = fig.add_subplot(panel_ny, panel_nx, i+1)
        # 
        ax.plot(lgMstar, lgPhiSFR1, c='k', ls='solid', lw=3.0, alpha=0.8, label=r'Davidzon+2017 SFG SMF'+'\n'+r' $\times$ %s MS'%(label_SF_MS))
        ax.plot(lgMstar, lgPhiSFR2, c='green', ls='dashed', lw=2.0, alpha=0.8, label=r'Ilbert+2013 SFG SMF'+'\n'+r' $\times$ %s MS'%(label_SF_MS))
        ax.plot(lgMstar, lgPhiSFR9, c='gold', ls='solid', lw=4.0, alpha=0.8, label=r'This work SFG SMF'+'\n'+r' $\times$ %s MS'%(label_SF_MS))
        # 
        if i < len(z_list)-1:
            ax.set_ylim([-6.8, 1])
            ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
            ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
            ax.text(0.95, 0.95, r'$z=%4g$'%(z), ha='right', va='top', fontsize=16, transform=ax.transAxes)
            
            if i == 0:
                plt.ylabel(r'$\log_{10} \, \Phi_{\mathrm{SFR}} \,/\, [\mathrm{M_{\odot} \, yr^{-1} \, Mpc^{-3} \, dex^{-1}}]$', fontsize=18)
                ax.yaxis.set_label_coords(-0.16, -0.10)
            
            if i % panel_nx != 0:
                ax.get_yaxis().set_ticklabels([])
            
            if i < len(z_list) - panel_nx:
                #ax.get_xaxis().set_ticklabels([])
                pass
            else:
                plt.xlabel(r'$\log_{10} (M_{\star}/\mathrm{M_{\odot}})$', fontsize=17)
        # 
        elif i == len(z_list)-1:
            ax.set_ylim([11, 17])
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            leg = ax.legend(fontsize=12, loc='center')
        # 
        #ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
        #ax.grid(True, ls='dotted', lw=0.25)
        #plt.xticks(fontsize=14)
        #plt.yticks(fontsize=14)
        #if i > 0:
        #    plt.yticks([])
    # 
    # savefig
    fig.savefig(output_name+'.pdf')
    print('Output to "%s"!' % (output_name+'.pdf') )
    os.system('open "%s"' % (output_name+'.pdf') )
    # 
    shutil.copy(output_name+'.pdf', output_name+'_with_'+re.sub(r'[^a-zA-Z0-9]', r'', label_SF_MS)+'_MS.pdf')






"""
Plot stellar mass function weighted cosmic specific-SFR per stellar mass bin distribution
"""
def plot_SMF_SSFRD():
    # 
    global calc_SFR_MS
    global label_SF_MS
    # 
    output_name = 'Plot_SMF_SSFRD'
    z_list = [0.05, 0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 9.5]
    panel_nx = 5
    panel_ny = 2
    fig = plt.figure(figsize=(15.0,5.6))
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.12, top=0.97, wspace=0.02, hspace=0.12)
    for i in range(len(z_list)):
        z = z_list[i]
        lgMstar = np.linspace(7.5, 13.0, num=100, endpoint=True)
        # 
        # Use SMF-1 and multiply SF MS SFR
        try:
            lgPhiMstar1 = calc_SMF_Davidzon2017(z, lgMstar)
            lgPhiSSFR1 = np.log10(10**lgPhiMstar1 * calc_SFR_MS(z, lgMstar)) - lgMstar + 9.0 # Gyr
        except ValueError as error:
            print(error)
            lgPhiMstar1 = lgMstar * 0 - 99
            lgPhiSSFR1 = lgMstar * 0 - 99
        # 
        # Use SMF-9 and multiply SF MS SFR
        lgPhiMstar9 = calc_SMF_dzliu2018(z, lgMstar, verbose=False)
        lgPhiSSFR9 = np.log10(10**lgPhiMstar9 * calc_SFR_MS(z, lgMstar)) - lgMstar + 9.0 # Gyr
        # 
        # ax
        ax = fig.add_subplot(panel_ny, panel_nx, i+1)
        # 
        ax.plot(lgMstar, lgPhiSSFR1, c='k', ls='solid', lw=3.0, alpha=0.8, label=r'Davidzon+2017 SFG SMF'+'\n'+r' $\times$ %s MS'%(label_SF_MS))
        ax.plot(lgMstar, lgPhiSSFR9, c='gold', ls='solid', lw=4.0, alpha=0.8, label=r'This work SFG SMF'+'\n'+r' $\times$ %s MS'%(label_SF_MS))
        # 
        if i < len(z_list)-1:
            ax.set_ylim([-6.8, 1])
            ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
            ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
            ax.text(0.95, 0.95, r'$z=%4g$'%(z), ha='right', va='top', fontsize=16, transform=ax.transAxes)
            
            if i == 0:
                plt.ylabel(r'$\log_{10} \, \Phi_{\mathrm{SSFR}} \,/\, [\mathrm{Gyr^{-1} \, Mpc^{-3} \, dex^{-1}}]$', fontsize=18)
                ax.yaxis.set_label_coords(-0.16, -0.10)
            
            if i % panel_nx != 0:
                ax.get_yaxis().set_ticklabels([])
            
            if i < len(z_list) - panel_nx:
                #ax.get_xaxis().set_ticklabels([])
                pass
            else:
                plt.xlabel(r'$\log_{10} (M_{\star}/\mathrm{M_{\odot}})$', fontsize=17)
        # 
        elif i == len(z_list)-1:
            ax.set_ylim([11, 17])
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            leg = ax.legend(fontsize=12, loc='center')
        # 
        #ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
        #ax.grid(True, ls='dotted', lw=0.25)
        #plt.xticks(fontsize=14)
        #plt.yticks(fontsize=14)
        #if i > 0:
        #    plt.yticks([])
    # 
    # savefig
    fig.savefig(output_name+'.pdf')
    print('Output to "%s"!' % (output_name+'.pdf') )
    os.system('open "%s"' % (output_name+'.pdf') )
    # 
    shutil.copy(output_name+'.pdf', output_name+'_with_'+re.sub(r'[^a-zA-Z0-9]', r'', label_SF_MS)+'_MS.pdf')






"""
Plot stellar mass function, i.e., number density of galaxies per stellar mass bin.
"""
def plot_SMF():
    output_name = 'Plot_SMF'
    z_list = [0.05, 0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 9.5]
    panel_nx = 5
    panel_ny = 2
    fig = plt.figure(figsize=(15.0,5.6))
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.12, top=0.97, wspace=0.02, hspace=0.12)
    for i in range(len(z_list)):
        z = z_list[i]
        lgMstar = np.linspace(7.5, 13.0, num=100, endpoint=True)
        # 
        # Use SMF-1
        try:
            lgPhiMstar1 = calc_SMF_Davidzon2017(z, lgMstar)
        except ValueError as error:
            print(error)
            lgPhiMstar1 = lgMstar * 0 - 99
        # 
        # Use SMF-2
        try:
            lgPhiMstar2 = calc_SMF_Ilbert2013(z, lgMstar)
        except ValueError as error:
            print(error)
            lgPhiMstar2 = lgMstar * 0 - 99
        # 
        # Use SMF-4
        try:
            lgPhiMstar4 = calc_SMF_Peng2010(z, lgMstar)
        except ValueError as error:
            print(error)
            lgPhiMstar4 = lgMstar * 0 - 99
        # 
        # Use SMF-2/3
        #lgPhiMstar2 = calc_SMF_Wright2018_single_component(z, lgMstar)
        #lgPhiMstar3 = calc_SMF_Wright2018_double_component(z, lgMstar)
        # 
        # Use SMF-9
        lgPhiMstar9 = calc_SMF_dzliu2018(z, lgMstar)
        # 
        # ax
        ax = fig.add_subplot(panel_ny, panel_nx, i+1)
        # 
        ax.plot(lgMstar, lgPhiMstar1, c='k', ls='solid', lw=3.0, alpha=0.8, label=r'Davidzon+2017, SFG SMF')
        ax.plot(lgMstar, lgPhiMstar2, c='green', ls='dashed', lw=2.0, alpha=0.8, label=r'Ilbert+2013, SFG SMF')
        ax.plot(lgMstar, lgPhiMstar4, c='k', ls='dashed', lw=3.0, alpha=0.8, label=r'Peng+2010, SFG SMF (z=0)')
        #ax.plot(lgMstar, lgPhiMstar2, c='green', ls=(0,(2,1)), lw=3.0, alpha=0.8, label=r'Wright+2018, single SMF(z)')
        #ax.plot(lgMstar, lgPhiMstar3, c='blue', ls=(0,(2,1)), lw=3.0, alpha=0.8, label=r'Wright+2018, double SMF(z)')
        ax.plot(lgMstar, lgPhiMstar9, c='gold', ls='solid', lw=4.0, alpha=0.8, label=r'This work, SFG SMF(z)')
        # 
        if i < len(z_list)-1:
            ax.set_ylim([-6.8, -1])
            ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
            ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
            ax.text(0.95, 0.95, r'$z=%4g$'%(z), ha='right', va='top', fontsize=16, transform=ax.transAxes)
            
            if i == 0:
                plt.ylabel(r'$\log_{10} \, \Phi_{M_{\star}} \,/\, [\mathrm{Mpc^{-3} \, dex^{-1}}]$', fontsize=18)
                ax.yaxis.set_label_coords(-0.16, -0.10)
            
            if i % panel_nx != 0:
                ax.get_yaxis().set_ticklabels([])
            
            if i < len(z_list) - panel_nx:
                #ax.get_xaxis().set_ticklabels([])
                pass
            else:
                plt.xlabel(r'$\log_{10} (M_{\star}/\mathrm{M_{\odot}})$', fontsize=17)
        # 
        elif i == len(z_list)-1:
            ax.set_ylim([11, 17])
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            leg = ax.legend(fontsize=12, loc='center')
        # 
        #ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
        #ax.grid(True, ls='dotted', lw=0.25)
        #plt.xticks(fontsize=14)
        #plt.yticks(fontsize=14)
        #if i > 0:
        #    plt.yticks([])
    # 
    # savefig
    fig.savefig(output_name+'.pdf')
    print('Output to "%s"!' % (output_name+'.pdf') )
    os.system('open "%s"' % (output_name+'.pdf') )














####################################
#               MAIN               #
####################################

if __name__ == '__main__':
    
    #plot_SMF()
    plot_SMF_SFRD()
    plot_SMF_SSFRD()
    
    #integrate_SMF_SFRD_for_CSFRD_evolution()
    #plot_CSFRD()
    
    #integrate_CSFRD_for_rho_Mstar() # output "datatable_integrating_CSFRD_for_rho_Mstar.txt"
    #plot_rho_Mstar() # output "Plot_CSMD.pdf"


















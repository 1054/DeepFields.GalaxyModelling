#!/usr/bin/env python
# 

from __future__ import print_function
import os, sys, re, json, shutil
import numpy as np
from astropy.table import Table
sys.path.append('/Users/dzliu/Cloud/Github/Crab.Toolkit.Python/lib/crab/crabgalaxy')
sys.path.append('/Users/dzliu/Cloud/Github/Crab.Toolkit.Python/lib/crab/crabtable')
sys.path.append('/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code')
from CrabGalaxy import CrabGalaxy
from CrabTable import CrabTableReadInfo
from calc_galaxy_stellar_mass_function import (calc_SMF_Davidzon2017, calc_SMF_Ilbert2013, calc_SMF_Peng2010, calc_SMF_Wright2018_single_component, calc_SMF_Wright2018_double_component, calc_SMF_dzliu2018)
from calc_galaxy_main_sequence import (calc_SFR_MS_Speagle2014, calc_SFR_MS_Sargent2014, calc_SFR_MS_Schreiber2015, calc_SFR_MS_Leslie20190710)
from calc_cosmic_star_formation_rate_density import (calc_CSFRD_Madau2014, calc_CSFRD_Liu2018, convert_age_to_z)
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from setup_matplotlib import setup_matplotlib; setup_matplotlib()
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=69.8, Om0=0.3, Tcmb0=2.725) # Freedman 2019ApJ...882...34F

#InfoDict = CrabTableReadInfo(InfoFile, verbose=0)

 
calc_SFR_MS = calc_SFR_MS_Leslie20190710 ; label_SF_MS = 'Leslie+2019'
#calc_SFR_MS = calc_SFR_MS_Speagle2014 ; label_SF_MS = 'Speagle+2014'
#calc_SFR_MS = calc_SFR_MS_Sargent2014 ; label_SF_MS = 'Sargent+2014'
#calc_SFR_MS = calc_SFR_MS_Schreiber2015 ; label_SF_MS = 'Schreiber+2015'






def generate_galaxy_number():
    pass


def model_galaxies():
    one_galaxy = CrabGalaxy(shape = 'sersic', size = {'major':1.0, 'minor':0.5, 'PA':20.0, 'n':0.5})







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


def plot_CSFRD():
    # 
    calc_CSFRD = calc_CSFRD_Madau2014
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
        ax1.set_ylim([-2.0, -0.2])
        # y minor ticks
        #ax1.set_yscale('log')
        #ax1.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=99))
        #ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=99))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(base=0.2))
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.01))
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
        ax1.plot(opz_list, np.log10(CSFRD_list), label='MD14')
    else:
        z_list = np.linspace(0.0, 10.0, num=100, endpoint=True)
        CSFRD_list = calc_CSFRD(z_list)
        ax1.plot(z_list, np.log10(CSFRD_list), label='MD14')
    # 
    # plot data point if it exists
    datatable_SMF_SFRD_file = 'datatable_SMF_SFRD' + '_with_' + re.sub(r'[^a-zA-Z0-9]', r'', label_SF_MS) + '_MS.txt'
    if os.path.isfile(datatable_SMF_SFRD_file):
        tb = Table.read(datatable_SMF_SFRD_file, format='ascii')
        if set_x_opz:
            ax1.plot(1.0+tb['z'], np.log10(tb['PhiSFR']), label='SMF SFRD '+label_SF_MS)
        else:
            ax1.plot(tb['z'], np.log10(tb['PhiSFR']), label='SMF SFRD '+label_SF_MS)
    # 
    ax1.legend(loc='lower left', ncol=2)
    # 
    # savefig
    fig.savefig(output_name+'.pdf', transparent=True)
    print('Output to "%s"' % (output_name+'.pdf'))
    os.system('open "%s"' % (output_name+'.pdf'))
    # 
    shutil.copy(output_name+'.pdf', output_name+'_with_'+re.sub(r'[^a-zA-Z0-9]', r'', label_SF_MS)+'_MS.pdf')



def Monte_Carlo_SMF_SFRD():
    pass



def integrate_SMF_SFRD():
    # 
    global calc_SFR_MS
    global label_SF_MS
    # 
    MS_SFR_center = 1.0 # 0.87 #<TODO># see Bethermin+2017 Sect. 2.6, The distribution of main-sequence galaxies is centered on 0.87 SFR_MS and 5.3 SFR_MS for the starbursts (Schreiber et al. 2015).
    # 
    lgMstar_min = 9.0
    lgMstar_max = 13.0
    # 
    output_name = 'datatable_SMF_SFRD' + '_with_' + re.sub(r'[^a-zA-Z0-9]', r'', label_SF_MS) + '_MS'
    # 
    opz_list = np.logspace(np.log10(1.0+0.0), np.log10(1.0+10.0), num=25, endpoint=True)
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
        lgPhiMstar = calc_SMF_dzliu2018(z, lgMstar, verbose=False) # per dex
        SFR_MS = calc_SFR_MS(z, lgMstar)
        SSFR_MS = SFR_MS / Mstar
        PhiSFR_list[i] = np.sum(10**lgPhiMstar * SFR_MS * MS_SFR_center * (lgMstar[1]-lgMstar[0]) ) # Phi is already dex^{-1}, so directly integrate with d lg(M)
        #PhiSFR_list[i] = np.sum(10**lgPhiMstar * SFR_MS * (lgMstar[1]-lgMstar[0]) * np.log(10) ) # \int Phi d SFR = \int Phi SSFR d M = \int Phi SSFR M d lg(M) ln(10)
        # 
        # we need to loop stellar mass bins because simple integration does not work!
        #Mstar = np.linspace(10**lgMstar_min, 10**lgMstar_max, num=10000, endpoint=True)
        #lgMstar = np.log10(Mstar)
        #lgPhiMstar = calc_SMF_dzliu2018(z, lgMstar, verbose=False) # per dex
        #SFR_MS = calc_SFR_MS(z, lgMstar)
        #SSFR_MS = SFR_MS / Mstar
        #PhiSFR_list2[i] = np.sum(10**lgPhiMstar * SFR_MS / (Mstar[1]-Mstar[0]) * (Mstar * (10**(0.5) - 10**(-0.5))) ) # \int Phi d SFR = \int Phi SSFR d M ( = \int Phi SSFR M d lg(M) ln(10) )
    # 
    # output
    tbout = Table( {'z': z_list, 'PhiSFR': PhiSFR_list} )
    #tbout = Table( {'z': z_list, 'PhiSFR': PhiSFR_list, 'PhiSFR2': PhiSFR_list2} )
    tbout.write(output_name+'.txt', format='ascii.fixed_width', delimiter=' ', bookend=True, overwrite=True)
    with open(output_name+'.txt', 'r+') as fp:
        fp.seek(0); fp.write('#')
    print('Output to "%s"!' % (output_name+'.txt') )
    #os.system('open "%s"' % (output_name+'.txt') )





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








def check_integral_m_dm():
    lg_m = np.linspace(6.0, 13.0, num=10000, endpoint=True)
    f_m = lg_m * 0.0 + 1.0 # assuming flat number density distribution
    sum1 = np.sum(f_m * (10**lg_m) * (lg_m[1]-lg_m[0]) ) # \int f d ln(M) = \int f d M / ln(10)
    print('%e'%(sum1))
    # 
    m = np.linspace(10**6.0, 10**13.0, num=10000, endpoint=True)
    f_m = m * 0.0 + 1.0 # assuming flat number density distribution
    sum2 = np.sum(f_m * (m[1]-m[0]) ) # \int f * d M
    print('%e'%(sum2), sum2/sum1, np.log(10))


def check_integral_SFR_dSFR():
    z = 3.0
    calc_SFR_MS = calc_SFR_MS_Sargent2014
    # 
    lg_m = np.linspace(8.0, 13.0, num=10000, endpoint=True)
    m = 10**lg_m
    lg_phi_m = lg_m * (-0.0) + 1.0 # assuming flat number density distribution
    #lg_phi_m = calc_SMF_dzliu2018(z, lg_m, verbose=False) # per dex
    SFR_MS = calc_SFR_MS(z, lg_m)
    SSFR_MS = SFR_MS / m
    sum1 = np.sum(10**lg_phi_m * SSFR_MS * m * (lg_m[1]-lg_m[0]) ) # \int d SFR = \int SSFR(M) d M = \int SSFR(M) M d lg(M) ln(10)
    print('%e'%(sum1))
    # 
    m = np.linspace(10**8.0, 10**13.0, num=10000, endpoint=True)
    lg_m = np.log10(m)
    lg_phi_m = lg_m * (-0.0) + 1.0 # assuming flat number density distribution
    #lg_phi_m = calc_SMF_dzliu2018(z, lg_m, verbose=False) # per dex
    SFR_MS = calc_SFR_MS(z, lg_m)
    SSFR_MS = SFR_MS / m
    sum2 = np.sum(10**lg_phi_m * SSFR_MS * (m[1]-m[0]) ) # \int d SFR = \int SSFR(M) d M
    print('%e'%(sum2), sum2/sum1, np.log(10))
    # 
    print(m[2] * (lg_m[1]-lg_m[0]), (m[1]-m[0]))
    print(m[20] * (lg_m[1]-lg_m[0]), (m[1]-m[0]))


def check_integral_SFR_dSFR_exactly_the_same():
    if True:
        i = 0
        PhiSFR_list = [0.0]
        PhiSFR_list2 = [0.0]
        lgMstar_min = 8.5
        lgMstar_max = 13
        calc_SFR_MS = calc_SFR_MS_Sargent2014
        z = 3.0
        print('z = %s'%(z))
        lgMstar = np.linspace(lgMstar_min, lgMstar_max, num=10000, endpoint=True)
        Mstar = 10**lgMstar
        lgPhiMstar = calc_SMF_dzliu2018(z, lgMstar, verbose=False) # per dex
        SFR_MS = calc_SFR_MS(z, lgMstar)
        SSFR_MS = SFR_MS / Mstar
        PhiSFR_list[i] = np.sum(10**lgPhiMstar * SFR_MS * (lgMstar[1]-lgMstar[0]) * np.log(10) ) # \int Phi d SFR = \int Phi SSFR d M = \int Phi SSFR M d lg(M) ln(10)
        # 
        # we need to loop stellar mass bins because simple integration does not work!
        Mstar = np.linspace(10**lgMstar_min, 10**lgMstar_max, num=10000, endpoint=True)
        lgMstar = np.log10(Mstar)
        lgPhiMstar = calc_SMF_dzliu2018(z, lgMstar, verbose=False) # per dex
        SFR_MS = calc_SFR_MS(z, lgMstar)
        SSFR_MS = SFR_MS / Mstar
        PhiSFR_list2[i] = np.sum(10**lgPhiMstar * SFR_MS / (Mstar[1]-Mstar[0]) * (Mstar * (10**(0.5) - 10**(-0.5))) ) # \int Phi d SFR = \int Phi SSFR d M ( = \int Phi SSFR M d lg(M) ln(10) )
        # 
        # 
        print(Mstar[2] * (lgMstar[1]-lgMstar[0]), (Mstar[1]-Mstar[0]))
        print(Mstar[20] * (lgMstar[1]-lgMstar[0]), (Mstar[1]-Mstar[0]))
        # 
        # 
        print('%e'%(PhiSFR_list[i]))
        print('%e'%(PhiSFR_list2[i]), PhiSFR_list2[i]/PhiSFR_list[i], np.log(10))





####################################
#               MAIN               #
####################################

if __name__ == '__main__':
    #plot_SMF()
    #plot_SMF_SFRD()
    #plot_SMF_SSFRD()
    integrate_SMF_SFRD()
    plot_CSFRD()
    #check_integral_m_dm()
    #check_integral_SFR_dSFR()
    #check_integral_SFR_dSFR_exactly_the_same()









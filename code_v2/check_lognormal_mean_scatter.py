#!/usr/bin/env python
# 
from __future__ import print_function
import os, sys, re, json, shutil
import numpy as np
from astropy.table import Table
script_dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(os.path.join(script_dir, 'python', 'lib', 'crab', 'crabgalaxy'))
sys.path.append(os.getenv('HOME')+'/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code') #<TODO># replace by a3g
#from CrabGalaxy import CrabGalaxy
#from calc_galaxy_stellar_mass_function import (calc_SMF_Davidzon2017, calc_SMF_Ilbert2013, calc_SMF_Peng2010, calc_SMF_Wright2018_single_component, calc_SMF_Wright2018_double_component, calc_SMF_dzliu2018)
from calc_galaxy_main_sequence import (calc_SFR_MS_Speagle2014, calc_SFR_MS_Sargent2014, calc_SFR_MS_Schreiber2015, calc_SFR_MS_Leslie20190710)
#from calc_cosmic_star_formation_rate_density import (calc_CSFRD_Madau2014, calc_CSFRD_Liu2018, convert_age_to_z)
#from matplotlib import pyplot as plt
#from matplotlib import ticker as ticker
#from setup_matplotlib import setup_matplotlib; setup_matplotlib()


global calc_SFR_MS
#calc_SFR_MS = calc_SFR_MS_Leslie20190710 ; label_SF_MS = 'Leslie+2019'
calc_SFR_MS = calc_SFR_MS_Speagle2014 ; label_SF_MS = 'Speagle+2014'
#calc_SFR_MS = calc_SFR_MS_Sargent2014 ; label_SF_MS = 'Sargent+2014'
#calc_SFR_MS = calc_SFR_MS_Schreiber2015 ; label_SF_MS = 'Schreiber+2015'




"""
Function for SF_MS_scatter
"""
def SF_MS_scatter(lgMstar):
    # See https://arxiv.org/pdf/1812.07057.pdf Figure 5
    # It seems their SED fitting derived main sequence scatter is most stable at around 0.29 dex.
    #return np.interp(lgMstar, [6.0, 10.0, 11.0, 14.0], [0.3, 0.3, 0.4, 0.4])
    return np.interp(lgMstar, [6.0, 10.0, 11.0, 14.0], [0.29, 0.29, 0.29, 0.29])
    #return 0.20






def check_lognormal_mean_scatter():
    global calc_SFR_MS
    galaxy_z = 3.0
    galaxy_lgMstar = 10.7
    main_sequence_lgSFR = np.log10(calc_SFR_MS(galaxy_z, galaxy_lgMstar))  # this should be mean(lg(SFR)), 
    N = 3000
    galaxy_lgSFR = np.random.normal(main_sequence_lgSFR, SF_MS_scatter(galaxy_lgMstar), N)
    MS_SFR_meanlg_to_lgmean_correction = np.exp(-0.5*(SF_MS_scatter(galaxy_lgMstar)*np.log(10))**2) # this is lg(mean(x)) - mean(lg(x)), see also Padoan & Nordlund 2002; Leroy et al. 2017 (ApJ 835:217)
    lgmean = np.log10(np.mean(10**galaxy_lgSFR))
    meanlg = main_sequence_lgSFR
    print('sum of SFR:', np.sum(10**galaxy_lgSFR))
    print('mean(lg(SFR))*N:', 10**main_sequence_lgSFR * N)
    print('lg(mean(SFR)):', lgmean, 'mean(lg(SFR)):', meanlg, 'ratio(meanlg/lgmean):', 10**(meanlg-lgmean))
    print('MS_SFR_meanlg_to_lgmean_correction', MS_SFR_meanlg_to_lgmean_correction)
    # 
    # sum of SFR: 476583.56240495446
    # mean(lg(SFR))*N: 381591.52090311813
    # lg(mean(SFR)): 2.201017804524656 mean(lg(SFR)): 2.1044774612445245 ratio(meanlg/lgmean): 0.8006812466999833
    # MS_SFR_meanlg_to_lgmean_correction 0.8001590044142833
    # 









####################################
#               MAIN               #
####################################

if __name__ == '__main__':
    check_lognormal_mean_scatter()










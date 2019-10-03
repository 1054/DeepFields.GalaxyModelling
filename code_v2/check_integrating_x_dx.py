#!/usr/bin/env python
# 
from __future__ import print_function
import os, sys, re, json, shutil
import numpy as np
from astropy.table import Table
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(script_dir, 'python', 'lib', 'crab', 'crabgalaxy'))
sys.path.append('/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code') #<TODO># replace by a3g
from CrabGalaxy import CrabGalaxy
from calc_galaxy_stellar_mass_function import (calc_SMF_Davidzon2017, calc_SMF_Ilbert2013, calc_SMF_Peng2010, calc_SMF_Wright2018_single_component, calc_SMF_Wright2018_double_component, calc_SMF_dzliu2018)
from calc_galaxy_main_sequence import (calc_SFR_MS_Speagle2014, calc_SFR_MS_Sargent2014, calc_SFR_MS_Schreiber2015, calc_SFR_MS_Leslie20190710)
from calc_cosmic_star_formation_rate_density import (calc_CSFRD_Madau2014, calc_CSFRD_Liu2018, convert_age_to_z)
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from setup_matplotlib import setup_matplotlib; setup_matplotlib()
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=69.8, Om0=0.3, Tcmb0=2.725) # Freedman 2019ApJ...882...34F




def check_integral_m_dm():
    lg_m = np.linspace(6.0, 13.0, num=10000, endpoint=True) # actually logspace
    f_m = lg_m * 0.0 + 1.0 # assuming flat number density distribution
    sum1 = np.sum(f_m * (10**lg_m) * (lg_m[1]-lg_m[0]) ) # \int f(lg(M)) M d lg(M) = \int f(M) d M / ln(10)
    print('%e'%(sum1))
    # 
    m = np.linspace(10**6.0, 10**13.0, num=10000, endpoint=True)
    f_m = m * 0.0 + 1.0 # assuming flat number density distribution
    sum2 = np.sum(f_m * (m[1]-m[0]) ) # \int f(M) * d M
    print('%e'%(sum2), 'sum2/sum1', sum2/sum1, 'np.log(10)', np.log(10))
    # 
    lg_m = np.linspace(6.0, 13.0, num=10000, endpoint=True) # actually logspace
    f_m = lg_m * 0.0 + 1.0 # assuming flat number density distribution
    sum3 = np.sum((f_m[1:]+f_m[0:-1])/2.0 * (10**lg_m[1:]-10**lg_m[0:-1]) ) # \int f(lg(M)) M d lg(M) = \int f(M) d M / ln(10)
    print('%e'%(sum3)) # very close to sum2


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
    print('%e'%(sum2), 'sum2/sum1', sum2/sum1, 'np.log(10)', np.log(10))
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
    check_integral_m_dm()
    #check_integral_SFR_dSFR()
    #check_integral_SFR_dSFR_exactly_the_same()










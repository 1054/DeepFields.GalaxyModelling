#!/usr/bin/env python
# 
from __future__ import print_function
import os, sys, re, json, shutil, random, timeit
from tqdm import tqdm, trange
import numpy as np
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord, FK5
from astropy.stats import sigma_clip
#from astropy.modeling.functional_models import Gaussian2D
from astropy.convolution import Gaussian2DKernel, convolve
from scipy.spatial import KDTree, distance
#from scipy.stats import sigmaclip
global script_dir
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(script_dir, 'python', 'lib', 'crab', 'crabgalaxy'))
sys.path.append(os.path.join(script_dir, 'python', 'lib', 'crab', 'crabtable'))
sys.path.append('/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code') #<TODO># replace by a3g
from CrabGalaxy import CrabGalaxy, CrabGalaxyMorphology
from CrabTable import CrabTableReadInfo
from calc_galaxy_stellar_mass_function import (calc_SMF_Davidzon2017, calc_SMF_Ilbert2013, calc_SMF_Peng2010, calc_SMF_Wright2018_single_component, calc_SMF_Wright2018_double_component, calc_SMF_dzliu2018)
from calc_galaxy_main_sequence import (calc_SFR_MS_Speagle2014, calc_SFR_MS_Sargent2014, calc_SFR_MS_Schreiber2015, calc_SFR_MS_Leslie20190710)
from calc_cosmic_star_formation_rate_density import (calc_CSFRD_Madau2014, calc_CSFRD_Liu2018, convert_age_to_z)
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from setup_matplotlib import setup_matplotlib; setup_matplotlib()
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=69.8, Om0=0.3, Tcmb0=2.725) # Freedman 2019ApJ...882...34F
ln = np.log

#InfoDict = CrabTableReadInfo(InfoFile, verbose=0)


global calc_SFR_MS
#calc_SFR_MS = calc_SFR_MS_Leslie20190710 ; label_SF_MS = 'Leslie+2019'
calc_SFR_MS = calc_SFR_MS_Speagle2014 ; label_SF_MS = 'Speagle+2014'
#calc_SFR_MS = calc_SFR_MS_Sargent2014 ; label_SF_MS = 'Sargent+2014'
#calc_SFR_MS = calc_SFR_MS_Schreiber2015 ; label_SF_MS = 'Schreiber+2015'


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
    #return 0.20






"""
"""
def make_z_grid():
    #zArray_1 = np.arange(9.75, 0.75-0.5, -0.5)
    #zArray_2 = np.array([0.50, 0.25, 0.125, 0.0625, 0.03125, 0.0])
    opz_list = np.logspace(np.log10(1.0+0.0), np.log10(1.0+10.0), num=31, endpoint=True)
    z_list = opz_list - 1.0
    return z_list





"""
"""
def make_lgMstar_grid():
    #lgMstar_min = 9.0 # 
    lgMstar_min = 8.0 # 
    lgMstar_max = 13.0
    lgMstar_list = np.linspace(lgMstar_min, lgMstar_max, num=41, endpoint=True)
    return lgMstar_list





"""
"""
def calc_comoving_volume(z_lower, z_upper, area_arcmin2:float):
    obs_area = area_arcmin2 * u.arcmin**2 # 1.4*1.4*u.deg*u.deg #<TODO># 
    differntial_z_list = np.linspace(z_lower, z_upper, num=10, endpoint=True)
    comoving_volume = np.sum((cosmo.differential_comoving_volume(differntial_z_list[1:]) * np.diff(differntial_z_list) * obs_area.to(u.steradian)))
    #print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    return comoving_volume.value




"""
"""
def calc_Galaxy_Size_vanderWel2014(z, Mstar, galaxy_type = 'SFG', add_noise = False):
    # compute galaxy size (R_eff)
    # van del Wel et al. (2014) -- SFG -- Reff/kpc = 10**A * (Mstar/5e10)**N, scatter 0.16 - 0.19 dex. 
    if galaxy_type == 'SFG':
        cmorph_z = [0.00,   0.25, 0.75, 1.25, 1.75, 2.25, 2.75,    10.0] #<Note># the first and last items are added by dzliu, others are from Table 1 of van del Wel et al. (2014, doi:10.1088/0004-637X/788/1/28)
        cmorph_A = [0.86,   0.86, 0.78, 0.70, 0.65, 0.55, 0.51,   -0.05] #<Note># the first and last items are added by dzliu, others are from Table 1 of van del Wel et al. (2014, doi:10.1088/0004-637X/788/1/28)
        cmorph_N = [0.25,   0.25, 0.22, 0.22, 0.23, 0.22, 0.22,    0.22] #<Note># the first and last and second last items are added by dzliu, others are from Table 1 of van del Wel et al. (2014, doi:10.1088/0004-637X/788/1/28)
        cmorph_E = [0.16,   0.16, 0.16, 0.17, 0.18, 0.19, 0.19,    0.19] #<Note># the first and last items are added by dzliu, others are from Table 1 of van del Wel et al. (2014, doi:10.1088/0004-637X/788/1/28)
    elif galaxy_type == 'QG':
        cmorph_z = [0.25,   0.25, 0.75, 1.25, 1.75,  2.25,  2.75,    10.0] #<Note># the first and last items are added by dzliu, others are from Table 1 of van del Wel et al. (2014, doi:10.1088/0004-637X/788/1/28)
        cmorph_A = [0.60,   0.60, 0.42, 0.22, 0.09, -0.05, -0.06,   -0.06] #<Note># the first and last items are added by dzliu, others are from Table 1 of van del Wel et al. (2014, doi:10.1088/0004-637X/788/1/28)
        cmorph_N = [0.75,   0.75, 0.71, 0.76, 0.76,  0.76,  0.79,    0.79] #<Note># the first and last and second last items are added by dzliu, others are from Table 1 of van del Wel et al. (2014, doi:10.1088/0004-637X/788/1/28)
        cmorph_E = [0.10,   0.10, 0.11, 0.12, 0.14,  0.14,  0.14,    0.14] #<Note># the first and last items are added by dzliu, others are from Table 1 of van del Wel et al. (2014, doi:10.1088/0004-637X/788/1/28)
    A = np.interp(z, cmorph_z, cmorph_A)
    N = np.interp(z, cmorph_z, cmorph_N)
    E = np.interp(z, cmorph_z, cmorph_E)
    R_eff_kpc = (10**A * (Mstar/5e10)**N) * 2.0
    # 
    # add some noise 
    add_noise = False
    if add_noise:
        if np.isscalar(E):
            noise_sigma = E
        else:
            noise_sigma = np.mean(E)
        # add normally distributed noise
        R_eff_kpc = R_eff_kpc * (np.random.lognormal(0.0, noise_sigma, len(R_eff_kpc)))
    # 
    # return
    return R_eff_kpc
    




"""
"""
def generate_galaxy_numbers(area_arcmin2):
    #raise NotImplementedError()
    # 
    # check existing output file
    output_name = 'datatable_generated_galaxies'
    if os.path.isfile('%s.fits'%(output_name)):
        print('Found existing "%s.fits"! Will not overwrite it!'%(output_name))
        return
    # 
    # generate z, lgMstar grid
    z = make_z_grid()
    lgopz = np.log10(1.0+z)
    t_cosmic_age = cosmo.age(z).value
    lgMstar = make_lgMstar_grid()
    # 
    # calculate stellar mass function (only star-forming galaxies)
    lgPhiMstar = calc_SMF_dzliu2018((z[0:-1]+z[1:])/2.0, lgMstar, galaxy_type='SFG', tuning_params='D17-no-renorm;', verbose=False)
    # 
    # calculate main sequence scatter lg(mean(x)) mean(lg(x)) correction
    #MS_SFR_meanlg_to_lgmean_correction = np.exp(-0.5*(SF_MS_scatter(lgMstar)*np.log(10))**2) # see ~/⁨Work⁩/⁨DeepFields⁩/⁨Works⁩/⁨20180000_COSMOS_SB_CO_Excitation⁩/⁨Gas_Modeling⁩/a_dzliu_code_test_integrating_lognormal.py
    # 
    # prepare output array
    output_dict = {}
    output_dict['z'] = []
    output_dict['cosmoAge'] = []
    output_dict['lgMstar'] = []
    output_dict['lgSFR'] = []
    output_dict['SB'] = []
    output_dict['Maj_kpc'] = []
    output_dict['Min_kpc'] = []
    output_dict['PA'] = []
    # 
    # loop redshift bins
    for i in range(len(z)-1):
        # 
        # debug
        #if z[i] < 1.0:
        #    continue
        # 
        # calc galaxy number density per stellar mass bin
        PhiMstar = 10**(lgPhiMstar[:,i]) * (lgMstar[1]-lgMstar[0]) # per dex --> per lgMstar logspace grid
        # 
        # calc comoving volume
        comoving_volume = calc_comoving_volume(z[i], z[i+1], area_arcmin2)
        # 
        # loop each stellar mass bin and generate galaxies therein
        for j in range(len(lgMstar)-1):
            # 
            # debug
            #if lgMstar[j] < 10.0:
            #    continue
            # 
            # galaxy number density multiplying the comoving volume gives the absolute galaxy number
            galaxy_number = np.int64(np.round((PhiMstar[j+1]+PhiMstar[j])/2.0 * comoving_volume))
            if galaxy_number <= 0:
                continue
            # 
            # generate uniformly-randomly distributed stellar masses within this stellar mass bin
            galaxy_Mstar = np.random.rand(galaxy_number) * (10**lgMstar[j+1] - 10**lgMstar[j]) + 10**lgMstar[j]
            galaxy_lgMstar = np.log10(galaxy_Mstar)
            # 
            # generate uniformly-randomly distributed redshifts within this redshift bin
            #galaxy_t_cosmic_age = np.random.rand(galaxy_number) * (t_cosmic_age[i] - t_cosmic_age[i+1]) + t_cosmic_age[i+1]
            galaxy_lgopz = np.random.rand(galaxy_number) * (lgopz[i+1] - lgopz[i]) + lgopz[i]
            galaxy_z = 10**(galaxy_lgopz)-1.0
            galaxy_t_cosmic_age = cosmo.age(galaxy_z).value
            # 
            # calc main sequence SFR, base on which we will generate the real SFR of each galaxy following lognormal distribution, see below.
            main_sequence_lgSFR = np.log10(calc_SFR_MS(galaxy_z, galaxy_lgMstar))  # this should be mean(lg(SFR)), 
            # 
            # define a subset of galaxies being starbursts, while the remain majority are main sequence galaxies.
            # for starburst fraction, see 
            # -- Conselice 2014 ARAA -- https://ned.ipac.caltech.edu/level5/March14/Conselice/Conselice_contents.html -- https://arxiv.org/abs/1403.2783
            # -- Conselice 2008 -- see text before Section 5.4.2 -- http://adsabs.harvard.edu/abs/2008MNRAS.386..909C
            # -- Conselice 2008 -- comparing Fig. 14 -- http://adsabs.harvard.edu/abs/2008MNRAS.386..909C
            #merger_fraction = 0.15 #<TODO># 
            ##if z[i] > 2.0:
            ##    merger_fraction = 0.15 # * (Schechter_M/10**10.0)**(-0.2)
            ##else:
            ##    merger_fraction = 0.01 * (1+z[i])**2.465
            ##    if lgMstar[j] > 10.0:
            ##        merger_fraction = merger_fraction * 10**(-0.6*(lgMstar[j]-10.0))
            ##starburst_fraction = 0.20
            #merger_starburst_phase_probability = 1.0
            #starburst_fraction = merger_fraction * merger_starburst_phase_probability
            # 
            # define a subset of galaxies being starbursts, while the remain majority are main sequence galaxies.
            # for starburst fraction,  
            # method 2, 
            # I now exactly follow Bethermin2017 which follows Bethermin2012a
            merger_fraction = 0.03 if (z[i]>=1.0) else ((z[i]-0.0)/(1.0-0.0)*(0.030-0.015)+0.015) # 1.5% at z=0.0 to 3% at z=1.0, then being constant.
            starburst_fraction = 0.03 if (z[i]>=1.0) else ((z[i]-0.0)/(1.0-0.0)*(0.030-0.015)+0.015) # 1.5% at z=0.0 to 3% at z=1.0, then being constant.
            # 
            # test 3 -- no starburst at all
            starburst_fraction = 0.0
            # 
            # calculate starburst number and set starburst boosting factor (in dex)
            starburst_number = np.int64(np.round(galaxy_number * starburst_fraction))
            starburst_boost = 0.6 # dex in SFR at a given lgMstar, assuming no stellar mass dependency
            # 
            if starburst_number > 0:
                # mix of main sequence and starburst galaxies
                # random number array
                rand_number_array = np.arange(galaxy_number) ; np.random.shuffle(rand_number_array)
                mask_starburst = rand_number_array.argsort()[0:starburst_number] # non duplicated indices
                mask_nonstarburst = rand_number_array.argsort()[starburst_number:]
                # generate lgSFR for nonstarburst and starburst galaxies respectively
                galaxy_lgSFR = np.full(shape=galaxy_number, fill_value=np.nan, dtype=np.float64)
                galaxy_lgSFR[mask_nonstarburst] = np.random.normal(main_sequence_lgSFR[mask_nonstarburst], SF_MS_scatter(galaxy_lgMstar[mask_nonstarburst]), galaxy_number-starburst_number) # mean(lg()), sigma(lg()), N
                galaxy_lgSFR[mask_starburst] = np.random.normal(main_sequence_lgSFR[mask_starburst]+starburst_boost, SF_MS_scatter(galaxy_lgMstar[mask_starburst]), starburst_number) # mean(lg()), sigma(lg()), N
                galaxy_starbursty = np.full(shape=galaxy_number, fill_value=0, dtype=int)
                galaxy_starbursty[mask_starburst] = 1
                # generate galaxy sizes
                galaxy_Maj_kpc = np.full(shape=galaxy_number, fill_value=np.nan, dtype=np.float64)
                galaxy_Maj_kpc[mask_nonstarburst] = calc_Galaxy_Size_vanderWel2014(galaxy_z[mask_nonstarburst], galaxy_Mstar[mask_nonstarburst], galaxy_type = 'SFG', add_noise = True)
                galaxy_Maj_kpc[mask_starburst] = calc_Galaxy_Size_vanderWel2014(galaxy_z[mask_starburst], galaxy_Mstar[mask_starburst], galaxy_type = 'QG', add_noise = True) #<TODO># we assume starbursts are as small as Quiescent Galaxies (red galaxies)
                galaxy_Min_kpc = galaxy_Maj_kpc * (np.random.random(galaxy_number)*(0.9-0.2)+0.2 + np.random.normal(0.0,0.05,galaxy_number)) # assuming b/a axis ratio is uniformly distributed within 0.2 - 1.0, with a blurring effect. 
                galaxy_PA = np.random.random(galaxy_number) * 360.0
                # set disc Height and calculate projected Minor axis
                ##gal_IA = random() # inclination angle
                ##galaxy_Hei_kpc = 0.1 * gal_Maj_kpc # disc height, assuming 0.1 * Maj [kpc] # Eva, Elisabete, -- 500 pc -- ask Sharon!
                ##galaxy_Min_kpc = gal_Maj_kpc * sin(gal_IA) + gal_Hei_kpc * cos(gal_IA)
                # append to output dict
                output_dict['z'].extend(galaxy_z.tolist())
                output_dict['cosmoAge'].extend(galaxy_t_cosmic_age.tolist())
                output_dict['lgMstar'].extend(galaxy_lgMstar.tolist())
                output_dict['lgSFR'].extend(galaxy_lgSFR.tolist())
                output_dict['SB'].extend(galaxy_starbursty.tolist())
                output_dict['Maj_kpc'].extend(galaxy_Maj_kpc.tolist())
                output_dict['Min_kpc'].extend(galaxy_Min_kpc.tolist())
                output_dict['PA'].extend(galaxy_PA.tolist())
                
            elif starburst_number == galaxy_number:
                # all are starburst galaxies (which is not possible)
                raise Exception('All are starburst galaxies?')
                #galaxy_lgSFR = np.random.normal(main_sequence_lgSFR+starburst_boost, SF_MS_scatter(galaxy_lgMstar), galaxy_number) # mean(lg()), sigma(lg()), N
                #galaxy_starbursty = np.full(shape=galaxy_number, fill_value=1, dtype=int) # all are starburst
                
            else:
                # all are main sequence galaxies (should also be not possible for current starburst fraction definition)
                #raise Exception('All are main sequence galaxies?')
                # 
                # generate lgSFR for nonstarburst and starburst galaxies respectively
                galaxy_lgSFR = np.random.normal(main_sequence_lgSFR, SF_MS_scatter(galaxy_lgMstar), galaxy_number) # mean(lg()), sigma(lg()), N
                galaxy_starbursty = np.full(shape=galaxy_number, fill_value=0, dtype=int) # all are MS type, no starburst
                # generate galaxy sizes
                galaxy_Maj_kpc = calc_Galaxy_Size_vanderWel2014(galaxy_z, galaxy_Mstar, galaxy_type = 'SFG', add_noise = True)
                galaxy_Min_kpc = galaxy_Maj_kpc * (np.random.random(galaxy_number)*(0.9-0.2)+0.2 + np.random.normal(0.0,0.05,galaxy_number)) # assuming b/a axis ratio is uniformly distributed within 0.2 - 1.0, with a blurring effect. 
                check_lower_boundary = (galaxy_Min_kpc<0.1) # make sure no b/a axis ratio is less than 0.1.
                if np.count_nonzero(check_lower_boundary) > 0:
                    galaxy_Min_kpc[check_lower_boundary] = np.random.random(np.count_nonzero(check_lower_boundary))*(0.9-0.2)+0.2
                galaxy_PA = np.random.random(galaxy_number) * 360.0
                # convert to arcsec
                galaxy_dL = cosmo.luminosity_distance(galaxy_z).to('Mpc').value
                galaxy_dA = galaxy_dL / (1.0+galaxy_z)**2
                kpc2arcsec = 1e-3/(galaxy_dA)/np.pi*180.0*3600.0
                galaxy_Maj_arcsec = galaxy_Maj_kpc * kpc2arcsec
                galaxy_Min_arcsec = galaxy_Min_kpc * kpc2arcsec
                # append to output dict
                output_dict['z'].extend(galaxy_z.tolist())
                output_dict['cosmoAge'].extend(galaxy_t_cosmic_age.tolist())
                output_dict['lgMstar'].extend(galaxy_lgMstar.tolist())
                output_dict['lgSFR'].extend(galaxy_lgSFR.tolist())
                output_dict['SB'].extend(galaxy_starbursty.tolist())
                output_dict['Maj_kpc'].extend(galaxy_Maj_kpc.tolist())
                output_dict['Min_kpc'].extend(galaxy_Min_kpc.tolist())
                output_dict['Maj_arcsec'].extend(galaxy_Maj_arcsec.tolist())
                output_dict['Min_arcsec'].extend(galaxy_Min_arcsec.tolist())
                output_dict['PA'].extend(galaxy_PA.tolist())
                
            # 
            print('z = %5.3f - %5.3f, lgMstar = %4.2f - %4.2f, comoving_volume = %.3e Mpc3, galaxy_number = %d, starburst_number = %d, merger_fraction = %.2f%%'%(z[i], z[i+1], lgMstar[j], lgMstar[j+1], comoving_volume, galaxy_number, starburst_number, merger_fraction*100.0))
            # 
            # 
            # 
            # <TODO> SFR_UV SFR_IR
            # 
            # 
            # 
            #sys.exit()
            # 
            # debug
            #if z[i] > 3.0:
            #    sys.exit()
        # 
        # dump output dict
        if i % 3 == 1 or i == len(z)-2:
            output_table = Table(output_dict)
            output_table.write('%s.dump.fits'%(output_name), overwrite=True)
    # 
    # output
    output_table = Table(output_dict)
    output_table.write('%s.fits'%(output_name), overwrite=True)
    print('Output to "%s.fits"!'%(output_name))






"""
"""
def measure_nearest_distances_for_grid_pixels(nx, ny, xypairs, return_grid_xy = False):
    # 
    # first find all nearest points in xyparis2 for each element in xypairs1
    #t_KDTree = KDTree(xypairs2)
    #dists, indicies = t_KDTree.query(xypairs1) #--> TOO SLOW
    # 
    #dists = distance.cdist(xypairs1, xypairs2) # output shape = (xypairs1.shape[0], xypairs2.shape[0])
    # 
    # check input dtype
    #if xypairs1.dtype != int and xypairs1.dtype != np.int32 and xypairs1.dtype != np.int64:
    #    raise Exception('Error! The first input array of measure_nearest_distances_for_grid_pixels() must be a integer array generated by np.mgrid[XX:XX,XX:XX]!')
    # 
    #gridy, gridx = np.mgrid[0:ny, 0:nx]
    #gridxy = np.array([gridx.flatten(), gridy.flatten()]).T
    #print('gridx.shape', gridx.shape)
    #print('gridxy.shape', gridxy.shape, 'gridxy.dtype', gridxy.dtype)
    #xypairs1 = gridxy
    #xypairs2 = xypairs
    #print('xypairs.shape', xypairs.shape, 'xypairs.dtype', xypairs.dtype, 'type(xypairs)', type(xypairs))
    # 
    # prepare output array
    dists = np.full((ny,nx), np.nan)
    # 
    # make sky patches (chunks)
    xstep = 300
    ystep = 300
    xbuffer = 10
    ybuffer = 10
    for y0 in trange(0, ny, ystep, desc='1st loop', leave=True):
        y1 = y0+ystep if y0+ystep<ny else ny
        ymask = np.logical_and(xypairs[:,1]>=y0-ybuffer, xypairs[:,1]<=y1+ybuffer)
        for x0 in trange(0, nx, xstep, desc='2nd loop', leave=False):
            x1 = x0+xstep if x0+xstep<nx else nx
            gridy, gridx = np.mgrid[y0:y1, x0:x1]
            gridxy = np.array([gridx.flatten(), gridy.flatten()]).T
            xmask = np.logical_and(xypairs[:,0]>=x0-xbuffer, xypairs[:,0]<=x1+xbuffer)
            mask = np.logical_and(ymask, xmask)
            if np.count_nonzero(mask) > 0:
                #indices = np.argwhere(mask).flatten()
                #print('np.count_nonzero(mask)', np.count_nonzero(mask))
                #print('Running distance.cdist()', gridxy.shape, xypairs[mask,:].shape)
                time_start = timeit.timeit()
                dists_2D = distance.cdist(gridxy, xypairs[mask,:])
                dists[y0:y1, x0:x1] = np.min(dists_2D, axis=1).reshape(dists[y0:y1, x0:x1].shape)
                time_end = timeit.timeit()
                #print('Elapsed %s'%(time_end - time_start))
                
                #sys.exit()
                
            else:
                dists[y0:y1, x0:x1] = +np.inf
    # 
    return dists


def minmax(a, b):
    if a <= b:
        return a, b
    else:
        return b, a


def maxmin(a, b):
    if a <= b:
        return b, a
    else:
        return a, b



"""
"""
def generate_galaxy_coordinates(RA00:float, Dec00:float, RA11:float, Dec11:float):
    # 
    # check existing output file
    output_name = 'datatable_generated_galaxies_with_coordinates'
    if os.path.isfile('%s.fits'%(output_name)):
        print('Found existing "%s.fits"! Will not overwrite it!'%(output_name))
        return
    # 
    # determine RA Dec range
    #input_image_file = '/Volumes/GoogleDrive/Shared drives/DeepFields/Data/COSMOS_Photos/VLA_3GHz/vla_3ghz_msmf.fits.gz' # 'vla_3ghz_msmf.fits.gz'
    #print('Reading "%s"'%(input_image_file))
    #hdulist = fits.open(input_image_file)
    #hdu = hdulist[0]
    #data_image = hdu.data
    #data_header = hdu.header
    #data_wcs = WCS(data_header, naxis=2)
    #data_nx = data_header['NAXIS1']
    #data_ny = data_header['NAXIS2']
    #pixelscale = proj_plane_pixel_scales(data_wcs)[1] * 3600.0 # arcsec
    #RADec00 = data_wcs.wcs_pix2world([(0, 0)], 0)[0]
    #RADec11 = data_wcs.wcs_pix2world([(data_nx, data_ny)], 1)[0]
    #RA_range = (RADec11[0], RADec00[0])
    #Dec_range = (RADec00[1], RADec11[1])
    # 
    # determine RA Dec range from user input
    data_header = hdr = fits.Header()
    pixelscale = 0.2
    RA00, RA11 = maxmin(RA00, RA11)
    Dec00, Dec11 = minmax(Dec00, Dec11)
    SkyCoord00 = SkyCoord(RA00, Dec00, frame=FK5, unit=(u.deg, u.deg))
    SkyCoord11 = SkyCoord(RA11, Dec11, frame=FK5, unit=(u.deg, u.deg))
    SkyCoord01 = SkyCoord(RA00, Dec11, frame=FK5, unit=(u.deg, u.deg))
    SkyCoord10 = SkyCoord(RA11, Dec00, frame=FK5, unit=(u.deg, u.deg))
    RA_extent = SkyCoord00.separation(SkyCoord10)
    Dec_extent = SkyCoord00.separation(SkyCoord01)
    data_nx = np.int64(np.round(RA_extent.to(u.arcsec).value / pixelscale))
    data_ny = np.int64(np.round(Dec_extent.to(u.arcsec).value / pixelscale))
    data_header['NAXIS1'] = data_nx
    data_header['NAXIS2'] = data_ny
    data_header['CDELT1'] = -pixelscale/3600.0
    data_header['CDELT2'] = +pixelscale/3600.0
    data_header['CRPIX1'] = 0.5 #<TODO># set the pixel center 0.5,0.5 for the ref RA Dec position
    data_header['CRPIX2'] = 0.5 #<TODO># set the pixel center 0.5,0.5 for the ref RA Dec position
    data_header['CRVAL1'] = RA00
    data_header['CRVAL2'] = Dec00
    data_header['CTYPE1'] = 'RA---TAN'
    data_header['CTYPE2'] = 'DEC--TAN'
    #print(data_header)
    data_wcs = WCS(data_header, naxis=2)
    #sys.exit()
    # 
    # read real galaxy catalog
    if not os.path.isfile('dump_posxy_real_galaxies.fits'):
        input_real_galaxy_catalog_file = '/Users/dzliu/Work/AlmaCosmos/Catalogs/A3COSMOS_Master_Catalog_20181106/master_catalog_single_entry_v20181109.fits'
        print('Reading "%s"'%(input_real_galaxy_catalog_file))
        real_galaxy_catalog = Table.read(input_real_galaxy_catalog_file)
        RA_real_galaxies = real_galaxy_catalog['RA'].data
        Dec_real_galaxies = real_galaxy_catalog['Dec'].data
        print('Running data_wcs.wcs_world2pix to convert real galaxy RA Dec to pixel coordinates')
        posxy_real_galaxies = data_wcs.wcs_world2pix(np.column_stack([RA_real_galaxies, Dec_real_galaxies]), 0)
        print('Writing to "%s" ...'%('dump_posxy_real_galaxies.fits'))
        dump_data_table = Table( {'posxy':posxy_real_galaxies} )
        dump_data_table.write('dump_posxy_real_galaxies.fits')
        print('Written to "%s"'%('dump_posxy_real_galaxies.fits'))
    else:
        print('Reading "%s"'%('dump_posxy_real_galaxies.fits'))
        dump_data_table = Table.read('dump_posxy_real_galaxies.fits')
        posxy_real_galaxies = dump_data_table['posxy'].data
    # 
    # compute the nearest distance from each grid pixel to the sources in the real galaxy catalog
    if not os.path.isfile('dump_gridxy_nearest_distances.fits'):
        print('Running measure_nearest_distances_for_grid_pixels for each pixel grid cell to real galaxy catalog')
        dists = measure_nearest_distances_for_grid_pixels(data_nx, data_ny, posxy_real_galaxies, return_grid_xy = True)
        print('Writing to "%s" ...'%('dump_gridxy_nearest_distances.fits'))
        dump_data_table = Table( {'dists':dists.flatten()} )
        dump_data_table.write('dump_gridxy_nearest_distances.fits')
        print('Written to "%s"'%('dump_gridxy_nearest_distances.fits'))
    else:
        print('Reading "%s"'%('dump_gridxy_nearest_distances.fits'))
        dump_data_table = Table.read('dump_gridxy_nearest_distances.fits')
        dists = dump_data_table['dists'].data
        dists = dists.reshape((data_ny, data_nx))
    # 
    # mask the area where there are no real galaxies locating within 0.7 arcsec
    thresh_distance_arcsec = 0.7
    thresh_distance_pixel = thresh_distance_arcsec / pixelscale
    thresh_distance_pixel = thresh_distance_pixel + 1.0 # plus a 0.5 pixel buffer because later we will jiggle posxy_model_galaxies to make it not exactly at the center of each pixel
    data_mask = (dists >= thresh_distance_pixel)
    print('dists.shape', dists.shape)
    print('Masked area to total area ratio: %.2f%%'%(np.count_nonzero(data_mask)*100.0/(data_nx*data_ny)))
    # 
    # read model galaxy z,Mstar,SFR properties catalog
    print('Reading "datatable_generated_galaxies.fits"')
    model_galaxy_table = Table.read('datatable_generated_galaxies.fits')
    # 
    # genereate our model galaxies only in the masked area, and jiggle a fractional pixel coordinate
    print('Making grid np.mgrid[0:%d, 0:%d]'%(data_ny, data_nx))
    data_gridy, data_gridx = np.mgrid[0:data_ny, 0:data_nx]
    #print('Transposing gridxy')
    #data_gridxy = np.array([data_gridx.flatten(), data_gridy.flatten()]).T
    #print('data_gridxy.shape', data_gridxy.shape)
    print('data_mask.shape', data_mask.shape)
    print('data_gridx.shape', data_gridx.shape, 'data_gridy.shape', data_gridy.shape)
    print('Running np.random.choice')
    posx_model_galaxies = np.random.choice(data_gridx[data_mask].flatten(), len(model_galaxy_table)) + (np.random.random(len(model_galaxy_table))-0.5)
    posy_model_galaxies = np.random.choice(data_gridy[data_mask].flatten(), len(model_galaxy_table)) + (np.random.random(len(model_galaxy_table))-0.5)
    print('Combining posxy and transposing posxy')
    posxy_model_galaxies = np.array([posx_model_galaxies, posy_model_galaxies]).T
    RADec_model_galaxies = data_wcs.wcs_pix2world(posxy_model_galaxies, 0)
    model_galaxy_table['RA'] = RADec_model_galaxies[:,0]
    model_galaxy_table['Dec'] = RADec_model_galaxies[:,1]
    # 
    # save 
    model_galaxy_table.write('%s.fits'%(output_name), overwrite=True)
    print('Output to "%s.fits"!'%(output_name))







def calc_Bethermin2014_U(z, is_starburst):
    # 
    # calculate <U> according to the <U> evolution track of MS galaxies in Bethermin 2014 arXiv (2015 A&A)
    # input $1 is redshift
    # input $2 is starburst-ness: 0 for MS, 1 for SB. 
    # output $0 is <U>
    # 
    # see paper http://fr.arxiv.org/pdf/1409.5796v2
    #     Fig.7 caption: (3.0+-1.1)*(1+z)**(1.8+-0.4) for MS, 31+-3 for SB. 
    # 
    Umean = 3.0 * (1.0+z)**1.8
    Umean[(is_starburst==1)] = 31.0
    return Umean




def get_SED_template(Age = None, EBV = None, Umean = None, Mstar = None, LIR = None):
    sys.path.append('/Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin')
    from michi2_read_lib_SEDs import (lib_file_get_header, lib_file_get_data_lines)
    # 
    stellar_lib_file_path = 'lib.FSPS.CSP.tau.1Gyr.Padova.BaSeL.Z0.0190.EBV.SED'
    stellar_lib_header = lib_file_get_header(stellar_lib_file_path)
    print(stellar_lib_header)
    stellar_lib_data_first_line = lib_file_get_data_lines(stellar_lib_file_path, 1, Lib_header=stellar_lib_header) # get the first line of each model block
    print(stellar_lib_data_first_line)
    stellar_SED_params_dict = {}
    for i,t in enumerate(stellar_lib_header['colnames']): 
        stellar_SED_params_dict[t] = stellar_lib_data_first_line[0,:,i].astype(stellar_lib_header['coltypes'][i])
    print(stellar_SED_params_dict)
    # 
    warm_dust_lib_file_path = 'lib.DL07.2010.03.18.spec.2.0.HiExCom.SED'
    warm_dust_lib_header = lib_file_get_header(warm_dust_lib_file_path)
    warm_dust_lib_data_first_line = lib_file_get_data_lines(warm_dust_lib_file_path, 1, Lib_header=warm_dust_lib_header) # get the first line of each model block
    warm_dust_SED_params_dict = {}
    for i,t in enumerate(warm_dust_lib_header['colnames']): 
        warm_dust_SED_params_dict[t] = warm_dust_lib_data_first_line[0,:,i].astype(warm_dust_lib_header['coltypes'][i])
    print(warm_dust_SED_params_dict)
    # 
    cold_dust_lib_file_path = 'lib.DL07.2010.03.18.spec.2.0.LoExCom.SED'
    cold_dust_lib_header = lib_file_get_header(cold_dust_lib_file_path)
    cold_dust_lib_data_first_line = lib_file_get_data_lines(cold_dust_lib_file_path, 1, Lib_header=cold_dust_lib_header) # get the first line of each model block
    cold_dust_SED_params_dict = {}
    for i,t in enumerate(cold_dust_lib_header['colnames']): 
        cold_dust_SED_params_dict[t] = cold_dust_lib_data_first_line[0,:,i].astype(cold_dust_lib_header['coltypes'][i])
    print(cold_dust_SED_params_dict)
    # 
    
    # 



def get_SED_template_Magdis2012(Umean, LIR, SB, qIR, z, dL):
    # 
    # prepare table_z_Umean 
    table_i_Umean = [ [ 1,  2.2], 
                      [ 2,  3.3], 
                      [ 3,  4.9], 
                      [ 4,  6.1], 
                      [ 5,  9.7], 
                      [ 6, 12.1], 
                      [ 7, 14.5], 
                      [ 8,   18], 
                      [ 9,   25], 
                      [10,   30], 
                      [11,   35], 
                      [12,   40], 
                      [13,   45], 
                      [14,   50] ]
    table_i_Umean = np.array(table_i_Umean)
    #if Umean < np.min(table_i_Umean[:,1]) or Umean > np.max(table_i_Umean[:,1]):
    #    raise ValueError('The input Umean is out of the allowed range!')
    # 
    sys.path.append('/Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin')
    from michi2_read_lib_SEDs import (integrate_LIR, integrate_vLv, spline)
    # 
    SED_X = None
    dict_SED_Y = {}
    global script_dir
    for i in np.concatenate([table_i_Umean[:,0], [101]]):
        temp_table = Table.read(os.path.join(os.path.dirname(script_dir), 'data', 'SED_template_Magdis', 'sed_z%d_U%d_radio.txt'%(i, i)), format='ascii.no_header')
        SED_Y = temp_table[temp_table.colnames[1]]
        if SED_X is None:
            SED_X = temp_table[temp_table.colnames[0]]
        if len(SED_X) != len(SED_Y):
            raise Exception('Error! SED libraries are inconsistent!')
        SED_Y = SED_Y / (2.99792458e5/SED_X) # vLv --> Lv
        #SED_Y = SED_Y / (1.0/SED_X)  -  (((2.99792458e5/SED_X)/1.4)**(-0.8) / 10**old_qIR * 71.163) # remove radio part
        #print('integrate_LIR')
        SED_LIR = integrate_LIR(SED_X, SED_Y, 0.0, 1.0) # SED template at z = 0 and dL = 1 Mpc
        #print(i, SED_LIR)
        dict_SED_Y[i] = SED_Y / (SED_LIR / (4*np.pi*1.0**2)) # normalized so that LIR at z = 0 and dL = 1 Mpc is 1.0.
        # check radio_flux_at_1p4GHz
        #radio_flux_at_1p4GHz = integrate_vLv(SED_X, SED_Y, 0.0, 1.0) / (4*np.pi*1.0**2) / 93.75 / 10**2.5
        #radio_flux_at_1p4GHz_2 = spline(SED_X, SED_Y, 2.99792458e5/1.4, xlog=1, ylog=1)
        #print('radio_flux_at_1p4GHz', radio_flux_at_1p4GHz, radio_flux_at_1p4GHz_2) #--> basically consistent
        #sys.exit()
    # 
    # interpolate Umean
    i2 = np.searchsorted(table_i_Umean[:,1], Umean, side='right') # find all 'a[i-1] <= Umean < a[i]'
    i2[(i2 >= len(table_i_Umean[:,1]))] = len(table_i_Umean[:,1])-1
    i1 = i2 - 1
    i1[(i1 < 0)] = 0
    a2 = (np.log10(Umean) - np.log10(table_i_Umean[:,1][i1])) / (np.log10(table_i_Umean[:,1][i2]) - np.log10(table_i_Umean[:,1][i1]))
    a1 = (np.log10(table_i_Umean[:,1][i2]) - np.log10(Umean)) / (np.log10(table_i_Umean[:,1][i2]) - np.log10(table_i_Umean[:,1][i1]))
    #a1[k] + a2[k] = 1
    # 
    # loop and apply SB, z, dL
    SED_bands = [70.0, 100.0, 160.0, 250.0, 350.0, 450.0, 500.0, 850.0, 1100.0, 1200.0, 1300.0, 2000.0, 3000.0, 2.99792458e5/3.0, 2.99792458e5/1.4]
    SED_fluxes = np.full((len(i1), len(SED_bands)), np.nan)
    for k in tqdm(range(len(i1))):
        # 
        if SB[k] == 0:
            SED_Y = a1[k] * dict_SED_Y[table_i_Umean[:,0][i1[k]]] + a2[k] * dict_SED_Y[table_i_Umean[:,0][i2[k]]]
        else:
            SED_Y = dict_SED_Y[101]
        # 
        Obs_Y = SED_Y / (4*np.pi*dL[k]**2) * (1.0+z[k]) * LIR[k]
        Obs_X = SED_X * (1.0+z[k])
        Obs_LIR = integrate_LIR(Obs_X, Obs_Y, z[k], dL[k])
        # 
        # qIR
        #old_qIR = 2.5
        new_qIR = qIR[k]
        old_radio_flux_at_1p4GHz = spline(SED_X, Obs_Y, 2.99792458e5/1.4, xlog=1, ylog=1) # rest-frame 1.4 GHz
        new_radio_flux_at_1p4GHz = Obs_LIR / (4*np.pi*dL[k]**2) / 93.75 / 10**new_qIR
        old_radio_SED = ((2.99792458e5/SED_X)/1.4)**(-0.8) * old_radio_flux_at_1p4GHz
        new_radio_SED = ((2.99792458e5/SED_X)/1.4)**(-0.8) * new_radio_flux_at_1p4GHz
        Obs_Y = Obs_Y - old_radio_SED + new_radio_SED
        #print('LIR_SED %e vs LIR_galaxy %e, z = %.4f, dL = %.2f, SB = %d, Umean = %.1f (model %.1f - %.1f), qIR = %.2f'%(Obs_LIR, LIR[k], z[k], dL[k], SB[k], Umean[k], table_i_Umean[:,1][i1[k]], table_i_Umean[:,1][i2[k]], qIR[k]))
        #print(a1[k], a2[k])
        # 
        # dump
        #dump_data_table = Table( {'Obs_X':Obs_X, 'Obs_Y':Obs_Y, 'old_radio_SED':old_radio_SED, 'new_radio_SED':new_radio_SED} )
        #dump_data_table.write('dump_data_table_SED.fits', overwrite=True)
        #print('Dump to "dump_data_table_SED.fits"')
        # 
        # get SED fluxes
        SED_fluxes[k,:] = spline(SED_X, Obs_Y, SED_bands, xlog=1, ylog=1) # rest-frame 1.4 GHz
        # 
        # 
        #sys.exit()
        #if k > 1000:
        #    break
    # 
    return SED_bands, SED_fluxes



"""
"""
def generate_galaxy_SEDs():
    # 
    # check existing output file
    output_name = 'datatable_generated_galaxies_with_SED_params'
    if os.path.isfile('%s.fits'%(output_name)):
        print('Found existing "%s.fits"! Will not overwrite it!'%(output_name))
        return
    # 
    # read model galaxy table
    print('Reading "datatable_generated_galaxies_with_coordinates.fits"')
    model_galaxy_table = Table.read('datatable_generated_galaxies_with_coordinates.fits')
    z = model_galaxy_table['z'].data
    SB = model_galaxy_table['SB'].data
    Umean = calc_Bethermin2014_U(z, SB)
    model_galaxy_table['Umean'] = Umean
    model_galaxy_table['sSFR'] = 10**(model_galaxy_table['lgSFR']-model_galaxy_table['lgMstar']+9.0) # Gyr
    # 
    # calc 
    dL = cosmo.luminosity_distance(z).to(u.Mpc).value # Mpc
    sSFR = model_galaxy_table['sSFR'].data
    lgMstar = model_galaxy_table['lgMstar'].data
    Mstar = 10**lgMstar
    LIR = 10**(model_galaxy_table['lgSFR'].data) * 1e10
    zAge = model_galaxy_table['cosmoAge'].data
    #qIR = 2.35*(1.0+z)**(-0.12)+np.log10(1.91) # Magnelli 2015A%26A...573A..45M
    qIR = 2.85*(1.0+z)**(-0.22) # Delhaize 2017A&A...602A...4D
    # 
    # assuming galaxy age
    Age = 1.0/sSFR # set galaxy age to 1.0/sSFR (~5Gyr at z~0.1, ~1Gyr at z~1, to ~0.4Gyr at z~4)
    mask = (lgMstar<12)
    Age[mask] = Age[mask] * 10**((lgMstar[mask]-12)*0.12) # set younger galaxy age for lower mass galaxies -- 
    mask = (Age>zAge)
    Age[mask] = zAge[mask] # limit galaxy age to no younger than cosmic age
    mask = (Age>1.5)
    Age[mask] = 1.5 # limit galaxy age to no older than 1.5Gyr, because they are star-forming. 
    mask = (Age<0.3)
    Age[mask] = 0.3 # limit galaxy age to no younger than 300Myr. 
    # 
    EBV = np.full(len(LIR), 0.0)
    mask = (LIR>1e9)
    EBV[mask] = np.log10(LIR[mask]/1e9)*0.15 # E(B-V) increases with L_IR -- https://arxiv.org/pdf/1403.3615.pdf -- Figure 10
    mask = (Age>0.5)
    EBV[mask] = EBV[mask] * ((Age[mask]/0.5)**0.1) # E(B-V) increases with Age or Metallicity, and only when Age>500Myr. 
    # 
    model_galaxy_table['EBV'] = EBV
    model_galaxy_table['Age'] = Age
    model_galaxy_table['dL'] = dL
    # 
    ###dA = dL / (1.0+z)**2
    ###kpc2arcsec = 1e-3/(dA)/np.pi*180.0*3600.0
    ###model_galaxy_table['Maj_arcsec'] = model_galaxy_table['Maj_kpc'] * kpc2arcsec
    ###model_galaxy_table['Min_arcsec'] = model_galaxy_table['Min_kpc'] * kpc2arcsec
    # 
    # 
    # calc SED flux
    #SED = get_SED_template(EBV=EBV, Age=Age, Umean=Umean, Mstar=Mstar, LIR=LIR)
    # 
    # 
    # calc SED flux (method 2)
    #SED_bands, SED_fluxes = get_SED_template_Magdis2012(Umean, LIR, SB, qIR, z, dL)
    #for i in range(len(SED_bands)):
    #    model_galaxy_table['SED_%g'%(SED_bands[i])] = SED_fluxes[:,i]
    # 
    # 
    # calc SED flux (at rest-frame 1.4 GHz only)
    radio_flux_at_rest_frame_1p4GHz = LIR / (4*np.pi*dL**2) / 93.75 / 10**qIR
    radio_flux_at_obs_frame_3GHz = radio_flux_at_rest_frame_1p4GHz * (1.0+z) * (3.0*(1.0+z)/1.4)**(-0.8)
    model_galaxy_table['SED_3GHz'] = radio_flux_at_obs_frame_3GHz # mJy
    # 
    # 
    # save 
    model_galaxy_table.write('%s.fits'%(output_name), overwrite=True)
    print('Output to "%s.fits"!'%(output_name))
    # 
    #one_galaxy = CrabGalaxy(shape = 'sersic', size = {'major':1.0, 'minor':0.5, 'PA':20.0, 'n':0.5})












"""
"""
def make_noisy_image_with_noise_map(noise_map):
    # 
    image_data = np.full(noise_map.shape, np.nan)
    # 
    ny = noise_map.shape[-2]
    nx = noise_map.shape[-1]
    # 
    print('Making noisy image with the input noise map (%dx%d)'%(nx, ny))
    # 
    # make sky patches (chunks) and loop galaxies within each sky patch
    xstep = 300
    ystep = 300
    xbuffer = 0
    ybuffer = 0
    for y0 in trange(0, ny, ystep, desc='1st loop', leave=True):
    #for y0 in range(0, ny, ystep):
        y1 = y0+ystep if y0+ystep<ny else ny
        for x0 in trange(0, nx, xstep, desc='2nd loop', leave=False):
        #for x0 in range(0, nx, xstep):
            x1 = x0+xstep if x0+xstep<nx else nx
            #print('np.isnan(noise_map[%d:%d, %d:%d])'%(y0, y1, x0, x1))
            mask_nan = ~np.isnan(noise_map[y0:y1, x0:x1])
            if np.count_nonzero(mask_nan)>0:
                #print('np.mgrid[%d:%d, %d:%d]'%(y0, y1, x0, x1))
                gridy, gridx = np.mgrid[y0:y1, x0:x1]
                gridxy = np.array([gridx.flatten(), gridy.flatten()]).T
                noise_data = noise_map[y0:y1, x0:x1][mask_nan]
                #print('3-sigma clip median')
                noise_median = np.median(sigma_clip(noise_data, sigma=3, masked=False)) # use median to filter out 
                #noise_median = np.median(noise_data) # use median to filter out 
                #print('np.random.normal')
                image_data[y0:y1, x0:x1][mask_nan] = np.random.normal(0.0, noise_median, np.count_nonzero(mask_nan))
                # 
        #print('dump fits image')
        #dump_hdu = fits.PrimaryHDU(data=image_data)
        #dump_hdu.writeto('dump_noisy_image.dump.fits', overwrite=True)
        #print('Output to "dump_noisy_image.dump.fits"!')
    # 
    print('Now made the noisy image')
    # 
    return image_data









"""
"""
def generate_image(input_wavelength_um_or_band_name, input_noise_map_Jy_per_beam, input_beam_FWHM_arcsec):
    # 
    # check existing output file
    output_name = 'simulated_image'
    if os.path.isfile('%s.fits'%(output_name)):
        print('Found existing "%s.fits"! Will not overwrite it!'%(output_name))
        return
    # 
    # set debug mode
    debug_mode = False
    verbose_mode = False
    # 
    # read model galaxy table
    if debug_mode:
        print('Reading "datatable_generated_galaxies_with_SED_params.fits"')
        model_galaxy_table = Table.read('datatable_generated_galaxies_with_SED_params.fits')
    else:
        print('Reading "datatable_generated_galaxies_with_SED_params.fits"')
        model_galaxy_table = Table.read('datatable_generated_galaxies_with_SED_params.fits')
    z = model_galaxy_table['z'].data
    # 
    # find the input wavelength in the mock galaxy catalog
    if not ('SED_%s'%(input_wavelength_um_or_band_name) in model_galaxy_table.colnames):
        raise Exception('Error! The input wavelength %s is not in the mock galaxy catalog!'%(input_wavelength_um_or_band_name))
        sys.exit()
    # 
    galaxy_fluxes = model_galaxy_table['SED_%s'%(input_wavelength_um_or_band_name)].data
    galaxy_z = model_galaxy_table['z'].data
    galaxy_RA = model_galaxy_table['RA'].data
    galaxy_Dec = model_galaxy_table['Dec'].data
    galaxy_Maj_arcsec = model_galaxy_table['Maj_arcsec'].data
    galaxy_Min_arcsec = model_galaxy_table['Min_arcsec'].data
    galaxy_PA = model_galaxy_table['PA'].data
    z = galaxy_z
    # 
    # read the input noise map
    # create noisy data with the input noise information
    if os.path.isfile('dump_noisy_image.fits'):
        # read dump noisy image
        print('Reading "dump_noisy_image.fits"!')
        hdulist = fits.open('dump_noisy_image.fits')
        hdu = hdulist[0]
        data_header = hdu.header
        data_image = hdu.data
        data_wcs = WCS(data_header, naxis=2)
    else:
        # read input_noise_map_Jy_per_beam
        hdulist = fits.open(input_noise_map_Jy_per_beam)
        hdu = hdulist[0]
        data_header = hdu.header
        data_image = hdu.data
        data_wcs = WCS(data_header, naxis=2)
        # 
        # chop dimension
        while len(data_image.shape) > 2:
            data_image = data_image[0]
        # 
        # dump noisy image
        data_image = make_noisy_image_with_noise_map(data_image)
        print('Writing to "dump_noisy_image.fits" ...')
        dump_hdu = fits.PrimaryHDU(data = data_image, header = data_wcs.to_header())
        dump_hdu.writeto('dump_noisy_image.fits', overwrite=True)
        print('Written to "dump_noisy_image.fits"')
    # 
    # use dtype float32
    data_image = data_image.astype(np.float32, copy=False)
    # 
    # obtain pixscale
    pixscale = proj_plane_pixel_scales(data_wcs)[1] * 3600.0
    print('pixscale', pixscale)
    # 
    # generate pixel grid
    image_size_y, image_size_x = data_image.shape[-2:] # N_X, N_Y
    #mgrid_y, mgrid_x = np.mgrid[0:image_size_y, 0:image_size_x]
    print('image_size_x', image_size_x, data_header['NAXIS1'])
    print('image_size_y', image_size_y, data_header['NAXIS2'])
    # 
    # convert sky2xy
    if debug_mode:
        if os.path.isfile('dump_posxy_model_galaxies_when_generating_image.debug.fits'):
            print('Reading "%s"'%('dump_posxy_model_galaxies_when_generating_image.debug.fits'))
            dump_data_table = Table.read('dump_posxy_model_galaxies_when_generating_image.debug.fits')
            posxy = dump_data_table['posxy'].data
        else:
            print('Running wcs_world2pix for %d sources'%(len(galaxy_RA)))
            posxy = data_wcs.wcs_world2pix(np.column_stack((galaxy_RA, galaxy_Dec)), 0) # 0-based pixel coordinates
            print('Writing to "%s" ...'%('dump_posxy_model_galaxies_when_generating_image.debug.fits'))
            dump_data_table = Table( {'posxy':posxy} )
            dump_data_table.write('dump_posxy_model_galaxies_when_generating_image.debug.fits')
            print('Written to "%s"'%('dump_posxy_model_galaxies_when_generating_image.debug.fits'))
    else:
        if os.path.isfile('dump_posxy_model_galaxies_when_generating_image.fits'):
            print('Reading "%s"'%('dump_posxy_model_galaxies_when_generating_image.fits'))
            dump_data_table = Table.read('dump_posxy_model_galaxies_when_generating_image.fits')
            posxy = dump_data_table['posxy'].data
        else:
            print('Running wcs_world2pix for %d sources'%(len(galaxy_RA)))
            posxy = data_wcs.wcs_world2pix(np.column_stack((galaxy_RA, galaxy_Dec)), 0) # 0-based pixel coordinates
            print('Writing to "%s" ...'%('dump_posxy_model_galaxies_when_generating_image.fits'))
            dump_data_table = Table( {'posxy':posxy} )
            dump_data_table.write('dump_posxy_model_galaxies_when_generating_image.fits')
            print('Written to "%s"'%('dump_posxy_model_galaxies_when_generating_image.fits'))
    # 
    print('posxy.shape', posxy.shape)
    print('z.shape', z.shape)
    xypairs = posxy
    ny = image_size_y
    nx = image_size_x
    # 
    # set a redshift lower limit so that we do not have local galaxies
    z_lower_limit = 0.0016 # this make sure no source has flux above 25 mJy or so, but still have a few local sources have fluxes ~ 15 mJy, similar to the real data. 
    # 
    # set a flux upper limit
    #flux_upper_limit = np.nan #<TODO># 
    # 
    # resume from dumped iskypatch
    iskypatch_resume = -1
    if os.path.isfile('Output_simulated_image.dump.iskypatch.txt') and os.path.isfile('Output_simulated_image.dump.fits'):
        with open('Output_simulated_image.dump.iskypatch.txt', 'r') as fp:
            iskypatch_resume = np.int64(fp.readline())
        data_image = fits.getdata('Output_simulated_image.dump.fits')
        print('Resuming from iskypatch %d!'%(iskypatch_resume))
    # 
    # prepare to write skypatch boxes into a ds9 region file
    skypatch_region_file = open('ds9_regions_skypatch_boxes.reg', 'w')
    skypatch_region_file.write('# Region file format: DS9 version 4.1\n')
    skypatch_region_file.write('image\n')
    # 
    # make sky patches (chunks) and loop galaxies within each sky patch
    xstep = 800
    ystep = 800
    xbuffer = int(np.ceil(10.0/pixscale)) # 10 arcsec buffer at each side
    ybuffer = int(np.ceil(10.0/pixscale)) # 10 arcsec buffer at each side
    iskypatch = 0
    nskypatch = 0
    skypatch_ygrid, skypatch_xgrid = np.mgrid[0:ystep+ybuffer+ybuffer, 0:xstep+xbuffer+xbuffer] # skypatch mgrid with buffer, skypatch[0:,0:] is data[y0-ybuffer:y1+ybuffer,x0-xbuffer:x1+xbuffer]
    skypatch_image = np.full((ystep+ybuffer+ybuffer, xstep+xbuffer+xbuffer), 0.0)
    PSF_size_pixel = int(np.ceil(input_beam_FWHM_arcsec*10/pixscale)*2+1) # about 20 times the input FWHM, make it an odd number
    PSF_FWHM_pixel = input_beam_FWHM_arcsec/pixscale
    PSF_sigma_pixel = PSF_FWHM_pixel / (2.0*np.sqrt(2.0*np.log(2.0)))
    #PSF_kernel_func = Gaussian2D(x_mean=(PSF_FWHM_pixel-1)/2, y_mean=(PSF_FWHM_pixel-1)/2, x_fwhm=PSF_FWHM_pixel, y_fwhm=PSF_FWHM_pixel)
    #PSF_kernel_image = PSF_kernel_func.evaluate()
    PSF_kernel_image = Gaussian2DKernel(PSF_sigma_pixel)
    for y0 in trange(4000, ny, ystep, desc='1st loop', leave=True):
        y1 = y0+ystep if y0+ystep<ny else ny
        ymask = np.logical_and(xypairs[:,1]>=y0-ybuffer, xypairs[:,1]<y1+ybuffer)
        for x0 in trange(0, nx, xstep, desc='2nd loop', leave=False):
            x1 = x0+xstep if x0+xstep<nx else nx
            # 
            # check valid pixel within this skypatch
            mask_nan = ~np.isnan(data_image[y0:y1, x0:x1])
            # 
            # check whether we will resume from certain iskypatch or not
            skip_due_to_resume = False
            if iskypatch_resume > 0:
                if iskypatch <= iskypatch_resume:
                    skip_due_to_resume = True
            # 
            # if we have valid pixel within this skypatch, and we are not skipping this skypatch due to resuming
            if np.count_nonzero(mask_nan)>0 and skip_due_to_resume==False:
                xmask = np.logical_and(xypairs[:,0]>=x0-xbuffer, xypairs[:,0]<x1+xbuffer)
                mask = np.logical_and.reduce((ymask, xmask, z>z_lower_limit))
                if np.count_nonzero(mask) > 0:
                    if verbose_mode:
                        print('\rInjecting %d sources within sky patch %s-%s,%s-%s'%(np.count_nonzero(mask), x0, x1, y0, y1))
                        time_start = timeit.default_timer()
                    # 
                    buf_y0 = y0-ybuffer
                    buf_y1 = y1+ybuffer
                    buf_x0 = x0-xbuffer
                    buf_x1 = x1+xbuffer
                    buf_dy = buf_y1-buf_y0 # ystep+ybuffer+ybuffer
                    buf_dx = buf_x1-buf_x0 # xstep+xbuffer+xbuffer
                    skypatch_image = skypatch_image * 0.0
                    # 
                    skypatch_region_file.write('box(%.1f,%.1f,%.1f,%.1f,0) # text={box %d (%d galaxies)}\n'%((x0+x1)/2.0, (y0+y1)/2.0, (x1-x0), (y1-y0), iskypatch, np.count_nonzero(mask)))
                    # 
                    # 
                    idx = np.arange(xypairs.shape[0])[mask] # indicies of selected model galaxies
                    for k in idx:
                        # 
                        #one_galaxy = CrabGalaxy(shape = 'Gaussian', 
                        #                        size = {'major':galaxy_Maj_arcsec[k], 
                        #                                'minor':galaxy_Min_arcsec[k], 
                        #                                'PA':galaxy_PA[k] },
                        #                        #SED = {'w': 2.99792458e5/3.0, 
                        #                        #       'f': galaxy_fluxes[k] * 1e-3, 
                        #                        #       'f_err': 0.0, 
                        #                        #       'band_name': 'VLA_3GHz'} #-- not useful for now
                        #                        )
                        #                        # here major minor are FWHM [arcsec], but what modeled are R_eff. 
                        # 
                        #data_image[y0:y1, x0:x1] = one_galaxy.inject_into_image(data_image[y0:y1, x0:x1], pixscale, xypairs[k,0], xypairs[k,1], mgrid_x, mgrid_y)
                        #skypatch_image = one_galaxy.inject_into_image(skypatch_image, 
                        #                                              pixscale, 
                        #                                              flux = galaxy_fluxes[k] * 1e-3, # Jy
                        #                                              x = xypairs[k,0]-buf_x0, 
                        #                                              y = xypairs[k,1]-buf_y0, 
                        #                                              mgrid_x = skypatch_xgrid, 
                        #                                              mgrid_y = skypatch_ygrid, 
                        #                                             )
                        # 
                        morph = CrabGalaxyMorphology(shape = 'Gaussian', 
                                                     size = {'major':galaxy_Maj_arcsec[k]/pixscale, 
                                                             'minor':galaxy_Min_arcsec[k]/pixscale, 
                                                             'PA':galaxy_PA[k] }, 
                                                    )
                        flux = galaxy_fluxes[k] * 1e-3 # Jy
                        flux = flux * PSF_FWHM_pixel**2 # Jy/beam
                        x = xypairs[k,0]-buf_x0 # coordinate (0-based) inside the skypatch (with buffer)
                        y = xypairs[k,1]-buf_y0 # coordinate (0-based) inside the skypatch (with buffer)
                        if galaxy_Maj_arcsec[k]/pixscale < 1.0:
                            # if the galaxy is smaller than one pixel, we set it as a delta function (single-pixel point source) and align it onto the pixel center...
                            if verbose_mode:
                                print('\rSource %d is smaller than one pixel (%.3f arcsec), we set it as a point source and align its position %.3f,%.3f to the pixel center'%(k, galaxy_Maj_arcsec[k], x, y))
                            y, x = np.int64(np.round(y)), np.int64(np.round(x))
                            if y < buf_dy and x < buf_dx:
                                skypatch_image[y, x] += flux
                        else:
                            skypatch_image += flux * morph.shape.func((skypatch_xgrid-x), (skypatch_ygrid-y))
                        # 
                    # 
                    # convolve with beam
                    if verbose_mode:
                        print('\rConvolving with the PSF')
                    skypatch_image = convolve(skypatch_image, PSF_kernel_image)
                    data_image[y0:y1, x0:x1] += skypatch_image[ybuffer:-ybuffer, xbuffer:-xbuffer]
                    # 
                    # 
                    if verbose_mode:
                        time_end = timeit.default_timer()
                        print('\rElapsed %.3f seconds'%(time_end - time_start))
                    # 
                    nskypatch += 1
                    # 
                    # dump
                    if nskypatch > 0 and (nskypatch == 5 or nskypatch % 37 == 0):
                        if verbose_mode:
                            print('\rWriting to "Output_simulated_image.dump.fits" ...')
                        dump_hdu = fits.PrimaryHDU(data = data_image, header = data_wcs.to_header())
                        dump_hdu.writeto('Output_simulated_image.dump.fits', overwrite=True)
                        if verbose_mode:
                            print('\rWritten to "Output_simulated_image.dump.fits"')
                        # 
                        with open('Output_simulated_image.dump.iskypatch.txt', 'w') as fp:
                            fp.write('%d\n'%(iskypatch))
                        if verbose_mode:
                            print('\rWritten to "Output_simulated_image.dump.iskypatch.txt"')
                        # 
                        if verbose_mode:
                            print('\rWriting to "dump_PSF_kernel_image.fits" ...')
                        dump_hdu = fits.PrimaryHDU(data = PSF_kernel_image)
                        dump_hdu.writeto('dump_PSF_kernel_image.fits', overwrite=True)
                        if verbose_mode:
                            print('\rWritten to "dump_PSF_kernel_image.fits"')
                        # 
                        skypatch_region_file.flush()
                        # 
                        #sys.exit()
                    # 
                    if verbose_mode:
                        print('\r')
            # 
            iskypatch += 1
        # 
        #break
    # 
    skypatch_region_file.close()
    # 
    # 
    if debug_mode:
        print('Writing to "Output_simulated_image.debug.fits" ...')
        output_hdu = fits.PrimaryHDU(data = data_image, header = data_wcs.to_header())
        output_hdu.writeto('Output_simulated_image.debug.fits', overwrite=True)
        print('Output to "Output_simulated_image.debug.fits"!')
    else:
        print('Writing to "Output_simulated_image.fits" ...')
        output_hdu = fits.PrimaryHDU(data = data_image, header = data_wcs.to_header())
        output_hdu.writeto('Output_simulated_image.fits', overwrite=True)
        print('Output to "Output_simulated_image.fits"!')
    























####################################
#               MAIN               #
####################################

if __name__ == '__main__':
    # 
    RA00, RA11 = maxmin(149.53, 150.76)
    Dec00, Dec11 = minmax(1.615, 2.88)
    SkyCoord00 = SkyCoord(RA00, Dec00, frame=FK5, unit=(u.deg, u.deg))
    SkyCoord11 = SkyCoord(RA11, Dec11, frame=FK5, unit=(u.deg, u.deg))
    SkyCoord01 = SkyCoord(RA00, Dec11, frame=FK5, unit=(u.deg, u.deg))
    SkyCoord10 = SkyCoord(RA11, Dec00, frame=FK5, unit=(u.deg, u.deg))
    RA_extent_00 = SkyCoord00.separation(SkyCoord10)
    RA_extent_11 = SkyCoord01.separation(SkyCoord11)
    Dec_extent_00 = SkyCoord00.separation(SkyCoord01)
    Dec_extent_11 = SkyCoord10.separation(SkyCoord11)
    print('RA_extent_00 = %s deg'%(RA_extent_00.to(u.deg).value))
    print('RA_extent_11 = %s deg'%(RA_extent_11.to(u.deg).value))
    print('Dec_extent_00 = %s deg'%(Dec_extent_00.to(u.deg).value))
    print('Dec_extent_11 = %s deg'%(Dec_extent_11.to(u.deg).value))
    area = ((RA_extent_00 + RA_extent_11) / 2.0) * (Dec_extent_00)
    area_deg2 = area.to(u.deg**2).value
    print('area_deg2 = %s'%(area_deg2))
    area_arcmin2 = area.to(u.arcmin**2).value
    print('area_arcmin2 = %s'%(area_arcmin2))
    # 
    #generate_galaxy_numbers(area_arcmin2)
    # 
    #generate_galaxy_coordinates(RA00, Dec00, RA11, Dec11)
    # 
    generate_galaxy_SEDs()
    # 
    generate_image('3GHz', 'input_images/vla_3ghz_msmf.rms.fits', 0.75)








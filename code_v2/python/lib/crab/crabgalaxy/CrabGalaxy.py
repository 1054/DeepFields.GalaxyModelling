#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# see also '/Users/dzliu/Work/AlmaCosmos/Simulations/CASA_Sim/Sim_a_sample_of_circular_Gaussian_sources/run_stacking/a_dzliu_code_plot_radial_profile.py'
# Last update: 2020-01-17 11h33m CET
# 


import os
import sys
import re
#import six
#import math
import numpy as np
import scipy
import astropy
from astropy.table import Table
from scipy.special import gamma, gammainc, gammaincinv, binom
from scipy.interpolate import interp1d, interp2d
#from scipy.signal import resample
#from datetime import datetime
#from astropy.modeling.models import Gaussian2D
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve, convolve_fft
from astropy.io import fits
#from astropy.wcs import WCS
import matplotlib as mpl
import matplotlib.pyplot as plt





def func_interp1d_loglog(x, y, kind = 'linear', normalization = 1.0):
    if type(x) is list: x = np.array(x)
    if type(y) is list: y = np.array(y)
    mask = np.logical_and(x>0, y>0)
    func = interp1d(np.log10(x[mask]), np.log10(y[mask]), kind = kind, fill_value = (np.nan, np.nan), bounds_error = False)
    return lambda t: 10**(func(np.log10(t))) * normalization



class CrabGalaxySED(object):
    # 
    def __init__(self, w = None, f = None, f_err = None, f_unit = None, band_name = None):
        self.w = []
        self.f = []
        self.f_err = []
        self.f_unit = []
        self.band_name = []
        self.allowed_templates = ['magdis2012', 'dzliu', ]
        self.template = None # SED template func
        self.template_properties = {} # 'z', 'dL', etc.
        if not (w is None and f is None):
            self.set_data(w = w, f = f, f_err = f_err, f_unit = f_unit, band_name = band_name)
        #if not (template_name == ''):
        #    self.set_SED_template(template_name)
    # 
    def set_data(self, w = None, f = None, f_err = None, f_unit = None, band_name = None):
        if (w is None and f is not None) or (w is not None and f is None):
            raise ValueError('Error input of w and f! Both must be given!')
        else:
            if (np.isscalar(w) and not np.isscalar(f)) or (not np.isscalar(w) and np.isscalar(f)):
                # check type consistency
                raise ValueError('Error input of w and f! One is scalar while the other is not!')
            elif (np.isscalar(w) and np.isscalar(f)):
                # if both scalar, then append
                self.w.append(w)
                self.f.append(f)
                if f_err is None:
                    self.f_err.append(np.nan)
                else:
                    if not np.isscalar(f_err):
                        raise ValueError('Error input of f_err! Should be scalar!')
                    else:
                        self.f_err.append(f_err)
                if f_unit is None:
                    self.f_unit.append('Jy')
                else:
                    if not np.isscalar(f_unit):
                        raise ValueError('Error input of f_unit! Should be scalar!')
                    else:
                        self.f_unit.append(f_unit)
                if band_name is None:
                    self.band_name.append('N/A')
                else:
                    if not np.isscalar(band_name):
                        raise ValueError('Error input of band_name! Should be scalar!')
                    else:
                        self.band_name.append(band_name)
            else:
                # if both list, then extend
                if len(w) != len(f):
                    raise ValueError('Error input of w and f! Inconsistent length!')
                else:
                    self.w.extend(np.array(w).flatten().tolist())
                    self.f.extend(np.array(f).flatten().tolist())
                    if f_err is None:
                        self.f_err.extend(np.full(np.array(f).flatten().shape, np.nan))
                    else:
                        if len(w) != len(f_err):
                            raise ValueError('Error input of f_err! Inconsistent length!')
                        else:
                            self.f_err.extend(np.array(f_err).flatten().tolist())
                    if f_unit is None:
                        self.f_unit.extend(np.full(np.array(f).flatten().shape, 'Jy'))
                    else:
                        if len(w) != len(f_unit):
                            raise ValueError('Error input of f_unit! Inconsistent length!')
                        else:
                            self.f_unit.extend(np.array(f_unit).flatten().tolist())
                    if band_name is None:
                        self.band_name.extend(np.full(np.array(f).flatten().shape, 'N/A'))
                    else:
                        if len(w) != len(band_name):
                            raise ValueError('Error input of band_name! Inconsistent length!')
                        else:
                            self.band_name.extend(np.array(band_name).flatten().tolist())
    # 
    def set_SED_template(self, template_name:str, **kwargs):
        if template_name is None:
            raise ValueError('Error! Please input a template name to the called function CrabGalaxySED::set_SED_template(). It should be one of %s (case insensitive).'%(str(self.allowed_templates)))
        else:
            if not (template_name.lower() in self.allowed_templates):
                raise ValueError('Error! The input template name %s is not in the allowed list of %s (case insensitive)!'%(template_name, str(self.allowed_templates)))
            # 
            self.template_properties = {} # clear and reset
            for key in kwargs.keys(): # store user input SED template property keys
                if re.match('(redshift|z)', key, re.IGNORECASE):
                    self.template_properties['z'] = kwargs[key] # recognize redshift
                else:
                    self.template_properties[key] = kwargs[key]
            # 
            if template_name.lower() == 'magdis2012':
                if not ('z' in kwargs.keys()):
                    raise ValueError('Error! The SED template %s requires an input %s!'%(template_name, 'z'))
                self.set_SED_template_Magdis2012(kwargs['z'])
                #--> this will update self.template_properties['normalization']
            # 
            elif template_name.lower() == 'dzliu':
                for t in ['z', 'Mstar', 'SFR']:
                    if not (t in kwargs.keys()):
                        raise ValueError('Error! The SED template %s requires an input %s!'%(template_name, t))
                self.set_SED_template_dzliu(**kwargs)
            # 
            #elif <TODO> more templates
    # 
    def set_SED_template_Magdis2012(self, z:float):
        SED_template_meta_table_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'SED_template_Magdis2012', 'meta_table_dliu.txt')
        SED_template_meta_table = Table.read(SED_template_meta_table_file, format='ascii')
        temp_matched_idx = np.argwhere(np.logical_and(SED_template_meta_table['z1']<z, SED_template_meta_table['z2']>z)).flatten().tolist()
        if len(temp_matched_idx) <= 0:
            raise ValueError('Error! The input redshift %s is out of the allowed range %s - %s for SED template of Magdis2012!'%(z, np.min(SED_template_meta_table['z1']), np.max(SED_template_meta_table['z2'])))
        SED_template_data_table_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'SED_template_Magdis2012', 'sed_z%s_U%s_radio.txt'%(SED_template_meta_table['i'][temp_matched_idx[0]], SED_template_meta_table['i'][temp_matched_idx[0]]))
        SED_template_data_table = Table.read(SED_template_data_table_file, format='ascii.no_header', names=['x','y'])
        SED_template_data_table['y'] = SED_template_data_table['y'] * (SED_template_data_table['x'])
        self.template = func_interp1d_loglog(SED_template_data_table['x'].data, SED_template_data_table['y'].data)
        # 
        if self.has_data():
            self.fit_SED_template(z)
    # 
    def set_SED_template_dzliu(self, **kwargs):
        if not (os.getenv('HOME')+'/Cloud/Github/Crab.Toolkit.michi2/bin' in sys.path):
            sys.path.append(os.getenv('HOME')+'/Cloud/Github/Crab.Toolkit.michi2/bin')
        z = kwargs['z']
        from michi2_make_galaxy_SED import michi2_make_galaxy_SED
        SED_dict = michi2_make_galaxy_SED(**kwargs, Silent=True)
        self.template = func_interp1d_loglog(SED_dict['w_SED']/(1.0+z), SED_dict['f_SED'])
        # 
        if self.has_data():
            self.fit_SED_template(z)
    # 
    def get_SED_model(self, input_wavelengths):
        # input_wavelengths in units of um, obs-frame
        if self.has_template():
            opz = 1.0
            norm = 1.0
            if 'z' in self.template_properties:
                opz = 1.0+self.template_properties['z']
            if 'normalization' in self.template_properties:
                norm = self.template_properties['normalization']
            return self.template(input_wavelengths / opz) * norm
    # 
    def has_data(self):
        return (len(self.w)>0)
    # 
    def has_template(self):
        return (self.template is not None)
    # 
    def has_model(self):
        return (self.template is not None)
    # 
    def is_valid(self):
        return (self.has_data() | self.has_template())
    # 
    def fit_SED_template(self, z:float):
        if self.has_data() and self.has_template():
            data_x = np.array(self.w)
            data_y = np.array(self.f)
            error_y = np.array(self.f_err)
            model_y = np.full(data_y.shape, np.nan)
            for k in range(len(data_x)):
                if not np.isnan(data_y[k]):
                    # TODO: filter
                    filter_x = np.logspace(np.log10(data_x[k])-0.2, np.log10(data_x[k])+0.2, num=100, endpoint=True) / (1.0+z) # convert to rest-frame wavelength
                    filter_y = np.array(filter_x * 0.0 + 1.0)
                    filter_y = filter_y / np.sum(filter_y)
                    template_y = self.template(filter_x) # rest-frame
                    if np.count_nonzero(~np.isnan(template_y)) > 0:
                        model_y[k] = np.sum(template_y * filter_y)
                        if np.isnan(error_y[k]): error_y[k] = data_y[k] * 0.1
                        #<DEBUG>#
                        #print('CrabGalaxySED::fit_SED_template()', 'data_x', data_x[k], 'data_y', data_y[k], 'model_y', model_y[k], 'S/N', data_y[k]/error_y[k])
            mask = np.logical_and(~np.isnan(data_y), ~np.isnan(model_y))
            if np.count_nonzero(mask) > 0:
                #template_normalization = np.average(data_y[mask] / model_y[mask], weights = (data_y[mask] / error_y[mask])**2)
                template_normalization = 10**(np.average(np.log10(data_y[mask] / model_y[mask]), weights = (data_y[mask] / error_y[mask])**2)) # average in log10
                self.template_properties['normalization'] = template_normalization
                #<DEBUG>#
                #print('CrabGalaxySED::fit_SED_template()', 'template_normalization', template_normalization)



def rotate_matrix_2D(angle = None):
    # angle is angle in units of degree to rotate, and increases along counter-clockwise direction. 
    # rotate matrix see https://en.wikipedia.org/wiki/Rotation_matrix
    if angle is None:
        raise ValueError('Error input of angle to the rotate_matrix_2D() function!')
    rot_mat = rotate_matrix_3D(axis = (0.0, 0.0, 1.0), angle = angle)
    rot_mat = rot_mat[0:2, 0:2]
    #print('rot_mat', rot_mat)
    #print('rot_mat.T', rot_mat.T)
    return rot_mat



def rotate_matrix_3D(axis = (0.0, 0.0, 1.0), angle = None):
    # angle is angle in units of degree about the axis (right-handed orientation), and increases along counter-clockwise direction. 
    # axis in default is +Z axis (toward us). +X is to the right and +Y is to the up (in camera view). 
    # rotate matrix see https://en.wikipedia.org/wiki/Rotation_matrix
    if axis is None or len(axis) != 3 or angle is None:
        raise ValueError('Error input of axis and angle to the rotate_matrix_3D() function!')
    cosAng = np.cos(np.deg2rad(angle))
    sinAng = np.sin(np.deg2rad(angle))
    axis = np.array(axis) / np.sqrt(np.sum(np.square(axis)))
    ux, uy, uz = axis[0:3]
    uxx, uyy, uzz = np.square(axis[0:3])
    uxy, uyz, uzx = ux*uy, uy*uz, uz*ux
    uyx, uzy, uxz = uxy, uyz, uzx
    rot_mat = [ [uxx*(1-cosAng)+cosAng,    uxy*(1-cosAng)-uz*sinAng, uxz*(1-cosAng)+uy*sinAng ], \
                [uyx*(1-cosAng)+uz*sinAng, uyy*(1-cosAng)+cosAng,    uyz*(1-cosAng)-ux*sinAng ], \
                [uzx*(1-cosAng)-uy*sinAng, uzy*(1-cosAng)+ux*sinAng, uzz*(1-cosAng)+cosAng    ] ]
    return np.array(rot_mat)



class CrabGalaxyMorphology(object):
    # 
    def __init__(self, shape:str = None, size = None, wavelength:float = None):
        # 
        self.allowed_shapes = ['gaussian', 'sersic']
        self.shape = None
        self.wavelength = None
        self.set_shape_and_size(shape, size)
        self.set_wavelength(wavelength)
    # 
    class gaussian_shape_class:
        def __init__(self, major:float, minor:float, angle:float):
            # major and minor axis FWHM in units of arcsec, note that Gaussian sigma = FWHM / (2.0*np.sqrt(2.0*np.log(2.0)))
            # angle is position angle in units of degree, zero means +Y direction, and increases along counter-clockwise direction. 
            self.major = major
            self.minor = minor
            self.angle = angle
            self.name = 'Gaussian'
        def func(self, dx, dy):
            rot_mat = rotate_matrix_2D(-(self.angle+90.0)) # this rotates the input dx dy from sky coordinate to a "source" coordinate where +X is the major axis and +Y is the minor axis. 
            if np.isscalar(dx): dx = [dx]
            if np.isscalar(dy): dy = [dy]
            if type(dx) is list: dx = np.array(dx)
            if type(dy) is list: dy = np.array(dy)
            ddx = np.full(dx.shape, np.nan)
            ddy = np.full(dy.shape, np.nan)
            for k in range(len(dx)):
                ddx[k], ddy[k] = np.matmul(rot_mat, [dx[k], dy[k]])
            #<20191028><BUGGY># sigx = self.major/(np.sqrt(2.0*np.log(2.0))) # this is equivalent to 2.0 * Gaussian_sigma
            #<20191028><BUGGY># sigy = self.minor/(np.sqrt(2.0*np.log(2.0))) # this is equivalent to 2.0 * Gaussian_sigma
            #<20191028><BUGGY># # output array is normalized to have a total area of 1.0 / pixscale**2
            #<20191028><BUGGY># return np.exp( - ( (ddx/sigx)**2 + \
            #<20191028><BUGGY>#                    (ddy/sigy)**2 ) ) / (np.pi*sigx*sigy)
            # 
            sigx = self.major/(2.0*np.sqrt(2.0*np.log(2.0)))
            sigy = self.minor/(2.0*np.sqrt(2.0*np.log(2.0)))
            # output array is normalized to have a total area of 1.0 / pixscale**2
            return np.exp( - ( (ddx**2/(2.0*sigx**2)) + \
                               (ddy**2/(2.0*sigy**2)) ) ) / (2.0*np.pi*sigx*sigy)
        def __str__(self):
            return 'Shape: %s, major: %s [arcsec], minor: %s [arcsec], PA: %s [degree].'%(\
                        self.name, self.major, self.minor, self.angle)
    # 
    class sersic_shape_class:
        def __init__(self, major:float, minor:float, angle:float, sersic_index:float):
            # major and minor axis effective radius in units of arcsec, within that ellipse the integrated flux is half of the total flux.
            # angle is position angle in units of degree, zero means +Y direction, and increases along counter-clockwise direction. 
            self.major = major
            self.minor = minor
            self.angle = angle
            self.index = sersic_index
            self.name = 'Sersic'
        def func(self, dx, dy):
            rot_mat = rotate_matrix_2D(-(self.angle+90.0)) # this rotates the input dx dy from sky coordinate to a "source" coordinate where +X is the major axis and +Y is the minor axis. 
            if np.isscalar(dx): dx = [dx]
            if np.isscalar(dy): dy = [dy]
            if type(dx) is list: dx = np.array(dx)
            if type(dy) is list: dy = np.array(dy)
            ddx = np.full(dx.shape, np.nan)
            ddy = np.full(dy.shape, np.nan)
            for k in range(len(dx)):
                ddx[k], ddy[k] = np.matmul(rot_mat, [dx[k], dy[k]])
            #bn = gammaincinv(2.0*self.index, 0.5) # Returns bn such that gammainc(2n, bn) = 0.5, which represents half-(total-)light radius R_eff .
            bn = gammaincinv(2.0*self.index, 0.5*gamma(2*self.index)) # Returns bn such that gammainc(2n, bn) = Gamma(2n), which represents half-(total-)light radius R_eff .
            #print(np.exp(bn)) # exp(bn) is the peak value. 
            # 
            # Area of Sersic2D: \int 2 pi R e^(-bn ((R/Re)^(1/n)-1) ) dR 
            # let x = bn (R/Re)^(1/n)
            # so R = (x^n / bn^n) * Re, and dR = n x^(n-1) Re / bn^n
            # so Area = \int 2 pi (x^n / bn^n) * Re * e^(-x) * e^(bn) * (n x^(n-1) Re / bn^n) dx
            #         = 2 pi (Re^2 / bn^2n) * e^(bn) * n * \int x^(2n-1) e^(-x) dx
            #         = 2 pi (Re^2 / bn^2n) * e^(bn) * n * \Gamma(2n)
            # 
            # Value at R_eff: 
            a = self.major
            b = self.minor
            #perimeter = np.pi * ( 3*(a+b) - np.sqrt((3*a+b)*(a+3*b)) )
            b_a_ratio = b/a # so perimeter = np.pi * ((3*(1+b_a_ratio)) - np.sqrt((3+b_a_ratio)*(1+3*b_a_ratio))) * a
            #perimeter_factor = ((3*(1+b_a_ratio)) - np.sqrt((3+b_a_ratio)*(1+3*b_a_ratio)))
            #perimeter = # a better way -- see ellipse perimeter -- https://www.mathsisfun.com/geometry/ellipse-perimeter.html
            # 
            #h = (a-b)**2/(a+b)**2
            #perimeter_factor = (1+b_a_ratio) * (1 + binom(0.5,1)*h + binom(0.5,2)*h + binom(0.5,3)*h)
            #I_int = perimeter_factor * np.pi * (self.major)**2 / (bn**(2*self.index)) * np.exp(bn) * self.index * gamma(2.0*self.index) # http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
            #print('I_int', I_int)
            # 
            #raise NotImplementedError('TODO: normalize sersic profile')
            # 
            # output array is normalized to have a total area of 1.0 / pixscale**2
            output_array = np.exp( - bn * ( np.sqrt( (ddx / (self.major))**2 + \
                                                     (ddy / (self.minor))**2 )**(1.0 / self.index) - 1.0 \
                                          ) )
            #print('sum', np.sum(output_array))
            output_array = output_array / np.sum(output_array)
            return output_array
        def __str__(self):
            return 'Shape: %s, major: %s [arcsec], minor: %s [arcsec], PA: %s [degree], Sersic index: %s.'%(\
                        self.name, self.major, self.minor, self.angle, self.index)
    # 
    def set_shape_and_size(self, input_shape:str, input_size):
        # 
        if input_shape is None and input_size is None:
            return
        # 
        if input_shape is None or input_size is None:
            raise ValueError('Error input of shape and size to the function CrabGalaxyMorphology::set_shape_and_size()! One is None!')
        # 
        if not (input_shape.lower() in self.allowed_shapes):
            raise ValueError('Error input of shape to the function CrabGalaxyMorphology::set_shape_and_size()! Must be one of %s'%(str(self.allowed_shapes)))
        # 
        if input_shape.lower() == 'gaussian':
            # if input a dict
            if type(input_size) is dict:
                # convert keys to lower cases
                input_keys = list(input_size.keys())
                for t in input_keys:
                    if t != t.lower():
                        input_size[t.lower()] = input_size[t]
                # copy some keys
                if 'pa' in input_size and not ('angle' in input_size):
                    input_size['angle'] = input_size['pa']
                # try to fill in gaussian major minor FWHM sizes and position angle
                if 'major' in input_size and 'minor' in input_size and 'angle' in input_size:
                    self.shape = self.gaussian_shape_class(input_size['major'], input_size['minor'], input_size['angle'])
                elif 'fwhm' in input_size:
                    self.shape = self.gaussian_shape_class(input_size['fwhm'], input_size['fwhm'], 0.0)
                elif 'reff' in input_size:
                    self.shape = self.gaussian_shape_class(input_size['reff']*2.43, input_size['reff']*2.43, 0.0)
                elif 'r_eff' in input_size:
                    self.shape = self.gaussian_shape_class(input_size['r_eff']*2.43, input_size['r_eff']*2.43, 0.0)
                else:
                    raise ValueError('Error input of size to the function CrabGalaxyMorphology::set_shape_and_size()! For Gaussian it can be a dict with \'major\', \'minor\' and \'PA\', or simply one \'FWHM\' value.')
            # else if input only one float number as the size, then the shape is a circular gaussian
            elif np.isscalar(input_size):
                self.shape = self.gaussian_shape_class(input_size, input_size, 0)
            # otherwise try to fill in gaussian major minor FWHM sizes and position angle
            else:
                try:
                    self.shape = self.gaussian_shape_class(*input_size)
                except:
                    raise ValueError('Error input of size to the function CrabGalaxyMorphology::set_shape_and_size()! For Gaussian it should contain major axis FWHM, minor axis FWHM and position angle, in the order as stated.')
        # 
        elif input_shape.lower() == 'sersic':
            # if input a dict
            if type(input_size) is dict:
                # convert keys to lower cases
                input_keys = list(input_size.keys())
                for t in input_keys:
                    if t != t.lower():
                        input_size[t.lower()] = input_size[t]
                # copy some keys
                if 'pa' in input_size and not ('angle' in input_size):
                    input_size['angle'] = input_size['pa']
                # try to fill in sersic major and minor effective raidii and position angle and sersic index
                if 'major' in input_size and 'minor' in input_size and 'angle' in input_size and 'n' in input_size:
                    self.shape = self.sersic_shape_class(input_size['major'], input_size['minor'], input_size['angle'], input_size['n'])
                elif 'reff' in input_size and 'n' in input_size:
                    self.shape = self.sersic_shape_class(input_size['reff'], input_size['reff'], 0.0, input_size['n'])
                elif 'r_eff' in input_size and 'n' in input_size:
                    self.shape = self.sersic_shape_class(input_size['r_eff'], input_size['reff'], 0.0, input_size['n'])
                else:
                    raise ValueError('Error input of size to the function CrabGalaxyMorphology::set_shape_and_size()! For Sersic it can be a dict with \'major\', \'minor\', \'PA\' and \'n\', or simply \'reff\' and \'n\'.')
            # otherwise try to fill in sersic major and minor effective raidii and position angle and sersic index
            else:
                try:
                    self.shape = self.sersic_shape_class(*input_size)
                except:
                    raise ValueError('Error input of size to the function CrabGalaxyMorphology::set_shape_and_size()! For Sersic it should contain major axis effective radius, minor axis effective radius, position angle and sersic index, in the order as stated.')
    # 
    def set_wavelength(self, wavelength:float):
        if wavelength is not None:
            self.wavelength = wavelength
    # 
    def plot_cutout(self, cutout_size = 5.0, pixel_scale = 1.0):
        if self.shape is not None:
            #cutout_size = 5.0 # arcsec
            #pixel_scale = 1.0 # arcsec/pixel
            if np.isscalar(cutout_size): cutout_size = np.array([cutout_size, cutout_size])
            if len(cutout_size) < 2: raise ValueError('Error input of cutout_size to the function CrabGalaxyMorphology::plot_cutout()! It should be a list of two elements!')
            if type(cutout_size) is list: cutout_size = np.array(cutout_size)
            gsize = (np.ceil(cutout_size/pixel_scale/2.0)*2+1).astype(int)
            gy, gx = np.mgrid[0:gsize[1], 0:gsize[0]] #<NOTE><BUGGY># np.mgrid needs an order of Y,X !!
            dx = gx - (gsize[0]-1)/2
            dy = gy - (gsize[1]-1)/2
            data_2D = self.shape.func(dx * pixel_scale, dy * pixel_scale)
            print('sum(data_2D)', np.sum(data_2D)*pixel_scale*pixel_scale)
            print('max(data_2D)', np.max(data_2D))
            plt.imshow(data_2D, origin='lower', interpolation='nearest')
            plt.contour(data_2D, colors = 'white', alpha = 0.7, linewidths = 1.2) # levels=np.array([-3.,-2,2,3,4,5])*image_rms, 
            plt.show()
    # 
    def resample_mgrid(self, input_mgrid_y, input_mgrid_x, oversampling_factor):
        if len(input_mgrid_y.shape) != 2 or len(input_mgrid_x.shape) != 2:
            raise ValueError('Error! The input mgrid array of resample_mgrid() is not a 2D mgrid array!')
        #output_mgrid_y = np.full(np.array(list(input_mgrid_y.shape))*oversampling_factor, 0.0)
        #output_mgrid_x = np.full(np.array(list(input_mgrid_x.shape))*oversampling_factor, 0.0)
        # e.g., mgrid[351:354, 401:404]
        # mgrid_x = [ [401, 402, 403], 
        #             [401, 402, 403], 
        #             [401, 402, 403] ]
        # mgrid_y = [ [351, 351, 351], 
        #             [352, 352, 352], 
        #             [353, 353, 353] ]
        # after resampling by a factor of 2:
        # mgrid_x = [ [401.0, XXX.X, XXX.X, XXX.X, XXX.X, 403.0], 
        #             [401.0, XXX.X, XXX.X, XXX.X, XXX.X, 403.0], 
        #             [401.0, XXX.X, XXX.X, XXX.X, XXX.X, 403.0], 
        #             [401.0, XXX.X, XXX.X, XXX.X, XXX.X, 403.0], 
        #             [401.0, XXX.X, XXX.X, XXX.X, XXX.X, 403.0], 
        #             [401.0, XXX.X, XXX.X, XXX.X, XXX.X, 403.0] ]
        # mgrid_y = [ [351.0, 351.0, 351.0, 351.0, 351.0, 351.0], 
        #             [XXX.X, XXX.X, XXX.X, XXX.X, XXX.X, XXX.X], 
        #             [XXX.X, XXX.X, XXX.X, XXX.X, XXX.X, XXX.X], 
        #             [XXX.X, XXX.X, XXX.X, XXX.X, XXX.X, XXX.X], 
        #             [XXX.X, XXX.X, XXX.X, XXX.X, XXX.X, XXX.X], 
        #             [353.0, 353.0, 353.0, 353.0, 353.0, 353.0] ]
        # 
        #size_y, size_x = input_mgrid_x.shape
        #output_mgrid_x0 = np.interp(np.arange(size_x*oversampling_factor)/float(size_x*oversampling_factor-1)*float(size_x-1), 
        #                               np.arange(size_x), 
        #                               input_mgrid_x[0,:])
        #output_mgrid_x = np.repeat(output_mgrid_x0[np.newaxis,:], size_y*oversampling_factor, axis=0)
        ## 
        #output_mgrid_y0 = np.interp(np.arange(size_y*oversampling_factor)/float(size_y*oversampling_factor-1)*float(size_y-1), 
        #                               np.arange(size_y), 
        #                               input_mgrid_y[:,0])
        #output_mgrid_y = np.repeat(output_mgrid_y0[:,np.newaxis], size_x*oversampling_factor, axis=1)
        # 
        # ---- NONONO
        # ---- For pixel coordinates, what we need to do is to resample fractional pixels, so that [1,2,3] --> [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]. 
        #      Each mgrid coordinate (0-based) is the lower-left corner of each pixel, not the pixel center!
        # so if we have
        # e.g., mgrid[351:354, 401:404]
        # mgrid_x = [ [401, 402, 403], 
        #             [401, 402, 403], 
        #             [401, 402, 403] ]
        # mgrid_y = [ [351, 351, 351], 
        #             [352, 352, 352], 
        #             [353, 353, 353] ]
        # after resampling by a factor of 2:
        # mgrid_x = [ [401.0, XXX.X, 402.0, XXX.X, 403.0, XXX.X], 
        #             [401.0, XXX.X, 402.0, XXX.X, 403.0, XXX.X], 
        #             [401.0, XXX.X, 402.0, XXX.X, 403.0, XXX.X], 
        #             [401.0, XXX.X, 402.0, XXX.X, 403.0, XXX.X], 
        #             [401.0, XXX.X, 402.0, XXX.X, 403.0, XXX.X], 
        #             [401.0, XXX.X, 402.0, XXX.X, 403.0, XXX.X] ]
        # mgrid_y = [ [351.0, 351.0, 351.0, 351.0, 351.0, 351.0], 
        #             [XXX.X, XXX.X, XXX.X, XXX.X, XXX.X, XXX.X], 
        #             [352.0, 352.0, 352.0, 352.0, 352.0, 352.0], 
        #             [XXX.X, XXX.X, XXX.X, XXX.X, XXX.X, XXX.X], 
        #             [353.0, 353.0, 353.0, 353.0, 353.0, 353.0], 
        #             [XXX.X, XXX.X, XXX.X, XXX.X, XXX.X, XXX.X] ]
        size_y, size_x = input_mgrid_x.shape
        output_mgrid_y = np.full(np.array(list(input_mgrid_y.shape))*oversampling_factor, 0.0)
        output_mgrid_x = np.full(np.array(list(input_mgrid_x.shape))*oversampling_factor, 0.0)
        output_mgrid_x0 = np.interp(np.arange(size_x*oversampling_factor)/float(oversampling_factor), 
                                       np.arange(size_x+1), 
                                       np.concatenate([input_mgrid_x[0,:], [input_mgrid_x[0,-1]+1]]) )
        output_mgrid_x = np.repeat(output_mgrid_x0[np.newaxis,:], size_y*oversampling_factor, axis=0)
        # 
        output_mgrid_y0 = np.interp(np.arange(size_y*oversampling_factor)/float(oversampling_factor), 
                                       np.arange(size_y+1), 
                                       np.concatenate([input_mgrid_y[:,0], [input_mgrid_y[-1,0]+1]]) )
        output_mgrid_y = np.repeat(output_mgrid_y0[:,np.newaxis], size_x*oversampling_factor, axis=1)
        # 
        return output_mgrid_y, output_mgrid_x
    # 
    def inject_into_image(self, input_image_data, pixel_scale = 1.0, 
                                flux = None, band = None, 
                                x = None, y = None, 
                                mgrid_x = None, mgrid_y = None, 
                                convolving_image = None, 
                                convolving_beam = None, 
                                oversampling_factor = None, 
                                convolving_image_oversampled = None):
        if self.shape is not None:
            # 
            # user can supply mgrid_x and mgrid_y to speed up
            if mgrid_x is None or mgrid_y is None:
                image_size_y, image_size_x = input_image_data.shape # N_X, N_Y
                mgrid_y, mgrid_x = np.mgrid[0:image_size_y, 0:image_size_x]
            # 
            # user can supply a string for 'band' (representing the SED band name), or alternatively a float number for 'flux'
            if flux is None:
                if band is not None:
                    if self.SED.is_valid():
                        if band in self.SED.band_name:
                            flux = self.SED.f[self.SED.band_name.index(band)]
                        else:
                            raise ValueError('The input band %s was not found in the galaxy SED band names (%s) when calling inject_into_image()!'%(band, ', '.join(self.SED.band_name)))
                    else:
                        raise ValueError('The user has input band %s when calling inject_into_image(), but the galaxy has no valid SED!'%(band))
                else:
                    raise ValueError('Please input a flux or band when calling inject_into_image()!')
            # 
            # check user input x y
            if x is None or y is None:
                raise ValueError('Please input source pixel coordinate x and y (0-based) when calling inject_into_image()!') #<TODO># 
            # 
            # user can also supply a convolving_image so that we will do convolution
            if convolving_beam is None and convolving_image is None and convolving_image_oversampled is None:
                # no need to do convolution
                #print('[DEBUG] inject_into_image(): no convolution')
                # 
                if oversampling_factor is not None:
                    if oversampling_factor < 0:
                        print('Warning! inject_into_image(): the oversampling factor is %d! We will simply align the source at %s,%s to the nearest pixel center!'%(oversampling_factor, x, y))
                        # if the user has given a negative oversampling_factor, then we actually do not oversample but take the nearest pixel
                        x = np.round(x)
                        y = np.round(y)
                # 
                data_2D = self.shape.func((mgrid_x-x) * pixel_scale, (mgrid_y-y) * pixel_scale)
            else:
                # do convolution
                #print('[DEBUG] inject_into_image(): doing convolution')
                # 
                # check user input convolving_image, if it is a scalar, then we take it as the FWHM of the PSF
                if convolving_image is None:
                    if convolving_beam is not None:
                        PSF_FWHM_arcsec = float(convolving_beam)
                        PSF_FWHM_pixel = PSF_FWHM_arcsec / pixel_scale
                        PSF_sigma_pixel = PSF_FWHM_pixel / (2.0*np.sqrt(2.0*np.log(2.0)))
                        convolving_image = Gaussian2DKernel(PSF_sigma_pixel) # default size is 8.0 * sigma
                # 
                # check if we need to oversample the source image for fractional pixel or not
                if self.shape.minor/pixel_scale < 2.0:
                    # 
                    # if source size is smaller than 2.0 pixel, we oversample the image
                    if oversampling_factor is None:
                        oversampling_factor = int(np.ceil(2.0/(self.shape.minor/pixel_scale)))
                        # check whether user has input convolving_image_oversampled but not oversampling_factor
                        if convolving_image_oversampled is not None:
                            raise ValueError('Error! The user has input convolving_image_oversampled but not oversampling_factor when calling inject_into_image()!')
                    else:
                        oversampling_factor = int(oversampling_factor)
                    #print('[DEBUG] inject_into_image(): oversampling_factor', oversampling_factor)
                    # 
                    if oversampling_factor < 0:
                        print('Warning! inject_into_image(): the oversampling factor is %d! We will simply align the source at %s,%s to the nearest pixel center!'%(oversampling_factor, x, y))
                        # if the user has given a negative oversampling_factor, then we actually do not oversample but take the nearest pixel
                        x = np.round(x)
                        y = np.round(y)
                        data_2D = self.shape.func((mgrid_x-x) * pixel_scale, (mgrid_y-y) * pixel_scale)
                        data_2D = convolve(data_2D, convolving_image)
                    else:
                        if oversampling_factor > 10:
                            print('Warning! inject_into_image(): the oversampling factor of %d is large! Convolution will take a lot of time!'%(oversampling_factor))
                        # 
                        # user can also supply a convolving_image_oversampled so that we do not need to resample convolving_image
                        if convolving_image_oversampled is None:
                            if convolving_beam is not None:
                                PSF_FWHM_arcsec = float(convolving_beam)
                                PSF_FWHM_pixel = PSF_FWHM_arcsec / pixel_scale * oversampling_factor
                                PSF_sigma_pixel = PSF_FWHM_pixel / (2.0*np.sqrt(2.0*np.log(2.0)))
                                convolving_image_oversampled = Gaussian2DKernel(PSF_sigma_pixel) # default size is 8.0 * sigma
                                #print('sum(convolving_image)', np.sum(convolving_image.array))
                                #print('sum(convolving_image_oversampled)', np.sum(convolving_image_oversampled.array))
                            else:
                                convolving_image_oversampling_func = interp2d(np.linspace(0, 1, convolving_image.shape[1], endpoint=True), \
                                                                              np.linspace(0, 1, convolving_image.shape[0], endpoint=True), \
                                                                              convolving_image )
                                convolving_image_oversampled = convolving_image_oversampling_func(np.linspace(0, 1, (int(np.ceil(convolving_image.shape[1]*oversampling_factor/2))*2+1), endpoint=True), \
                                                                                                  np.linspace(0, 1, (int(np.ceil(convolving_image.shape[0]*oversampling_factor/2))*2+1), endpoint=True) )
                                convolving_image_oversampled = convolving_image_oversampled / np.sum(convolving_image_oversampled) * np.sum(convolving_image.array)
                                #print('type(convolving_image)', type(convolving_image))
                                #print('sum(convolving_image)', np.sum(convolving_image.array))
                                #print('sum(convolving_image_oversampled)', np.sum(convolving_image_oversampled))
                        # 
                        # debug
                        #print('\rWriting to "dump_convolving_image_oversampled.fits" ...')
                        #dump_hdu = fits.PrimaryHDU(convolving_image_oversampled)
                        #dump_hdu.writeto('dump_convolving_image_oversampled.fits', overwrite=True)
                        #print('\rWritten to "dump_convolving_image_oversampled.fits"')
                        # debug
                        #print('\rWriting to "dump_convolving_image.fits" ...')
                        #dump_hdu = fits.PrimaryHDU(convolving_image)
                        #dump_hdu.writeto('dump_convolving_image.fits', overwrite=True)
                        #print('\rWritten to "dump_convolving_image.fits"')
                        # 
                        # 
                        # first oversample the grid, then generate source shape, then convolve, then resample back. 
                        data_2D_oversampled = np.full(np.array(list(input_image_data.shape))*oversampling_factor, 0.0)
                        mgrid_y_oversampled, mgrid_x_oversampled = self.resample_mgrid(mgrid_y, mgrid_x, oversampling_factor)
                        data_2D_oversampled = self.shape.func((mgrid_x_oversampled-x) * pixel_scale, (mgrid_y_oversampled-y) * pixel_scale)
                        #print('np.sum(data_2D_oversampled)', np.sum(data_2D_oversampled))
                        data_2D_oversampled = convolve(data_2D_oversampled, convolving_image_oversampled)
                        #print('np.sum(data_2D_oversampled)', np.sum(data_2D_oversampled))
                        data_2D = data_2D_oversampled[::oversampling_factor, ::oversampling_factor]
                        #print('np.sum(data_2D)', np.sum(data_2D))
                        data_2D = data_2D * float(oversampling_factor)**2
                        #print('np.sum(data_2D)', np.sum(data_2D))
                else:
                    # else if the source size is large enough, we do not need oversampling. 
                    #data_2D = np.full([input_image_data.shape[0], input_image_data.shape[1]], 0.0)
                    data_2D = self.shape.func((mgrid_x-x) * pixel_scale, (mgrid_y-y) * pixel_scale)
                    data_2D = convolve(data_2D, convolving_image)
                # 
            #print('flux', flux, 'sum(data_2D)', np.sum(data_2D))
            return input_image_data + flux * data_2D
        # 
        return None
    # 
    def is_valid(self):
        return (self.shape is not None)




class CrabGalaxy(object):
    # 
    def __init__(self, z:float = None, dL:float = None, SED = None, shape = None, size = None):
        self.z = None
        self.dL = None
        self.SED = CrabGalaxySED()
        self.Morph = CrabGalaxyMorphology()
        self.set_z(z)
        self.set_dL(dL)
        self.set_SED(SED)
        self.set_Morph(shape, size)
    # 
    def set_z(self, z):
        if z is None:
            pass
        else:
            self.z = float(z)
    # 
    def set_dL(self, dL):
        if dL is None:
            pass
        else:
            self.dL = float(dL)
    # 
    def set_SED(self, SED):
        if SED is None:
            pass
        else:
            if type(SED) is CrabGalaxySED:
                self.SED = SED
            # 
            elif type(SED) is dict:
                if 'w' in SED and 'f' in SED:
                    temp_w = SED['w']
                    temp_f = SED['f']
                    temp_f_err = None
                    temp_f_unit = None
                    temp_band_name = None
                    if 'f_err' in SED:
                        temp_f_err = SED['f_err']
                    if 'f_unit' in SED:
                        temp_f_unit = SED['f_unit']
                    if 'band_name' in SED:
                        temp_band_name = SED['band_name']
                    self.SED.set_data(w = temp_w, f = temp_f, f_err = temp_f_err, f_unit = temp_f_unit, band_name = temp_band_name)
                else:
                    raise ValueError('Error input of SED! If it is a dict, it should contain at least \'w\' and \'f\', and optionally \'f_err\' and \'band_name\'.')
            # 
            elif type(SED) is Table:
                if 'w' in SED.colnames and 'f' in SED.colnames:
                    temp_w = SED['w'].data.tolist()
                    temp_f = SED['f'].data.tolist()
                    temp_f_err = None
                    temp_band_name = None
                    if 'f_err' in SED:
                        temp_f_err = SED['f_err'].data.tolist()
                    if 'f_unit' in SED:
                        temp_f_unit = SED['f_unit'].data.tolist()
                    if 'band_name' in SED:
                        temp_band_name = SED['band_name'].data.tolist()
                    self.SED.set_data(w = temp_w, f = temp_f, f_err = temp_f_err, f_unit = temp_f_unit, band_name = temp_band_name)
                else:
                    raise ValueError('Error input of SED! If it is a dict, it should contain at least \'w\' and \'f\', and optionally \'f_err\' and \'band_name\'.')
            # 
            else:
                try:
                    temp_w = SED[0]
                    temp_f = SED[1]
                    temp_f_err = None
                    temp_f_unit = None
                    temp_band_name = None
                    try: 
                        temp_f_err = SED[2]
                    except:
                        pass
                    try: 
                        temp_f_unit = SED[3]
                    except:
                        pass
                    try: 
                        temp_band_name = SED[4]
                    except:
                        pass
                    self.SED.set_data(w = temp_w, f = temp_f, f_err = temp_f_err, f_unit = temp_f_unit, band_name = temp_band_name)
                except:
                    raise ValueError('Error input of SED! It seems the input is a list but it does not contain enough data? It is better to input a astropy.table.Table or dict with at least \'w\' and \'f\' columns, and optionally \'f_err\' and \'band_name\'.')
    # 
    def set_Morph(self, input_shape:str, input_size):
        if input_shape is None and input_size is None:
            pass
        else:
            self.Morph.set_shape_and_size(input_shape, input_size)
    # 
    def set_SED_template(self, input_name:str, **kwargs):
        if input_name is not None:
            self.SED.set_SED_template(input_name, **kwargs)
                
    # 
    def plot_cutout(self, cutout_size = 5.0, pixel_scale = 1.0):
        if self.Morph.is_valid():
            self.Morph.plot_cutout(cutout_size, pixel_scale)
    # 
    def print_galaxy_size(self):
        if self.Morph.is_valid():
            print(self.Morph.shape)
    # 
    def print_galaxy_shape(self):
        if self.Morph.is_valid():
            print(self.Morph.shape)
    # 
    def get_galaxy_size(self):
        if self.Morph.is_valid():
            return self.Morph.shape.major, self.Morph.shape.minor
    # 
    def plot_SED(self):
        if self.SED.is_valid():
            ax = plt.gca()
            if self.SED.has_data():
                ax.plot(self.SED.w, self.SED.f, ls='none', marker='o', markersize=3)
            if self.SED.has_model():
                gx = np.logspace(np.log10(0.1), np.log10(3e6), num=1000, endpoint=True)
                SED_x = gx
                if self.z is not None:
                    SED_x = gx / (1.0+self.z)
                SED_y = self.SED.get_SED_model(SED_x)
                ax.plot(gx, SED_y, ls='solid', lw=1.2)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim([0.07, 3e5])
            ax.set_ylim([1e-6, 1e3])
            ax.set_xlabel(r'Wavelength [$\mu$m]')
            ax.set_ylabel(r'Flux Density [mJy]') #<TODO># f_unit
            # 
            plt.show()
    # 
    def inject_into_image(self, input_image_data, pixel_scale, 
                                flux = None, band = None, 
                                x = None, y = None, 
                                mgrid_x = None, mgrid_y = None, 
                                convolving_beam = None, 
                                convolving_image = None, 
                                oversampling_factor = None, 
                                convolving_image_oversampled = None):
        if self.Morph.is_valid():
            return self.Morph.inject_into_image(input_image_data, pixel_scale, 
                                                flux = flux, band = band, 
                                                x = x, y = y, 
                                                mgrid_x = mgrid_x, mgrid_y = mgrid_y, 
                                                convolving_beam = convolving_beam, 
                                                convolving_image = convolving_image, 
                                                oversampling_factor = oversampling_factor, 
                                                convolving_image_oversampled = convolving_image_oversampled)
        return None








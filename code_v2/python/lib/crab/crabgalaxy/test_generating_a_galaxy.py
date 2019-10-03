#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# 

import os, sys, re
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from CrabGalaxy import CrabGalaxy

#my_galaxy = CrabGalaxy(shape = 'gaussian', size = [1.0, 0.5, 20.0])
#my_galaxy = CrabGalaxy(shape = 'gaussian', size = {'major':1.0, 'minor':0.5, 'PA':20.0})
#my_galaxy = CrabGalaxy(shape = 'sersic', size = {'major':1.0, 'minor':0.5, 'PA':20.0, 'n':0.5})
my_galaxy = CrabGalaxy(shape = 'sersic', size = {'major':1.0, 'minor':0.5, 'PA':20.0, 'n':0.5}, 
                       SED = {'w':[3.6, 4.5, 500], 'f':[1e-3, 1.2e-3, 5]})

#my_galaxy.plot_cutout(cutout_size = 3.0, pixel_scale = 0.1)

#my_galaxy.set_SED_template('Magdis2012', z = 3.0)
my_galaxy.set_SED_template('dzliu', z = 3.0, Mstar = 1e10, SFR = 300, EBV = 0.5)
#my_galaxy.set_SED_template('dzliu', z = 3.0, Mstar = 1e10, SFR = 300, EBV = 0.5, Age = 500e6, tau = 1)

#my_galaxy.plot_SED()

my_galaxy.print_galaxy_shape()

major, minor = my_galaxy.get_galaxy_size()





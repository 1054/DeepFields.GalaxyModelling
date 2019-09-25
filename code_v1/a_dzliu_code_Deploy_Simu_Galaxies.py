#!/usr/bin/env python2.7
# 

import os, sys

sys.path.append(os.path.expanduser('~')+'Cloud/Github/Crab.Toolkit.CAAP/bin')

from caap_google_drive_operator.py import CAAP_Google_Drive_Operator



# 
# Main Code
# 

gdo = CAAP_Google_Drive_Operator()
files = gdo.get_file_by_name('Simu*/Cosmological_Galaxy_Modelling_for_COSMOS/*/result_simu_galaxies.tar.gz')
if len(files)>0:
    gdo.download_files(files[0])





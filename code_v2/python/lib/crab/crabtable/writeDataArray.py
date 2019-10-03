#!/usr/local/bin/python
# -*- coding: utf-8 -*-
# 
# 
# This script will write data array to file
# 
# Usage:
#   execfile('/home/dzliu/Software/ipy/setup.py')
#   TempDataArray = asciitable.read('Arp193_FromHIPE_SLW.CSV')
#   crabtable.writeDataArray(TempDataArray, "output.txt")
# 
# Last update:
#   2015-11-17 created
# 


def writeDataArray(DataArray, FilePath, ColumnNames=[]):
    if len(DataArray)>0 and len(FilePath)>0:
        if len(ColumnNames)>0: 
            asciitable.write(DataArray[ColumnNames], FilePath, Writer=asciitable.FixedWidthTwoLine, delimiter_pad=' ', position_char=' ')
        else:
            asciitable.write(DataArray, FilePath, Writer=asciitable.FixedWidthTwoLine, delimiter_pad=' ', position_char=' ')
        # add the '#' mark
        os.system("sed -i -e '1s/^ /#/g' %s"%(FilePath))
        os.system("sed -i -e '2s/^ /#/g' %s"%(FilePath))
        print("crabtable.writeDataArray: Written to %s!"%(FilePath))
    pass



# the dependent importing should be put at the end
# http://effbot.org/zone/import-confusion.htm
try:
    import astropy.io.ascii as asciitable
except ImportError: 
    print "Error! Could not import astropy.io.ascii as asciitable!"
try:
    import glob
except ImportError: 
    print "Error! Could not import glob!"
try:
    import re
except ImportError: 
    print "Error! Could not import re!"
try:
    import numpy
except ImportError: 
    print "Error! Could not import numpy!"
try:
    import numpy.lib.recfunctions as rec
except ImportError: 
    print "Error! Could not import numpy.lib.recfunctions!"
try:
    import os
except ImportError: 
    print "Error! Could not import os!"




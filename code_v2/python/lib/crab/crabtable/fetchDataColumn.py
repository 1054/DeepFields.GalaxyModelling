#!/usr/local/bin/python
# -*- coding: utf-8 -*-
# 
# 
# This script will read data column(s) from data table(s)
# 
# input a list of data table file path
# input a list of data table header name
# output an array of data columns
# 
# Usage:
#   execfile('/Users/dliu/Software/ipy/setup.py')
#   data = crabtable.fetchDataColumn('Arp193_FromHIPE_SLW.CSV','LineName')
#   
# Example: 
#   # get all LineName=[CI](3P_2-3P_1) 
# 
# Last update:
#   2015-03-20 created
#  


def fetchDataColumn(DataTableFilePathList,DataTableHeaderNameList,Verbose=True):
    # the list of data table file path
    TempDataTableFilePathList = []
    if type(DataTableFilePathList) is str:
        TempDataTableFilePathList = glob.glob(DataTableFilePathList)
    elif type(DataTableFilePathList) is list:
        for TempDataTableFilePath in DataTableFilePathList:
            if len(TempDataTableFilePathList)>0:
                TempDataTableFilePathList.extend(sorted(glob.glob(TempDataTableFilePath)))
            else:
                TempDataTableFilePathList = sorted(glob.glob(TempDataTableFilePath))
        # TempDataTableFilePathList = DataTableFilePathList
    
    # the data table header name to fetch (must be exactly the same) (supports utf-8 characters)
    TempTableHeaderNameList = []
    if type(DataTableHeaderNameList) is str:
        TempTableHeaderNameList = [DataTableHeaderNameList]
    elif type(DataTableHeaderNameList) is list:
        TempTableHeaderNameList = DataTableHeaderNameList
    
    # fecth data array
    FetchedDataArray = []
    for TempDataFile in TempDataTableFilePathList:
        if Verbose: print "fetchDataColumn: Reading %s"%(TempDataFile)
        TempDataArray = asciitable.read(TempDataFile)
        if Verbose: print "fetchDataColumn: Read %d rows"%(len(TempDataArray))
        if len(FetchedDataArray) == 0:
            FetchedDataArray = TempDataArray[TempTableHeaderNameList]
        else:
            # FetchedDataArray.append(TempDataArray[TempTableHeaderNameList])
            FetchedDataArray = rec.stack_arrays((FetchedDataArray,TempDataArray[TempTableHeaderNameList]),usemask=False,asrecarray=True,autoconvert=True)
        # for np recarray, just use [names] to select columns
        # http://stackoverflow.com/questions/15575878/how-do-you-remove-a-column-from-a-structured-numpy-array
    return FetchedDataArray



# the dependent importing should be put at the end?
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




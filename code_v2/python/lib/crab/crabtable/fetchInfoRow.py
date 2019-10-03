#!/usr/local/bin/python
# -*- coding: utf-8 -*-
# 
# 
# This script will read info row(s) from data table(s)
# 
# input a list of data table file path
# input a list of data table row string pattern
# output an array of info row(s)
# 
# Usage:
#   execfile('/Users/dliu/Software/ipy/setup.py')
#   data = crabtable.fetchInfoRow('/Users/dliu/Software/ipy/crabtable/test/3C084_FromHIPE_S?W.CSV',
#                                 '^#.*(CO[(][0-9-]*[)]) *at ([0-9.]*) GHz.* ([0-9.]*) arcsec.* ([0-9.]*) sigma.* [=] ([0-9.]*) Jy or flux [=] ([0-9.]*) Jy .*',
#                                 OutputHeaderList=['LineName','LineFreqObs','BeamSize','LineSigma','LinePeak','LineFPrm'])
#   
#   TempDataArray = asciitable.read('Arp193_FromHIPE_SLW.CSV')
# 
# Last update:
#   2015-03-20 created
# 


def fetchInfoRow(DataTableFilePathList,MatchStringPatternList,OutputHeaderList=[],OutputDataTablePattern=[],Verbose=True,Debug=False):
    # the list of data table file path
    TempDataTableFilePathList = []
    if type(DataTableFilePathList) is str:
        TempDataTableFilePathList = sorted(glob.glob(DataTableFilePathList))
    elif type(DataTableFilePathList) is list:
        for TempDataTableFilePath in DataTableFilePathList:
            if len(TempDataTableFilePathList)>0:
                TempDataTableFilePathList.extend(sorted(glob.glob(TempDataTableFilePath)))
            else:
                TempDataTableFilePathList = sorted(glob.glob(TempDataTableFilePath))
            # TempDataTableFilePathList = DataTableFilePathList #<DONE># foreach glob.glob(DataTableFilePathList)
    # <DEBUG> print TempDataTableFilePathList
    # <DEBUG> return []
    
    # the string pattern to search using regular expression (regex)
    TempStringPatternList = []
    if type(MatchStringPatternList) is str:
        TempStringPatternList = [MatchStringPatternList]
    elif type(MatchStringPatternList) is list:
        TempStringPatternList = MatchStringPatternList
    
    # convert string pattern to regex pattern variables
    TempRegexPatternList = []
    for TempStringPattern in TempStringPatternList:
        TempRegexPattern = re.compile(TempStringPattern)
        TempRegexPatternList.append(TempRegexPattern)
    
    # the string pattern to apply when output DataTable
    OutputDataTableStringPatternList = []
    if type(OutputDataTablePattern) is str:
        OutputDataTableStringPatternList = [OutputDataTablePattern]
    elif type(OutputDataTablePattern) is list:
        OutputDataTableStringPatternList = OutputDataTablePattern
    
    # convert string pattern to apply when output DataTable
    OutputDataTableRegexPatternList = []
    for OutputDataTableStringPattern in OutputDataTableStringPatternList:
        OutputDataTableRegexPattern = re.compile(OutputDataTableStringPattern)
        OutputDataTableRegexPatternList.append(OutputDataTableRegexPattern)
    
    # fecth info lines (info lines is every text things including commented things)
    FetchedDataArray = []
    for TempDataFile in TempDataTableFilePathList:
        if Verbose: print "fetchInfoRow: Reading %s\n"%(TempDataFile)
        with open(TempDataFile) as TempDataFPTR:
            TempInfoLines = TempDataFPTR.readlines()
        for TempInfoLine in TempInfoLines:
            TempRowMatched = 0
            for TempRegexPattern in TempRegexPatternList:
                if TempRegexPattern.match(TempInfoLine) != None:
                    if Verbose: print "fetchInfoRow: Matched %s to Row %s in %s\n"%(TempRegexPattern.pattern,TempInfoLine,TempDataFile)
                    TempRowMatched = 1
                    TempRegexSearch = TempRegexPattern.search(TempInfoLine)
                    TempRowContent = TempRegexSearch.groups()
                    # if len(TempRegexSearch.groups())>0:
                    #     TempRowContent = TempRegexSearch.groups()
                    # else:
                    #     #<TODO># TempRowContent = TempRegexSearch.group() #<fixed><20150413><dzliu>#
                    if Debug:
                        print("TempRowContent")
                        print( TempRowContent )
                    break
            if TempRowMatched == 1:
                try:
                    if len(OutputDataTableRegexPatternList)>0:
                        for OutputDataTableRegexPattern in OutputDataTableRegexPatternList:
                            TempRegexSearch = OutputDataTableRegexPattern.search(TempDataFile)
                            TempDataFile = '-'.join(TempRegexSearch.groups()).replace('_','-') #<TODO># join file name by '-'
                    TempNewContent = rec.append_fields(TempRowContent,'DataTable',[TempDataFile],usemask=False,asrecarray=True)
                    # <TODO> rec.append_fields(TempRowContent,'DataTable',replicate(TempDataFile,len(TempRowContent),,)
                except:
                    TempNewContent = TempRowContent
                
                if len(FetchedDataArray) == 0:
                    FetchedDataArray = TempNewContent
                    del(TempNewContent)
                else:
                    TempOldContent = FetchedDataArray
                    FetchedDataArray = rec.stack_arrays((TempOldContent,TempNewContent),usemask=False,asrecarray=True,autoconvert=True)
                    del(TempOldContent)
                    del(TempNewContent)
        # for np recarray, just use [names] to select columns
        # http://stackoverflow.com/questions/15575878/how-do-you-remove-a-column-from-a-structured-numpy-array
    # 
    if len(FetchedDataArray)>0:
        # now try to assign header list
        if type(OutputHeaderList) is list:
            NewHeaderList = []
            for TempIndex in range(len(FetchedDataArray.dtype.names)):
                if TempIndex<len(OutputHeaderList):
                    NewHeaderList.append(OutputHeaderList[TempIndex])
                else:
                    NewHeaderList.append(FetchedDataArray.dtype.names[TempIndex])
            FetchedDataArray.dtype.names = NewHeaderList
    if len(FetchedDataArray)>1:
        if Verbose:
            print "fetchInfoRow: Read %d rows\nDone!\n\n"%(len(FetchedDataArray))
    else:
        if Verbose:
            print "fetchInfoRow: Read %d row\nDone!\n\n"%(len(FetchedDataArray))
    # 
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




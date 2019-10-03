#!/usr/local/bin/python
# -*- coding: utf-8 -*-
# 
# 
# This script will read data row(s) from data table(s)
# 
# input a list of data table file path
# input a list of data table row string pattern
# output an array of data row(s)
# 
# Usage:
#   execfile('/Users/dliu/Software/ipy/setup.py')
#   data = crabtable.fetchDataRow('Arp193_FromHIPE_S?W.CSV','.*CI.*2.*1.*')
#   data = crabtable.fetchDataRow('/data/SpireLines/DataRelease/2014-05-31/ListOfSourcesPart?/Arp*/sfit/SPIRE_FTS_SOF1_SincGauss/Arp*_FromHIPE_S?W.CSV','.*CI.*2.*1.*',[],'[^ ]*/([^/]*)/sfit/[^ ]*')
#   asciitable.write(data[['DataTable','LineName','FluxPrime','σFluxPrime']],"FetchedDataTable.DAT",Writer=asciitable.FixedWidth,delimiter=' ',delimiter_pad=' ',bookend=False)
#   
#   TempDataArray = asciitable.read('Arp193_FromHIPE_SLW.CSV')
# 
# Last update:
#   2015-03-20 created
# 


def fetchDataRow(DataTableFilePathList,MatchStringPatternList,MatchHeaderList=[],OutputHeaderList=[],OutputDataTablePattern=[],Verbose=True):
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
    
    # fecth data array
    FetchedDataArray = []
    for TempDataFile in TempDataTableFilePathList:
        if Verbose: print "fetchDataRow: Reading %s\n"%(TempDataFile)
        TempDataArray = asciitable.read(TempDataFile)
        for TempRowContent in TempDataArray:
            TempRowMatched = 0
            for TempIndex in range(len(TempRowContent.dtype)):
                if len(MatchHeaderList)>0:
                    if not TempRowContent.dtype.names[TempIndex] in MatchHeaderList:
                        continue
                try:
                    if TempRowContent.dtype[TempIndex].name.index("string")>=0:
                        # TempRowContent.dtype[TempIndex].name is the dtype type name e.g. float64 string etc
                        # TempRowContent.dtype.names[TempIndex] is the dtype field name e.g. Header1 Header2 etc
                        for TempRegexPattern in TempRegexPatternList:
                            if TempRegexPattern.match(TempRowContent[TempRowContent.dtype.names[TempIndex]]) != None:
                                if Verbose: print "fetchDataRow: Matched %s to Column %s of Row %s in %s\n"%(TempRegexPattern.pattern,TempRowContent.dtype.names[TempIndex],TempRowContent,TempDataFile)
                                TempRowMatched = 1
                                break
                        if TempRowMatched == 1:
                            break
                except:
                    pass
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
                
                # <fixed><20151226><dzliu> headers
                if len(TempDataArray.dtype.names)>0 and len(TempNewContent)>0:
                    NewHeaderList = []
                    for TempIndex in range(len(TempNewContent.dtype.names)):
                        if TempIndex<len(TempDataArray.dtype.names):
                            NewHeaderList.append(TempDataArray.dtype.names[TempIndex])
                        else:
                            NewHeaderList.append(TempNewContent.dtype.names[TempIndex])
                    TempNewContent.dtype.names = NewHeaderList
                
                # append new row to FetchedDataArray
                if len(FetchedDataArray) == 0:
                    FetchedDataArray = TempNewContent
                    del(TempNewContent)
                else:
                    TempOldContent = FetchedDataArray
                    FetchedDataArray = rec.stack_arrays((TempOldContent,TempNewContent),usemask=False,asrecarray=True,autoconvert=True)
                    del(TempOldContent)
                    del(TempNewContent)
        ### <fixed><20151226><dzliu> headers
        # # <added><201512010><dzliu> headers
        # if len(TempDataArray.dtype.names)>0 and len(FetchedDataArray)>0:
        #     NewHeaderList = []
        #     for TempIndex in range(len(FetchedDataArray.dtype.names)):
        #         if TempIndex<len(TempDataArray.dtype.names):
        #             NewHeaderList.append(TempDataArray.dtype.names[TempIndex])
        #         else:
        #             NewHeaderList.append(FetchedDataArray.dtype.names[TempIndex])
        #     FetchedDataArray.dtype.names = NewHeaderList
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
    # 
    if len(FetchedDataArray) > 1:
        if Verbose:
            print "fetchDataRow: Read %d rows\nfetchDataRow: Done!\n\n"%(len(FetchedDataArray))
    else:
        if Verbose:
            print "fetchDataRow: Read %d row\nfetchDataRow: Done!\n\n"%(len(FetchedDataArray))
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




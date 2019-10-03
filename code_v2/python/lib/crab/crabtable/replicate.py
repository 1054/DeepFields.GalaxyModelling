#!/usr/local/bin/python
# -*- coding: utf-8 -*-
# 
# 
# This script will replicate one value to make an array
# 
# Last update:
#   2015-03-27 created
# 


def replicate(Value,Dimension):
    NewList = []
    # make Value a list
    if type(Value) is not list:
        TempValue = numpy.array([Value])
    else:
        TempValue = numpy.array(Value)
    # loop
    for TempCounter in range(Dimension):
        if TempCounter == 0:
            NewList = TempValue
        else:
            NewList = numpy.append(NewList,TempValue)
    return NewList




# the dependent importing should be put at the end
# http://effbot.org/zone/import-confusion.htm
try:
    import numpy
except ImportError: 
    print "Error! Could not import numpy!"


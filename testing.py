# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 16:54:59 2016

@author: rink
"""
from pyfast import CalcLookupTable,ReadFASTFile

#TurbName = 'WP_0.75MW'
#ModlDir  = 'C:\\Users\\rink\\Dropbox\\contracts\\' + \
#                'nrel-windpact\\fast_input_files\\v7\\'+TurbName
#
#WriteLookupTable(TurbName,ModlDir,WindDir=None,
#                     version=7,WindSpeeds=None,TMax=140.,Tss=80.,
#                     overwrite=0)

FilePath = 'C:\\Users\\rink\\Dropbox\\contracts\\nrel-windpact\\' + \
                'fast_input_files\\v8\\WP_0.75MW\\WP_0.75MW.out'

df = ReadFASTFile(FilePath)

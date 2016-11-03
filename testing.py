# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 16:54:59 2016

@author: rink
"""
from pyfast import CalcLookupTable,ReadFASTFile
import os

TurbName = 'WP_0.75MW'
version  = 8
#ModlDir  = 'C:\\Users\\rink\\Dropbox\\contracts\\' + \
#                'nrel-windpact\\fast_input_files\\v7\\'+TurbName
ModlDir  = os.path.join('C:\\Users\\jrinker\\Dropbox\\contracts\\' + \
                'nrel-windpact\\fast_input_files','v'+str(version),TurbName)

FastExe = 'FAST_v{}.exe'.format(version)
WindSpeeds = [5]
CalcLookupTable(TurbName,ModlDir,WindDir=ModlDir,
                    WindSpeeds=WindSpeeds,TMax=140.,Tss=80.,
                    FastExe=FastExe,clean=0,overwrite=1)

#FilePath = 'C:\\Users\\rink\\Dropbox\\contracts\\nrel-windpact\\' + \
#                'fast_input_files\\v8\\WP_0.75MW\\WP_0.75MW.out'
#
#df = ReadFASTFile(FilePath)

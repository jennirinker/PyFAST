"""
PyFAST: a collection of tools for working with FAST,
        NREL's open-source wind turbine simulator
        
Author: Jennifer M. Rinker

"""
import os, scipy.io as scio

# =============================================================================
#                  Reading/writing FAST input files
# =============================================================================

def WriteFastAD(TurbName,WindPath,ModlDir,
                FastDir=None,FastName=None,
                version=7,verbose=0,lin=0,
                **kwargs):
    """ Write FAST and AeroDyn input files for a specified wind file
    
        - Possible keyword arguments inclue any FAST initial conditions,
          TMax, and TStart.
        - If model does not have a pre-calculated steady-state
          look-up table, default values for initial conditions will be used.
        - AeroDyn and FAST template files must be present in subdirectory 
          "templates/" with the naming convention "<TurbName>_AD_template.ipt"
          and "<TurbName>_template.fst".
        - Written FAST/AeroDyn files will be named "<FastName>.fst" and 
          "<FastName>_AD.ipt". Default value for <FastName> is <TurbName>.
    
        Args:
            TurbName (str)  : turbine name
            WindPath (str)  : path to wind file
            ModlDir  (str)  : directory with wind-independent files (e.g.,
                              Blade, Tower, Pitch files)
            FastDir  (str)  : directory to write FAST & AeroDyn files to [opt]
            FastName (str)  : name for .fst file [opt]
            version  (int)  : FAST version (7 or 8) [opt]
            verbose  (int)  : flag to suppress print statements [opt]
            lin      (int)  : flag to indicate linearization analysis [opt]
            kwargs  (dict)  : keyword arguments to WriteFastADOne [opt]
                 
    """
    
    # check that the templates subfolder exists
    TmplDir    = os.path.join(ModlDir,'templates')       # template dir
    if not os.path.isdir(TmplDir):
        raise ValueError('Template directory for ')
    
    # print status update if verbose
    if verbose:
        print('\nWriting FAST v{:.0f} files for \"{:s}\" '.format(version,TurbName) + \
                'with wind file {:s}...'.format(WindPath))
    
    # set FastName and FastDir if not given
    if FastName is None:
        FastName = TurbName     # .fst name will be turbine name
    if FastDir is None:
        FastDir  = ModlDir      # .fst written to folder with blade, etc. files
    
    # define dictionary of default wind-dependent parameters
    WindDict = {'BlPitch(1)':0.,'BlPitch(2)':0.,'BlPitch(3)':0.,
                'OoPDefl':0.,'IPDefl':0.,'TeetDefl':0.,'Azimuth':0.,
                'RotSpeed':0.,'NacYaw':0.,'TTDspFA':0.,'TTDspSS':0.,
                'TMax':630.,'TStart':30.,'WindFile':WindPath,
                'StallMod':'BEDDOES','AnalMode':1,'PCMode':1,
                'Gravity':9.80665,'GenDOF':'True'}
#TODO : finish separating linearization vs. not inputs
    if lin:
        WindDict['LinFile']  = 'unused'
        WindDict['CompAero'] = 'False'
    else:
        WindDict['LinFile']  = 'unused'
        WindDict['CompAero'] = 'True'
            
    # set values of passed-in keyword arguments
    for key in kwargs:
        if key in WindDict.keys():
            WindDict[key] = kwargs[key]
    
    # check if steady-state look-up table exists
    LUTPath = os.path.join(ModlDir,'steady_state',
                           TurbName+'_SS.mat')
                                      
    #    if LUT exists, load unspecific initial condition values
    if os.path.exists(LUTPath):
        if verbose:
            print('  Interpolating unspecified IC values from ' + \
                    'look-up table {:s}'.format(LUTPath))
                
        # get list of keys corresponding to initial conditions and LUT keys
        IC_keys  = GetICKeys()
        mdict    = scio.loadmat(LUTPath,squeeze_me=True)
        LUT_keys = [str(s).strip() for s in mdict['Fields']]
        LUT      = mdict['SS']
        
        # loop through IC keys
        for i_key in range(len(IC_keys)):
            IC_key = IC_keys[i_key]
            
            # if it's a blade pitch, different FAST key than LUT
            if ('BlPitch' in IC_key):
                IC_key = 'BldPitch'                    
                
            # see if IC key is in LUT, get LUT key if applicable
            try:
                LUT_key = [s for s in LUT_keys if IC_key in s][0]
            except IndexError:
                LUT_key = []
            
            # if IC was not passed in and key is in LUT
            if ((not [key for key in kwargs if IC_key in kwargs]) and \
                LUT_key):
                    
                # get grid-averaged first wind speed
                u0 = jr_wind.GetFirstWind(WindDict['WindFile'])
                
                # linearly interpolate initial condition
                IC = np.interp(u0,LUT[:,LUT_keys.index('WindVxi')],
                                    LUT[:,LUT_keys.index(LUT_key)])
                                    
                # save initial condition
                WindDict[IC_keys[i_key]] = IC

    # create and save filenames
    ADTempPath = os.path.join(TmplDir,TurbName+'_AD_template.ipt')
    ADName   = '{:s}_AD.ipt'.format(FastName)
    ADPath     = os.path.join(FastDir,ADName)
    FastName = '{:s}.fst'.format(FastName)
    FastTempPath = os.path.join(TmplDir,TurbName+'_template.fst')
    FastPath  = os.path.join(FastDir,FastName)
    
    WindDict['ADFile'] = os.path.join(FastDir,ADPath)       # version 7 key
    WindDict['AeroFile'] = os.path.join(FastDir,ADPath)     # version 8 key
    
    # write AeroDyn file
    with open(ADTempPath,'r') as f_temp:
        with open(ADPath,'w') as f_write:
            for line in f_temp:
                if ('{:' in line):
                    field = line.split()[1]
                    f_write.write(line.format(WindDict[field]))
                else:
                    f_write.write(line)
                
    # write FAST file
    with open(FastTempPath,'r') as f_temp:
        with open(FastPath,'w') as f_write:
            for line in f_temp:
                if ('{:' in line):
                    field = line.split()[1]
                    f_write.write(line.format(WindDict[field]))
                else:
                    f_write.write(line)
        
    # if FAST version is 8, also write ElastoDyn file
    if (version == 8):
        
        # create and save filepaths
        EDTempPath = os.path.join(TmplDir,TurbName+'_ElastoDyn_template.dat')
        EDName   = '{:s}_ElastoDyn.ipt'.format(FastName)
        EDPath     = os.path.join(FastDir,EDName)
              
        # write ElastoDyn file
        with open(EDTempPath,'r') as f_temp:
            with open(EDPath,'w') as f_write:
                for line in f_temp:
                    if ('{:' in line):
                        field = line.split()[1]
                        f_write.write(line.format(WindDict[field]))
                    else:
                        f_write.write(line)
                             
    return
    


# =============================================================================
#                  Miscellaneous utilities
# =============================================================================

def GetICKeys():
    """ List of keys in .fst that are initial conditions
    
        Returns:
            IC_keys (list): list of FAST keys that are initial conditions
    """
    
    IC_keys = ['BlPitch(1)','BlPitch(2)','BlPitch(3)',
                    'OoPDefl','IPDefl','TeetDefl','Azimuth',
                    'RotSpeed','TTDspFA','TTDspSS']
    
    return IC_keys
	
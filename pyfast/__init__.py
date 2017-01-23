"""
PyFAST : A collection of tools for working with FAST,
         NREL's open-source wind turbine simulator
        
Author : Jennifer M. Rinker

Date   : 20-oct-2016
"""
import numpy as np
import os
import pandas as pd
from   PyTurbSimLite import readModel
import scipy.io as scio
import shutil 
    

def CalcLookupTable(TurbName,ModlDir,
                    WindSpeeds=None,FastExe='FAST.exe',
                    TMax=140.,Tss=80.,overwrite=0,
                    **kwargs):
    """ Calculate and save steady-state look-up table
    
        Runs a series of steady-state simulations using the previous
        steady-state values as initial conditions. Saves all steady-state
        values into a table and saves it.
        
        Note: the FAST executable must either be on the system path or in the
        model directory.
        
        Args:
            TurbName    (str) : name of turbine model
            ModlDir     (str) : path to model directory (no trailing slash)
            WindSpeeds (iter) : iterable of wind speeds at which to calculate
                                steady-state values
            FastExe     (str) : name of FAST executable to simulate turbine
            TMax        (flt) : total simulation time
            Tss         (flt) : averaging time for calculating steady-state
            overwrite   (int) : force directory overwrite without asking
            
        Returns:
            LUT        (dict) : look-up table 
    """
    
    # check if model directory exists
    if not os.path.isdir(ModlDir):
        raise ValueError('Model directory {:s} does not exist.'.format(ModlDir))
    
    # check if steady-state directory exists, overwrite if user allows
    SSDir = os.path.join(ModlDir,'steady-state')
    if os.path.isdir(SSDir) and not overwrite:
        try:                    # bind raw_input to input for Python 2
            input = raw_input
        except NameError:
            pass
        UserResp = input('Steady-state directory already exists. Overwrite? [y/n] ')
        if UserResp in ['y','Y',1]:
            shutil.rmtree(SSDir)
        elif UserResp in ['n','N',0]:
            return None
        else:
            raise ValueError('Unknown response {}'.format(UserResp))
    elif os.path.isdir(SSDir) and overwrite:
        shutil.rmtree(SSDir)
    os.mkdir(SSDir)
    
    # define default wind speeds, ensure monotonically increasing
    if WindSpeeds is None:
        WindSpeeds = np.arange(3,25,0.5)
    else:
        WindSpeeds = np.sort(np.array(WindSpeeds))
    
    # define LUT name and path
    SSName = TurbName + '_SS.mat'
    SSPath = os.path.join(SSDir,SSName)
    
    # change directory to steady-state directory to run FAST
    os.chdir(ModlDir)
    
    # define initial dictionary of default wind-dependent parameters
    WindDict = {'BlPitch(1)':0.,'BlPitch(2)':0.,'BlPitch(3)':0.,
                'OoPDefl':0.,'IPDefl':0.,'TeetDefl':0.,'Azimuth':0.,
                'RotSpeed':0.,'NacYaw':0.,'TTDspFA':0.,'TTDspSS':0.,
                'TMax':TMax,'TStart':0.}
                
    # set initial values of passed-in keyword arguments
    for key in kwargs:
        if key in WindDict.keys():
            WindDict[key] = kwargs[key]
    
    # loop through wind speeds
    NumIDs = int(np.ceil(np.log10(len(WindSpeeds))))
    
    for iWS in range(len(WindSpeeds)):
        WindSpeed = WindSpeeds[iWS]
        
        print('Processing wind speed {} of {}...'.format(iWS+1,len(WindSpeeds)))
        
        # define file ID for FAST run
        fileID     = '{:.0f}'.format(iWS).zfill(NumIDs)
        
        # create wind filename
        WindName = 'NoShr_'+'{:2.1f}'.format(WindSpeed).zfill(4)+'.wnd'
        
        # check if wind file exists, make it if not
        WindPath = os.path.join(SSDir,WindName)
        if not os.path.exists(WindPath):
            WriteSteadyWind(WindSpeed,WindPath)
                
        # write FAST input files
        FastName = TurbName + '_' + fileID
        FastPath = os.path.join(ModlDir,FastName)
        SSFastPath = os.path.join(SSDir,FastName)
        WriteFastAD(TurbName,WindPath,ModlDir,
                    FastDir=ModlDir,FastName=FastName,
                    **WindDict)
                  
        # run FAST
        command  = ' '.join([FastExe,FastName + '.fst' ])
        ExitCode = os.system(command)
        if ExitCode:
            raise ValueError('FAST did not complete successfully ' + \
                    '(Wind speed {:.1f}, exit code {:.0f})'.format(WindSpeed,ExitCode))
                      
        # load FAST files
        FASTdf = ReadFASTFile(FastName + '.out')
        Fields = [s for s in FASTdf.columns]
        
        # initialize LUT if it doesn't exist
        if iWS == 0: LUT = np.empty((len(WindSpeeds),len(Fields)))
        
        # loop through and save steady-state values
        n_t = int(FASTdf['Time'].size*Tss/TMax)
        for i_parm in range(len(Fields)):
            
            # get data
            parm = Fields[i_parm]
            x = FASTdf[parm]
                          
            # calculate and save last value
            x_SS = np.mean(x[-n_t:])
            LUT[iWS,i_parm] = x_SS
                                  
        # set initial conditions for next round
        for key in WindDict.keys():
            if key in Fields:
                WindDict[key] = LUT[iWS,Fields.index(key)]
            elif key+'1' in Fields:
                WindDict[key] = LUT[iWS,Fields.index(key+'1')]

        # delete fst, ipt files, move .out
        os.system('del {}'.format(WindPath))
        os.system('del {}.fst'.format(FastPath))
        os.system('del {}_AD.ipt'.format(FastPath))
        os.system('del {}.sum'.format(FastPath))
        os.system('move {}.out {}.out'.format(FastPath,
                                                SSFastPath))
    
    print('Simulations completed.')    
       
    # rearrange to increasing wind speed
    LUT = LUT[LUT[:,Fields.index('WindVxi')].argsort()]
       
    # save LUT   
    LUTdict = {}
    LUTdict['SS'] = LUT   
    LUTdict['Fields'] = Fields 
    scio.savemat(SSPath,LUTdict)
    print('Look-up table saved.')
    
    return LUT 
    

def GetDefaultOutputs():
    """ Defaults outputs for FAST time-marching analysis
    
        Returns:
            FastOutputs (list) : default outputs for time-marching analyses
    """
    FastOutputs = ['"WindVxi,WindVyi,WindVzi"          		- Wind-speed components\n',
    '"OoPDefl1,IPDefl1,TipDzb1,TwrClrnc1"	- Blade 1 tip motions\n',
    '"OoPDefl2,IPDefl2,TipDzb2,TwrClrnc2"	- Blade 2 tip motions\n',
    '"OoPDefl3,IPDefl3,TipDzb3,TwrClrnc3"	- Blade 3 tip motions\n',
    '"BldPitch1,BldPitch2,BldPitch3"         - Blade pitch motions\n',
    '"Azimuth,RotSpeed,RotAccel"		        - LSS motion\n',
    '"GenSpeed,GenAccel"				    	- HSS motion\n',
    '"TTDspFA,TTDspSS,TTDspAx"				- Towertop motions\n',
    '"YawBrTAxp,YawBrTAyp,YawBrTAzp"			- Towertop accelerations\n',
    '"RootFzb1,RootMIP1,RootMOoP1,RootMzb1"	- Blade 1 root loads (1/2)\n',
    '"RootMEdg1,RootMFlp1"					- Blade 1 root loads (2/2)\n',
    '"RootFzb2,RootMIP2,RootMOoP2,RootMzb2"	- Blade 2 root loads (1/2)\n',
    '"RootMEdg2,RootMFlp2"					- Blade 2 root loads (2/2)\n',
    '"RootFzb3,RootMIP3,RootMOoP3,RootMzb3"	- Blade 3 root loads (1/2)\n',
    '"RootMEdg3,RootMFlp3"					- Blade 3 root loads (2/2)\n',
    '"RotThrust,LSSGagFya,LSSGagFza"			- Hub and rotor loads (1/3)\n',
    '"LSSGagFys,LSSGagFzs,RotTorq"			- Hub and rotor loads (2/3)\n',
    '"LSShftPwr"								- Hub and rotor loads (3/3)\n',
    '"HSShftTq,HSShftPwr"					- Generator and HSS loads\n',
    '"GenPwr"								- Generator power\n',
    '"YawBrFzp,YawBrMzp"						- Tower top yaw-bearing loads\n',
    '"TwrBsFxt,TwrBsFyt,TwrBsFzt"			- Tower base loads (1/2)\n',
    '"TwrBsMxt,TwrBsMyt,TwrBsMzt"			- Tower base loads (2/2)\n']
    
    return FastOutputs
	
 
def GetFirstWind(WindPath):
    """ First wind speed from wind file
    
        Args:
            WindPath (string): path to wind file
            
        Returns:
            u0 (float): initial wind value
    """
    
    # if it's a .wnd file
    if WindPath.endswith('.wnd'):
        
        # try to read it as a text file
        try:
            with open(WindPath,'r') as f:
                f.readline()
                f.readline()
                f.readline()
                first_line = f.readline().split()
                u0 = float(first_line[1])
                
        # if error, try to read as binary file
        except:
            turb  = readModel(WindPath)         # read file
            u0    = turb[0,:,:,0].mean()        # take mean of u(t0)
    
    # if it's a .bts file
    elif WindPath.endswith('.bts'):
        
        turb  = readModel(WindPath)             # read file
        u0    = turb[0,:,:,0].mean()            # take mean of u(t0)
        
    # if it's a .bl file
    elif WindPath.endswith('.bl'):
        
        turb  = readModel(WindPath)             # read file
        u0    = turb[0,:,:,0].mean()            # take mean of u(t0)

    else:
        errStr = 'Uncoded file extension ' + \
                        '\"{:s}\"'.format(WindPath.endswith())
        raise ValueError(errStr)
        
    return u0  
    

def GetICKeys():
    """ List of keys in .fst that are initial conditions
    
        Returns:
            IC_keys (list): list of FAST keys that are initial conditions
    """
    
    IC_keys = ['BlPitch(1)','BlPitch(2)','BlPitch(3)',
                'OoPDefl','IPDefl','TeetDefl','Azimuth',
                'RotSpeed','TTDspFA','TTDspSS']
    
    return IC_keys


def GetTemplateDict():
    """ Dictionary of number formats for creating template file
                 
    """
    
    TmplDict = {'ADFile':'\"{:s}\"   	   ',
                'AeroFile':'\"{:s}\"        ',
                'AnalMode':'   {:.0f}        ',
                'Azimuth':'{:6.1f}     ',
                'BlPitch(1)':'{:6.1f}      ',
                'BlPitch(2)':'{:6.1f}      ',
                'BlPitch(3)':'{:6.1f}      ',
                'CompAero':'{:<5s}       ',
                'FastOutPuts':'{:s} ',
                'GenDOF':'{:<5s}       ',
                'Gravity':'   {:7.5f}  ',
                'IPDefl':'{:6.1f}     ',
                'LinFile':'\"{:s}\"    ',
                'OoPDefl':'{:6.1f}     ',
                'NacYaw':'{:6.1f}     ',
                'PCMode':'   {:.0f}        ',
                'RotSpeed':'{:6.1f}     ',
                'StallMod':'{:<7s}                                ',
                'TeetDefl':'{:6.1f}     ',
                'TStart':'{:5.1f}       ',
                'TTDspFA':'{:6.1f}     ',
                'TTDspSS':'{:6.1f}     ',
                'TMax':'{:8.1f}    ',
                'WindFile':'\"{:s}\"			   '}
                      
    return TmplDict


def ReadFASTFile(FilePath):
    """ Read FAST data into pandas dataframe
    
        Args:
            FilePath (string): path to .fst file
    
        Returns:
            FASTDict (dictionary): dictionary with FAST data
    """
    
    # number of lines to skip at the beginning of the file
    SkipRows = 8
    
    # get field names and units
    with open(FilePath,'r') as f:
        for i in range(6):
            f.readline()
        Fields = f.readline().strip('\n').split()
        Units  = f.readline().strip('\n').split()
    Fields = [str(s) for s in Fields]  # replace weird character
    Units  = [s.replace('\xb7','-') for s in Units]  # replace weird character
        
    # read as pandas dataframe
    FASTdf = pd.read_csv(FilePath,
                         delim_whitespace=True,header=0,names=Fields,
                         skiprows=SkipRows)
    FASTdf.units = Units
    
    return FASTdf
    

def ReplaceVal(Line,WriteFile,WindDict):
    """ Replace the value in the line with value in WindDict
        
        Args:
            Line       (str) : line with value to replace
            WriteFile (file) : file object to be written to
            WindDict  (dict) : dictionary with input file values
    """
    
    # determine which field line corresponds to
    Field = Line.split()[1]
    
    # if it's the output line, write all outputs for time-marching analysis or
    #    none for a linearization analysis
    if Field == 'FastOutputs':
        FastOutputs = WindDict['FastOutputs']
        if FastOutputs is None:
            WriteFile.write('')
        else:
            for OutLine in FastOutputs:
                WriteFile.write(OutLine)
                
    # for all other fields, just substitue the field into the line
    else:
        NewLine = Field.join([Line.split(Field)[0].format(WindDict[Field]),
                              Line.split(Field)[1]])
        WriteFile.write(NewLine)
    
    return


def WriteFastAD(TurbName,WindPath,ModlDir,
                FastDir=None,FastName=None,version=7,verbose=0,lin=0,
                **kwargs):
    """ Write FAST and AeroDyn input files for a specified wind file
    
        - Possible keyword arguments inclue any FAST initial conditions,
          TMax, and TStart.
        - If model does not have a pre-calculated steady-state
          look-up table, zero values for initial conditions will be used.
        - AeroDyn and FAST template files must be present in subdirectory 
          "<ModlDir>/templates/" with the naming convention 
          "<TurbName>_AD_template.ipt" and "<TurbName>_template.fst", 
          respectively.
        - Written FAST/AeroDyn files will be named "<FastName>.fst" and 
          "<FastName>_AD.ipt". Default value for <FastName> is <TurbName>.
    
        Args:
            TurbName (str)  : turbine name
            WindPath (str)  : path to wind file
            ModlDir  (str)  : directory with wind-independent files (e.g.,
                              Blade, Tower, Pitch files)
            FastDir  (str)  : directory to write FAST & AeroDyn files to [opt]
            FastName (str)  : name for written .fst file, no extension [opt]
            version  (int)  : FAST version (7 or 8) [opt]
            verbose  (int)  : flag to suppress print statements [opt]
            lin      (int)  : flag to indicate linearization analysis [opt]
            kwargs  (dict)  : keyword arguments to WriteFastADOne [opt]
                 
    """
    
    # check that the templates subfolder exists
    TmplDir    = os.path.join(ModlDir,'templates')       # template directory
    if not os.path.isdir(TmplDir):
        raise ValueError('Template directory not found in ' + ModlDir)
    
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
                'WindFile':WindPath}
             
    # define default parameters specific to linearization
    if lin:
        WindDict['LinFile']     = os.path.join(FastDir,TurbName+'_Linear.dat')
        WindDict['CompAero']    = 'False'
        WindDict['PCMode']      = 0
        WindDict['TMax']        = 30.
        WindDict['TStart']      = 0.
        WindDict['StallMod']    = 'STEADY'
        WindDict['AnalMode']    = 2
        WindDict['Gravity']     = 0.
        WindDict['GenDOF']      = 'False'
        WindDict['FastOutputs'] = None
    else:
        WindDict['LinFile']     = 'unused'
        WindDict['CompAero']    = 'True'
        WindDict['PCMode']      = 1
        WindDict['TMax']        = 630.
        WindDict['TStart']      = 30.
        WindDict['StallMod']    = 'BEDDOES'
        WindDict['AnalMode']    = 1
        WindDict['Gravity']     = 9.80665
        WindDict['GenDOF']      = 'True'
        WindDict['FastOutputs'] = GetDefaultOutputs()
            
    # set values of passed-in keyword arguments
    for key in kwargs:
        if key in WindDict.keys():
            WindDict[key] = kwargs[key]
    
    # check if steady-state look-up table exists
    LUTPath = os.path.join(ModlDir,'steady-state',
                           TurbName+'_SS.mat')
                                      
    # print message if verbose and LUT doesn't exist
    if not os.path.exists(LUTPath):
        if verbose:
            print('  WARNING: look-up table {:s} '.format(LUTPath) + \
                    'not found. Using defaults.')
        
    # if LUT exists, load unspecific initial condition values
    else:
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
                u0 = GetFirstWind(WindDict['WindFile'])
                
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
    with open(ADTempPath,'r') as TmplFile:
        with open(ADPath,'w') as WriteFile:
            for Line in TmplFile:
                if ('{:' in Line):
                    ReplaceVal(Line,WriteFile,WindDict)
                else:
                    WriteFile.write(Line)
                
    # write FAST file
    with open(FastTempPath,'r') as TmplFile:
        with open(FastPath,'w') as WriteFile:
            for Line in TmplFile:
                if ('{:' in Line):
                    ReplaceVal(Line,WriteFile,WindDict)
                else:
                    WriteFile.write(Line)
        
    # if FAST version is 8, also write ElastoDyn file
    if (version == 8):
        
        # create and save filepaths
        EDTempPath = os.path.join(TmplDir,TurbName+'_ElastoDyn_template.dat')
        EDName   = '{:s}_ElastoDyn.ipt'.format(FastName)
        EDPath     = os.path.join(FastDir,EDName)
              
        # write ElastoDyn file
        with open(EDTempPath,'r') as TmplFile:
            with open(EDPath,'w') as WriteFile:
                for Line in TmplFile:
                    if ('{:' in Line):
                        ReplaceVal(Line,WriteFile,WindDict)
                    else:
                        WriteFile.write(Line)
    
    # print status update if verbose
    if verbose:
        print('   done.')
                      
    return


def WriteFastADAll(TurbName,ModlDir,WindDir,
                   FastDir=None,version=7,
                   **kwargs):
    """ Write FAST/AeroDyn input files for all wind files in chosen directory
    
        - Written files will be labeled "<TurbName>_<WindFileID>.fst" and 
          "<TurbName>_<WindFileID>_AD.ipt", respectively, where <WindFileID>
          corresponds to the ordered list of wind files analyzed.
    
        Args:
            TurbName   (str) : turbine name
            ModlDir    (str) : directory with wind-independent files (e.g.,
                               Blade, Tower, Pitch files)
            WindDir    (str) : directory with wind files
            FastDir    (str) : directory to write FAST & AeroDyn files to
            version    (int) : FAST version (7 or 8) [opt]
            kwargs    (dict) : keyword arguments to WriteFastADOne [opt]
        Returns:
            WindFiles (list) : ordered list of wind files
    """
    
    # check that the templates subfolder exists
    TmplDir    = os.path.join(ModlDir,'templates')       # template directory
    if not os.path.isdir(TmplDir):
        raise ValueError('Template directory not found in ' + ModlDir)
    
    # possible wind file endings
    wind_ends = ('.bts','.wnd','.bl')
            
    # get list of wind files from directory
    WindNames     = [f for f in os.listdir(WindDir) if f.endswith(wind_ends)]
    NumWinds      = len(WindNames)
    LenWindFileID = int(np.ceil(np.log10(NumWinds)))

    # loop through wind files
    WindPaths = []
    for iWind in range(NumWinds):
        WindName   = WindNames[iWind]
        WindPath   = os.path.join(WindDir,WindName)
        WindPaths.append(WindPath)
        
        # set file ID number and define name for .fst file to be written
        WindFileID = str(iWind).zfill(LenWindFileID)
        FastName   = '_'.join([TurbName,WindFileID])
        
        # write FAST/AD files for the wind file
        WriteFastAD(TurbName,WindPath,ModlDir,
                    FastDir=FastDir,FastName=FastName,
                    version=version,verbose=0,lin=0,
                    **kwargs)
    
    return WindPaths


def WriteSteadyWind(WindSpeed,WindPath):
    """ Write steady-state wind file
    
        Args:
            WindSpeed (float): wind value
            WindPath (string): optional path to directory with wind files
    """
    
    with open(WindPath,'w') as f:
        f.write('! Wind file for steady {:.1f} m/s wind without shear.\n'.format(WindSpeed))
        f.write('! Time	Wind	Wind	Vert.	Horiz.	Vert.	LinV	Gust\n')
        f.write('!	Speed	Dir	Speed	Shear	Shear	Shear	Speed\n')
        f.write('   0.0\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(WindSpeed))
        f.write('   0.1\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(WindSpeed))
        f.write('9999.9\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(WindSpeed))
        
    return


def WriteTemplate(TurbName,FastDir,
                  version=7,**kwargs):
    """ Create template file from FAST model
                 
    """
    
    # create and save filenames
    TmplDir = os.path.join(FastDir,'templates')
    ADTmplPath = os.path.join(TmplDir,TurbName+'_AD_template.ipt')
    ADName   = '{:s}_AD.ipt'.format(TurbName)
    ADPath     = os.path.join(FastDir,ADName)
    FastName = '{:s}.fst'.format(TurbName)
    FastTmplPath = os.path.join(TmplDir,TurbName+'_template.fst')
    FastPath  = os.path.join(FastDir,FastName)
    
    # get dictionary of template keys
    TmplDict = GetTemplateDict()
    
    # write AeroDyn template file
    with open(ADPath,'r') as OrigFile:
        with open(ADTmplPath,'w') as TmplFile:
            for Line in OrigFile:
                NewLine = Line
                SplitLine = Line.split()
                if len(SplitLine) > 1:
                    if SplitLine[1] in TmplDict.keys():
                        Key = SplitLine[1]
                        Fmt = TmplDict[Key]
                        NewLine = Fmt + Key + Line.split(Key)[1]
                TmplFile.write(NewLine)
                
    # write FAST file
    OutputLine = 0
    with open(FastPath,'r') as OrigFile:
        with open(FastTmplPath,'w') as TmplFile:
            for Line in OrigFile:
                NewLine = Line
                SplitLine = Line.split()
                if 'END' in SplitLine[0]:
                    OutputLine = 0
                elif OutputLine:
                    NewLine = ''
                elif len(SplitLine) > 1:
                    if 'OutList' in SplitLine[0]:
                        OutputLine = 1
                        NewLine += '{:s} FastOutputs\n'
                    elif SplitLine[1] in TmplDict.keys():
                        Key = SplitLine[1]
                        Fmt = TmplDict[Key]
                        NewLine = Fmt + Key + Line.split(Key)[1]
                TmplFile.write(NewLine)
        
    # if FAST version is 8, also write ElastoDyn file
    if (version == 8):
        # TODO: CONTINUE SORTING OUT format keys
        # create and save filepaths
        EDTmplPath = os.path.join(TmplDir,TurbName+'_ElastoDyn_template.dat')
        EDName   = '{:s}_ElastoDyn.dat'.format(TurbName)
        EDPath     = os.path.join(FastDir,EDName)
              
        # write ElastoDyn file
        OutputLine = 0
        with open(EDPath,'r') as OrigFile:
            with open(EDTmplPath,'w') as TmplFile:
                for Line in OrigFile:
                    NewLine = Line
                    SplitLine = Line.split()
                    if 'END' in SplitLine[0]:
                        OutputLine = 0
                    elif OutputLine:
                        NewLine = ''
                    elif len(SplitLine) > 1:
                        if 'OutList' in SplitLine[0]:
                            OutputLine = 1
                            NewLine += '{:s} FastOutputs\n'
                        elif SplitLine[1] in TmplDict.keys():
                            Key = SplitLine[1]
                            Fmt = TmplDict[Key]
                            NewLine = Fmt + Key + Line.split(Key)[1]
                    TmplFile.write(NewLine)
    
                      
    return


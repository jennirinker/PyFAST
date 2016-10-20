"""
PyFAST : A collection of tools for working with FAST,
         NREL's open-source wind turbine simulator
        
Author : Jennifer M. Rinker

Date   : 20-oct-2016
"""
import os
import scipy.io as scio
import numpy as np
from struct import unpack
from warnings import warn

# =============================================================================
#                  Reading/writing FAST input files
# =============================================================================

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
        WindDict['LinFile']     = 'unused'
        WindDict['CompAero']    = 'False'
        WindDict['PCMode']      = 0
        WindDict['TMax']        = 10.
        WindDict['TStart']      = 0.
        WindDict['StallMod']    = 'STEADY'
        WindDict['AnalMode']    = 2
        WindDict['Gravity']     = 0.
        WindDict['GenDOF']      = 'True'
        WindDict['LinFile']     = os.path.join(FastDir,TurbName+'_Linear.dat')
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
        WindDict['FastOutputs'] = DefaultOutputs()
            
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
    
# =============================================================================
#                  Miscellaneous utilities
# =============================================================================

def DefaultOutputs():
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

def GetICKeys():
    """ List of keys in .fst that are initial conditions
    
        Returns:
            IC_keys (list): list of FAST keys that are initial conditions
    """
    
    IC_keys = ['BlPitch(1)','BlPitch(2)','BlPitch(3)',
                'OoPDefl','IPDefl','TeetDefl','Azimuth',
                'RotSpeed','TTDspFA','TTDspSS']
    
    return IC_keys
	
 
def GetFirstWind(WindPath):
    """ First wind speed from file
    
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
    

# =============================================================================
#                  Utilities copied from PyTurbSim
# =============================================================================
# Subsequent code copied from PyTurbSim functions to load turbsim output.
# Modified only slightly to be placed in this module.
#
# Code author:    Levi Kilcher
# TurbSim GitHub: http://lkilcher.github.io/pyTurbSim/

def readModel(fname, ):
    """
    Read a TurbSim data and input file and return a
    :class:`tsdata <pyts.main.tsdata>` data object.

    Parameters
    ----------
    fname : str
            The filename to load.
            If the file ends in:
              .bl or .wnd,  the file is assumed to be a bladed-format file.
              .bts, the file is assumed to be a TurbSim-format file.
    Returns
    -------
    turb : :class:`numpy.ndarray`
             [3 x n_z x n_y x n_t] array of wind velocity values
    """
    
    if (fname.endswith('wnd')):
        return bladed(fname,)
    elif (fname.endswith('bl')):
        return bladed(fname,)
    elif (fname.endswith('bts')):
        return turbsim(fname,)

    # Otherwise try reading it as a .wnd file.
    bladed(fname)  # This will raise an error if it doesn't work.    
    
    
def bladed(fname,):
    """
    Read Bladed format (.wnd, .bl) full-field time-series binary data files.

    Parameters
    ----------
    fname : str
            The filename from which to read the data.

    Returns
    -------
    turb : :class:`numpy.ndarray`
             [3 x n_z x n_y x n_t] array of wind velocity values

    """
    e = '<'         # define "endianness"
    fname = checkname(fname, ['.wnd', '.bl'])
    with file(fname, 'rb') as fl:
        junk, nffc, ncomp, lat, z0, center = unpack(e + '2hl3f', fl.read(20))
        if junk != -99 or nffc != 4:
            raise IOError("The file %s does not appear to be a valid 'bladed (.bts)' format file."
                          % fname)
        ti = np.array(unpack(e + '3f', fl.read(12))) / 100
        dz, dy, dx, n_f, uhub = unpack(e + '3flf', fl.read(20))
        n_t = int(2 * n_f)
        fl.seek(12, 1)  # Unused bytes
        clockwise, randseed, n_z, n_y = unpack(e + '4l', fl.read(16))
        fl.seek(24, 1)  # Unused bytes
        nbt = ncomp * n_y * n_z * n_t
        turb = np.rollaxis(np.fromstring(fl.read(2 * nbt), dtype=np.int16)
                          .astype(np.float32).reshape([ncomp,
                                                       n_y,
                                                       n_z,
                                                       n_t], order='F'),
                          2, 1)
    turb[0] += 1000.0 / ti[0]
    turb /= 1000. / (uhub * ti[:, None, None, None])
    # Create the grid object:
    dt = dx / uhub
    # Determine the clockwise value.
    if clockwise == 0:
        try:
            d = sum_scan(convname(fname, '.sum'))
            clockwise = d['clockwise']
        except IOError:
            warn("Value of 'CLOCKWISE' not specified in binary file, "
                 "and no .sum file found. Assuming CLOCKWISE = True.")
            clockwise = True
        except KeyError:
            warn("Value of 'CLOCKWISE' not specified in binary file, "
                 "and %s has no line containing 'clockwise'. Assuming "
                 "CLOCKWISE = True." % convname(fname, '.sum'))
            clockwise = True
    else:
        clockwise = bool(clockwise - 1)
    if clockwise:
        # flip the data back
        turb = turb[:, :, ::-1, :]

    return turb, dt


def turbsim(fname):
    """
    Read TurbSim format (.bts) full-field time-series binary
    data files.

    Parameters
    ----------
    fname : str
            The filename from which to read the data.

    Returns
    -------
    turb : :class:`numpy.ndarray`
             [3 x n_z x n_y x n_t] array of wind velocity values

    """
    e = '<'         # define "endianness"
    fname = checkname(fname, ['.bts'])
    u_scl = np.zeros(3, np.float32)
    u_off = np.zeros(3, np.float32)
    fl = file(fname, 'rb')
    (junk,
     n_z,
     n_y,
     n_tower,
     n_t,
     dz,
     dy,
     dt,
     uhub,
     zhub,
     z0,
     u_scl[0],
     u_off[0],
     u_scl[1],
     u_off[1],
     u_scl[2],
     u_off[2],
     strlen) = unpack(e + 'h4l12fl', fl.read(70))    
    desc_str = fl.read(strlen)  # skip these bytes.
    # load turbulent field
    nbt = 3 * n_y * n_z * n_t
    turb = np.rollaxis(np.fromstring(fl.read(2 * nbt), dtype=np.int16).astype(
        np.float32).reshape([3, n_y, n_z, n_t], order='F'), 2, 1)
    turb -= u_off[:, None, None, None]
    turb /= u_scl[:, None, None, None]
    return turb

def sum_scan(filename,):
    """
    Scan a sum file for specific variables.

    Parameters
    ----------
    filename : string
        The file to scan.

    Returns
    -------
    out : dict
        A dictionary of values identified.
    """
    # Currently this routine only searches for 'clockwise'.
    out = dict()
    with open(checkname(filename, ['.sum', '.SUM']), 'r') as infl:
        for ln in infl:
            ln = ln.lower()
            if 'clockwise' in ln.lower():
                v = ln.split()[0]
                if v in ['t', 'y']:
                    out['clockwise'] = True
                else:
                    out['clockwise'] = False
    return out
    

def convname(fname, extension=None):
    """
    Change the file extension.
    """
    if extension is None:
        return fname
    if extension != '' and not extension.startswith('.'):
        extension = '.' + extension
    return fname.rsplit('.', 1)[0] + extension

def checkname(fname, extensions=[]):
    """Test whether fname exists.

    If it does not, change the file extension in the list of
    extensions until a file is found. If no file is found this
    function raises IOError.
    """
    if os.path.isfile(fname):
        return fname
    if isinstance(extensions, basestring):
        # If extensions is a string make it a single-element list.
        extensions = [extensions]
    for e in extensions:
        fnm = convname(fname, e)
        if os.path.isfile(fnm):
            return fnm
    raise IOError("No such file or directory: '%s', and no "
                  "files found with specified extensions." % fname)
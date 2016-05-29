"""
nenzfr_input_utils.py -- Functions to check the nenzfr input dictionary that
    is pulled in from the config file and add default vales where required.

.. Author: Chris James

   Modified: Luke Doherty 03-Oct-2014
             Added checks/defaults for 'reactionSchemeFile',
             'energySchemeFile' and nitrogen gas name (can be
             either 'n2', 'N2', 'nitrogen'). 
"""
import os

#---------------------------------------
# Utility Functions

def check_reaction_scheme(cfg):
    """Takes the config dictionary and checks that the reaction scheme file
    is in the current directory. If not, and if we recognise the gas-name, we
    will copy the relevant scheme from the repository.

    """
    if not os.path.exists(cfg['reactionSchemeFile']):
        # Set directory path to where reaction-scheme files are stored
        basePath = '$HOME/cfcfd3/lib/gas/reaction-schemes/'
        
        # If the reaction-scheme file is one that we know, we copy it from the repository
        # while also checking that some of the 'hard-coded' settings are correct. If it 
        # is not one that we know, we simply copy it without modification.
        if cfg['reactionSchemeFile'] in ['nitrogen-5sp-6r.lua']:
            fp = open(os.path.expandvars(basePath+'nitrogen/nitrogen-5sp-6r.lua'),'r')
            fout = open('./nitrogen-5sp-6r.lua','w')
            for line in fp.readlines():
                items = line.strip().split()
                if (len(items)>0 and items[0] == 'WITH_IONIZATION' and items[-1] is not 'false'):
                    fout.write('WITH_IONIZATION = false\n')
                else:
                    fout.write('%s' % line)
            fout.close()

        elif cfg['reactionSchemeFile'] in ['gupta_etal_air_reactions.lua']:
            fp = open(os.path.expandvars(basePath+'air/gupta_etal_air_reactions.lua'),'r')
            fout = open('./gupta_etal_air_reactions.lua','w')
            for line in fp.readlines():
                items = line.strip().split()
                if (len(items)>0 and items[0] == 'NO_SPECIES' and items[-1] is not 5):
                    fout.write('NO_SPECIES = 5\n')
                else:
                    fout.write('%s' % line)
            fout.close()
        
        else: 
            # gasName:folderName
            folderNames = {'air5species': 'air', 'n2': 'nitrogen', cfg['gasName']:cfg['gasName']}
            fp = open(os.path.expandvars(basePath+folderNames[cfg['gasName']]+'/'+cfg['reactionSchemeFile']))
            fout = open('./'+cfg['reactionSchemeFile'],'w')
            for line in fp.readlines():
                fout.write('%s' % line)
            fout.close()
    return

#------------------------------------------------
# Main function 

def input_checker(cfg):
    """Takes the config dictionary and checks that everything is there. Will 
    bail out if important variables are missing, and set default values for 
    those that are not crucial.
       
    Returns the config dictionary after it has been checked over and will
    tell nenzfr to bail out if it found an issue.
    
    """
    
    # Define default the default nozzles
    default_nozzles = {'t4-m4':{'contourFileName':'Bezier-control-pts-t4-m4.data',
                                'include_throat_block':True, 'Lthr':12.5e-3,
                                'truncate_nozzle':False, 'Lnoz':5.118100e-01,
                                'nni':600, 'nnj':40, 'nbi':60, 'nbj':1, 'bx':1.10, 'by':1.02,
                                'fully_contoured_nozzle':False},
                       #
                       't4-m4b':{'contourFileName':'Bezier-control-pts-t4-m4b.data',
                                'include_throat_block':True, 'Lthr':12.5e-3,
                                'truncate_nozzle':False, 'Lnoz':5.118100e-01,
                                'nni':600, 'nnj':40, 'nbi':30, 'nbj':2, 'bx':1.10, 'by':1.02,
                                'fully_contoured_nozzle':True},
                       #
                       't4-m6':{'contourFileName':'Bezier-control-pts-t4-m6.data',
                                'include_throat_block':True, 'Lthr':12.5e-3,
                                'truncate_nozzle':False, 'Lnoz':7.970280e-01,
                                'nni':600, 'nnj':40, 'nbi':60, 'nbj':1, 'bx':1.10, 'by':1.02,
                                'fully_contoured_nozzle':False},
                       #
                       't4-m7':{'contourFileName':'Bezier-control-pts-t4-m7.data',
                                'include_throat_block':True, 'Lthr':20.0e-3,
                                'truncate_nozzle':False, 'Lnoz':1.00,
                                'nni':600, 'nnj':60, 'nbi':60, 'nbj':4, 'bx':1.10, 'by':1.02,
                                'fully_contoured_nozzle':True},
                       #
                       't4-m8b':{'contourFileName':'Bezier-control-pts-t4-m8.data',
                                'include_throat_block':True, 'Lthr':8.2e-3,
                                'truncate_nozzle':False, 'Lnoz':1.11,
                                'nni':1440, 'nnj':80, 'nbi':144, 'nbj':4, 'bx':1.05, 'by':1.002,
                                'fully_contoured_nozzle':True},
                       #
                       't4-m10b':{'contourFileName':'Bezier-control-pts-t4-m10.data',
                                'include_throat_block':True, 'Lthr':10.0e-3,
                                'truncate_nozzle':True, 'Lnoz':1.633197,
                                'nni':1800, 'nnj':100, 'nbi':180, 'nbj':4, 'bx':1.05, 'by':1.002,
                                'fully_contoured_nozzle':True},
                       #
                       'x2-m10':{'contourFileName':'Bezier-control-pts-x2-m10.data',
                                'include_throat_block':True, 'Lthr':42.5e-3,
                                'truncate_nozzle':False, 'Lnoz':1.40,
                                'nni':1800, 'nnj':100, 'nbi':180, 'nbj':4, 'bx':1.05, 'by':1.002,
                                'fully_contoured_nozzle':True},
                       #
                       'x3-m10':{'contourFileName':'Bezier-control-pts-x3-m10.data',
                                'include_throat_block':True, 'Lthr':91.3e-3,
                                'truncate_nozzle':False, 'Lnoz':2.59,
                                'nni':1800, 'nnj':100, 'nbi':180, 'nbj':4, 'bx':1.05, 'by':1.002,
                                'fully_contoured_nozzle':True},
                       #
                       'x3r-m7':{'contourFileName':'Bezier-control-pts-x3r-m7.data',
                                'include_throat_block':True, 'Lthr':60.0e-3,
                                'truncate_nozzle':False, 'Lnoz':2.857143,
                                'nni':1800, 'nnj':100, 'nbi':180, 'nbj':4, 'bx':1.05, 'by':1.002,
                                'fully_contoured_nozzle':True}
                       }
    
    print "Checking nenzfr inputs."
    
    cfg['bad_input'] = False
    
    if 'facility' not in cfg:
        print "Need to specify a type of facility to simulate."
        cfg['bad_input'] = True
    
    if cfg['facility'] == 'reflected-shock-tunnel':
        if 'p1' not in cfg:
            print "Need to supply a float value for p1."
            cfg['bad_input'] = True
        if 'T1' not in cfg:
            print "Need to supply a float value for T1."
            cfg['bad_input'] = True
        if 'Vs' not in cfg:
            print "Need to supply a float value for Vs."
            cfg['bad_input'] = True
        if 'pe' not in cfg:
            print "Need to supply a float value for pe."
            cfg['bad_input'] = True
        #now need to check that the inputs are actually numbers
        if not isinstance(cfg['p1'], (int, float)):
            print "Need to supply a float value for p1."
            cfg['bad_input'] = True
        if not isinstance(cfg['T1'], (int, float)):
            print "Need to supply a float value for T1."
            cfg['bad_input'] = True 
        if not isinstance(cfg['Vs'], (int, float)):
            print "Need to supply a float value for Vs."
            cfg['bad_input'] = True
        if not isinstance(cfg['pe'], (int, float)):
            print "Need to supply a float value for pe."
            cfg['bad_input'] = True
            
    elif cfg['facility'] == 'expansion-tube' and 'pitot_input_file' not in cfg:
        if 'p7' not in cfg:
            print "Need to supply a float value for p7."
            cfg['bad_input'] = True
        if 'T7' not in cfg:
            print "Need to supply a float value for T7."
            cfg['bad_input'] = True
        if 'V7' not in cfg:
            print "Need to supply a float value for V7."
            cfg['bad_input'] = True    
        #now need to check that the inputs are actually numbers
        if not isinstance(cfg['p7'], (int, float)):
            print "Need to supply a float value for p7."
            cfg['bad_input'] = True
        if not isinstance(cfg['T7'], (int, float)):
            print "Need to supply a float value for T7."
            cfg['bad_input'] = True 
        if not isinstance(cfg['V7'], (int, float)):
            print "Need to supply a float value for V7."
            cfg['bad_input'] = True
                  
    elif cfg['facility'] == 'gun-tunnel':
        if 'p0' not in cfg:
            print "Need to supply a float value for p0."
            cfg['bad_input'] = True
        if 'T0' not in cfg:
            print "Need to supply a float value for T0."
            cfg['bad_input'] = True
        #now need to check that the inputs are actually numbers
        if not isinstance(cfg['p0'], (int, float)):
            print "Need to supply a float value for p0."
            cfg['bad_input'] = True
        if not isinstance(cfg['T0'], (int, float)):
            print "Need to supply a float value for T0."
            cfg['bad_input'] = True 
              
    if 'gasName' not in cfg:
        cfg['gasName'] = 'air5species'
        print "No gas model chosen."
        print "    Setting it to default value of '{0}'.".format(cfg['gasName'])

    if cfg['gasName'] in ['n2','N2','nitrogen']:
        cfg['gasName'] = 'n2'

    if 'chemModel' not in cfg:
        cfg['chemModel'] = 'neq'
        print "No chemistry model chosen. Setting it to default value of '{0}'."\
        .format(cfg['chemModel'])
    
    if 'reactionSchemeFile' not in cfg:
        if cfg['chemModel'] in ['neq','tc-neq']:
            if cfg['gasName'] in ['air']:
                # I have deprecated this combination
                print ("The combination of gasName='air' and chemModel='neq' ('tc-neq')\n"
                      "    is not supported due to a lack of suitable reaction scheme.")
                cfg['reactionSchemeFile'] = 'Invalid-inputs'
                cfg['bad_input'] = True
            elif cfg['gasName'] in ['air5species']:
                cfg['reactionSchemeFile'] = 'gupta_etal_air_reactions.lua'
                # check that the file is present in working directory
                check_reaction_scheme(cfg)
            elif cfg['gasName'] in ['n2']:
                cfg['reactionSchemeFile'] = 'nitrogen-5sp-6r.lua'
                # check that the file is present in the working directory
                check_reaction_scheme(cfg)
            else:
                print ("The specified gasName is not one of the defaults, you MUST "
                      "therefore define a 'reactionSchemeFile'")
                cfg['bad_input'] = True
        else:
            cfg['reactionSchemeFile'] = 'Not-needed'
        #
        #print "No reaction scheme specified." 
        #print "    Setting it to '{0}'.".format(cfg['reactionSchemeFile'])
    print "Reaction scheme is '{0}'.".format(cfg['reactionSchemeFile'])
    
    if 'thermalSchemeFile' not in cfg:
        if cfg['chemModel'] in ['tc-neq']:
            print ("No energy-exchange scheme file specified.\n"
                  "    The default energy-schemes have yet to be defined.\n"
                  "    This is coming...")
            cfg['bad_input'] = True
        else:
            cfg['thermalSchemeFile'] = 'Not-needed'
    #
    print "Energy-exchange scheme is '{0}'.".format(cfg['thermalSchemeFile'])
            
    if cfg['facility'] == 'reflected-shock-tunnel' and 'areaRatio' not in cfg:
        cfg['areaRatio'] = 1581.165
        print "No area ratio selected for estcj calculation."
        print "    Setting area ratio to default value of {0}".format(cfg['areaRatio'])
        
    if 'jobName' not in cfg:
        cfg['jobName'] = 'nozzle'
        print "No job name set. Setting default job name of '{0}'."\
        .format(cfg['jobName'])
         
    if 'exitSliceFileName' not in cfg:
        cfg['exitSliceFileName'] = 'nozzle-exit.data'
        print "No filename specified for file that holds nozzle-exit data."
        print "    Setting it to default value of '{0}'.".format(cfg['exitSliceFileName'])
        
    if 'justStats' not in cfg:
        cfg['justStats'] = False
        print "Switch to skip detailed calculations not specified."
        print "    Setting it to default value of justStats = {0}".format(cfg['justStats'])

    if 'noStats' not in cfg:
        cfg['noStats'] = False
        print "Switch to skip flow statistic calculations not specified."
        print "    Setting it to default value of noStats = {0}".format(cfg['noStats'])
    
    if 'blockMarching' not in cfg:
        cfg['blockMarching'] = True
        print "Switch to use or not use Block Marching mode not specified."
        print "    Setting it to default value of blockMarching = {0}".format(cfg['blockMarching'])
     
    #config_backwards_compatability = {'Bezier-control-pts-t4-m4.data' : 't4-m4',
    #                                  'Bezier-control-pts-t4-m6.data' : 't4-m6',
    #                                  'Bezier-control-pts-t4-m7.data' : 't4-m7',
    #                                  'Bezier-control-pts-t4-m8.data' : 't4-m8b',
    #                                  'Bezier-control-pts-t4-m10.data': 't4-m10b',
    #                                  'Bezier-control-pts-x2-m10.data': 'x2-m10',
    #                                  'Bezier-control-pts-x3-m10.data': 'x3-m10'}   
    #if 'nozzle' not in cfg and cfg['contourFileName'] in config_backwards_compatability:
    #    cfg['nozzle'] = nozzle_backwards_compatability[cfg['contourFileName']]
    #
    if 'nozzle' not in cfg:
        cfg['bad_input'] = True
        print "Your must define a nozzle (either by way of using a default nozzle," 
        print "    defining your own or by specifying an external grid file)."
    else:
        nozzle = cfg['nozzle']
        if nozzle in default_nozzles:
            # Use a default nozzle but perhaps with altered grid dimensions
            if 'nni' not in cfg:
                cfg['nni'] = default_nozzles[nozzle]['nni']
                print "Number of axial cells not set."
                print "    Setting it to default value of {0}.".format(cfg['nni'])
        
            if 'nnj' not in cfg:
                cfg['nnj'] = default_nozzles[nozzle]['nnj']
                print "Number of radial cells not set."
                print "    Setting it to default value of {0}.".format(cfg['nnj'])
     
            if 'nbi' not in cfg:
                cfg['nbi'] = default_nozzles[nozzle]['nbi']
                print "Number of axial blocks for the divergence section (nozzle_blk) not set"
                print "    Setting it to default value of {0}.".format(cfg['nbi']) 
            
            if 'nbj' not in cfg:
                cfg['nbj'] = default_nozzles[nozzle]['nbj']
                print "Number of radial blocks not set."
                print "    Setting it to default value of {0}.".format(cfg['nbj'])
            
            if 'bx' not in cfg:
                cfg['bx'] = default_nozzles[nozzle]['bx']
                print "Clustering in the axial direction not set."
                print "    Setting it to default value of {0}.".format(cfg['bx'])
            
            if 'by' not in cfg:
                cfg['by'] = default_nozzles[nozzle]['by']
                print "Clustering in the radial direction not set."
                print "    Setting it to default value of {0}.".format(cfg['by'])
            
            # Assign the rest of the required inputs. We don't allow users access to these when 
            # using a default nozzle. If they want to change them, they should use a custom nozzle
            # name
            cfg['contourFileName'] = default_nozzles[nozzle]['contourFileName']
            cfg['gridFileName'] = 'None'
            cfg['include_throat_block'] = default_nozzles[nozzle]['include_throat_block']
            cfg['Lthr'] = default_nozzles[nozzle]['Lthr']
            cfg['truncate_nozzle'] = default_nozzles[nozzle]['truncate_nozzle']
            cfg['Lnoz'] = default_nozzles[nozzle]['Lnoz']
            cfg['fully_contoured_nozzle'] = default_nozzles[nozzle]['fully_contoured_nozzle']
            
        else:
            # Use a non-standard grid. We must either provide a grid-file to load in or define
            # everything we need to build the grid within nenzfr/eilmer
            if 'contourFileName' not in cfg:
                if 'gridFileName' not in cfg:
                    cfg['bad_input'] = True
                    print "You appear to want to use a non-default nozzle but have not defined"
                    print "    either a contourFileName or a gridFileName. Check your input."
                else:
                    if cfg['gridFileName'] in ['None']:
                        cfg['bad_input'] = True
                        print "You appear to want to use a non-default nozzle and import a grid"
                        print "    but have not defined a grid file. Check your input."
                    else:
                        # The user will be importing a predefined grid so we don't need various
                        # inputs.
                        cfg['nni'] = None
                        cfg['nnj'] = None
                        cfg['nbi'] = None
                        cfg['nbj'] = None
                        cfg['bx'] = None
                        cfg['by'] = None
                        cfg['include_throat_block'] = None
                        cfg['Lthr'] = None
                        cfg['trucate_nozzle'] = None
                        cfg['Lnoz'] = None
                        cfg['fully_contoured_nozzle'] = None
            else:
                if cfg['contourFileName'] in ['None']:
                    if 'gridFileName' not in cfg:
                        cfg['bad_input'] = True
                        print "You appear to want to use a non-default nozzle but have not defined"
                        print "    either the contourFileName or gridFileName appropriately. Check your input."
                    else:
                        if cfg['gridFileName'] in ['None']:
                            cfg['bad_input'] = True
                            print "You appear to want to use a non-default nozzle but have not defined"
                            print ("    either the contourFileName or gridFilename appropriately. "
                                  "Check your input.")
                        else:
                            # The user will be importing a predefined grid so we don't need various
                            # inputs.
                            cfg['nni'] = None
                            cfg['nnj'] = None
                            cfg['nbi'] = None
                            cfg['nbj'] = None
                            cfg['bx'] = None
                            cfg['by'] = None
                            cfg['include_throat_block'] = None
                            cfg['Lthr'] = None
                            cfg['trucate_nozzle'] = None
                            cfg['Lnoz'] = None
                            cfg['fully_contoured_nozzle'] = None    
                else:
                    # A contour file has been specified so now we check that the rest of the 
                    # grid-related inputs are present
                    cfg['gridFileName'] = 'None'
                 
                    if 'nni' not in cfg:
                        cfg['bad_input'] = True
                        print "To build a custom nozzle grid you must specify the"
                        print "    number of axial cells (nni). Check your input."
            
                    if 'nnj' not in cfg:
                        cfg['bad_input'] = True
                        print "To build a custom nozzle grid you must specify the"
                        print "    number of radial cells (nnj). Check your input."
            
                    if 'nbi' not in cfg:
                        cfg['bad_input'] = True
                        print "To build a custom nozzle grid you must specify the"
                        print "    number of axial blocks (nbi). Check your input."
                
                    if 'nbj' not in cfg:
                        cfg['bad_input'] = True
                        print "To build a custom nozzle grid you must specify the"
                        print "    number of radial blocks (nbj). Check your input."
                
                    if 'bx' not in cfg:
                        cfg['bad_input'] = True
                        print "To build a custom nozzle grid you must specify the"
                        print "    axial cell clustering (bx). Check your input."
            
                    if 'by' not in cfg:
                        cfg['bad_input'] = True
                        print "To build a custom nozzle grid you must specify the"
                        print "    radial cell clustering (by). Check your input."
                
                    if 'include_throat_block' not in cfg:
                        cfg['bad_input'] = True
                        print "To build a custom nozzle grid you must specify whether or not"
                        print "    to include a throat block (include_throat_block=True/False)."
                        print "    Check your input."
                    else:
                        if cfg['include_throat_block']:
                            if 'Lthr' not in cfg:
                                cfg['bad_input'] = True
                                print "Length of the throat block (Lthr) not set."
                                print "    Check your input."
                        else:
                            if 'Lthr' not in cfg:
                                cfg['Lthr'] = 0.0
                    
                    if 'truncate_nozzle' not in cfg:
                        cfg['bad_input'] = True
                        print "To build a custom nozzle grid you must specify whether or not"
                        print "    nozzle is trucated wrt the Bezier-contour (truncate_nozzle=True/False)."
                        print "    Check your input."
                    else:
                        if cfg['truncate_nozzle']:
                            if 'Lnoz' not in cfg:
                                cfg['bad_input'] = True
                                print "Length of truncated nozzle (Lnoz) not set."
                                print "    Check your input."
                        else:
                            if 'Lnoz' not in cfg:
                                cfg['Lnoz'] = 1.0
                
                    if 'fully_contoured_nozzle' not in cfg:
                        cfg['bad_input'] = True
                        print "To build a custom nozzle grid you must specify if the nozzle is"
                        print "    fully contoured (or not). Check your input."
                
    if 'max_time' not in cfg:
        cfg['max_time'] = 6.0e-3
        print "Overall simulation time for nozzle flow not set."
        print "    Setting it to default value of {0}.".format(cfg['max_time'])

    if 'max_step' not in cfg:
        cfg['max_step'] = 800000
        print "Maximum simulation steps allowed not set."
        print "    Setting to to default value of {0}.".format(cfg['max_step'])
        
    if 'Tw' not in cfg:
        cfg['Tw'] = 300.0
        print "Nozzle wall temperature not set. Setting it to default value of {0} K."\
        .format(cfg['Tw'])
        
    if 'BLTrans' not in cfg:
        cfg['BLTrans'] = "x_c[-1]*1.1"
        print "Transition location for the boundary layer (used to define the"
        print "    turbulent portion of the nozzle) not set."
        print "    Setting to to default value of > nozzle length (ie. laminar nozzle)."
        
    if 'TurbVisRatio' not in cfg:
        cfg['TurbVisRatio'] = 100.0
        print "Turbulent to Laminar Viscoscity ratio not set."
        print "    Setting it to default value of {0}.".format(cfg['TurbVisRatio'])
        
    if 'TurbInten' not in cfg:
        cfg['TurbInten'] = 0.05
        print "Turbulence intensity at the throat not set."
        print "    Setting it to default value of {0}.".format(cfg['TurbInten'])
        
    if 'coreRfraction' not in cfg:
        cfg['coreRfraction'] = 2.0/3.0
        print "Radius of core flow as a fraction of the nozzle exit radius not set."
        print "    Setting it to default value of {0}.".format(cfg['coreRfraction'])
    
    print "Done checking nenzfr inputs."
        
    return cfg
    
# Function below is used by nenzfr_perturbed.py
    
def nenzfr_perturbed_input_checker(cfg):
    """Takes the config dictionary and checks that all of the extra
    nenzfr perturbed inputs are in there. Will bail out if important variables 
    are missing, and set default values for those that are not crucial.
       
    Returns the config dictionary after it has been checked over and will
    tell nenzfr perturbed to bail out if it finds an issue.
    
    """
    
    print "Checking specific nenzfr perturbed inputs."
    
    if 'Cluster' not in cfg:
        cfg['Cluster'] = 'Mango'
        print "No Cluster slected. Setting to to default value of '{0}."\
        .format(cfg['Cluster'])
    
    if 'runCMD' not in cfg and cfg['Cluster'] == 'Mango':
        cfg['runCMD'] = './'
        print "Run command for the cluster not specified."
        print "Setting to a default Mango run command ('{0}')"\
        .format(cfg['runCMD'])
    elif 'runCMD' not in cfg and cfg['Cluster'] == 'Barrine':
        cfg['runCMD'] = 'qsub '
        print "Run command for the cluster not specified."
        print "Setting to a default Barrine run command ('{0}')"\
        .format(cfg['runCMD'])        
                
    if 'levels' not in cfg:
        cfg['levels'] = 3
        print "Levels to use not specified. Setting it to default value of {0}"\
        .format(cfg['levels'])
        
    if 'createRSA' not in cfg:
        print "Response surface approximation variable not set."
        print "Setting it to False."
        
    # Convert to integer. NB. I don't specify the type as int in the option
    # as this negates the ability to have choices.
    if cfg['levels'] == '3-reduced':
        if cfg['createRSA']:
            cfg['levels'] = 2.5
        else:
            print "setting levels to '3-reduced' is only valid when creating a"
            print "response surface. Changing this to levels = 3."
            cfg['levels'] = 3
    else:
       cfg['levels'] = int(cfg['levels'])
       
    if 'perturbedVariables' not in cfg and not 'createRSA':
        # will not bail out here if they want to create an RSA (as it will make that list for them)
        print "You have not specified a perturbed variables list."
        print "Bailing out."
        cfg['bad_input'] = True
       
    # Re-set the perturbed variables list is the user wants to do a response surface calculation
    if cfg['createRSA']:
        print "Create Response Surface (RS) option selected. Perturbing only Vs and pe"
        print "NOTE: This will override other selected perturbed variables."
        perturbedVariables = ['Vs','pe']
        
    if cfg['BLTrans'] == "x_c[-1]*1.1": # Laminar nozzle. Don't perturb BLTrans
        # For fully laminar nozzles the inflow turbulent parameters are forced
        # to zero so we shouldn't perturb them.
        print "Nozzle is assumed laminar (based on input BLTrans), therefore\n"+\
                  "BLTrans, TurbVisRatio and TurbInten are not perturbed"
        print "Deleting any of these variables from the perturbed variables list."
        if 'BLTrans' in cfg['perturbedVariables']:
            del cfg['perturbedVariables'][cfg['perturbedVariables'].index('BLTrans')]
        if 'TurbVisRatio' in cfg['perturbedVariables']:
            del cfg['perturbedVariables'][cfg['perturbedVariables'].index('TurbVisRatio')]
        if 'TurbInten' in cfg['perturbedVariables']:
            del cfg['perturbedVariables'][cfg['perturbedVariables'].index('TurbInten')]
    elif isinstance(cfg['BLTrans'], float) or isinstance(cfg['BLTrans'], int):
	pass
    else:
        BLTransValues = cfg['BLTrans'].strip(']').strip('[').split(',')
        if float(BLTransValues[0]) == 0.0: # Full turbulent nozzle. Don't perturb BLTrans
            print "Nozzle is assumed fully turbulent, therefore BLTrans\n"+\
                      "will not be perturbed"
        if 'BLTrans' in cfg['perturbedVariables']:
            del cfg['perturbedVariables'][cfg['perturbedVariables'].index('BLTrans')]
            
    # Now go through and check our perturbed variables, and assign defaults 
    # where asked to.
    
    cfg['perturbedDict'] = {} #this will get populated below
    
    # This will make a dictionary for each variable to be perturbed that will go
    # inside the perturbed variables dictionary
    
    for variable in cfg['perturbedVariables']:
        # Check the cfg file for what we need and set defaults if we don't find anything
        if variable+'_perturbation' not in cfg:
            cfg[variable+'_perturbation'] = 'relative'
            print "You have not specified the type of perturbation for {0}.".format(variable)
            print "Setting perturbation to default value of '{0}'."\
            .format(cfg[variable+'_perturbation'])
        if variable+'_deltas' not in cfg:
            cfg[variable+'_deltas'] = [cfg['defaultPerturbations'][variable]]
            print "You have not specified any perturbations for {0}.".format(variable)
            print "Setting pertubation to default value for this variable ({0})."\
            .format(cfg[variable+'_deltas'])
            
        # if a single delta is specified not as a list, make it one
        if not isinstance(cfg[variable+'_deltas'], list):
            print "Delta value for {0} not a list. Making it one.".format(variable)
            cfg[variable+'_deltas'] = [cfg[variable+'_deltas']]  
            
        bad_delta = False    
        #check the defaults are numbers, change them all if one isn't                
        for delta in cfg[variable+'_deltas']:
            if not isinstance(delta, (int, float)):
                print "One of the deltas for {0} is not a number.".format(variable)
                print "Making them all floats if possible."
                bad_delta = True
            if bad_delta:
                new_delta_list = []
                for delta in cfg[variable+'_deltas']:
                    new_delta_list.append(float(delta))
      
        # Time to start actually doing thigs.
        
        # Build value in the perturbed variable dictionary
        cfg['perturbedDict'][variable] = {}
        
        #now start to populate it
        cfg['perturbedDict'][variable]['nominal'] = cfg[variable]
        cfg['perturbedDict'][variable]['pertubation'] = cfg[variable+'_perturbation']
        cfg['perturbedDict'][variable]['deltas'] = cfg[variable+'_deltas']
        
        #now need to run this through the set_perturbed_values function
        
        from nenzfr_perturbed import set_perturbed_values
        
        # The current perturbed dictionary value is then overwritten by the 
        # perturbed variable list returned by set_perturbed_values function.
        
        cfg['perturbedDict'][variable] = set_perturbed_values(variable,\
            [cfg['perturbedDict'][variable]['nominal']] + cfg['perturbedDict'][variable]['deltas'],
            cfg['defaultPerturbations'], cfg['perturbedDict'][variable]['pertubation'],
            cfg['levels'])
            
        print "Done with perturbed variable {0}.".format(variable)

    print "Done checking specific nenzfr perturbed inputs."
               
    return cfg
    
#Function below is used for nenzfr_sensitivity.py
        
def nenzfr_sensitivity_input_checker(cfg):
    """Takes the config dictionary and checks that all of the extra
    nenzfr perturbed inputs are in there. Will bail out if important variables 
    are missing, and set default values for those that are not crucial.
       
    Returns the config dictionary after it has been checked over and will
    tell nenzfr perturbed to bail out if it finds an issue.
    
    """
    
    print "Checking specific nenzfr sensitivity inputs."

    # Make empty input uncertainties dictionary that we will populate
    cfg['inputUncertainties'] = {}        
        
    for variable in cfg['perturbedVariables']:
        
        # Set a default relative uncertainty if nothing is specified
        if variable+'_rel_uncertainty' not in cfg:
            cfg[variable+'_rel_uncertainty'] = [cfg['default_rel_uncertainties'][variable]]
            print "You have not specified a relative uncertainty for {0}.".format(variable)
            print "Setting relative uncertainty to default value for this variable ({0})."\
            .format( cfg[variable+'_rel_uncertainty'])
        
        #Then add it to the dictionary.
        
        cfg['inputUncertainties'][variable] = cfg[variable+'_rel_uncertainty']
        
        print "Done with perturbed variable {0}.".format(variable)
        
    print "Done checking specific nenzfr sensitivity inputs."
        
    return cfg
        
        

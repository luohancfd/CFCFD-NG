"""
nenzfr_input_utils.py -- Functions to check the nenzfr input dictionary that
    is pulled in from the config file and add default vales where required.

.. Authors: Chris James
"""

def input_checker(cfg):
    """Takes the config dictionary and checks that everything is there. Will 
    bail out if important variables are missing, and set default values for 
    those that are not crucial.
       
    Returns the config dictionary after it has been checked over and will
    tell nenzfr to bail out if it found an issue.
    
    """
    
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
        print "No gas model chosen. Setting it to default value of '{0}'."\
        .format(cfg['gasName'])

    if 'chemModel' not in cfg:
        cfg['chemModel'] = 'neq'
        print "No chemistry model chosen. Setting it to default value of '{0}'."\
        .format(cfg['chemModel'])
        
    if cfg['facility'] == 'reflected-shock-tunnel' and 'areaRatio' not in cfg:
        cfg['areaRatio'] = 1581.165
        print "No area ratio selected for estcj calculation."
        print "Setting area ratio to default value of {0}".format(cfg['areaRatio'])
        
    if 'jobName' not in cfg:
        cfg['jobName'] = 'nozzle'
        print "No job name set. Setting default job name of '{0}'."\
        .format(cfg['jobName'])
        
    if 'contourFileName' not in cfg:
        cfg['contourFileName'] = 'Bezier-control-pts-t4-m10.data'
        print "No nozzle contour file selected. Choosing '{0}'."\
        .format(cfg['contourFileName'])
        
    if 'exitSliceFileName' not in cfg:
        cfg['exitSliceFileName'] = 'nozzle-exit.data'
        print "No filename specified for file that holds nozzle-exit data."
        print "Setting it to default value of '{0}'.".format(cfg['exitSliceFileName'])
        
    if 'justStats' not in cfg:
        cfg['justStats'] = False
        print "Switch to skip detailed calculations not specified."
        print "Setting it to default value of justStats = {0}".format(cfg['justStats'])

    if 'noStats' not in cfg:
        cfg['noStats'] = False
        print "Switch to skip flow statistic calculations not specified."
        print "Setting it to default value of noStats = {0}".format(cfg['noStats'])
    
    if 'blockMarching' not in cfg:
        cfg['blockMarching'] = True
        print "Switch to use or not use Block Marching mode not specified."
        print "Setting it to default value of blockMarching = {0}".format(cfg['blockMarching'])
        
    #these default values below are based on Luke's Mach 10 calculations.
    
    if 'nni' not in cfg:
        cfg['nni'] = 1800
        print "Number of axial cells not set. Setting it to default value of {0}."\
        .format(cfg['nni'])
        
    if 'nnj' not in cfg:
        cfg['nnj'] = 100
        print "Number of radial cells not set. Setting it to default value of {0}."\
        .format(cfg['nnj'])
    
    if 'nbi' not in cfg:
        cfg['nbi'] = 180
        print "Number of axial blocks for the divergence section (nozzle_blk) not set"
        print "Setting it to default value of {0}.".format(cfg['nbi']) 
        
    if 'nbj' not in cfg:
        cfg['nbj'] = 1
        print "Number of radial blocks not set. Setting it to default value of {0}."\
        .format(cfg['nbj'])
        
    if 'bx' not in cfg:
        cfg['bx'] = 1.10
        print "Clustering in the axial direction not set. Setting it to default value of {0}."\
        .format(cfg['bx'])
    
    if 'by' not in cfg:
        cfg['by'] = 1.002
        print "Clustering in the radial direction not set. Setting it to default value of {0}."\
        .format(cfg['by'])

    if 'max_time' not in cfg:
        cfg['max_time'] = 6.0e-3
        print "Overall simulation time for nozzle flow not set."
        print "Setting it to default value of {0}.".format(cfg['max_time'])

    if 'max_step' not in cfg:
        cfg['max_step'] = 800000
        print "Maximum simulation steps allowed not set."
        print "Setting to to default value of {0}.".format(cfg['max_step'])
        
    if 'Tw' not in cfg:
        cfg['Tw'] = 300.0
        print "Nozzle wall temperature not set. Setting it to default value of {0} K."\
        .format(cfg['Tw'])
        
    if 'BLTrans' not in cfg:
        cfg['BLTrans'] = "x_c[-1]*1.1"
        print "Transition location for the boundary layer (used to define the"
        print "turbulent portion of the nozzle) not set."
        print "Setting to to default value of > nozzle length (ie. laminar nozzle)."
        
    if 'TurbVisRatio' not in cfg:
        cfg['TurbVisRatio'] = 100.0
        print "Turbulent to Laminar Viscoscity ratio not set."
        print "Setting it to default value of {0}.".format(cfg['TurbVisRatio'])
        
    if 'TurbInten' not in cfg:
        cfg['TurbInten'] = 0.05
        print "Turbulence intensity at the throat not set."
        print "Setting it to default value of {0}.".format(cfg['TurbInten'])
        
    if 'coreRfraction' not in cfg:
        cfg['coreRfraction'] = 2.0/3.0
        print "Radius of core flow as a fraction of the nozzle exit radius not set."
        print "Setting it to default value of {0}.".format(cfg['coreRfraction'])
    
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
    
    if 'runCMD' not in cfg:
        cfg['runCMD'] = './'
        print "Run command for the cluster not specified."
        print "Setting to a default value of the Mango run command ('{0}')"\
        .format(cfg['runCMD'])
        print "This may or may not be correct."
        
    if 'Cluster' not in cfg:
        cfg['Cluster'] = 'Mango'
        print "No Cluster slected. Setting to to default value of '{0}."\
        .format(cfg['Cluster'])
        
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
        
        

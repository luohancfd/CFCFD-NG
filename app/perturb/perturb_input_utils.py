def perturb_input_checker(cfg):
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

    if 'VariablesToPerturb' not in cfg:
        print "You have not specified a list of which Variables are to be perturbed."
        print "Bailing out."
        cfg['bad_input'] = True
    
    # Re-set the perturbed variables list is the user wants to do a response surface calculation
    if cfg['createRSA']:
        print "Create Response Surface (RS) option selected."
        #print "NOTE: This will override other selected perturbed variables."
        #perturbedVariables = ['Vs','pe']

    #if cfg['BLTrans'] == "x_c[-1]*1.1": # Laminar nozzle. Don't perturb BLTrans
    #    # For fully laminar nozzles the inflow turbulent parameters are forced
    #    # to zero so we shouldn't perturb them.
    #    print "Nozzle is assumed laminar (based on input BLTrans), therefore\n"+\
    #              "BLTrans, TurbVisRatio and TurbInten are not perturbed"
    #    print "Deleting any of these variables from the perturbed variables list."
    #    if 'BLTrans' in cfg['perturbedVariables']:
    #        del cfg['perturbedVariables'][cfg['perturbedVariables'].index('BLTrans')]
    #    if 'TurbVisRatio' in cfg['perturbedVariables']:
    #        del cfg['perturbedVariables'][cfg['perturbedVariables'].index('TurbVisRatio')]
    #    if 'TurbInten' in cfg['perturbedVariables']:
    #        del cfg['perturbedVariables'][cfg['perturbedVariables'].index('TurbInten')]
    #elif isinstance(cfg['BLTrans'], float) or isinstance(cfg['BLTrans'], int):
    #    pass
    #else:
    #    BLTransValues = cfg['BLTrans'].strip(']').strip('[').split(',')
    #if float(BLTransValues[0]) == 0.0: # Full turbulent nozzle. Don't perturb BLTrans
    #        print "Nozzle is assumed fully turbulent, therefore BLTrans\n"+\
    #                  "will not be perturbed"
    #    if 'BLTrans' in cfg['perturbedVariables']:
    #        del cfg['perturbedVariables'][cfg['perturbedVariables'].index('BLTrans')]

    # Now go through and check our perturbed variables, and assign defaults 
    # where asked to.

    cfg['perturbedDict'] = {} #this will get populated below

    # This will make a dictionary for each variable to be perturbed that will go
    # inside the perturbed variables dictionary
    
    #for variable in cfg['perturbedVariables']:
    #    # Check the cfg file for what we need and set defaults if we don't find anything
    #    if variable+'_perturbation' not in cfg:
    #        cfg[variable+'_perturbation'] = 'relative'
    #        print "You have not specified the type of perturbation for {0}.".format(variable)
    #        print "Setting perturbation to default value of '{0}'."\
    #        .format(cfg[variable+'_perturbation'])
    #    if variable+'_deltas' not in cfg:
    #        cfg[variable+'_deltas'] = [cfg['defaultPerturbations'][variable]]
    #        print "You have not specified any perturbations for {0}.".format(variable)
    #        print "Setting pertubation to default value for this variable ({0})."\
    #        .format(cfg[variable+'_deltas'])

        # if a single delta is specified not as a list, make it one
    #    if not isinstance(cfg[variable+'_deltas'], list):
    #        print "Delta value for {0} not a list. Making it one.".format(variable)
    #        cfg[variable+'_deltas'] = [cfg[variable+'_deltas']]

    #    bad_delta = False
    #    #check the defaults are numbers, change them all if one isn't                
    #    for delta in cfg[variable+'_deltas']:
    #        if not isinstance(delta, (int, float)):
    #            print "One of the deltas for {0} is not a number.".format(variable)
    #            print "Making them all floats if possible."
    #            bad_delta = True
    #        if bad_delta:
    #            new_delta_list = []
    #            for delta in cfg[variable+'_deltas']:
    #                new_delta_list.append(float(delta))
    
    for variable in cfg['VariablesToPerturb']: 
        # Time to start actually doing thigs.

        # Build value in the perturbed variable dictionary
        #cfg['perturbedDict'][variable] = {}

        #now start to populate it
        #cfg['perturbedDict'][variable]['nominal'] = cfg[variable]
        #cfg['perturbedDict'][variable]['pertubation'] = cfg[variable+'_perturbation']
        #cfg['perturbedDict'][variable]['deltas'] = cfg[variable+'_deltas']

        #now need to run this through the set_perturbed_values function

        from perturb_utils import set_perturbed_values

        # The current perturbed dictionary value is then overwritten by the 
        # perturbed variable list returned by set_perturbed_values function.
        
        cfg['perturbedDict'][variable] = set_perturbed_values(cfg['NominalValues'][variable], cfg['PerturbationMagnitudes'][variable], cfg['TypeOfPerturbation'][variable], cfg['levels'])

        #cfg['perturbedDict'][variable] = set_perturbed_values(variable,\
        #    [cfg['perturbedDict'][variable]['nominal']] + cfg['perturbedDict'][variable]['deltas'],
        #    cfg['defaultPerturbations'], cfg['perturbedDict'][variable]['pertubation'],
        #    cfg['levels'])

        print "Done with perturbed variable {0}.".format(variable)

    print "Done checking specific nenzfr perturbed inputs."

    return cfg


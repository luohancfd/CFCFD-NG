def sensitivity_input_checker(cfg):
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


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
            
    elif cfg['facility'] == 'expansion-tube':
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
        
    if 'contourFilename' not in cfg:
        cfg['contourFilename'] = 'Bezier-control-pts-t4-m10.data'
        print "No nozzle contour file selected. Choosing '{0}'."\
        .format(cfg['contourFilename'])
        
    if 'exitSliceFileName' not in cfg:
        cfg['exitSliceFileName'] = 'nozzle-exit.data'
        print "No filename specified for file that holds nozzle-exit data."
        print "Setting it to default value of '{0}'.".format(cfg['exitSliceFileName'])
        
    if 'justStats' not in cfg:
        cfg['justStats'] = False
        print "Switch to skip detailed calculations not specified."
        print "Setting it to default value of justStats = {0}".format(cfg['justStats'])
    
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
            
    return cfg
        
#! /usr/bin/env python
"""
pitot_condition_builder.py: pitot condition builder

This file will take a normal pitot input file with a few
extra variables and do a series of pitot runs checking how different 
facility parameters will affect the condition.

Chris James (c.james4@uq.edu.au) - 29/12/13 

"""

VERSION_STRING = "1-Oct-2014"


import sys

from pitot import run_pitot
from pitot_input_utils import *

def check_new_inputs(cfg):
    """Takes the input file and checks that the extra inputs required for the
       condition builder are working..
    
    Returns the checked over input file and will tell the bigger program to 
    bail out if it finds an issue.
    
    """
    
    print "Starting check of condition builder specific inputs."
    
    if 'driver_condition_list' not in cfg:
        print "'driver_condition_list' variable not specified. Bailing out."
        print "Variable needs to be a list containing valid driver conditions."
        cfg['bad_input'] = True
    driver_condition_check_list = ['He:0.80,Ar:0.20', 'He:0.90,Ar:0.10', 'He:1.0']
    for driver_condition in cfg['driver_condition_list']:
        if driver_condition not in driver_condition_check_list:
            print "Invalid driver condition found in 'driver_condition_list'. Bailing out."
            cfg['bad_input'] = True
            
    if 'secondary_list' not in cfg:
        print "'secondary_list' variable not specified. Bailing out."
        print "Variable needs to be a list containing boolean statements."
        cfg['bad_input'] = True
    for secondary_value in cfg['secondary_list']:
        if not isinstance(secondary_value, bool):
            print "Value in 'secondary_list' is not a boolean statement. Bailing out."
            cfg['bad_input'] = True
    if len(cfg['secondary_list']) == 2 and not True in cfg['secondary_list'] or \
        len(cfg['secondary_list']) == 2 and not False in cfg['secondary_list']:
        print "'secondary_list' has two values but it does not contain True and False"
        print "Rethink your input. Bailing out."
        cfg['bad_input'] = True
    if len(cfg['secondary_list']) > 2:
        print "'secondary_list' has more than two values. This is not possible."
        print "Rethink your input. Bailing out."
        cfg['bad_input'] = True
            
    if True in cfg['secondary_list'] and 'psd1_list' not in cfg:
        print "You have chosen to use the secondary driver and not specified a 'psd1_list'. Bailing out."
        cfg['bad_input'] = True
        
    if 'p1_list' not in cfg:
        print "You have not specified a 'p1_list'. Bailing out."
        cfg['bad_input'] = True
        
    if cfg['tunnel_mode'] == 'expansion-tube' and 'p5_list' not in cfg:
        print "You have not specified a 'p5_list'. Bailing out."
        cfg['bad_input'] = True
        
    if 'store_electron_concentration' not in cfg:
        cfg['store_electron_concentration'] = False
        
    if cfg['bad_input']: #bail out here if you end up having issues with your input
        print "Config failed check. Bailing out now."
        print '-'*60
        sys.exit(1)
        
    if not cfg['bad_input']:
        print "Extra input check completed. Test will now run."
        
    return cfg
    
def calculate_number_of_test_runs(cfg):
    """Function that uses a simple function to calculate the amount of test
       runs that the program has to perform.
    """
    if cfg['tunnel_mode'] == 'expansion-tube':
        if True in cfg['secondary_list']:
            total_with_sd = len(cfg['driver_condition_list'])*\
            len(cfg['psd1_list'])*len(cfg['p1_list'])*len(cfg['p5_list'])
            total = total_with_sd
        
        if False in cfg['secondary_list']:
            total_without_sd = len(cfg['driver_condition_list'])*len(cfg['p1_list'])\
            *len(cfg['p5_list'])
            total = total_without_sd
            
        if len(cfg['secondary_list']) == 2:
            total = total_with_sd + total_without_sd
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel':
        if True in cfg['secondary_list']:
            total_with_sd = len(cfg['driver_condition_list'])*\
            len(cfg['psd1_list'])*len(cfg['p1_list'])
            total = total_with_sd
        
        if False in cfg['secondary_list']:
            total_without_sd = len(cfg['driver_condition_list'])*len(cfg['p1_list'])
            total = total_without_sd
            
        if len(cfg['secondary_list']) == 2:
            total = total_with_sd + total_without_sd
    
    return total
    
def condition_builder_test_run(cfg, condition_builder_output, results):
    """Function that takes the fully built config dictionary
       and the text file that is being used for the program output
       and does the test run then adds a line to the output file.
    """
    
    condition_status = True #This will be turned to False if the condition fails
    
    cfg['filename'] = cfg['original_filename'] + '-test-{0}'.format(cfg['test_number'])
    
    print '-'*60
    print "Running test {0} of {1}.".format(cfg['test_number'], cfg['number_of_test_runs'])
    try:
        cfg, states, V, M = run_pitot(cfg = cfg)
    except Exception:          
        print "Test {0} failed. Result will not be printed to csv output.".format(cfg['test_number'])
        condition_status = False
    if cfg['secondary'] and cfg['Vsd'] > cfg['Vs1']:
        print "Vsd is faster than Vs1, condition cannot be simulated by Pitot properly."
        print "Test {0} is considered failed, and result will not be printed to csv output.".format(cfg['test_number'])
        condition_status = False
    if condition_status:
        results = add_new_result_to_results_dict(cfg, states, V, M, results)
        string_output = output_builder (cfg, states, V, M, results)
        condition_builder_output.write(string_output + '\n')
        # need to remove Vs values from the dictionary or it will bail out
        # on the next run
        cfg.pop('Vsd'); cfg.pop('Vs1'); cfg.pop('Vs2')
    else:
        # need to remove Vs values from the dictionary or it will bail out
        # on the next run         
        if cfg['secondary']: cfg.pop('Vsd') 
        cfg.pop('Vs1'); cfg.pop('Vs2') 
            
    return condition_status, results
    
def output_builder(cfg, states, V, M, results):
    """Function that takes the four dictionaries from the completed
       Pitot run and builds a string output that can be added to the
       csv output file.
    """
    
    if not cfg['secondary']: #need to fill psd and Vsd with string values
        cfg['psd1'] = 'Not used'
        cfg['Vsd'] = 'N/A'
    
    # needed to change these as the extra comma was screwing up the csv
    if cfg['driver_gas'] == 'He:1.0':
        driver_gas = cfg['driver_gas']
    elif cfg['driver_gas'] == 'He:0.90,Ar:0.10':
        driver_gas = 'He:0.9 Ar:0.1'
    elif cfg['driver_gas'] == 'He:0.80,Ar:0.20':
        driver_gas = 'He:0.8 Ar:0.2'
    
    # Now make the basic string
    
    # something to add a line if the stagnation enthalpy did not calculate
    # properly in the last run
    # (it will be set to None)
    
    if cfg['stagnation_enthalpy'] == None:
        cfg['stagnation_enthalpy'] = 0.0
    
    basic = "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18}"\
            .format(cfg['test_number'], driver_gas, cfg['psd1'], 
                    cfg['p1'], cfg['p5'], cfg['Vsd'], cfg['Vs1'], cfg['Vs2'],
                    states['s2'].p, states['s2'].T, states['s2'].rho, 
                    V['s2'], M['s2'],
                    states['s7'].p, states['s7'].T, states['s7'].rho, 
                    V['s7'], M['s7'],cfg['stagnation_enthalpy']/10**6)
                    
    # then make other strings that are needed and add what is required.               
    
    if cfg['nozzle']:
        nozzle = ",{0},{1},{2},{3},{4},{5}"\
                 .format(cfg['area_ratio'], states['s8'].p, states['s8'].T, 
                         states['s8'].rho, V['s8'], M['s8'])        
    if cfg['conehead']:
        conehead = ",{0},{1},{2},{3}".format(states['s10c'].p, states['s10c'].T, 
                                             states['s10c'].rho, V['s10c'])
    if cfg['shock_over_model']:
        shock_over_model = ",{0},{1},{2},{3},{4},{5},{6},{7}"\
                           .format(states['s10f'].p, states['s10f'].T, 
                                   states['s10f'].rho, V['s10f'],
                                   states['s10e'].p, states['s10e'].T, 
                                   states['s10e'].rho, V['s10e'])
    
    if cfg['store_electron_concentration']:
        # I decided to start just pulling these from the results dict lists 
        # instead of refinding them, that seems lazy 
        store_electron_concentration = ",{0},{1},{2},{3}"\
                                       .format(results['s2ec'][-1],
                                               results['s7ec'][-1],
                                               results['s8ec'][-1],
                                               results['s10ec'][-1])        
                                   
    # now put it all together
    
    if cfg['nozzle']:
        basic = basic + nozzle
    
    if cfg['conehead'] and not cfg['shock_over_model'] and not cfg['store_electron_concentration']:
        string_output = basic + conehead
    elif cfg['shock_over_model'] and not cfg['conehead'] and not cfg['store_electron_concentration']:
        string_output = basic + shock_over_model
    elif cfg['shock_over_model'] and cfg['conehead'] and not cfg['store_electron_concentration']:
        string_output = basic + conehead + shock_over_model   
    elif cfg['conehead'] and not cfg['shock_over_model'] and cfg['store_electron_concentration']:
        string_output = basic + conehead + store_electron_concentration
    elif cfg['shock_over_model'] and not cfg['conehead'] and cfg['store_electron_concentration']:
        string_output = basic + shock_over_model + store_electron_concentration
    elif cfg['shock_over_model'] and cfg['conehead'] and cfg['store_electron_concentration']:
        string_output = basic + conehead + shock_over_model + store_electron_concentration
    else:
        string_output = basic
    
    return string_output
    
def add_new_result_to_results_dict(cfg, states, V, M, results):
    """Function that takes a completed test run and adds the tunnel
       configuration and results to the results dictionary.
    """ 
    
    if not cfg['secondary']: #need to fill psd and Vsd with string values
        cfg['psd1'] = 'Not used'
        cfg['Vsd'] = 'N/A'
    if cfg['tunnel_mode'] == 'nr-shock-tunnel':
        cfg['p5'] = 'Not used'
        cfg['Vs2'] = 'N/A'
    
    results['test number'].append(cfg['test_number'])
    results['driver condition'].append(cfg['driver_gas'])
    results['psd1'].append(cfg['psd1'])
    results['p1'].append(cfg['p1']) 
    results['p5'].append(cfg['p5'])
    results['Vsd'].append(cfg['Vsd'])
    results['Vs1'].append(cfg['Vs1'])
    results['Vs2'].append(cfg['Vs2'])
    results['p2'].append(states['s2'].p)
    results['T2'].append(states['s2'].T)
    results['rho2'].append(states['s2'].rho)
    results['V2'].append(V['s2'])
    results['M2'].append(M['s2'])
    results['p7'].append(states['s7'].p)
    results['T7'].append(states['s7'].T)
    results['rho7'].append(states['s7'].rho)
    results['V7'].append(V['s7'])
    results['M7'].append(M['s7'])
    if cfg['stagnation_enthalpy']:
        results['Ht'].append(cfg['stagnation_enthalpy']/10**6)
    else:
        results['Ht'].append('did not solve')
    
    if cfg['nozzle']:
        results['arearatio'].append(cfg['area_ratio'])
        results['p8'].append(states['s8'].p)
        results['T8'].append(states['s8'].T)
        results['rho8'].append(states['s8'].rho)
        results['V8'].append(V['s8'])
        results['M8'].append(M['s8'])

    if cfg['conehead']:
        results['p10c'].append(states['s10c'].p)
        results['T10c'].append(states['s10c'].T)
        results['rho10c'].append(states['s10c'].rho)
        results['V10c'].append(V['s10c'])
        
    if cfg['shock_over_model']:
        results['p10f'].append(states['s10f'].p)
        results['T10f'].append(states['s10f'].T)
        results['rho10f'].append(states['s10f'].rho)
        results['V10f'].append(V['s10f'])
        results['p10e'].append(states['s10e'].p)
        results['T10e'].append(states['s10e'].T)
        results['rho10e'].append(states['s10e'].rho)
        results['V10e'].append(V['s10e'])
        
    if cfg['store_electron_concentration']:
        if 'e-' in states['s2'].species.keys():
            results['s2ec'].append(states['s2'].species['e-'])
        else:
            results['s2ec'].append(0.0)
            
        if cfg['tunnel_mode'] == 'expansion-tube':
            if 'e-' in states['s7'].species.keys():
                results['s7ec'].append(states['s7'].species['e-'])
            else:
                results['s7ec'].append(0.0)   
        else:
            results['s7ec'].append('N/A')
            
        if cfg['nozzle']:
            if 'e-' in states['s8'].species.keys():
                results['s8ec'].append(states['s8'].species['e-'])
            else:
                results['s8ec'].append(0.0) 
        if cfg['shock_over_model']:
            if 's10e' in states.keys():
                if 'e-' in states['s10e'].species.keys():
                    results['s10ec'].append(states['s10e'].species['e-'])
                else:
                    results['s10ec'].append(0.0)
            else:
                results['s10ec'].append('did not solve')
                     
    return results
    
def condition_builder_summary_builder(cfg, results, condition_builder_summary_file):
    """Function that takes the config dictionary and results dictionary 
       made throughout the running of the program and prints a summary of 
       the run to the screen and to a summary text file
    """
    
    print '-'*60
    print "Printing summary to screen and to a text document."
    
    summary_line_3 = "{0} tests ran. {1} ({2:.1f}%) were successful."\
        .format(cfg['number_of_test_runs'], len(results['test number']),
        float(len(results['test number']))/float(cfg['number_of_test_runs'])*100.0)
    print summary_line_3
    condition_builder_summary_file.write(summary_line_3 + '\n')        
    
    if len(cfg['driver_condition_list']) > 1:
        summary_line_4 = "{0} driver conditions were tested ({1})."\
        .format(len(cfg['driver_condition_list']), cfg['driver_condition_list'])
    else:
        summary_line_4 = "Only the {0} driver condition was tested."\
        .format(cfg['driver_condition_list'][0])
    print summary_line_4
    condition_builder_summary_file.write(summary_line_4 + '\n')

    if len(cfg['secondary_list']) == 2: 
        summary_line_5 = "Calculations performed with and without a secondary driver."
    elif cfg['secondary_list'][0]:
        summary_line_5 = "Calculations only performed with a secondary driver."
    elif not cfg['secondary_list'][0]:
        summary_line_5 = "Calculations only performed WITHOUT a secondary driver."
    print summary_line_5
    condition_builder_summary_file.write(summary_line_5 + '\n')

    for variable in results['full_list']:
        # first check it's not a variable that doesn't need to be summarised
        if variable not in ['test number', 'driver condition', 'arearatio']:
            # now check the list and remove any string values
            data_list = []
            for value in results[variable]:
                if isinstance(value, (float,int)): data_list.append(value)
            if len(data_list) > 0: #ie. don't bother summarising if there is no numerical data there
                max_value = max(data_list)
                min_value = min(data_list)
                if variable[0] == 'p':
                    summary_line = "Variable {0} varies from {1:.2f} - {2:.2f} Pa."\
                    .format(variable, min_value, max_value)
                elif variable[0] == 'V':
                    summary_line = "Variable {0} varies from {1:.2f} - {2:.2f} m/s."\
                    .format(variable, min_value, max_value)
                elif variable[0] == 'T':
                    summary_line = "Variable {0} varies from {1:.2f} - {2:.2f} K."\
                    .format(variable, min_value, max_value)
                elif variable[0] == 'r':
                    summary_line = "Variable {0} varies from {1:.7f} - {2:.7f} kg/m**3."\
                    .format(variable, min_value, max_value)
                elif variable[0] == 'H':
                    summary_line = "Variable {0} varies from {1:.7f} - {2:.7f} MJ/kg."\
                    .format(variable, min_value, max_value)
                elif variable[-2:] == 'ec':
                    summary_line = "Variable {0} varies from {1:.7f} - {2:.7f} (mole fraction)."\
                    .format(variable, min_value, max_value)
                else:
                    summary_line = "Variable {0} varies from {1:.1f} - {2:.1f}."\
                    .format(variable, min_value, max_value)
                print summary_line
                condition_builder_summary_file.write(summary_line + '\n')
                
    return
            
def run_pitot_condition_builder(cfg = {}, config_file = None):
    """
    
    Chris James (c.james4@uq.edu.au) 27/12/13
    
    run_pitot_condition_builder(dict) - > depends
    
    """
    
    import time
    
    #---------------------- get the inputs set up --------------------------
    
    if config_file:
        cfg = config_loader(config_file)
        
    #----------------- check inputs ----------------------------------------
    
    cfg = input_checker(cfg)
    
    cfg = check_new_inputs(cfg)
    
    # make a counter so we can work out what test we're running
    # also make one to store how many runs are successful
    
    cfg['number_of_test_runs'] = calculate_number_of_test_runs(cfg)
    
    import copy
    cfg['original_filename'] = copy.copy(cfg['filename'])
    
    counter = 0
    good_counter = 0
    
    # print how many tests we're going to run, and the ranges.
    
    print '-'*60    
    print "{0} tests will be run.".format(cfg['number_of_test_runs'])
    
    if True in cfg['secondary_list']:
        if len(cfg['psd1_list']) > 1:
            print 'psd1 will be changed from from {0} - {1} Pa in increments of {2} Pa.'\
            .format(cfg['psd1_list'][0], cfg['psd1_list'][-1], cfg['psd1_list'][1] - cfg['psd1_list'][0])
        else:
            print 'psd1 is kept at {0} Pa.'.format(cfg['psd1_list'][0])
            
    if len(cfg['p1_list']) > 1:
        print 'p1 will be changed from from {0} - {1} Pa in increments of {2} Pa.'\
        .format(cfg['p1_list'][0], cfg['p1_list'][-1], cfg['p1_list'][1] - cfg['p1_list'][0])
    else:
        print 'p1 is kept at {0} Pa.'.format(cfg['p1_list'][0])
        
    if cfg['tunnel_mode'] == 'expansion-tube':
        if len(cfg['p5_list']) > 1:
            print 'p5 will be changed from from {0} - {1} Pa in increments of {2} Pa.'\
            .format(cfg['p5_list'][0], cfg['p5_list'][-1], cfg['p5_list'][1] - cfg['p5_list'][0])
        else:
            print 'p5 is kept at {0} Pa.'.format(cfg['p5_list'][0])
        
    #open our csv file ready to go
    
    # open a file to start saving results
    condition_builder_output = open(cfg['filename']+'-condition-builder.csv',"w")  #csv_output file creation
    # print a line explaining the results
    intro_line_1 = "# Output of pitot area condition building program."
    condition_builder_output.write(intro_line_1 + '\n')
    basic = "#test number,driver condition,psd1,p1,p5,Vsd,Vs1,Vs2,p2,T2,rho2,V2,M2,p7,T7,rho7,V7,M7,Ht"
    nozzle = ",arearatio,p8,T8,rho8,V8,M8"
    if cfg['nozzle']:
        basic = basic + nozzle
    conehead = ",p10c,T10c,rho10c,V10c"
    shock_over_model = ",p10f,T10f,rho10f,V10f,p10e,T10e,rho10f,V10e"
    store_electron_concentration = ",s2ec, s7ec, s8ec, s10ec"
    if cfg['conehead'] and not cfg['shock_over_model'] and not cfg['store_electron_concentration']:
        intro_line_2 = basic + conehead
    elif cfg['shock_over_model'] and not cfg['conehead'] and not cfg['store_electron_concentration']:
        intro_line_2 = basic + shock_over_model
    elif cfg['shock_over_model'] and cfg['conehead'] and not cfg['store_electron_concentration']:
        intro_line_2 = basic + conehead + shock_over_model   
    if cfg['conehead'] and not cfg['shock_over_model'] and cfg['store_electron_concentration']:
        intro_line_2 = basic + conehead + store_electron_concentration
    elif cfg['shock_over_model'] and not cfg['conehead'] and cfg['store_electron_concentration']:
        intro_line_2 = basic + shock_over_model + store_electron_concentration
    elif cfg['shock_over_model'] and cfg['conehead'] and cfg['store_electron_concentration']:
        intro_line_2 = basic + conehead + shock_over_model + store_electron_concentration
    elif cfg['store_electron_concentration'] and not cfg['shock_over_model'] and not cfg['conehead']:
        intro_line_2 = basic + store_electron_concentration
    else:
        intro_line_2 = basic
    condition_builder_output.write(intro_line_2 + '\n')
    
    # then make a dictionary of lists to store results in the Python memory
        
    # need to make a list to create a series of empty lists in the results
    # dictionary to store the data. the list is tailored to the test condition
    basic_list = ['test number','driver condition','psd1','p1','p5','Vsd','Vs1',
                  'Vs2','p2','T2','rho2','V2','M2','p7','T7','rho7','V7','M7','Ht']
    nozzle_list = ['arearatio','p8','T8','rho8','V8','M8']
    if cfg['nozzle']:
        basic_list = basic_list + nozzle_list
    conehead_list = ['p10c','T10c','rho10c','V10c']
    shock_over_model_list = ['p10f','T10f','rho10f','V10f','p10e','T10e','rho10e','V10e']
    store_electron_concentration_list = ['s2ec','s7ec','s8ec','s10ec']
    if cfg['conehead'] and not cfg['shock_over_model'] and not cfg['store_electron_concentration']:
        full_list = basic_list + conehead_list
    elif cfg['shock_over_model'] and not cfg['conehead'] and not cfg['store_electron_concentration']:
        full_list = basic_list + shock_over_model_list
    elif cfg['shock_over_model'] and cfg['conehead'] and not cfg['store_electron_concentration']:
        full_list = basic_list + conehead_list + shock_over_model_list
    elif cfg['conehead'] and not cfg['shock_over_model'] and cfg['store_electron_concentration']:
        full_list = basic_list + conehead_list + store_electron_concentration_list
    elif cfg['shock_over_model'] and not cfg['conehead'] and cfg['store_electron_concentration']:
        full_list = basic_list + shock_over_model_list + store_electron_concentration_list
    elif cfg['shock_over_model'] and cfg['conehead'] and cfg['store_electron_concentration']:
        full_list = basic_list + conehead_list + shock_over_model_list + store_electron_concentration_list 
    elif cfg['store_electron_concentration'] and not cfg['conehead'] and not cfg['shock_over_model']:
        full_list = basic_list + store_electron_concentration_list 
    else:
        full_list = basic_list
    
    results = {title : [] for title in full_list}
    
    # add the list of titles in case we want to use it in future
    
    results['full_list'] = full_list
    
    have_checked_time = False
    
    #now start up the for loops and get running

    if cfg['tunnel_mode'] == 'expansion-tube':   
        for driver_condition in cfg['driver_condition_list']:
            cfg['driver_gas'] = driver_condition
            for secondary_value in cfg['secondary_list']:
                cfg['secondary'] = secondary_value
                if cfg['secondary']:
                    for psd1 in cfg['psd1_list']:
                        cfg['psd1'] = psd1
                        for p1 in cfg['p1_list']:
                            cfg['p1'] = p1
                            for p5 in cfg['p5_list']:
                                cfg['p5'] = p5
                                counter += 1
                                cfg['test_number'] = counter
                                if not have_checked_time:
                                    start_time = time.time()
                                run_status, results = condition_builder_test_run(cfg, condition_builder_output, results) 
                                if run_status:
                                    good_counter += 1
                                    if not have_checked_time:
                                        test_time = time.time() - start_time
                                        print '-'*60
                                        print "Time to complete first test was {0:.2f} seconds."\
                                        .format(test_time)
                                        print "If every test takes this long. It will take roughly {0:.2f} hours to perform all {1} tests."\
                                        .format(test_time*cfg['number_of_test_runs']/3600.0, cfg['number_of_test_runs'])
                                        have_checked_time = True
                else:
                    for p1 in cfg['p1_list']:
                        cfg['p1'] = p1
                        for p5 in cfg['p5_list']:
                            cfg['p5'] = p5
                            counter += 1
                            cfg['test_number'] = counter
                            if not have_checked_time:
                                start_time = time.time()
                            run_status, results = condition_builder_test_run(cfg, condition_builder_output, results) 
                            if run_status:
                                good_counter += 1
                                if not have_checked_time:
                                    test_time = time.time() - start_time
                                    print '-'*60
                                    print "Time to complete first test was {0:.2f} seconds."\
                                    .format(test_time)
                                    print "If every test takes this long. It will take roughly {0:.2f} hours to perform all {1} tests."\
                                    .format(test_time*cfg['number_of_test_runs']/3600.0, cfg['number_of_test_runs'])
                                    have_checked_time = True

    elif cfg['tunnel_mode'] == 'nr-shock-tunnel':   
        for driver_condition in cfg['driver_condition_list']:
            cfg['driver_gas'] = driver_condition
            for secondary_value in cfg['secondary_list']:
                cfg['secondary'] = secondary_value
                if cfg['secondary']:
                    for psd1 in cfg['psd1_list']:
                        cfg['psd1'] = psd1
                        for p1 in cfg['p1_list']:
                            cfg['p1'] = p1
                            counter += 1
                            cfg['test_number'] = counter
                            if not have_checked_time:
                                start_time = time.time()
                            run_status, results = condition_builder_test_run(cfg, condition_builder_output, results) 
                            if run_status:
                                good_counter += 1
                                if not have_checked_time:
                                    test_time = time.time() - start_time
                                    print '-'*60
                                    print "Time to complete first test was {0:.2f} seconds."\
                                    .format(test_time)
                                    print "If every test takes this long. It will take roughly {0:.2f} hours to perform all {1} tests."\
                                    .format(test_time*cfg['number_of_test_runs']/3600.0, cfg['number_of_test_runs'])
                                    have_checked_time = True
                else:
                    for p1 in cfg['p1_list']:
                        cfg['p1'] = p1
                        counter += 1
                        cfg['test_number'] = counter
                        if not have_checked_time:
                            start_time = time.time()
                        run_status, results = condition_builder_test_run(cfg, condition_builder_output, results) 
                        if run_status:
                            good_counter += 1
                            if not have_checked_time:
                                test_time = time.time() - start_time
                                print '-'*60
                                print "Time to complete first test was {0:.2f} seconds."\
                                .format(test_time)
                                print "If every test takes this long. It will take roughly {0:.2f} hours to perform all {1} tests."\
                                .format(test_time*cfg['number_of_test_runs']/3600.0, cfg['number_of_test_runs'])
                                have_checked_time = True       
                        
    condition_builder_output.close()
    
    # now analyse results dictionary and print some results to the screen
    # and another external file
    
    condition_builder_summary_file = open(cfg['original_filename']+'-condition-builder-summary.txt',"w")
    # print a line explaining the results
    summary_line_1 = "# Summary of pitot condition building program output."
    condition_builder_summary_file.write(summary_line_1 + '\n')
    summary_line_2 = "# Summary performed using Version {0} of the condition building program.".format(VERSION_STRING)
    condition_builder_summary_file.write(summary_line_2 + '\n')
    
    condition_builder_summary_builder(cfg, results, condition_builder_summary_file)
    
    condition_builder_summary_file.close()    
    
    return
                                
#----------------------------------------------------------------------------

def main():
    
    import optparse  
    op = optparse.OptionParser(version=VERSION_STRING)   
    op.add_option('-c', '--config_file', dest='config_file',
                  help=("filename where the configuration file is located"))    

    opt, args = op.parse_args()
    config_file = opt.config_file
           
    run_pitot_condition_builder(cfg = {}, config_file = config_file)
    
    return
    
#----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "pitot_condition_builder.py - Pitot Equilibrium expansion tube simulator condition building tool"
        print "start with --help for help with inputs"
        
    else:
        main()

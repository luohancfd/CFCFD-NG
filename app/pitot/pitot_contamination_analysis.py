#! /usr/bin/env python
"""
pitot_contamination_analysis.py: pitot air contamination analysing tool

This file with take a normal pitot input file with two variables added, one 
list or array telling it what range of air contamination to analyse,
and another telling it whether those values are by volume (mole fraction)
or by mass (mass fraction). If the second variable is not supplied, it will
default to volume as pitot generally uses everything by volume.

This is basically a cut down version of the pitot condition building program
(pitot_condition_builder.py) and it prints data to the screen and a text file
in a similar way. 

Chris James (c.james4@uq.edu.au) - 12/09/14

"""

VERSION_STRING = "15-Dec-2016"

from pitot_condition_builder import stream_tee, cleanup_old_files, pickle_result_data, zip_result_and_log_files, results_csv_builder, normalised_results_csv_builder, pickle_intermediate_data

import sys, os

from pitot import run_pitot
from pitot_input_utils import *

def check_new_inputs(cfg):
    """Takes the input file and checks that the extra inputs required for the
       condition builder are working..
    
    Returns the checked over input file and will tell the bigger program to 
    bail out if it finds an issue.
    
    """
    
    print "Starting check of condition builder specific inputs."
        
    if 'air_contamination_list' not in cfg:
        print "You have not specified an 'air_contamination_list'. Bailing out."
        cfg['bad_input'] = True
        
    if 'air_contamination_inputUnits' not in cfg:
        print "'air_contamination_inputUnits' not in cfg."
        print "Setting it to 'moles'."
        cfg['air_contamination_inputUnits']
        
    if 'normalise_results_by' not in cfg:
        print "'normalise_results_by' not in cfg."
        print "Setting it to 'first value'."
        cfg['normalise_results_by'] = 'first value'   
        
    # add an extra check to make sure that the inputUnits match between the contamination
    # and the test gas if the test gas is custom. I think it's pretty stupid
    # and confusing if they aren't the same
    
    if cfg['test_gas'] == 'custom' and cfg['air_contamination_inputUnits'] != cfg['test_gas_inputUnits']:
        print "'air_contamination_inputUnits' = {0} and 'test_gas_inputUnits' = {1}."\
        .format(cfg['air_contamination_inputUnits'], cfg['test_gas_inputUnits'])
        print "This isn't the smartest idea. Bailing out."
        cfg['bad_input'] = True
        
    if 'store_electron_concentration' not in cfg:
        cfg['store_electron_concentration'] = False
        
    if 'cleanup_old_files' not in cfg:
        print "'cleanup_old_files' variable not set. Setting to default value of 'False'"
        cfg['cleanup_old_files'] = False
                
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
    return len(cfg['air_contamination_list'])
    
def start_message(cfg):
    """Function that takes the cfg file and prints a start message for the
       program, detailing what it will do.
    """
    
    # print how many tests we're going to run, and the ranges.
    
    print '-'*60    
    print 'Running Pitot Air Contamination Analysis Version: {0}.'.format(VERSION_STRING)    
    print "{0} tests will be run.".format(cfg['number_of_test_runs'])
    
    if cfg['air_contamination_inputUnits'] == 'moles':
        print "Air contamination will be tested from {0} - {1} % in increments of {2} as a mole fraction."\
        .format(cfg['air_contamination_list'][0], cfg['air_contamination_list'][-1], 
                cfg['air_contamination_list'][1] - cfg['air_contamination_list'][0])
    elif cfg['air_contamination_inputUnits'] == 'massf':
        print "Air contamination will be tested from {0} - {1} % in increments of {2} as a massf fraction."\
        .format(cfg['air_contamination_list'][0], cfg['air_contamination_list'][-1], 
                cfg['air_contamination_list'][1] - cfg['air_contamination_list'][0])
            
    return cfg
    
def build_results_dict(cfg):
    """Function that looks at the cfg dictionary and works out what data needs
       to be stored for the type of test that we're running. Then it populates
       a dictionary called 'results' with empty lists for the data to be stored in.
       The dictionary is then returned.
    """
    
    # need to make a list to create a series of empty lists in the results
    # dictionary to store the data. the list is tailored to the test we're running
    
    full_list = []
    if cfg['secondary']:
        basic_list = ['test number','air contamination','psd1','p1','p5','Vsd',
                      'Vs1', 'Vs2', 'Ht','h','u_eq', 'rho1', 'gamma1', 'R1', 'MW1',
                      'p2','T2','rho2','V2','M2', 'a2', 'gamma2', 'R2', 'Vs1 - V2', 'Ht2',
                      'p6','T6','rho6','V6','M6','p7','T7','rho7','V7','M7']
    else:
       basic_list = ['test number','air contamination','p1','p5','Vs1', 'Vs2', 
                     'Ht','h','u_eq', 'rho1', 'gamma1', 'R1', 'MW1',
                     'p2','T2','rho2','V2','M2', 'a2', 'gamma2', 'R2', 'Vs1 - V2', 'Ht2',
                     'p6','T6','rho6','V6','M6','p7','T7','rho7','V7','M7']        
    full_list += basic_list
    
    if cfg['nozzle']:
        nozzle_list = ['arearatio','p8','T8','rho8','V8','M8']
        full_list += nozzle_list
    if cfg['conehead']:
         conehead_list = ['p10c','T10c','rho10c','V10c']
         full_list += conehead_list
    if cfg['shock_over_model']:
        shock_over_model_list = ['p10f','T10f','rho10f','V10f','p10e','T10e','rho10e','V10e']
        full_list += shock_over_model_list
    if cfg['store_electron_concentration']:     
        store_electron_concentration_list = ['s2ec','s7ec','s8ec','s10ec']
        full_list += store_electron_concentration_list

    # now populate the dictionary with a bunch of empty lists based on that list

    results = {title : [] for title in full_list}
    
    # add the list of titles in case we want to use it in future
    
    results['full_list'] = full_list
    
    #add a list where we can store unsuccesful run numbers for analysis
    results['unsuccessful_runs'] = []
    
    print '-'*60
    print "The full list of variables to be added to the output are:"
    print full_list
    
    return results
    
def build_test_condition_input_details_dictionary(cfg):
    """Function that builds a dictionary that tells the condition
       builder what tests to run. It will include a list called 'test_names'
       that will be cycled through by the main condition program.
       
    """
    
    test_condition_input_details = {}
    test_condition_input_details['test_names'] = []
    
    test_name = 0 # we will add to this as we go, first test_name will be 1
    
    print '-'*60
    print "Building the test condition details dictionary containing a dictionary for each simulation." 
    
    for contamination_percentage in cfg['air_contamination_list']:
        # change the test name
        test_name += 1
        # store that test name
        test_condition_input_details['test_names'].append(test_name)
        # now make a dictionary in the test_condition_input_details dictionary for this simulation
        # using the diluent percantage
        # this is simple here as diluent_percentage is the only variable...
        # (for the normal condition builder it was a bit more complicated...)
        input_dictionary = {'contamination_percentage': contamination_percentage}
        test_condition_input_details[test_name] = input_dictionary
                
    print "The test_names list for this simulation is:"
    print test_condition_input_details['test_names']
                
    return test_condition_input_details    
    
def contamination_analysis_test_run(cfg, results):
    """Function that takes the fully built config dictionary
       and the text file that is being used for the program output
       and does the test run then adds a line to the output file.
    """
    
    condition_status = True #This will be turned to False if the condition fails
    
    cfg['filename'] = cfg['original_filename'] + '-test-{0}-result'.format(cfg['test_number'])
    
    # some code here to make a copy of the stdout printouts for each test and store it
    
    import sys
    
    test_log = open(cfg['original_filename'] + '-test-{0}-log.txt'.format(cfg['test_number']),"w")
    sys.stdout = stream_tee(sys.stdout, test_log)
    
    print '-'*60
    print "Running test {0} of {1}.".format(cfg['test_number'], cfg['number_of_test_runs'])
    print "Current level of air contamination is {0} % (by {1})."\
         .format(cfg['contamination_percentage'], cfg['air_contamination_inputUnits'])
    try:
        cfg, states, V, M = run_pitot(cfg = cfg)
    except Exception as e:  
        print "Error {0}".format(str(e))
        print "Test {0} failed. Result will not be printed to csv output.".format(cfg['test_number'])
        condition_status = False

    if cfg['secondary'] and cfg['Vsd'] > cfg['Vs1']:
        print "Vsd is faster than Vs1, condition cannot be simulated by Pitot properly."
        print "Test {0} is considered failed, and result will not be printed to csv output.".format(cfg['test_number'])
        condition_status = False
    if condition_status:
        results = add_new_result_to_results_dict(cfg, states, V, M, results)
        # need to remove Vs values from the dictionary or it will bail out
        # on the next run
        if cfg['secondary']: cfg.pop('Vsd') 
        cfg.pop('Vs1'); cfg.pop('Vs2')
    else:
        # need to remove Vs values from the dictionary or it will bail out
        # on the next run         
        if cfg['secondary'] and 'Vsd' in cfg: cfg.pop('Vsd')
        if 'Vs1' in cfg: cfg.pop('Vs1') 
        if 'Vs2' in cfg: cfg.pop('Vs2')
        cfg['last_run_successful'] = False
        results['unsuccessful_runs'].append(cfg['test_number'])
        
    # now we end the stream teeing here by pulling out the original stdout object
    # and overwriting the stream tee with that, then closing the log file
    sys.stdout = sys.stdout.stream1   
    test_log.close()
            
    return condition_status, results
    
def add_new_result_to_results_dict(cfg, states, V, M, results):
    """Function that takes a completed test run and adds the tunnel
       configuration and results to the results dictionary.
    """ 
    
    results['test number'].append(cfg['test_number'])
    results['air contamination'].append(cfg['contamination_percentage'])
    if cfg['secondary']:
        results['psd1'].append(cfg['psd1'])
    results['p1'].append(cfg['p1']) 
    results['p5'].append(cfg['p5'])
    if cfg['secondary']:
        results['Vsd'].append(cfg['Vsd'])
    results['Vs1'].append(cfg['Vs1'])
    results['Vs2'].append(cfg['Vs2'])
    
    if cfg['stagnation_enthalpy']:
        results['Ht'].append(cfg['stagnation_enthalpy']/10**6)
    else:
        results['Ht'].append('did not solve')
    results['h'].append(cfg['freestream_enthalpy']/10**6)
    if cfg['u_eq']:
        results['u_eq'].append(cfg['u_eq'])
    else:
        results['u_eq'].append('did not solve')
    
    results['rho1'].append(states['s1'].rho)
    results['gamma1'].append(states['s1'].gam)
    results['R1'].append(states['s1'].R)
    results['MW1'].append(states['s1'].Mmass)
    
    results['p2'].append(states['s2'].p)
    results['T2'].append(states['s2'].T)
    results['rho2'].append(states['s2'].rho)
    results['V2'].append(V['s2'])
    results['M2'].append(M['s2'])
    results['a2'].append(states['s2'].a)
    results['gamma2'].append(states['s2'].gam)
    results['R2'].append(states['s2'].R)
    results['Vs1 - V2'].append(cfg['Vs1'] - V['s2'])
    
    if cfg['Ht2']:
        results['Ht2'].append(cfg['Ht2']/10**6)
    else:
        results['Ht2'].append('did not solve')
        
    if cfg['tunnel_mode'] == 'expansion-tube':
        results['p6'].append(states['s6'].p)
        results['T6'].append(states['s6'].T)
        results['rho6'].append(states['s6'].rho)
        results['V6'].append(V['s6'])
        results['M6'].append(M['s6'])
        
        results['p7'].append(states['s7'].p)
        results['T7'].append(states['s7'].T)
        results['rho7'].append(states['s7'].rho)
        results['V7'].append(V['s7'])
        results['M7'].append(M['s7'])
  
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
        if 's10f' in states.keys():
            results['p10f'].append(states['s10f'].p)
            results['T10f'].append(states['s10f'].T)
            results['rho10f'].append(states['s10f'].rho)
            results['V10f'].append(V['s10f'])
        else:
            results['p10f'].append('did not solve')
            results['T10f'].append('did not solve')
            results['rho10f'].append('did not solve')
            results['V10f'].append('did not solve')
        if 's10e' in states.keys():
            results['p10e'].append(states['s10e'].p)
            results['T10e'].append(states['s10e'].T)
            results['rho10e'].append(states['s10e'].rho)
            results['V10e'].append(V['s10e'])
        else:
            results['p10e'].append('did not solve')
            results['T10e'].append('did not solve')
            results['rho10e'].append('did not solve')
            results['V10e'].append('did not solve')  
        
    return results
    
def contamination_analysis_summary(cfg, results):
    """Function that takes the config dictionary and results dictionary 
       made throughout the running of the program and prints a summary of 
       the run to the screen and to a summary text file
    """
    
    print '-'*60
    print "Printing summary to screen and to a text document."
    
    contamination_analysis_summary_file = open(cfg['original_filename']+'-contamination-analysis-summary.txt',"w")
    # print a line explaining the results
    summary_line_1 = "# Summary of pitot contamination analysis program output."
    contamination_analysis_summary_file.write(summary_line_1 + '\n')
    print summary_line_1
    summary_line_2 = "# Summary performed using Version {0} of the contamination analysis program.".format(VERSION_STRING)
    contamination_analysis_summary_file.write(summary_line_2 + '\n')
    print summary_line_2
    
    summary_line_3 = "{0} tests ran. {1} ({2:.1f}%) were successful."\
        .format(cfg['number_of_test_runs'], len(results['test number']),
        float(len(results['test number']))/float(cfg['number_of_test_runs'])*100.0)
    print summary_line_3
    contamination_analysis_summary_file.write(summary_line_3 + '\n')  

    if results['unsuccessful_runs']: 
        summary_line_4 = "Unsucessful runs were run numbers {0}.".format(results['unsuccessful_runs'])
        print summary_line_4
        contamination_analysis_summary_file.write(summary_line_4 + '\n')         

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
                elif variable[0] == 'H' or variable[0] == 'h':
                    summary_line = "Variable {0} varies from {1:.7f} - {2:.7f} MJ/kg."\
                    .format(variable, min_value, max_value)                    
                else:
                    summary_line = "Variable {0} varies from {1:.1f} - {2:.1f}."\
                    .format(variable, min_value, max_value)
                print summary_line
                contamination_analysis_summary_file.write(summary_line + '\n')
                
    contamination_analysis_summary_file.close()    
                
    return
            
def run_pitot_contamination_analysis(cfg = {}, config_file = None, force_restart = False):
    """
    
    Chris James (c.james4@uq.edu.au) 12/09/14
    
    run_pitot_contamination_analysis(dict) - > depends
    
    """
    
    import time, os
    
    #---------------------- get the inputs set up --------------------------
    
    if config_file:
        cfg = config_loader(config_file)
        
    #----------------- check inputs ----------------------------------------
    
    cfg = input_checker(cfg)
    
    cfg = check_new_inputs(cfg)
    
    intermediate_filename = cfg['filename']+'-contamination-analysis-intermediate-result-pickle.dat'
    
    # now check if we have attempted an old run before or not....
    if not os.path.isfile(intermediate_filename) or force_restart: 
        # if not, we set up a new one...    
    
        # clean up any old files if the user has asked for it
        
        if cfg['cleanup_old_files']:
            # build a list of auxiliary files and run the cleanup_old_files function frim pitot_condition_builder.py
            auxiliary_file_list = ['-contamination-analysis-log-and-result-files.zip',
                                   '-contamination-analysis-final-result-pickle.dat',
                                   '-contamination-analysis-normalised.csv',
                                   '-contamination-analysis-summary.txt',
                                   '-contamination-analysis.csv']
            cleanup_old_files(auxiliary_file_list)
        
        # make a counter so we can work out what test we're running
        # also make one to store how many runs are successful
        
        cfg['number_of_test_runs'] = calculate_number_of_test_runs(cfg)
        
        import copy
        cfg['original_filename'] = copy.copy(cfg['filename'])
        cfg['original_test_gas'] = copy.copy(cfg['test_gas'])
                
        # print a start message to get us going
        
        cfg = start_message(cfg)
        
        # work out what we need in our results dictionary and make the dictionary
        
        results = build_results_dict(cfg)
        
        # build the dictionary with the details of all the tests we want to run...
        
        test_condition_input_details = build_test_condition_input_details_dictionary(cfg)
        
        # so we can make the first run check the time...                   
        cfg['have_checked_time'] = False
        
        # counter to count the experiments that don't fail
        cfg['good_counter'] = 0
        # list that will be filled by the experiment numbers as they finish...
        cfg['finished_simulations'] = []   
        
    else:
        # we can load an unfinished simulation and finishing it...!
        print '-'*60
        print "It appears that an unfinished simulation was found in this folder"
        print "We are now going to attempt to finish this simulation."
        print "If this is not what you want, please delete the file '{0}' and run the condition builder again.".format(intermediate_filename)
    
        import pickle
        with open(intermediate_filename,"rU") as data_file:
            intermediate_result = pickle.load(data_file)
            cfg = intermediate_result['cfg']
            results = intermediate_result['results']
            test_condition_input_details = intermediate_result['test_condition_input_details']
            data_file.close()
    
    
     #now start the main for loop for the simulation...
    
    for test_name in test_condition_input_details['test_names']:
        # this is for a re-loaded simulation, to make sure we don't re-run what has already been run
        if test_name in cfg['finished_simulations']:
            print '-'*60
            print "test_name '{0}' is already in the 'finished_simulations' list.".format(test_name)
            print "Therefore it has probably already been run and will be skipped..."
            continue # code to skip iteration
        
        # first set the test number variable...
        cfg['test_number'] = test_name
                
        # then set the variables that the simulation will need set
        cfg['contamination_percentage'] = test_condition_input_details[test_name]['contamination_percentage']
        percentage_not_air = 100.0 - cfg['contamination_percentage']
        
        # this is easy if the original test gas was already custom
        if cfg['contamination_percentage'] > 0.0: #obviously don't do anything special when the percentage is 0
            if cfg['original_test_gas'] == 'custom':
                if cfg['test_gas_inputUnits'] == 'moles' and cfg['air_contamination_inputUnits'] == 'moles' or \
                cfg['test_gas_inputUnits'] == 'massf' and cfg['air_contamination_inputUnits'] == 'massf':
                    for species in cfg['test_gas_composition'].keys():
                        cfg['test_gas_composition'][species] = cfg['test_gas_composition'][species]*percentage_not_air
                    # if air was already in the mix, just increase the amount
                    if 'Air' in cfg['test_gas_composition'].keys():
                        cfg['test_gas_composition']['Air'] = cfg['test_gas_composition']['Air'] + cfg['contamination_percentage'] / 100.0
                    # if not, add the air
                    else:
                        cfg['test_gas_composition']['Air'] = cfg['contamination_percentage'] / 100.0              
            else: #noncustom test gas
                from pitot_input_utils import make_test_gas
                if cfg['air_contamination_inputUnits'] == 'moles':
                    original_test_gas, not_needed = make_test_gas(cfg['original_test_gas'], outputUnits='moles')
                elif cfg['air_contamination_inputUnits'] == 'massf':
                    original_test_gas, not_needed = make_test_gas(cfg['original_test_gas'], outputUnits='massf')
                # need to set a gas state here to make it output the species at room temperature and pressure
                # (we could have just used the input values but then we douldn't get the
                # change from moles to massf if we want that)
                original_test_gas.set_pT(101300.0, 300.0)
                # now set our current test gas to custom and populate a reactants
                # dictionary for the custom test gas
                cfg['test_gas'] = 'custom'
                cfg['test_gas_inputUnits'] = cfg['air_contamination_inputUnits']
                cfg['test_gas_composition'] = {}
                # as we set the output units of our custom gas to the units of our
                # air contamination no annoying conversion has to be done here
                for species in original_test_gas.species.keys():
                    cfg['test_gas_composition'][species] = original_test_gas.species[species]*percentage_not_air
                # if air was already in the mix, just increase the amount
                if 'Air' in original_test_gas.species.keys():
                    cfg['test_gas_composition']['Air'] = original_test_gas.species['Air'] + cfg['contamination_percentage'] 
                # if not, add the air
                else:
                    cfg['test_gas_composition']['Air'] = cfg['contamination_percentage']  

        if not cfg['have_checked_time']:
            start_time = time.time()
        run_status, results = contamination_analysis_test_run(cfg, results) 
        if run_status:
            cfg['good_counter'] += 1
            if not cfg['have_checked_time']:
                test_time = time.time() - start_time
                print '-'*60
                print "Time to complete first test was {0:.2f} seconds."\
                .format(test_time)
                print "If every test takes this long. It will take roughly {0:.2f} hours to perform all {1} tests."\
                .format(test_time*cfg['number_of_test_runs']/3600.0, cfg['number_of_test_runs'])
                have_checked_time = True
        
        # add this to finished simulations list, regardless of whether it finished correctly or not...
        cfg['finished_simulations'].append(test_name)
        
        # now pickle the intermediate result so we can result the simulation if needed...
        pickle_intermediate_data(cfg, results, test_condition_input_details, filename = intermediate_filename)

    # now that we're done we can dump the results to the results csv 
    intro_line = "Output of pitot contamination analysis program Version {0}.".format(VERSION_STRING)     
    results_csv_builder(results, test_name = cfg['original_filename'],  
                        intro_line = intro_line, filename = cfg['original_filename'] + '-contamination-analysis.csv')
                        
    #and a normalised csv also
    normalised_results_csv_builder(results, cfg, test_name = cfg['original_filename'],  
                                   intro_line = intro_line, normalised_by = cfg['normalise_results_by'],
                                   filename = cfg['original_filename']+'-contamination-analysis-normalised.csv',
                                   extra_normalise_exceptions = ['test number', 'air contamination'])  

    # and use the pickle function from pitot_condition_builder.py
    # so we can pick the results and config dictionaries                  
    pickle_result_data(cfg, results, filename = cfg['original_filename']+'-contamination-analysis-final-result-pickle.dat')
    
    # now analyse results dictionary and print some results to the screen
    # and another external file
    
    contamination_analysis_summary(cfg, results)
    
    # now delete the intermediate pickle that we made during the simulation...
    print "Removing the final intermediate pickle file."
    
    if os.path.isfile(intermediate_filename): 
        os.remove(intermediate_filename)
    
    #zip up the final result using the function from pitot_condition_builder.py
    zip_result_and_log_files(cfg, output_filename = cfg['original_filename'] + '-contamination-analysis-log-and-result-files.zip')
    
    return
                                
#----------------------------------------------------------------------------

def main():
    
    import optparse  
    op = optparse.OptionParser(version=VERSION_STRING)   
    op.add_option('-c', '--config_file', '--config-file', dest='config_file',
                  help=("filename where the configuration file is located"))    
    op.add_option('-f', '--force_restart', action = "store_true", dest='force_restart',
                  help=("flag that can be used to force the simulation to restart"
                        "It stops it looking for an unfinished simulation.")) 
    opt, args = op.parse_args()
    config_file = opt.config_file
    force_restart = opt.force_restart
           
    run_pitot_contamination_analysis(cfg = {}, config_file = config_file, force_restart = force_restart)
    
    return
    
#----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "pitot_contamination_analysis.py - Pitot Equilibrium expansion tube simulator air contamination analysis tool"
        print "start with --help for help with inputs"
        
    else:
        main()

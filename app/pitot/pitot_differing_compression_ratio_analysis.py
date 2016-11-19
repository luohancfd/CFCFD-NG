#! /usr/bin/env python
"""
pitot_differing_compression_ratio_analysis.py: pitot differing compression ratio analysis tool

I had to do some analysis with different compression ratios, so I decided that I
should make this. Requires the user to specify a compression ratio list, and then whether
we want to specify driver fill pressure (and expand up to each compression ratio),
or specify the burst pressure and start from the pressure required...

This is basically a cut down version of the pitot condition building program
(pitot_condition_builder.py) and it prints data to the screen and a text file
in a similar way. 

Chris James (c.james4@uq.edu.au) - 25/09/15

"""

VERSION_STRING = "19-Nov-2016"

from pitot_condition_builder import stream_tee, pickle_result_data, pickle_intermediate_data, results_csv_builder, normalised_results_csv_builder, zip_result_and_log_files, cleanup_old_files

import sys, os

from pitot import run_pitot
from pitot_input_utils import *

def check_new_inputs(cfg):
    """Takes the input file and checks that the extra inputs required for the
       compression ratio analysis diluent analysis are working..
    
       Returns the checked over input file and will tell the bigger program to 
        bail out if it finds an issue.
    
    """
    
    print "Starting check of differing compression ratio analysis specific inputs."
    
    if 'compression_ratio_list' not in cfg:
        print "You have not specified a 'compression_ratio_list'. Bailing out."
        cfg['bad_input'] = True
        
    if 'specified_pressure' not in cfg:
        print "You have not specified a 'specified_pressure' value."
        print "This must either be set to 'p4' or 'driver_p'."
        print "Bailing out."
        cfg['bad_input'] = True
    else:
        if cfg['specified_pressure'] not in ['p4', 'driver_p']:
            print "'specified_pressure' value must be set to either 'p4' or 'driver_p'"
            print "Bailing out."
            cfg['bad_input'] = True
        if cfg['specified_pressure'] == 'p4' and 'p4' not in cfg:
            print "'specified_pressure' is set to 'p4' but 'p4' cannot be found in cfg."
            print "Bailing out."
            cfg['bad_input'] = True
        if cfg['specified_pressure'] == 'driver_p' and 'driver_p' not in cfg:
            print "'specified_pressure' is set to 'driver_p' but 'driver_p' cannot be found in cfg."
            print "Bailing out."
            cfg['bad_input'] = True
        if 'driver_T' not in cfg:
            cfg['driver_T'] = 298.15 # K
              
    if 'store_electron_concentration' not in cfg:
        cfg['store_electron_concentration'] = False
        
    if 'calculate_modified_bsp' not in cfg:
        cfg['calculate_modified_bsp'] = False    
        
    if 'normalise_results_by' not in cfg:
        print "'normalise_results_by' not in cfg."
        print "Setting it to 'first value'."
        cfg['normalise_results_by'] = 'first value'  
        
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
    return len(cfg['compression_ratio_list'])
    
def start_message(cfg):
    """Function that takes the cfg file and prints a start message for the
       program, detailing what it will do.
    """
    
    # print how many tests we're going to run, and the ranges.
    
    print '-'*60   
    print 'Running Pitot Gas Giant Differing Compression Ratio Analysis Version: {0}.'.format(VERSION_STRING)    
    print "{0} tests will be run.".format(cfg['number_of_test_runs'])
    
    print "Differing compression ratios will be tested from {0} - {1} in increments of {2}."\
          .format(cfg['compression_ratio_list'][0], cfg['compression_ratio_list'][-1],
                  cfg['compression_ratio_list'][0] - cfg['compression_ratio_list'][1])
    
    if cfg['specified_pressure'] == 'driver_p':
        print "Using an initial driver fill pressure of {0} Pa.".format(cfg['driver_p'])
    elif cfg['specified_pressure'] == 'p4':
        print "Using a driver burst pressure (p4) of {0} Pa.".format(cfg['p4'])        
                   
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
        basic_list = ['test number','compression_ratio','driver_p', 'driver_T', 
                      'p4', 'T4','gamma4', 'R4', 'psd1','p1','p5','Vsd',
                      'Vs1', 'Vs2', 'Ht','h','u_eq', 'rho1', 'gamma1', 'R1', 'MW1',
                      'p2','T2','rho2','V2','M2', 'a2', 'gamma2', 'R2', 'Ht2',
                      'p6','T6','rho6','V6','M6','p7','T7','rho7','V7','M7']
    else:
        basic_list = ['test number','compression_ratio', 'driver_p', 'driver_T', 
                      'p4', 'T4','gamma4', 'R4','p1','p5',
                      'Vs1', 'Vs2', 'Ht','h','u_eq','rho1', 'gamma1', 'R1', 'MW1',
                      'p2','T2','rho2','V2','M2', 'a2', 'gamma2', 'R2', 'Ht2',
                      'p6','T6','rho6','V6','M6','p7','T7','rho7','V7','M7']
    full_list += basic_list
    
    if cfg['nozzle']:
        nozzle_list = ['arearatio','p8','T8','rho8','V8','M8']
        full_list += nozzle_list
    if cfg['conehead']:
         conehead_list = ['p10c','T10c','rho10c','V10c']
         full_list += conehead_list
    if cfg['shock_over_model']:
        shock_over_model_list = ['p10f','T10f','rho10f','V10f',
                                 'p10e','T10e','rho10e','V10e']
        full_list += shock_over_model_list
    if cfg['store_electron_concentration']:     
        store_electron_concentration_list = ['s2ec','s7ec','s8ec','s10ec']
        full_list += store_electron_concentration_list
    if cfg ['calculate_modified_bsp']:
        full_list += ['modified_bsp']

    # now populate the dictionary with a bunch of empty lists based on that list

    results = {title : [] for title in full_list}
    
    # add the list of titles in case we want to use it in future
    
    results['full_list'] = full_list
    
    # and I would like to also store the diluent so we can access it when 
    # outputting results later on
       
    #add a list where we can store unsuccesful run numbers for analysis
    results['unsuccessful_runs'] = []
    
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
    
    for compression_ratio in cfg['compression_ratio_list']:
        # change the test name
        test_name += 1
        # store that test name
        test_condition_input_details['test_names'].append(test_name)
        # now make a dictionary in the test_condition_input_details dictionary for this simulation
        # using the compression ratio
        # this is simple here as compression_ratio is the only variable...
        # (for the normal condition builder it was a bit more complicated...)
        input_dictionary = {'compression_ratio': compression_ratio}
        test_condition_input_details[test_name] = input_dictionary
                
    print "The test_names list for this simulation is:"
    print test_condition_input_details['test_names']
                
    return test_condition_input_details
    
def differing_compression_ratio_analysis_test_run(cfg, results):
    """Function that takes the fully built config dictionary
       and the text file that is being used for the program output
       and does the test run then adds a line to the output file.
    """
    
    condition_status = True #This will be turned to False if the condition fails
    
    cfg['filename'] = cfg['original_filename'] + '-test-{0}-result'.format(cfg['test_number'])
    
    if not cfg['have_checked_time']:
        # we check the amount of time the first run takes and then tell the user...
        import time
        start_time = time.time()
    
    # some code here to make a copy of the stdout printouts for each test and store it
    
    import sys
    
    test_log = open(cfg['original_filename'] + '-test-{0}-log.txt'.format(cfg['test_number']),"w")
    sys.stdout = stream_tee(sys.stdout, test_log)
    
    print '-'*60
    print "Running test {0} of {1}.".format(cfg['test_number'], cfg['number_of_test_runs'])
    print "Current compression ratio is {0}.".format(cfg['compression_ratio'],)
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
    else:
        cfg['last_run_successful'] = False
        results['unsuccessful_runs'].append(cfg['test_number'])
        
    # now we end the stream teeing here by pulling out the original stdout object
    # and overwriting the stream tee with that, then closing the log file
    sys.stdout = sys.stdout.stream1   
    test_log.close()
    
    # need to remove Vs values from the dictionary or it will bail out
    # on the next run         
    if cfg['secondary'] and 'Vsd' in cfg: cfg.pop('Vsd')
    if 'Vs1' in cfg: cfg.pop('Vs1') 
    if 'Vs2' in cfg: cfg.pop('Vs2')
    # for the condition builder stuff also need to remove 'T4' so the next run works fine
    if 'T4' in cfg: cfg.pop('T4')
        
    # we now need to go through and purge the guesses and bounds of the secant solvers
    # if they were not custom, so then new guesses can be made next time
        
    secant_solver_variables = ['Vsd_lower', 'Vsd_upper', 'Vsd_guess_1', 'Vsd_guess_2', 
                               'Vs1_lower', 'Vs1_upper', 'Vs1_guess_1', 'Vs1_guess_2',
                               'Vs2_lower', 'Vs2_upper', 'Vs2_guess_1', 'Vs2_guess_2']
    
    for variable in secant_solver_variables:
        # if the variables are not in the original cfg, they are not custom inputs
        # and were added by the last run, so we remove them
        if variable not in cfg['cfg_original']:
            if variable in cfg: cfg.pop(variable)
                
    if not cfg['have_checked_time']:
        test_time = time.time() - start_time
        print '-'*60
        print "Time to complete first test was {0:.2f} seconds."\
        .format(test_time)
        print "If every test takes this long. It will take roughly {0:.2f} hours to perform all {1} tests."\
        .format(test_time*cfg['number_of_test_runs']/3600.0, cfg['number_of_test_runs'])
        cfg['have_checked_time'] = True    
            
    return condition_status, results
           
def add_new_result_to_results_dict(cfg, states, V, M, results):
    """Function that takes a completed test run and adds the tunnel
       configuration and results to the results dictionary.
    """ 
    
    results['test number'].append(cfg['test_number'])
    results['compression_ratio'].append(cfg['compression_ratio'])
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

    results['driver_p'].append(cfg['driver_p'])
    results['driver_T'].append(cfg['driver_T'])
    results['p4'].append(states['s4'].p)
    results['T4'].append(states['s4'].T)
    results['gamma4'].append(states['s4'].gam)
    results['R4'].append(states['s4'].R)
    
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
    
    if cfg['calculate_modified_bsp']:
        # if the user has asked for it, here we calculate, omega, the modified
        # binary scaling paramater from Stalker and Edward's 1998 Paper
        # Hypersonic Blunt-Body Flows in Hydrogen-Neon Mixtures
        # JOURNAL OF SPACECRAFT AND ROCKETS 
        # Vol. 35, No. 6, November-December 1998
    
        r = cfg['diluent_percentage'] / 100.0
        p = states['s10f'].p
        if cfg['mode'] == 'expansion-tube':
            if cfg['nozzle']:
                epsilon = states['s8'].rho / states['s10f'].rho
                U = V['s8']
            else:
                epsilon = states['s7'].rho / states['s10f'].rho
                U = V['s7']
            
        else:
            if cfg['nozzle']:
                epsilon = states['s8'].rho / states['s10f'].rho
                U = V['s8']
            else:
                epsilon = states['s2'].rho / states['s10f'].rho  
                U = V['s2']
        
        modified_bsp = (r*p*(epsilon*(1.0 - epsilon))**-0.5)/U
        
        results['modified_bsp'].append(modified_bsp)
        
    return results
    
def differing_compression_ratio_analysis_summary(cfg, results):
    """Function that takes the config dictionary and results dictionary 
       made throughout the running of the program and prints a summary of 
       the run to the screen and to a summary text file
    """
    
    print '-'*60
    print "Printing summary to screen and to a text document."
    
    p_d_cr_analysis_summary_file = open(cfg['original_filename']+'-differing-compression-ratio-analysis-summary.txt',"w")
    
    # print lines explaining the results
    summary_line_1 = "# Summary of pitot different compression ratio analysis program output."
    p_d_cr_analysis_summary_file.write(summary_line_1 + '\n')
    summary_line_2 = "# Summary performed using Version {0} of the different compression ratio analysis program.".format(VERSION_STRING)
    p_d_cr_analysis_summary_file.write(summary_line_2 + '\n')
    
    summary_line_3 = "{0} tests ran. {1} ({2:.1f}%) were successful."\
        .format(cfg['number_of_test_runs'], len(results['test number']),
        float(len(results['test number']))/float(cfg['number_of_test_runs'])*100.0)
    print summary_line_3
    p_d_cr_analysis_summary_file.write(summary_line_3 + '\n')

    if results['unsuccessful_runs']: 
        summary_line_4 = "Unsucessful runs were run numbers {0}.".format(results['unsuccessful_runs'])
        print summary_line_4
        p_d_cr_analysis_summary_file.write(summary_line_4 + '\n')          

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
                p_d_cr_analysis_summary_file.write(summary_line + '\n')
                
    p_d_cr_analysis_summary_file.close()   
                
    return
            
def run_pitot_differing_compression_ratio_analysis(cfg = {}, config_file = None, force_restart = False):
    """
    
    Chris James (c.james4@uq.edu.au) 23/12/14
    
    run_pitot_differing_compression_ratio_analysis(dict) - > depends
    
    force_restart can be used to force the simulation to start again instead of 
    looking for an unfinished simulation that may be there...
    
    """
    
    #---------------------- get the inputs set up --------------------------
    
    if config_file:
        cfg = config_loader(config_file)
    
    # some dummy stuff here to pass the initial pitot input test if values are not there
    
    cfg['facility'] = 'custom'
    cfg['compression_ratio'] = 10.0
        
    #----------------- check inputs ----------------------------------------
    
    cfg = input_checker(cfg)
    
    cfg = check_new_inputs(cfg)
    
    intermediate_filename = cfg['filename']+'-differing-compression-ratio-analysis-intermediate-result-pickle.dat'
    
    # now check if we have attempted an old run before or not....
    if not os.path.isfile(intermediate_filename) or force_restart: 
        # if not, we set up a new one...
        
        # clean up any old files if the user has asked for it
        
        if cfg['cleanup_old_files']:
            # build a list of auxiliary files and run the cleanup_old_files function frim pitot_condition_builder.py
            auxiliary_file_list = ['-differing-compression-ratio-analysis-log-and-result-files.zip',
                                   '-differing-compression-ratio-analysis-final-result-pickle.dat',
                                   '-differing-compression-ratio-analysis-normalised.csv',
                                   '-differing-compression-ratio-analysis-summary.txt',
                                   '-differing-compression-ratio-analysis.csv']
            cleanup_old_files(auxiliary_file_list)
            cleanup_old_files()
        
        # make a counter so we can work out what test we're running
        # also make one to store how many runs are successful
        
        cfg['number_of_test_runs'] = calculate_number_of_test_runs(cfg)
        
        # we use this variable to try to speed some things up later on
        cfg['last_run_successful'] = None
        
        import copy
        cfg['original_filename'] = copy.copy(cfg['filename'])
        
        # I think install the original cfg too
    
        cfg['cfg_original'] = copy.copy(cfg)
            
        # print the start message
           
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
        # then set the details that every simulation will have
        cfg['compression_ratio'] = test_condition_input_details[test_name]['compression_ratio']
        
        # if 'specified_pressure' is 'driver_p' we don't need to do anything here,
        # if it's 'p4' we need to get the related 'driver_p'
        
        if cfg['specified_pressure'] == 'p4':
            from cfpylib.gasdyn.cea2_gas import Gas
            primary_driver = Gas(cfg['driver_composition'],inputUnits=cfg['driver_inputUnits'],
                                 outputUnits=cfg['driver_inputUnits'])
            primary_driver.set_pT(101300.0, 300.0)
            pressure_ratio = cfg['compression_ratio']**primary_driver.gam #pressure ratio is compression ratio to the power of gamma
            cfg['driver_p'] = cfg['p4']/pressure_ratio

        run_status, results = differing_compression_ratio_analysis_test_run(cfg, results) 
        if run_status:
            cfg['good_counter'] += 1
        
        # add this to finished simulations list, regardless of whether it finished correctly or not...
        cfg['finished_simulations'].append(test_name)
        
        # now pickle the intermediate result so we can result the simulation if needed...
        pickle_intermediate_data(cfg, results, test_condition_input_details, filename = intermediate_filename)
    
    # now that we're done we can dump the results in various ways using tools from pitot condition builder...
    intro_line = "Output of pitot differing compression ratio analysis program Version {0}."\
                 .format(VERSION_STRING)            
    results_csv_builder(results, test_name = cfg['original_filename'],  
                        intro_line = intro_line, filename = cfg['original_filename']+'-differing-compression-ratio-analysis.csv')
                        
    #and a normalised csv also
    normalised_results_csv_builder(results, test_name = cfg['original_filename'],  
                        intro_line = intro_line, 
                        normalised_by = cfg['normalise_results_by'],
                        filename = cfg['original_filename']+'-differing-compression-ratio-analysis-normalised.csv',
                        extra_normalise_exceptions = ['compression_ratio'])  
    
                       
    pickle_result_data(cfg, results, filename = cfg['original_filename']+'-differing-compression-ratio-analysis-final-result-pickle.dat')
    
    # now delete the intermediate pickle that we made during the simulation...
    print "Removing the final intermediate pickle file."
    
    if os.path.isfile(intermediate_filename): 
        os.remove(intermediate_filename)
    
    # now analyse results dictionary and print some results to the screen
    # and another external file
    
    differing_compression_ratio_analysis_summary(cfg, results)
    
    #zip up the final result using the function from pitot_condition_builder.py
    zip_result_and_log_files(cfg, output_filename = cfg['original_filename'] + '-differing-compression-ratio-analysis-log-and-result-files.zip')    
    
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
           
    run_pitot_differing_compression_ratio_analysis(cfg = {}, config_file = config_file, force_restart = force_restart)
    
    return
    
#----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "pitot_differing_compression_ratio_analysis.py - Pitot Equilibrium expansion tube simulator differing compression ratio analysis"
        print "start with --help for help with inputs"
        
    else:
        main()

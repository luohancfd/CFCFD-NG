#! /usr/bin/env python
"""
pitot_condition_builder.py: pitot condition builder

This file will take a normal pitot input file with a few
extra variables and do a series of pitot runs checking how different 
facility parameters will affect the condition.

Chris James (c.james4@uq.edu.au) - 29/12/13 

"""

VERSION_STRING = "14-June-2016"

import sys

from pitot import run_pitot
from pitot_input_utils import *

class stream_tee(object):
    # Based on https://gist.github.com/327585 by Anand Kunal
    def __init__(self, stream1, stream2):
        self.stream1 = stream1
        self.stream2 = stream2
        self.__missing_method_name = None # Hack!
 
    def __getattribute__(self, name):
        return object.__getattribute__(self, name)
 
    def __getattr__(self, name):
        self.__missing_method_name = name # Could also be a property
        return getattr(self, '__methodmissing__')
 
    def __methodmissing__(self, *args, **kwargs):
            # Emit method call to the log copy
            callable2 = getattr(self.stream2, self.__missing_method_name)
            callable2(*args, **kwargs)
 
            # Emit method call to stdout (stream 1)
            callable1 = getattr(self.stream1, self.__missing_method_name)
            return callable1(*args, **kwargs)

def check_new_inputs(cfg):
    """Takes the input file and checks that the extra inputs required for the
       condition builder are working..
    
    Returns the checked over input file and will tell the bigger program to 
    bail out if it finds an issue.
    
    """
    
    cfg['bad_input'] = False
    
    print "Starting check of condition builder specific inputs."
    
    if 'driver_condition_list' not in cfg:
        print "'driver_condition_list' variable not specified. Bailing out."
        print "Variable needs to be a list containing valid driver conditions."
        cfg['bad_input'] = True
        
    driver_condition_check_list = ['He:0.80,Ar:0.20', 'He:0.90,Ar:0.10', 'He:1.0', 
                                   'He:0.60,Ar:0.40','custom']
                                   
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
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel' or cfg['tunnel_mode'] == 'reflected-shock-tunnel':
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
    
def start_message(cfg):
    """Function that takes the cfg file and prints a start message for the
       program, detailing what it will do.
    """
    
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
            
    return cfg
    
def build_results_dict(cfg):
    """Function that looks at the cfg dictionary and works out what data needs
       to be stored for the type of test that we're running. Then it populates
       a dictionary called 'results' with empty lists for the data to be stored in.
       The dictionary is then returned.
    """
    
    # need to make a list to create a series of empty lists in the results
    # dictionary to store the data. the list is tailored to the test we're running
    
    print '-'*60
    print "Building the results dictionary."
    
    full_list = ['test number','driver condition']
       
    if cfg['tunnel_mode'] == 'expansion-tube':
        if True in cfg['secondary_list']:
            pressure_shock_list = ['psd1','p1','p5','Vsd','Vs1', 'Vs2']
        else:
            pressure_shock_list = ['p1','p5','Vs1', 'Vs2']
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel':
        if True in cfg['secondary_list']:
            pressure_shock_list = ['psd1','p1','Vsd','Vs1']
        else:
            pressure_shock_list = ['p1','Vs1']
    elif cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        if True in cfg['secondary_list']:
            pressure_shock_list = ['psd1','p1','Vsd','Vs1','Vr','Mr','Vr_d','Mr_d']
        else:
            pressure_shock_list = ['p1','Vs1','Vr','Mr','Vr_d','Mr_d']            
            
    full_list += pressure_shock_list
    
    if True in cfg['secondary_list']:
        statesd2_list = ['psd2','Tsd2','rhosd2','Vsd2','Msd2', 'asd2', 'gammasd2', 'Rsd2']
        full_list += statesd2_list
    if cfg['rs_out_of_sd']:
        statesd2_reflected_list = ['Vr-sd','Mr-sd','psd2r','Tsd2r','rhosd2r','Vsd2r','Msd2r', 'asd2r', 'gammasd2r', 'Rsd2r']
        full_list += statesd2_reflected_list

    state2_list = ['p2','T2','rho2','V2','M2', 'a2', 'gamma2', 'R2', 'Vs1 - V2', 'Ht2']
    full_list += state2_list

    if cfg['rs_out_of_st']:
        state2_reflected_list = ['Vr-st','Mr-st','p2r','T2r','rho2r','V2r','M2r', 'a2r', 'gamma2r', 'R2r', 'Ht2r']
        full_list += state2_reflected_list
    
    if cfg['tunnel_mode'] == 'expansion-tube':
        state7_list = ['p7','T7','rho7','V7','M7']
        full_list += state7_list
        
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        state5_list = ['p5','T5','rho5','V5','M5', 'p5_d','T5_d','rho5_d','V5_d','M5_d']
        full_list += state5_list        

    enthalpy_list = ['Ht','h','u_eq']
    full_list += enthalpy_list 
    
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
        
    if cfg['calculate_scaling_information']:
        calculate_scaling_information_freestream_list = ['freestream_mu', 'freestream_rhoL', 'freestream_pL', 'freestream_Re', 'freestream_Kn']
        full_list += calculate_scaling_information_freestream_list
        
        if cfg['shock_over_model']:
            calculate_scaling_information_normal_shock_list = ['s10e_mu', 's10e_rhoL', 's10e_pL', 's10e_Re', 's10e_Kn']
            full_list += calculate_scaling_information_normal_shock_list            
            
    # now populate the dictionary with a bunch of empty lists based on that list

    results = {title : [] for title in full_list}
    
    # add the list of titles in case we want to use it in future
    
    results['full_list'] = full_list
    
    #add a list where we can store unsuccesful run numbers for analysis
    results['unsuccessful_runs'] = []
    
    print "The full list of variables to be added to the output are:"
    print full_list
    
    return results
    
def build_test_condition_input_details_dictionary(cfg):
    """Function that builds a dictionary that tells the condition
       builder what tests to run. It will include a list called 'test_names'
       that will be cycled through by the main condition program.
       
       iterools.product() is used to perform a cartesian product of all
       of the input variables (ie. to produce every combination from the right)
    """
    
    from itertools import product # bring in our cartesian product function
    
    test_condition_input_details = {}
    test_condition_input_details['test_names'] = []
    
    test_name = 0 # we will add to this as we go, first test_name will be 1
    
    print '-'*60
    print "Building the test condition details dictionary containing a dictionary for each simulation."    
    
    if cfg['tunnel_mode'] == 'expansion-tube':
        if True in cfg['secondary_list']:
            for test_settings in product(cfg['driver_condition_list'], cfg['secondary_list'],
                                         cfg['psd1_list'], cfg['p1_list'], cfg['p5_list']):
                # change the test name
                test_name += 1
                # store that test name
                test_condition_input_details['test_names'].append(test_name)
                # now make a dictionary in the test_condition_input_details dictionary for this simulation
                # using the current test_settings combination
                input_dictionary = {'driver_condition': test_settings[0],
                                    'secondary_value': test_settings[1],
                                    'psd1': test_settings[2], 'p1': test_settings[3],
                                    'p5': test_settings[4]}
                test_condition_input_details[test_name] = input_dictionary
                
        else:
            for test_settings in product(cfg['driver_condition_list'], cfg['secondary_list'],
                                         cfg['p1_list'], cfg['p5_list']):
                # change the test name
                test_name += 1
                # store that test name
                test_condition_input_details['test_names'].append(test_name)
                # now make a dictionary in the test_condition_input_details dictionary for this simulation
                # using the current test_settings combination
                input_dictionary = {'driver_condition': test_settings[0],
                                    'secondary_value': test_settings[1],
                                    'p1': test_settings[2],'p5': test_settings[3]}
                test_condition_input_details[test_name] = input_dictionary
    
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel' or cfg['tunnel_mode'] == 'reflected-shock-tunnel':   
        if True in cfg['secondary_list']:
            for test_settings in product(cfg['driver_condition_list'], cfg['secondary_list'],
                                         cfg['psd1_list'], cfg['p1_list']):
                # change the test name
                test_name += 1
                # store that test name
                test_condition_input_details['test_names'].append(test_name)
                # now make a dictionary in the test_condition_input_details dictionary for this simulation
                # using the current test_settings combination
                input_dictionary = {'driver_condition': test_settings[0],
                                    'secondary_value': test_settings[1],
                                    'psd1': test_settings[2], 'p1': test_settings[3]}
                test_condition_input_details[test_name] = input_dictionary
        else:
            for test_settings in product(cfg['driver_condition_list'], cfg['secondary_list'], cfg['p1_list']):
                # change the test name
                test_name += 1
                # store that test name
                test_condition_input_details['test_names'].append(test_name)
                # now make a dictionary in the test_condition_input_details dictionary for this simulation
                # using the current test_settings combination
                input_dictionary = {'driver_condition': test_settings[0],
                                    'secondary_value': test_settings[1],
                                    'p1': test_settings[2]}
                test_condition_input_details[test_name] = input_dictionary
                
    print "The test_names list for this simulation is:"
    print test_condition_input_details['test_names']
                
    return test_condition_input_details
    
def condition_builder_test_run(cfg, results):
    """Function that takes the fully built config dictionary
       and the text file that is being used for the program output
       and does the test run then adds a line to the output file.
    """
    
    condition_status = True #This will be turned to False if the condition fails
    
    if not cfg['have_checked_time']:
        # we check the amount of time the first run takes and then tell the user...
        import time
        start_time = time.time()
    
    # first we check if we should slightly modify our guesses based on the last
    # successful run to speed up the code.
    
    if cfg['last_run_successful']:
        cfg, results = guess_modifier(cfg, results)
    
    cfg['filename'] = cfg['original_filename'] + '-test-{0}'.format(cfg['test_number'])
    
    # some code here to make a copy of the stdout printouts for each test and store it
    
    import sys
    
    test_log = open(cfg['filename']+'-log.txt',"w")
    sys.stdout = stream_tee(sys.stdout, test_log)
    
    print '-'*60
    print "Running test {0} of {1}.".format(cfg['test_number'], cfg['number_of_test_runs'])
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
        print "Simulation {0} was successful. Result will be added to the output.".format(cfg['test_number'])
        results = add_new_result_to_results_dict(cfg, states, V, M, results)
        cfg['last_run_successful'] = True
    else:
        print "Simulation {0} was unsuccessful. Result will NOT be added to the output.".format(cfg['test_number'])
        cfg['last_run_successful'] = False
        results['unsuccessful_runs'].append(cfg['test_number'])
    
    # need to remove Vs values from the dictionary or it will bail out
    # on the next run                 
    if cfg['secondary'] and 'Vsd' in cfg: cfg.pop('Vsd')
    if 'Vs1' in cfg: cfg.pop('Vs1') 
    if 'Vs2' in cfg: cfg.pop('Vs2')
    
    # now we end the stream teeing here by pulling out the original stdout object
    # and overwriting the stream tee with that, then closing the log file
    sys.stdout = sys.stdout.stream1   
    test_log.close()
    
    if not cfg['have_checked_time']:
        test_time = time.time() - start_time
        print '-'*60
        print "Time to complete first test was {0:.2f} seconds."\
        .format(test_time)
        print "If every test takes this long. It will take roughly {0:.2f} hours to perform all {1} tests."\
        .format(test_time*cfg['number_of_test_runs']/3600.0, cfg['number_of_test_runs'])
        cfg['have_checked_time'] = True    
            
    return condition_status, results
    
def guess_modifier(cfg, results):
    """Function that checks the results dictionary and sees how similar the last test
       was to the one about to run, and will modify the starting guesses for the
       current run to try to speed it up if they are similar.
    """
    
    print "-"*60
    print "Checking if we can modify any starting guesses to speed up the program."
    
    # start by checking they both used the same driver condition
    if results['driver condition'][-1] == cfg['driver_gas']:
        if cfg['secondary']:
            # check the fill pressures were the same
            if results['psd1'][-1] == cfg['psd1']:
                # then make the starting guesses basically the right answer
                cfg['Vsd_guess_1'] = results['Vsd'][-1] - 1.0
                cfg['Vsd_guess_2'] = results['Vsd'][-1] + 1.0
                print "Changed guesses for Vsd secant solver."
                print "('Vsd_guess_1' = {0} m/s and 'Vsd_guess_2' = {1} m/s)".\
                  format(cfg['Vsd_guess_1'], cfg['Vsd_guess_2'])
                # it's important that we do this part in here, or otherwise we may
                # be setting new p1 guesses which are not suitable if psd1 is changing
                if results['p1'][-1] == cfg['p1']:
                    cfg['Vs1_guess_1'] = results['Vs1'][-1] - 1.0
                    cfg['Vs1_guess_2'] = results['Vs1'][-1] + 1.0
                    print "Changed guesses for Vs1 secant solver."
                    print "('Vs1_guess_1' = {0} m/s and 'Vs1_guess_2' = {1} m/s)".\
                  format(cfg['Vs1_guess_1'], cfg['Vs1_guess_2'])
        else:
            if results['p1'][-1] == cfg['p1']:
                cfg['Vs1_guess_1'] = results['Vs1'][-1] - 1.0
                cfg['Vs1_guess_2'] = results['Vs1'][-1] + 1.0
                print "Changed guesses for Vs1 secant solver."
                print "('Vs1_guess_1' = {0} m/s and 'Vs1_guess_2' = {1} m/s)".\
                      format(cfg['Vs1_guess_1'], cfg['Vs1_guess_2'])
                
    return cfg, results
       
def results_csv_builder(results, test_name = 'pitot_run',  intro_line = None, filename = None):
    """Function that takes the final results dictionary (which must include a 
       list called 'full_list' that tells this function what to print and in 
       what order) and then outputs a results csv. It will also add an intro line
       if a string with that is added. The name of the test is also required.
    """
    
#   # make a file name if nothing is specified...    
    
    if not filename:
        filename = test_name + '-condition-builder.csv'
       
    with open(filename,"w") as condition_builder_output:
    
        # print a line explaining the results if the user gives it
        if intro_line:
            intro_line_optional = "# " + intro_line
            condition_builder_output.write(intro_line_optional + '\n')
        
        #now we'll make the code build us the second intro line
        intro_line = '#'
        for value in results['full_list']:
            if value != results['full_list'][-1]:
                intro_line += "{0},".format(value)
            else: #don't put the comma if it's the last value
                intro_line += "{0}".format(value)
            
        condition_builder_output.write(intro_line + '\n')
        
        # now we need to go through every test run and print the data.
        # we'll use 'full_list' to guide our way through
        
        # get the number of the test runs from the length of the first data list mentioned
        # in 'full_list'. need to assume the user hasn't screwed up and got lists of
        # different lengths
        number_of_test_runs = len(results[results['full_list'][0]])
        
        for i in range(0, number_of_test_runs, 1):
            output_line = ''
            for value in results['full_list']:
                if value != results['full_list'][-1]:
                    output_line += "{0},".format(results[value][i])
                else: #don't put the comma if it's the last value in the csv
                    output_line += "{0}".format(results[value][i])
            
            condition_builder_output.write(output_line + '\n')  
    
        condition_builder_output.close()              
                                  
    return 
    
def normalised_results_csv_builder(results, test_name = 'pitot_run',  
                                   intro_line = None, normalised_by = 'first value'):
    """Function that takes the final results dictionary (which must include a 
       list called 'full_list' that tells this function what to print and in 
       what order) and then outputs a normalised version of the results csv.
       You can tell it to normalise by other values, but 'first value' is default.
       
       It will also add an intro line if a string with that is added. 
       The name of the test is also required.
    """
    
    # open a file to start saving results
    with open(test_name + '-condition-builder-normalised.csv',"w") as condition_builder_output:
    
        # print a line explaining the results if the user gives it
        if intro_line:
            intro_line_optional = "# " + intro_line
            condition_builder_output.write(intro_line_optional + '\n')
            
        normalised_intro_line = "# all variables normalised by {0}".format(normalised_by)
        condition_builder_output.write(normalised_intro_line + '\n')
        
        # 'test number' and 'diluent percentage' and the species concentrations
        # will not be normalised
        normalise_exceptions = ['test number', 'driver condition', 'psd1','Vsd']
        
        # we need to make sure that we don't try to normalise these values if they are all empty
        if 'Vsd2r' in results:
            all_zeros = True
            for value in results['Vsd2r']:
                if bool(value):
                    all_zeros = False
            if all_zeros:
                normalise_exceptions += ['Vsd2r', 'Msd2r']
       
        if 'V2r' in results:
            all_zeros = True
            for value in results['V2r']:
                if bool(value):
                    all_zeros = False
            if all_zeros:
                normalise_exceptions += ['V2r', 'M2r']
        
        #now we'll make the code build us the second intro line
        intro_line = '#'
        for value in results['full_list']:
            if value != results['full_list'][-1]:
                if value in normalise_exceptions:
                    intro_line += "{0},".format(value)
                else:
                    intro_line += "{0} normalised,".format(value)
            else: #don't put the comma if it's the last value
                if value in normalise_exceptions:
                    intro_line += "{0}".format(value)
                else:
                    intro_line += "{0} normalised".format(value)
            
        condition_builder_output.write(intro_line + '\n')
        
        # now we need to go through every test run and print the data.
        # we'll use 'full_list' to guide our way through
        
        # get the number of the test runs from the length of the first data list mentioned
        # in 'full_list'. need to assume the user hasn't screwed up and got lists of
        # different lengths
        number_of_test_runs = len(results[results['full_list'][0]])
        
        # build a dictionary to store all of our normalisation values
        normalising_value_dict = {}
        
        for value in results['full_list']:
            if normalised_by == 'first value':
                normalising_value_dict[value] = results[value][0]
            elif normalised_by == 'maximum value':
                normalising_value_dict[value] = max(results[value])
            elif normalised_by == 'last value':
                normalising_value_dict[value] = results[value][-1] 
        
        for i in range(0, number_of_test_runs, 1):
            output_line = ''
            for value in results['full_list']:
                if value != results['full_list'][-1]:
                    # don't normalise selected exceptions, or values that are not numbers
                    # or a value that is not a number
                    if value in normalise_exceptions or \
                    not isinstance(results[value][i], (int, float)):
                        output_line += "{0},".format(results[value][i])
                    else:
                        output_line += "{0},".format(results[value][i]/normalising_value_dict[value])
                else: #don't put the comma if it's the last value in the csv
                    if value in normalise_exceptions or \
                    not isinstance(results[value][i], (int, float)):
                        output_line += "{0},".format(results[value][i])
                    else:  # only normalise if the value is a number
                        output_line += "{0}".format(results[value][i]/normalising_value_dict[value])
            
            condition_builder_output.write(output_line + '\n')  
    
        condition_builder_output.close()              
                                  
    return     
    
def add_new_result_to_results_dict(cfg, states, V, M, results):
    """Function that takes a completed test run and adds the tunnel
       configuration and results to the results dictionary.
    """ 
           
    # needed to change these as the extra comma was screwing up the csv
    if cfg['driver_gas'] == 'He:1.0':
        driver_condition = cfg['driver_gas']
    elif cfg['driver_gas'] == 'He:0.90,Ar:0.10':
        driver_condition = 'He:0.9 Ar:0.1'
    elif cfg['driver_gas'] == 'He:0.80,Ar:0.20':
        driver_condition = 'He:0.8 Ar:0.2'
    elif cfg['driver_gas'] == 'He:0.60,Ar:0.40':
        driver_condition = 'He:0.6 Ar:0.4'
    elif cfg['driver_gas'] == 'custom':
        driver_condition = 'custom'
    
    results['test number'].append(cfg['test_number'])
    results['driver condition'].append(driver_condition)
    if True in cfg['secondary_list']:
        if cfg['secondary']:
            results['psd1'].append(cfg['psd1'])
        else:
            results['psd1'].append('Not used')
    results['p1'].append(cfg['p1'])
    if cfg['tunnel_mode'] == 'expansion-tube':
        results['p5'].append(cfg['p5'])
            
    if True in cfg['secondary_list']:
        if cfg['secondary']:
            results['Vsd'].append(cfg['Vsd'])
        else:
            results['Vsd'].append('N/A')
    results['Vs1'].append(cfg['Vs1'])
    if cfg['tunnel_mode'] == 'expansion-tube':
        results['Vs2'].append(cfg['Vs2'])
        
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        results['Vr'].append(cfg['Vr'])
        results['Mr'].append(cfg['Mr'])   
        results['Vr_d'].append(cfg['Vr_d'])
        results['Mr_d'].append(cfg['Mr_d'])   
      
    if True in cfg['secondary_list']: 
        if cfg['secondary']:        
            results['psd2'].append(states['sd2'].p)
            results['Tsd2'].append(states['sd2'].T)
            results['rhosd2'].append(states['sd2'].rho)
            results['Vsd2'].append(V['sd2'])
            results['Msd2'].append(M['sd2'])
            results['asd2'].append(states['sd2'].a)
            results['gammasd2'].append(states['sd2'].gam)
            results['Rsd2'].append(states['sd2'].R)
            if cfg['rs_out_of_sd']:
                results['Vr-sd'].append(cfg['Vr-sd'])
                results['Mr-sd'].append(cfg['Mr-sd'])
                results['psd2r'].append(states['sd2r'].p)
                results['Tsd2r'].append(states['sd2r'].T)
                results['rhosd2r'].append(states['sd2r'].rho)
                results['Vsd2r'].append(V['sd2r'])
                results['Msd2r'].append(M['sd2r'])
                results['asd2r'].append(states['sd2r'].a)
                results['gammasd2r'].append(states['sd2r'].gam)
                results['Rsd2r'].append(states['sd2r'].R)
        else:
            results['psd2'].append('Not used')
            results['Tsd2'].append('Not used')
            results['rhosd2'].append('Not used')
            results['Vsd2'].append('Not used')
            results['Msd2'].append('Not used')
            results['asd2'].append('Not used')
            results['gammasd2'].append('Not used')
            results['Rsd2'].append('Not used') 
            if cfg['rs_out_of_sd']:
                results['Vr-sd'].append('Not used')
                results['Mr-sd'].append('Not used')
                results['psd2r'].append('Not used')
                results['Tsd2r'].append('Not used')
                results['rhosd2r'].append('Not used')
                results['Vsd2r'].append('Not used')
                results['Msd2r'].append('Not used')
                results['asd2r'].append('Not used')
                results['gammasd2r'].append('Not used')
                results['Rsd2r'].append('Not used')
            
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
        
    if cfg['rs_out_of_st']:
        results['Vr-st'].append(cfg['Vr-st'])
        results['Mr-st'].append(cfg['Mr-st'])
        results['p2r'].append(states['s2r'].p)
        results['T2r'].append(states['s2r'].T)
        results['rho2r'].append(states['s2r'].rho)
        results['V2r'].append(V['s2r'])
        results['M2r'].append(M['s2r'])
        results['a2r'].append(states['s2r'].a)
        results['gamma2r'].append(states['s2r'].gam)
        results['R2r'].append(states['s2r'].R)
        if cfg['Ht2r']:
            results['Ht2r'].append(cfg['Ht2r']/10**6)
        else:
            results['Ht2r'].append('did not solve')
  
    if cfg['tunnel_mode'] == 'expansion-tube':
        results['p7'].append(states['s7'].p)
        results['T7'].append(states['s7'].T)
        results['rho7'].append(states['s7'].rho)
        results['V7'].append(V['s7'])
        results['M7'].append(M['s7'])
        
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        results['p5'].append(states['s5'].p)
        results['T5'].append(states['s5'].T)
        results['rho5'].append(states['s5'].rho)
        results['V5'].append(V['s5'])
        results['M5'].append(M['s5'])
        
        results['p5_d'].append(states['s5_d'].p)
        results['T5_d'].append(states['s5_d'].T)
        results['rho5_d'].append(states['s5_d'].rho)
        results['V5_d'].append(V['s5_d'])
        results['M5_d'].append(M['s5_d'])
        
    if cfg['stagnation_enthalpy']:
        results['Ht'].append(cfg['stagnation_enthalpy']/10**6)
    else:
        results['Ht'].append('did not solve')
    results['h'].append(cfg['freestream_enthalpy']/10**6)    
    if cfg['u_eq']:
        results['u_eq'].append(cfg['u_eq'])
    else:
        results['u_eq'].append('did not solve')
    
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
                
    if cfg['calculate_scaling_information']:
        results['freestream_mu'].append(states[cfg['test_section_state']].mu)        
        results['freestream_rhoL'].append(cfg['rho_l_product_freestream'])
        results['freestream_pL'].append(cfg['pressure_l_product_freestream'])  
        results['freestream_Re'].append(cfg['reynolds_number_freestream']) 
        results['freestream_Kn'].append(cfg['knudsen_number_freestream'])          
        
        if cfg['shock_over_model']:
            results['s10e_mu'].append(states['s10e'].mu) 
            results['s10e_rhoL'].append(cfg['rho_l_product_state10e'])
            results['s10e_pL'].append(cfg['pressure_l_product_state10e'])  
            results['s10e_Re'].append(cfg['reynolds_number_state10e']) 
            results['s10e_Kn'].append(cfg['knudsen_number_state10e'])             
                     
    return results
    
def condition_builder_summary_builder(cfg, results):
    """Function that takes the config dictionary and results dictionary 
       made throughout the running of the program and prints a summary of 
       the run to the screen and to a summary text file.
    """
    
    condition_builder_summary_file = open(cfg['original_filename']+'-condition-builder-summary.txt',"w")
    # print a line explaining the results
    summary_line_1 = "# Summary of pitot condition building program output."
    condition_builder_summary_file.write(summary_line_1 + '\n')
    summary_line_2 = "# Summary performed using Version {0} of the condition building program.".format(VERSION_STRING)
    condition_builder_summary_file.write(summary_line_2 + '\n')
    
    print '-'*60
    print "Printing summary to screen and to a text document."
    
    summary_line_3 = "{0} tests ran. {1} ({2:.1f}%) were successful."\
        .format(cfg['number_of_test_runs'], len(results['test number']),
        float(len(results['test number']))/float(cfg['number_of_test_runs'])*100.0)
    print summary_line_3
    condition_builder_summary_file.write(summary_line_3 + '\n')  

    if results['unsuccessful_runs']: 
        summary_line_4 = "Unsucessful runs were run numbers {0}.".format(results['unsuccessful_runs'])
        print summary_line_4
        condition_builder_summary_file.write(summary_line_4 + '\n')  

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
                elif variable[0] == 'V' or variable[0] == 'u':
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
                elif variable[-2:] == 'ec':
                    summary_line = "Variable {0} varies from {1:.7f} - {2:.7f} (mole fraction)."\
                    .format(variable, min_value, max_value)
                elif 'mu' in variable:
                    summary_line = "Variable {0} varies from {1:.4e} - {2:.4e} Pa.s."\
                    .format(variable, min_value, max_value)
                elif 'rhoL' in variable:
                    summary_line = "Variable {0} varies from {1:.7f} - {2:.7f} kg/m**2."\
                    .format(variable, min_value, max_value)      
                elif 'pL' in variable:
                    summary_line = "Variable {0} varies from {1:.2f} - {2:.2f} Pa.m."\
                    .format(variable, min_value, max_value)  
                elif 'Re' in variable:
                    summary_line = "Variable {0} varies from {1:.2f} - {2:.2f}."\
                    .format(variable, min_value, max_value)   
                else:
                    summary_line = "Variable {0} varies from {1:.1f} - {2:.1f}."\
                    .format(variable, min_value, max_value)
                print summary_line
                condition_builder_summary_file.write(summary_line + '\n')
    
    condition_builder_summary_file.close()   
            
    return
    
def pickle_result_data(cfg, results):
    """Function that takes the config and results dictionaries 
       made throughout the running of the program and dumps them in another
       dictionary in a pickle object. Basically, this means the dictionaries can
       be "unpickled" and analysed by the user directly without needing to data 
       import the csv.
       
       The file can then be opened like this:
       
       import pickle
       data_file = open('file_location')
       cfg_and_results = pickle.load(data_file)
       data_file.close()
    """
    
    import pickle
    
    print '-'*60
    print "Pickling final cfg and results dictionaries to store the end result."
    
    pickle_file = open(cfg['original_filename']+'-condition-builder-final-result-pickle.dat',"w")
    
    cfg_and_results = {'cfg':cfg, 'results':results}
    
    pickle.dump(cfg_and_results, pickle_file)
    pickle_file.close()
   
    return
    
def pickle_intermediate_data(cfg, results, test_condition_input_details):
    """This function mainly exists for the code so that it can be restart itself if necessary.
        At the end of each run, the code pickles the current cfg, results, and 
        test_condition_input_details dictionaries so that the code can then
        pick them back up and keep going if the code ever bails out.
    """
    
    import pickle
    
    print '-'*60
    print "Pickling cfg, results, and test_condition_input_details of the intermediate result."
    print "This functionality allows the simulation to be restarted if it is stopped."
    
    with open(cfg['original_filename']+'-condition-builder-intermediate-result-pickle.dat',"w") as pickle_file:
    
        intermediate_result = {'cfg':cfg, 'results':results, 'test_condition_input_details': test_condition_input_details}
        
        pickle.dump( intermediate_result, pickle_file)
        pickle_file.close()
   
    return
    
def cleanup_old_files():
    """Function that will remove any files in the current directory ending with
       .txt, .csv, or .dat. This is handy because even if you think the code will
       overwrite any old runs with new ones, you may have a situation where the new run
       fails and the old run files are kept instead. Better to start clean.
    """
    
    print "Cleaning up any old condition builder files in the current folder."
    
    import os
    
    # start by getting our current working directory
    cwd = os.getcwd()
    
    # now get a list of files and folders in the current working directory
    
    file_list = os.listdir(cwd)

    # now loop through and remove anything ending in '.csv', '.txt', or '.dat'

    for filename in file_list:
        if filename[-4:] in ['.csv', '.txt', '.dat']:
            if os.path.isfile(filename): os.remove(filename)
    
    return
            
def run_pitot_condition_builder(cfg = {}, config_file = None, force_restart = False):
    """
    
    Chris James (c.james4@uq.edu.au) 27/12/13
    
    run_pitot_condition_builder(dict) - > depends
    
    force_restart can be used to force the simulation to start again instead of 
    looking for an unfinished simulation that may be there...
    
    """
    
    import os       
        
    #---------------------- get the inputs set up --------------------------
    
    if config_file:
        cfg = config_loader(config_file)
        
    #----------------- check inputs ----------------------------------------
    
    # need to get secondary into cfg
    
    if True in cfg['secondary_list']:
        cfg['secondary'] = True
    else:
        cfg['secondary'] = False
    
    cfg = input_checker(cfg, condition_builder = True)
    
    cfg = check_new_inputs(cfg)
    
    # now check if we have attempted an old run before or not....
    if not os.path.isfile(cfg['filename']+'-condition-builder-intermediate-result-pickle.dat') or force_restart: 
        # if not, we set up a new one...
        
        # clean up any old files if the user has asked for it
        
        if cfg['cleanup_old_files']:
            cleanup_old_files()
        
        # make a counter so we can work out what test we're running
        # also make one to store how many runs are successful
        
        cfg['number_of_test_runs'] = calculate_number_of_test_runs(cfg)
        
        # we use this variable to try to speed some things up later on
        cfg['last_run_successful'] = None
        
        import copy
        cfg['original_filename'] = copy.copy(cfg['filename'])
            
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
        print "If this is not what you want, please delete the file '{0}' and run the condition builder again.".format(cfg['filename']+'-condition-builder-intermediate-result-pickle.dat')
    
        import pickle
        with open(cfg['filename']+'-condition-builder-intermediate-result-pickle.dat',"rU") as data_file:
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
        # first set the details that every simulation will have
        cfg['driver_gas'] = test_condition_input_details[test_name]['driver_condition']
        cfg['secondary'] = test_condition_input_details[test_name]['secondary_value']
        cfg['p1'] = test_condition_input_details[test_name]['p1']
        # now the optional parameters...
        if 'psd1' in test_condition_input_details[test_name]:
            cfg['psd1'] = test_condition_input_details[test_name]['psd1']  
        if 'p5' in test_condition_input_details[test_name]:
            cfg['p5'] = test_condition_input_details[test_name]['p5']  

        run_status, results = condition_builder_test_run(cfg, results) 
        if run_status:
            cfg['good_counter'] += 1
        
        # add this to finished simulations list, regardless of whether it finished correctly or not...
        cfg['finished_simulations'].append(test_name)
        
        # now pickle the intermediate result so we can result the simulation if needed...
        pickle_intermediate_data(cfg, results, test_condition_input_details)

    # now analyse results dictionary and print some results to the screen
    # and another external file
    
    condition_builder_summary_builder(cfg, results) 
    
    # and dump the results to a csv
    intro_line = "Output of pitot condition building program Version {0}.".format(VERSION_STRING)            
    results_csv_builder(results, test_name = cfg['original_filename'],  
                        intro_line = intro_line)
    normalised_results_csv_builder(results, test_name = cfg['original_filename'],  
                        intro_line = intro_line)
                        
    # and a to pickled object the user can load with pickle
    # (this allows the cfg and results dictionaries to be loaded directly)
    # it just pickles the dictionaries to pitot should not be needed to load
    # this data
    pickle_result_data(cfg, results)
    
    # now delete the intermediate pickle that we made during the simulation...
    print "Removing the final intermediate pickle file."
    import os
    
    if os.path.isfile(cfg['original_filename']+'-condition-builder-intermediate-result-pickle.dat'): 
        os.remove(cfg['original_filename']+'-condition-builder-intermediate-result-pickle.dat')
        
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
           
    run_pitot_condition_builder(cfg = {}, config_file = config_file, force_restart = force_restart)
    
    return
    
#----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "pitot_condition_builder.py - Pitot Equilibrium expansion tube simulator condition building tool"
        print "start with --help for help with inputs"
        
    else:
        main()

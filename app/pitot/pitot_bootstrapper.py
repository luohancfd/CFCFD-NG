#! /usr/bin/env python
"""
pitot_bootstrapper.py: pitot condition builder

This is a file to do a statistical bootstrapping analysis with PITOT.

Chris James (c.james4@uq.edu.au) - 07/12/16

"""

VERSION_STRING = "15-Dec-2016"

import sys, os

from pitot import run_pitot
from pitot_input_utils import *

from pitot_condition_builder import stream_tee, pickle_result_data, pickle_intermediate_data, results_csv_builder, normalised_results_csv_builder, cleanup_old_files, zip_result_and_log_files, build_results_dict, add_new_result_to_results_dict 

def check_pitot_bootstrapper_inputs(cfg):
    """Takes the input file and checks that the extra inputs required for the
       gas giant diluent analysis are working..
    
       Returns the checked over input file and will tell the bigger program to 
      bail out if it finds an issue.
    
    """
    
    print "Starting check of pitot bootstrapper specific inputs."
                
    if 'store_electron_concentration' not in cfg:
        cfg['store_electron_concentration'] = False
        
    if 'calculate_modified_bsp' not in cfg:
        cfg['calculate_modified_bsp'] = False    
        
    if 'cleanup_old_files' not in cfg:
        print "'cleanup_old_files' variable not set. Setting to default value of 'False'"
        cfg['cleanup_old_files'] = False
        
    if 'number_of_test_runs' not in cfg:
        print "'number_of_test_runs' not in cfg"
        print "Pitot bootstrapper cannot be ran without a number of test runs. Please provide one."
        raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() 'number_of_test_runs' not in cfg"
        
    if not isinstance(cfg['number_of_test_runs'], int):
        print "'number_of_test_runs' is not an int."
        raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() 'number_of_test_runs' is not an int."
                
    if 'variable_list' not in cfg:
        print "'variable_list' not in cfg"
        print "Pitot bootstrapper cannot be ran without a variable list. Please provide one."
        raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() 'variable_list' not in cfg"
        
    if 'store_mole_fractions' not in cfg:
        print "'store_mole_fractions' variable not set. Setting to default value of 'False'"
        cfg['store_mole_fractions'] = False   
        
    if not isinstance(cfg['variable_list'], list):
        print "'variable_list' is not actually a list!"
        raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() 'variable_list' is not actually a list!"
    
    distribution_list = ['uniform', 'normal','lognormal','gaussian']
    
    for variable in cfg['variable_list']:
        if not isinstance(variable, str):
             print "variable {0} in 'variable_list' is not a string!".format(variable)
             raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() 'variable_list' contains a non string variable ({0})".format(variable)
        if variable not in cfg:
             print "variable {0} in 'variable_list' is not actually in the cfg!".format(variable)
             raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() variable {0} in 'variable_list' is not actually in the cfg!".format(variable)
        if variable + '_distribution' not in cfg:
             print "variable {0} in 'variable_list' does not have a related {0}_distribution value.".format(variable)
             raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() variable {0} in 'variable_list' does not have a related {0}_distribution value.".format(variable)
        if not isinstance(cfg[variable + '_distribution'], str):
            print 'input for {0} is not a string value {1}.'.format(variable + '_distribution', cfg[variable + '_distribution'])
            raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() input for {0} is not a string value".format(variable + '_distribution')
        if cfg[variable + '_distribution'] not in distribution_list:
            print "Distribution value for this variable {0} ({1}) not in valid distribution list.".format(variable, cfg[variable + '_distribution'])
            print "Valid distributions are: ", distribution_list
            raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() Distribution value for this variable {0} ({1}) not in valid distribution list.".format(variable, cfg[variable + '_distribution'])

        # now we move onto things needed for each distribution
        if cfg[variable + '_distribution'] == 'uniform':
            # the main value is the mean, and we either need to provide a delta around that
            # or a range. 
            if variable + '_delta' in cfg and variable + '_range' in cfg:
                print "Variable {0} has both an _delta and an _range value. You must pick only one.".format(variable)
                raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() Variable {0} has both an _delta and an _range value. You must pick only one.".format(variable)
                
            if variable + '_delta' not in cfg and variable + '_range' not in cfg:
                print "Variable {0} has neither an _delta nor an _range value. You must have one of these.".format(variable)
                raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() Variable {0} has neither an _delta nor an _range value. You must have one of these.".format(variable)
            
            if variable + '_delta' in cfg and not isinstance(cfg[variable + '_delta'], (float,int)):
                print "The _delta value for variable {0} is not either a float or an int.".format(variable)
                raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() The _delta value for variable {0} is not either a float or a int.".format(variable)
            if variable + '_range' in cfg and not isinstance(cfg[variable + '_range'], list):
                print "The _range value for variable {0} is not either a list.".format(variable)
                raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() The _range value for variable {0} is not either a list.".format(variable)
            if variable + '_range' in cfg:
                for value in cfg[variable + '_range']:
                    if not isinstance(value, (float,int)):
                        print "A value in the _range for variable {0} is not either a float or an int ({1}).".format(variable, value)
                        raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() value in the _range for variable {0} is not either a float or an int ({1}).".format(variable, value)
         
        elif cfg[variable + '_distribution'] in ['normal','lognormal','gaussian']:
            if variable + '_std_dev' not in cfg:
                print "variable {0} in 'variable_list' does not have a related {0}_std_dev value.".format(variable)
                raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() variable {0} in 'variable_list' does not have a related {0}_std_dev value.".format(variable)
            if not isinstance(cfg[variable + '_std_dev'], (float,int)):
                print "The _std_dev value for variable {0} is not either a float or an int.".format(variable)
                raise Exception, "pitot_bootstrapper.check_pitot_boostrapper_inputs() The _std_dev value for variable {0} is not either a float or a int.".format(variable)

                 
    return cfg
    
def build_pitot_bootstrapper_condition_input_details_dictionary(cfg):
    """Function that builds a dictionary that tells the condition
       builder what tests to run. It will include a list called 'test_names'
       that will be cycled through by the main condition program.
       
    """
    
    import random
    
    test_condition_input_details = {}
    test_condition_input_details['test_names'] = []
    
    test_name = 0 # we will add to this as we go, first test_name will be 1
    
    print '-'*60
    print "Building the test condition details dictionary containing a dictionary for each simulation." 
    
    # for this we just need to randomly mess around with all of the input variables until we have
    # the amount of tests as the 'number_of_test_runs' variable    
    
    for i in range(0, cfg['number_of_test_runs']):
        # change the test name
        test_name += 1
        # store that test name
        test_condition_input_details['test_names'].append(test_name)
        # now go through each variable in the variable list and assign values...
        # for now we just want to randomise nominal or above or below
        # I will use random.randint from -1 to 1
        input_dictionary = {}
        
        for variable in cfg['variable_list']:
            
            if cfg[variable + '_distribution'] == 'uniform':
            # we will either have a random distribution between a range, or with deltas...
                if variable + '_delta' in cfg:
                    range_lower = cfg[variable] - cfg[variable + '_delta']
                    range_upper = cfg[variable] + cfg[variable + '_delta']
                elif variable + '_range' in cfg:
                    range_lower = cfg[variable + '_range'][0]
                    range_upper = cfg[variable + '_range'][1]
                 # now our value will be a random number in between these bounds...   
                input_dictionary[variable] = random.uniform(range_lower, range_upper)
            elif cfg[variable + '_distribution'] == 'normal':
                # our nominal value will be the mean and the input must include a standard deviation value
                mean_value = cfg[variable]
                std_dev_value = cfg[variable + '_std_dev']
                # now our value will be a random normal variate number using these...  
                input_dictionary[variable] = random.normalvariate(mean_value, std_dev_value)
            elif cfg[variable + '_distribution'] == 'lognormal':
                # our nominal value will be the mean and the input must include a standard deviation value
                mean_value = cfg[variable]
                std_dev_value = cfg[variable + '_std_dev']
                # now our value will be a random lognormal variate number using these...  
                input_dictionary[variable] = random.lognormvariate(mean_value, std_dev_value)
            elif cfg[variable + '_distribution'] == 'gaussian':
                # our nominal value will be the mean and the input must include a standard deviation value
                mean_value = cfg[variable]
                std_dev_value = cfg[variable + '_std_dev']
                # now our value will be a random gaussian number using these...  
                input_dictionary[variable] = random.gauss(mean_value, std_dev_value)
                
        test_condition_input_details[test_name] = input_dictionary
                
    print "The test_names list for this simulation is:"
    print test_condition_input_details['test_names']
                
    return test_condition_input_details

def pitot_bootstrapper_test_run(current_cfg, cfg, results):
    """Function that takes the fully built config dictionary
       and the text file that is being used for the program output
       and does the test run then adds a line to the output file.
       
        current_cfg: cfg for pitot for the current run
        cfg: overall cfg for the full simulation set...
        results: result_dict
       
    """
    
    condition_status = True #This will be turned to False if the condition fails
    
    current_cfg['filename'] = cfg['original_filename'] + '-test-{0}-result'.format(cfg['test_number'])
    
    # some code here to make a copy of the stdout printouts for each test and store it
    
    import sys
    
    test_log = open(cfg['original_filename'] + '-test-{0}-log.txt'.format(cfg['test_number']),"w")
    sys.stdout = stream_tee(sys.stdout, test_log)
    
    print '-'*60
    print "Running test {0} of {1}.".format(cfg['test_number'], cfg['number_of_test_runs'])
    
    try:
        current_cfg, states, V, M = run_pitot(cfg = current_cfg)
    except Exception as e:
         print "Error {0}".format(str(e))
         print "Test {0} failed. Result will not be printed to csv output.".format(cfg['test_number'])
         condition_status = False
    if current_cfg['secondary'] and current_cfg['Vsd'] > current_cfg['Vs1']:
        print "Vsd is faster than Vs1, condition cannot be simulated by Pitot properly."
        print "Test {0} is considered failed, and result will not be printed to csv output.".format(cfg['test_number'])
        condition_status = False
    if condition_status:
        # this is to trick the function from pitot_condition_builder to work more generally for now...
        if current_cfg['secondary']:
            current_cfg['secondary_list'] = [True]
        else:
            current_cfg['secondary_list'] = [False]
        results = add_new_result_to_results_dict(current_cfg, states, V, M, results)
    else:
        cfg['last_run_successful'] = False
        results['unsuccessful_runs'].append(cfg['test_number'])
        
    # now we end the stream teeing here by pulling out the original stdout object
    # and overwriting the stream tee with that, then closing the log file
    sys.stdout = sys.stdout.stream1   
    test_log.close()
            
    return condition_status, results  
    
def pitot_bootstrapper_summary(cfg, results):
    """Function that takes the config dictionary and results dictionary 
       made throughout the running of the program and prints a summary of 
       the run to the screen and to a summary text file
    """
    
    from numpy import mean, sqrt, std, array
    
    print '-'*60
    print "Printing summary to screen and to a text document."
    
    
    with open(cfg['original_filename']+'-bootstrapper-summary.txt',"w") as summary_file:
    
        # print lines explaining the results
        summary_line_1 = "# Summary of pitot boostrapper using Version {0} of the program.".format(VERSION_STRING)
        summary_file.write(summary_line_1 + '\n')
        
        summary_line_2 = "{0} tests ran. {1} ({2:.1f}%) were successful."\
            .format(cfg['number_of_test_runs'], len(results['test number']),
            float(len(results['test number']))/float(cfg['number_of_test_runs'])*100.0)
        print summary_line_2
        summary_file.write(summary_line_2 + '\n')
    
        if results['unsuccessful_runs']: 
            summary_line_3 = "Unsucessful runs were run numbers {0}.".format(results['unsuccessful_runs'])
            print summary_line_3
            summary_file.write(summary_line_3 + '\n')          
    
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
                    
                    # get mean and mean error
                    mean_value = mean(data_list)
                    results[variable + '_mean'] = mean_value

                    std_dev = std(data_list, ddof=1)
                    c_i_95 = std_dev * 1.96
                    
                    results[variable + '_std_dev'] = std_dev
                    results[variable + '_95_pc_CI'] = c_i_95
                    
                    if variable[0] == 'p':
                        summary_line = "Variable {0} varies from {1:.2f} - {2:.2f} Pa. Mean is {3:.2f} Pa. Std dev is {4:.2f} Pa. 95pc CI is {5:.2f} Pa."\
                        .format(variable, min_value, max_value, mean_value, std_dev, c_i_95)
                    elif variable[0] == 'V':
                        summary_line = "Variable {0} varies from {1:.2f} - {2:.2f} m/s. Mean is {3:.2f} m/s. Std dev is {4:.2f} m/s. 95pc CI is {5:.2f} m/s."\
                        .format(variable, min_value, max_value, mean_value, std_dev, c_i_95)
                    elif variable[0] == 'T':
                        summary_line = "Variable {0} varies from {1:.2f} - {2:.2f} K. Mean is {3:.2f} K. Std dev is {4:.2f} K. 95pc CI is {5:.2f} K."\
                        .format(variable, min_value, max_value, mean_value, std_dev, c_i_95)
                    elif variable[0] == 'r':
                        summary_line = "Variable {0} varies from {1:.7f} - {2:.7f} kg/m**3. Mean is {3:.7f} kg/m**3. Std dev is {4:.7f} kg/m**3. 95pc CI is {5:.7f} kg/m**3."\
                        .format(variable, min_value, max_value, mean_value, std_dev, c_i_95)
                    elif variable[0] == 'H' or variable[0] == 'h':
                        summary_line = "Variable {0} varies from {1:.7f} - {2:.7f} MJ/kg. Mean is {3:.7f} MJ/kg. Std dev is {4:.7f} MJ/kg. 95pc CI is {5:.7f} MJ/kg."\
                        .format(variable, min_value, max_value, mean_value, std_dev, c_i_95)  
                    elif '%' in variable or 'y' in variable:
                        summary_line = "Variable {0} varies from {1} - {2}. Mean is {3}. Std dev is {4}. 95pc CI is {5}."\
                        .format(variable, min_value, max_value, mean_value, std_dev, c_i_95)
                    else:
                        summary_line = "Variable {0} varies from {1:.1f} - {2:.1f}. Mean is {3:.1f}. Std dev is {4:.1f}. 95pc CI is {5:.1f}."\
                        .format(variable, min_value, max_value, mean_value, std_dev, c_i_95)
                    print summary_line
                    summary_file.write(summary_line + '\n')
                    
        summary_file.close()   
        
    return results

def run_pitot_bootstrapper(cfg = {}, config_file = None, force_restart = None):
    """
    
    Chris James (c.james4@uq.edu.au) 23/12/14
    
    run_pitot_gg_differing_diluent_analysis(dict) - > depends
    
    force_restart can be used to force the simulation to start again instead of 
    looking for an unfinished simulation that may be there...
    
    """
    
    import time, copy
    
    #---------------------- get the inputs set up --------------------------
    
    if config_file:
        cfg = config_loader(config_file)
    
    # set our test gas to custom here and give it a dummy test gas
    # inputUnits and test gas input to pass the input test    
        
    #----------------- check inputs ----------------------------------------
    
    cfg = input_checker(cfg)
    
    cfg = check_pitot_bootstrapper_inputs(cfg)
    
    intermediate_filename = cfg['filename']+'-bootstrapper-intermediate-result-pickle.dat'
    
    # now check if we have attempted an old run before or not....
    if not os.path.isfile(intermediate_filename) or force_restart: 
        # if not, we set up a new one...
        
        # clean up any old files if the user has asked for it
        
        if cfg['cleanup_old_files']:
            # build a list of auxiliary files and run the cleanup_old_files function frim pitot_condition_builder.py
            auxiliary_file_list = ['-bootstrapper-log-and-result-files.zip',
                                   '-bootstrapper-final-result-pickle.dat',
                                   '-bootstrapper-summary.txt',
                                   '-bootstrapper.csv']
            cleanup_old_files(auxiliary_file_list)
                
        # we use this variable to try to speed some things up later on
        cfg['last_run_successful'] = None
        
        cfg['original_filename'] = copy.copy(cfg['filename'])
        
        # I think store the original cfg too
        # do a deepcopy so we don't just have references to internal stuff...
        cfg['cfg_original'] = copy.deepcopy(cfg)
            
        # work out what we need in our results dictionary and make the dictionary
        # this is just the normal pitot_condition_builder one...
        results = build_results_dict(cfg, extra_variable_list = cfg['variable_list'])
        
        # build the dictionary with the details of all the tests we want to run...
        
        test_condition_input_details = build_pitot_bootstrapper_condition_input_details_dictionary(cfg)
        
        # print the start message
           
        print '-'*60
        print 'Running Pitot Bootstrapper Version: {0}.'.format(VERSION_STRING)    
        print "{0} tests will be run.".format(cfg['number_of_test_runs'])
        print "{0} variables will be used in the analysis.".format(len(cfg['variable_list']))
        print "Those variables are:"
        print cfg['variable_list']
        
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
        print "If this is not what you want, please delete the file '{0}' and run the bootstrapper again.".format(intermediate_filename)
    
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
        
        # grab the original cfg
        current_cfg = copy.deepcopy(cfg['cfg_original'])   
        # use update to add any of the changes...
        current_cfg.update(test_condition_input_details[test_name])
        
        current_cfg['test_number'] = test_name
          
        run_status, results = pitot_bootstrapper_test_run(current_cfg, cfg, results) 
        if run_status:
            cfg['good_counter'] += 1
        
        # add this to finished simulations list, regardless of whether it finished correctly or not...
        cfg['finished_simulations'].append(test_name)
        
        # now pickle the intermediate result so we can result the simulation if needed...
        pickle_intermediate_data(cfg, results, test_condition_input_details, filename = intermediate_filename)
    
    # now that we're done we can dump the results to the results csv 
    intro_line = "Output of pitot bootstrapper Version {0}.".format(VERSION_STRING)            
    results_csv_builder(results, test_name = cfg['original_filename'],  
                        intro_line = intro_line, filename = cfg['original_filename'] + '-bootstrapper.csv')
    
    # print the summary
    
    results = pitot_bootstrapper_summary(cfg, results)
                    
    # and pull in the pickle function from pitot_condition_builder.py
    # so we can pick the results and config dictionaries                  
    
    pickle_result_data(cfg, results, filename = cfg['original_filename']+'-bootstrapper-final-result-pickle.dat')
    
    # now delete the intermediate pickle that we made during the simulation...
    print "Removing the final intermediate pickle file."
    
    if os.path.isfile(intermediate_filename): 
        os.remove(intermediate_filename)
    
    # now analyse results dictionary and print some results to the screen
    # and another external file
    
    # need to add a summary!
    #gg_differing_diluent_analysis_summary(cfg, results) 
    
    #zip up the final result using the function from pitot_condition_builder.py
    zip_result_and_log_files(cfg, output_filename = cfg['original_filename'] + '-bootstrapper-log-and-result-files.zip')
    
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
           
    run_pitot_bootstrapper(cfg = {}, config_file = config_file, force_restart = force_restart)
    
    return
    
#----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "pitot_bootstrapper.py - Pitot bootstrapping utility"
        print "start with --help for help with inputs"
        
    else:
        main()

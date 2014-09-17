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

VERSION_STRING = "12-Sep-2014"


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
        
    if 'air_contamination_list' not in cfg:
        print "You have not specified an 'air_contamination_list'. Bailing out."
        cfg['bad_input'] = True
        
    if 'air_contamination_inputUnits' not in cfg:
        print "'air_contamination_inputUnits' not in cfg."
        print "Setting it to 'moles'."
        cfg['air_contamination_inputUnits']
        
    # add an extra check to make sure that the inputUnits match between the contamination
    # and the test gas if the test gas is custom. I think it's pretty stupid
    # and confusing if they aren't the same
    
    if cfg['test_gas'] == 'custom' and cfg['air_contamination_inputUnits'] != cfg['test_gas_inputUnits']:
        print "'air_contamination_inputUnits' = {0} and 'test_gas_inputUnits' = {1}."\
        .format(cfg['air_contamination_inputUnits'], cfg['test_gas_inputUnits'])
        print "This isn't the smartest idea. Bailing out."
        cfg['bad_input'] = True
                
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
    
def contamination_analysis_test_run(cfg, contamination_analysis_output, results):
    """Function that takes the fully built config dictionary
       and the text file that is being used for the program output
       and does the test run then adds a line to the output file.
    """
    
    condition_status = True #This will be turned to False if the condition fails
    
    cfg['filename'] = cfg['original_filename'] + '-test-{0}'.format(cfg['test_number'])
    
    print '-'*60
    print "Running test {0} of {1}.".format(cfg['test_number'], cfg['number_of_test_runs'])
    print "Current level of air contamination is {0} % (by {1})."\
         .format(cfg['contamination_percentage'], cfg['air_contamination_inputUnits'])
    try:
        cfg, states, V, M = run_pitot(cfg = cfg)
    except Exception:
        cfg['state7_no_ions'] = True
        # need to remove Vs values from the dictionary or it will bail out
        # on the next run            
        if cfg['secondary']: cfg.pop('Vsd') 
        cfg.pop('Vs1'); cfg.pop('Vs2')        
        print "Original test failed, trying again with 'state7_no_ions' turned on."
        try:
            cfg, states, V, M = run_pitot(cfg = cfg)
        except Exception:
            # need to remove Vs values from the dictionary or it will bail out
            # on the next run            
            print "Test {0} failed. Result will not be printed to csv output.".format(cfg['test_number'])
            condition_status = False
    if cfg['secondary'] and cfg['Vsd'] > cfg['Vs1']:
        print "Vsd is faster than Vs1, condition cannot be simulated by Pitot properly."
        print "Test {0} is considered failed, and result will not be printed to csv output.".format(cfg['test_number'])
        condition_status = False
    if condition_status:
        string_output = output_builder (cfg, states, V, M)
        contamination_analysis_output.write(string_output + '\n')
        results = add_new_result_to_results_dict(cfg, states, V, M, results)
        # need to remove Vs values from the dictionary or it will bail out
        # on the next run
        cfg.pop('Vsd'); cfg.pop('Vs1'); cfg.pop('Vs2')
    else:
        # need to remove Vs values from the dictionary or it will bail out
        # on the next run         
        if cfg['secondary']: cfg.pop('Vsd') 
        cfg.pop('Vs1'); cfg.pop('Vs2') 
            
    return condition_status, results
    
def output_builder(cfg, states, V, M):
    """Function that takes the four dictionaries from the completed
       Pitot run and builds a string output that can be added to the
       csv output file.
    """
    
    if not cfg['secondary']: #need to fill psd and Vsd with string values
        cfg['psd1'] = 'Not used'
        cfg['Vsd'] = 'N/A'
        
    # Now make the basic string   
    
    basic = "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18}"\
            .format(cfg['test_number'], cfg['contamination_percentage'], cfg['psd1'], 
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
                                   
    # now put it all together
    
    if cfg['nozzle']:
        basic = basic + nozzle
    
    if cfg['conehead'] and not cfg['shock_over_model']:
        string_output = basic + conehead
    elif cfg['shock_over_model'] and not cfg['conehead']:
        string_output = basic + shock_over_model
    elif cfg['shock_over_model'] and cfg['conehead']:
        string_output = basic + conehead + shock_over_model        
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
    
    results['test number'].append(cfg['test_number'])
    results['air contamination'].append(cfg['contamination_percentage'])
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
    results['Ht'].append(cfg['stagnation_enthalpy']/10**6)
    
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
        
    return results
    
def contamination_analysis_summary_builder(cfg, results, contamination_analysis_summary_file):
    """Function that takes the config dictionary and results dictionary 
       made throughout the running of the program and prints a summary of 
       the run to the screen and to a summary text file
    """
    
    print '-'*60
    print "Printing summary to screen and to a text document."
    
    summary_line_1 = "{0} tests ran. {1} ({2:.1f}%) were successful."\
        .format(cfg['number_of_test_runs'], len(results['test number']),
        float(len(results['test number']))/float(cfg['number_of_test_runs'])*100.0)
    print summary_line_1
    contamination_analysis_summary_file.write(summary_line_1 + '\n')        

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
                else:
                    summary_line = "Variable {0} varies from {1:.1f} - {2:.1f}."\
                    .format(variable, min_value, max_value)
                print summary_line
                contamination_analysis_summary_file.write(summary_line + '\n')
                
    return
            
def run_pitot_contamination_analysis(cfg = {}, config_file = None):
    """
    
    Chris James (c.james4@uq.edu.au) 12/09/14
    
    run_pitot_contamination_analysis(dict) - > depends
    
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
    cfg['original_test_gas'] = copy.copy(cfg['test_gas'])
    
    counter = 0
    good_counter = 0
    
    # print how many tests we're going to run, and the ranges.
    
    print '-'*60    
    print "{0} tests will be run.".format(cfg['number_of_test_runs'])
    
    if cfg['air_contamination_inputUnits'] == 'moles':
        print "Air contamination will be tested from {0} - {1} % in increments of {2} as a mole fraction."\
        .format(cfg['air_contamination_list'][0], cfg['air_contamination_list'][-1], 
                cfg['air_contamination_list'][1] - cfg['air_contamination_list'][0])
    elif cfg['air_contamination_inputUnits'] == 'massf':
        print "Air contamination will be tested from {0} - {1} % in increments of {2} as a massf fraction."\
        .format(cfg['air_contamination_list'][0], cfg['air_contamination_list'][-1], 
                cfg['air_contamination_list'][1] - cfg['air_contamination_list'][0])
            
    #open our csv file ready to go
    
    # open a file to start saving results
    contamination_analysis_output = open(cfg['filename']+'-contamination-analysis.csv',"w")  #csv_output file creation
    # print a line explaining the results
    intro_line_1 = "# Output of pitot air contamination program."
    contamination_analysis_output.write(intro_line_1 + '\n')
    basic = "#test number,air contamination percentage ({0}),psd1,p1,p5,Vsd,Vs1,Vs2,p2,T2,rho2,V2,M2,p7,T7,rho7,V7,M7,Ht"\
    .format(cfg['air_contamination_inputUnits'])
    nozzle = ",arearatio,p8,T8,rho8,V8,M8"
    if cfg['nozzle']:
        basic = basic + nozzle
    conehead = ",p10c,T10c,rho10c,V10c"
    shock_over_model = ",p10f,T10f,rho10f,V10f,p10e,T10e,rho10f,V10e"
    if cfg['conehead'] and not cfg['shock_over_model']:
        intro_line_2 = basic + conehead
    elif cfg['shock_over_model'] and not cfg['conehead']:
        intro_line_2 = basic + shock_over_model
    elif cfg['shock_over_model'] and cfg['conehead']:
        intro_line_2 = basic + conehead + shock_over_model        
    else:
        intro_line_2 = basic
    contamination_analysis_output.write(intro_line_2 + '\n')
    
    # then make a dictionary of lists to store results in the Python memory
        
    # need to make a list to create a series of empty lists in the results
    # dictionary to store the data. the list is tailored to the test condition
    basic_list = ['test number','air contamination','psd1','p1','p5','Vsd','Vs1',
                  'Vs2','p2','T2','rho2','V2','M2','p7','T7','rho7','V7','M7','Ht']
    nozzle_list = ['arearatio','p8','T8','rho8','V8','M8']
    if cfg['nozzle']:
        basic_list = basic_list + nozzle_list
    conehead_list = ['p10c','T10c','rho10c','V10c']
    shock_over_model_list = ['p10f','T10f','rho10f','V10f','p10e','T10e','rho10e','V10e']
    if cfg['conehead'] and not cfg['shock_over_model']:
        full_list = basic_list + conehead_list
    elif cfg['shock_over_model'] and not cfg['conehead']:
        full_list = basic_list + shock_over_model_list
    elif cfg['shock_over_model'] and cfg['conehead']:
        full_list = basic_list + conehead_list + shock_over_model_list        
    else:
        full_list = basic_list
    
    results = {title : [] for title in full_list}
    
    # add the list of titles in case we want to use it in future
    
    results['full_list'] = full_list
    
    have_checked_time = False
    
    #now start up the for loops and get running    
    
    for contamination_percentage in cfg['air_contamination_list']:
        # need to extract our original test gas here and then we
        # need to set a custom test gas and then work out how much
        # we need to adjust the original values by to fit in with
        # the amount of contamination
        
        cfg['contamination_percentage'] = contamination_percentage
        percentage_not_air = 100.0 - contamination_percentage
        
        # this is easy if the original test gas was already custom
        if contamination_percentage > 0.0: #obviously don't do anything special when the percentage is 0
            if cfg['original_test_gas'] == 'custom':
                if cfg['test_gas_inputUnits'] == 'moles' and cfg['air_contamination_inputUnits'] == 'moles' or \
                cfg['test_gas_inputUnits'] == 'massf' and cfg['air_contamination_inputUnits'] == 'massf':
                    for species in cfg['test_gas_composition'].keys():
                        cfg['test_gas_composition'][species] = cfg['test_gas_composition'][species]*percentage_not_air
                    # if air was already in the mix, just increase the amount
                    if 'Air' in cfg['test_gas_composition'].keys():
                        cfg['test_gas_composition']['Air'] = cfg['test_gas_composition']['Air'] + contamination_percentage / 100.0
                    # if not, add the air
                    else:
                        cfg['test_gas_composition']['Air'] = contamination_percentage / 100.0              
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
                    cfg['test_gas_composition']['Air'] = original_test_gas.species['Air'] + contamination_percentage 
                # if not, add the air
                else:
                    cfg['test_gas_composition']['Air'] = contamination_percentage   

        counter += 1
        cfg['test_number'] = counter
        if not have_checked_time:
            start_time = time.time()
        run_status, results = contamination_analysis_test_run(cfg, contamination_analysis_output, results) 
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

                        
    contamination_analysis_output.close()
    
    # now analyse results dictionary and print some results to the screen
    # and another external file
    
    contamination_analysis_summary_file = open(cfg['original_filename']+'-contamination-analysis-summary.txt',"w")
    # print a line explaining the results
    summary_line_1 = "# Summary of pitot contamination analysis program output."
    contamination_analysis_summary_file.write(summary_line_1 + '\n')
    summary_line_2 = "# Summary performed using Version {0} of the contamination analysis program.".format(VERSION_STRING)
    contamination_analysis_summary_file.write(summary_line_2 + '\n')
    
    contamination_analysis_summary_builder(cfg, results, contamination_analysis_summary_file)
    
    contamination_analysis_summary_file.close()    
    
    return
                                
#----------------------------------------------------------------------------

def main():
    
    import optparse  
    op = optparse.OptionParser(version=VERSION_STRING)   
    op.add_option('-c', '--config_file', dest='config_file',
                  help=("filename where the configuration file is located"))    

    opt, args = op.parse_args()
    config_file = opt.config_file
           
    run_pitot_contamination_analysis(cfg = {}, config_file = config_file)
    
    return
    
#----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "pitot_contamination_analysis.py - Pitot Equilibrium expansion tube simulator air contamination analysis tool"
        print "start with --help for help with inputs"
        
    else:
        main()

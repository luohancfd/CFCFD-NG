#! /usr/bin/env python
"""
pitot_area_ratio_check.py: pitot area ratio check

This file includes functions that run through multiple specified
area ratios after a normal pitot run is completed.

Chris James (c.james4@uq.edu.au) - 20/08/13 

"""

import sys, copy

from pitot_flow_functions import nozzle_expansion, conehead_calculation, shock_over_model_calculation, calculate_scaling_information

def area_ratio_check_results_dict_builder(cfg):
    """Function that takes the config file and based on the details of the run,
       will build a results dictionary that the area ratio check will use.
    """
    
    #first pull out the nozzle entry state number for the nozzle ratios...
    # (it will be a single digit so just pull out last value...)
    nozzle_entry_number = cfg['nozzle_entry_state'][-1]  
    
    full_list = ['test number','area ratio','p8','T8','rho8','V8','M8',
                 'p8/p{0}'.format(nozzle_entry_number), 'T8/T{0}'.format(nozzle_entry_number),
                 'rho8/rho{0}'.format(nozzle_entry_number), 'V8/V{0}'.format(nozzle_entry_number),
                 'M8/M{0}'.format(nozzle_entry_number),]
                 
    if cfg['calculate_scaling_information']:
        calculate_scaling_information_freestream_list = ['s8_mu', 's8_rhoL', 's8_pL', 's8_Re', 's8_Kn']
        full_list += calculate_scaling_information_freestream_list
    
    if cfg['conehead']:
         conehead_list = ['p10c','T10c','rho10c','V10c']
         full_list += conehead_list
    if cfg['shock_over_model']:
        shock_over_model_list = ['p10f','T10f','rho10f','V10f','p10e','T10e','rho10e','V10e']
        full_list += shock_over_model_list
        if cfg['calculate_scaling_information']:
            calculate_scaling_information_normal_shock_list = ['s10e_mu', 's10e_rhoL', 's10e_pL', 's10e_Re', 's10e_Kn']
            full_list += calculate_scaling_information_normal_shock_list  

    # now populate the dictionary with a bunch of empty lists based on that list

    results = {title : [] for title in full_list}
    
    # add the list of titles in case we want to use it in future
    
    results['full_list'] = full_list
    
    print '-'*60
    print "The full list of variables to be added to the area ratio check output are:"
    print full_list
    
    #add a list where we can store unsuccesful run numbers for analysis
    results['unsuccessful_runs'] = []
      
    return results
    
def area_ratio_add_new_result_to_results_dict(cfg, states, V, M, results):
    """Function that takes a completed test run and adds the tunnel
       configuration and results to the results dictionary.
    """ 
           
    results['test number'].append(cfg['counter'])
        
    results['area ratio'].append(cfg['area_ratio'])
    results['p8'].append(states['s8'].p)
    results['T8'].append(states['s8'].T)
    results['rho8'].append(states['s8'].rho)
    results['V8'].append(V['s8'])
    results['M8'].append(M['s8'])
    
    #first pull out the nozzle entry state number to store the ratios (it will be a single digit so just pull out last value...)
    nozzle_entry_number = cfg['nozzle_entry_state'][-1]  
    
    results['p8/p{0}'.format(nozzle_entry_number)].append(cfg['nozzle_pressure_ratio'])
    results['T8/T{0}'.format(nozzle_entry_number)].append(cfg['nozzle_temperature_ratio'])
    results['rho8/rho{0}'.format(nozzle_entry_number)].append(cfg['nozzle_density_ratio']) 
    results['V8/V{0}'.format(nozzle_entry_number)].append(cfg['nozzle_velocity_ratio'])
    results['M8/M{0}'.format(nozzle_entry_number)].append(cfg['nozzle_mach_number_ratio'])    
    
    if cfg['calculate_scaling_information']:
        results['s8_mu'].append(states[cfg['test_section_state']].mu)        
        results['s8_rhoL'].append(cfg['rho_l_product_freestream'])
        results['s8_pL'].append(cfg['pressure_l_product_freestream'])  
        results['s8_Re'].append(cfg['reynolds_number_freestream']) 
        results['s8_Kn'].append(cfg['knudsen_number_freestream'])  

    if cfg['conehead'] and 's10c' in states:
        results['p10c'].append(states['s10c'].p)
        results['T10c'].append(states['s10c'].T)
        results['rho10c'].append(states['s10c'].rho)
        results['V10c'].append(V['s10c'])
    elif cfg['conehead'] and 's10c' not in states:
        results['p10c'].append('did not solve')
        results['T10c'].append('did not solve')
        results['rho10c'].append('did not solve')
        results['V10c'].append('did not solve')        
        
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
            
        if cfg['calculate_scaling_information']:
            results['s10e_mu'].append(states['s10e'].mu) 
            results['s10e_rhoL'].append(cfg['rho_l_product_state10e'])
            results['s10e_pL'].append(cfg['pressure_l_product_state10e'])  
            results['s10e_Re'].append(cfg['reynolds_number_state10e']) 
            results['s10e_Kn'].append(cfg['knudsen_number_state10e']) 
                             
    return results
    
def area_ratio_check_results_csv_builder(results, test_name = 'pitot_run',  intro_line = None):
    """Function that takes the final results dictionary (which must include a 
       list called 'full_list' that tells this function what to print and in 
       what order) and then outputs a results csv. It will also add an intro line
       if a string with that is added. The name of the test is also required.
    """
    
    # open a file to start saving results
    area_ratio_check_output = open(test_name + '-area-ratio-check.csv',"w")  #csv_output file creation
    
    # print a line explaining the results if the user gives it
    if intro_line:
        intro_line_optional = "# " + intro_line
        area_ratio_check_output.write(intro_line_optional + '\n')
    
    #now we'll make the code build us the second intro line
    intro_line = '#'
    for value in results['full_list']:
        if value != results['full_list'][-1]:
            intro_line += "{0},".format(value)
        else: #don't put the comma if it's the last value
            intro_line += "{0}".format(value)
        
    area_ratio_check_output.write(intro_line + '\n')
    
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
        
        area_ratio_check_output.write(output_line + '\n')  

    area_ratio_check_output.close()              
                                  
    return 
    
def area_ratio_check_normalised_results_csv_builder(results, test_name = 'pitot_run',  
                                   intro_line = None, normalised_by = 'first value',
                                   original_value_result = None):
    """Function that takes the final results dictionary (which must include a 
       list called 'full_list' that tells this function what to print and in 
       what order) and then outputs a normalised version of the results csv.
       You can tell it to normalise by other values, but 'first value' is default.
       
       If the user wants to normalise by the result for the inputted original result
       (the area ratio set in PITOT), they must provide a results dictionary for that 
       original result too...
       
       It will also add an intro line if a string with that is added. 
       The name of the test is also required.
    """
    
    # open a file to start saving results
    area_ratio_check_output = open(test_name + '-area-ratio-check-normalised.csv',"w")  #csv_output file creation
    
    # print a line explaining the results if the user gives it
    if intro_line:
        intro_line_optional = "# " + intro_line
        area_ratio_check_output.write(intro_line_optional + '\n')
        
    normalised_intro_line = "# all variables normalised by {0}".format(normalised_by)
    area_ratio_check_output.write(normalised_intro_line + '\n')
    
    # 'test number' and 'diluent percentage' and the species concentrations
    # will not be normalised
    normalise_exceptions = ['test number', 'area ratio']
    
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
        
    area_ratio_check_output.write(intro_line + '\n')
    
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
        elif normalised_by == 'original value':
            normalising_value_dict[value] = original_value_result[value][-1]
    
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
        
        area_ratio_check_output.write(output_line + '\n')  

    area_ratio_check_output.close()              
                                  
    return     
    
def area_ratio_check_pickle_data(cfg, results, test_name = 'pitot_run'):
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
    print "Pickling cfg and results dictionaries."
    
    pickle_file = open(test_name + '-area-ratio-check-pickle.dat',"w")
    
    cfg_and_results = {'cfg':cfg, 'results':results}
    
    pickle.dump(cfg_and_results, pickle_file)
    pickle_file.close()
   
    return

def area_ratio_check(cfg, states, V, M):
    """Overarching area ratio check function."""
       
    # start by storing old area ratio so it can be retained later

    old_area_ratio = copy.copy(cfg['area_ratio'])
    
    # now build the results dictionary that we will place the results in
    
    results = area_ratio_check_results_dict_builder(cfg)
    
    print 60*"-"    
    print "Performing area ratio check by running through {0} different area ratios."\
    .format(len(cfg['area_ratio_check_list']))
    
    cfg['counter'] = 0 #counter used to tell user how far through the calculations we are
    
    # now do everything...           
    for area_ratio in cfg['area_ratio_check_list']:
        # add the current area ratio
        cfg['counter'] += 1
        cfg['area_ratio'] = area_ratio
        print 60*"-"
        print "Test {0} of {1} (Current area ratio = {2})."\
        .format(cfg['counter'], len(cfg['area_ratio_check_list']), area_ratio)   
        # run the nozzle expansion
        cfg, states, V, M = nozzle_expansion(cfg, states, V, M)
        if cfg['conehead']: #do the conehead calculation if required
            cfg, states, V, M = conehead_calculation(cfg, states, V, M)
        if cfg['shock_over_model']: #do the shock over model calc if required
            cfg, states, V, M = shock_over_model_calculation(cfg, states, V, M)
        if cfg['calculate_scaling_information']:
            cfg = calculate_scaling_information(cfg, states, V, M)
        # print some stuff to the screen
        #I commented all of this out as the code generally does this as it runs through now,
        # but I kept it here incase we ever want to re-add it...
#        print "V8 = {0} m/s, M8 = {1}.".format(V['s8'], M['s8'])
#        print "State 8 (freestream at the nozzle exit):"
#        states['s8'].write_state(sys.stdout)
#        if cfg['conehead'] and cfg['conehead_completed']:
#            print '-'*60
#            print "V10c = {0} m/s.".format(V['s10c'])
#            print "State 10c (surface of a 15 degree conehead):"
#            states['s10c'].write_state(sys.stdout)
#        elif cfg['conehead'] and not cfg['conehead_completed']:
#            print "Conehead calculation failed so result is not being printed."
#        if cfg['shock_over_model']:
#            print '-'*60
#            print "V10f = {0} m/s.".format(V['s10f'])
#            print "State 10f (frozen normal shock over the test model):"
#            states['s10f'].write_state(sys.stdout)    
#            print '-'*60
#            print "V10e = {0} m/s.".format(V['s10e'])
#            print "State 10e (equilibrium normal shock over the test model):"
#            states['s10e'].write_state(sys.stdout)
            
        # now add the result to the results dictionary
        results = area_ratio_add_new_result_to_results_dict(cfg, states, V, M, results)

    #return the original area ratio and values when we leave
    print 60*"-"
    print "Now returning original area ratio and values..."
    cfg['area_ratio'] = old_area_ratio
    # run the nozzle expansion
    cfg, states, V, M = nozzle_expansion(cfg, states, V, M)
    if cfg['conehead']: #do the conehead calculation if required
        cfg, states, V, M = conehead_calculation(cfg, states, V, M)
    if cfg['shock_over_model']: #do the shock over model calc if required
        cfg, states, V, M = shock_over_model_calculation(cfg, states, V, M)

    # now that the run is finished, make the cfg output
    intro_line = "Output of pitot area ratio checking program performed using Pitot version {0}.".format(cfg['VERSION_STRING'])
    area_ratio_check_results_csv_builder(results, test_name = cfg['filename'],  intro_line = intro_line)    
    
    # and the normalised output
    if 'normalise_results_by' in cfg:
        if cfg['normalise_results_by'] == 'original value':
            # we need to make a dummy results dictionary that we put the resulta data
            # of the original run in if we want to normalise by the original value
            original_value_result = area_ratio_check_results_dict_builder(cfg)
            original_value_result = area_ratio_add_new_result_to_results_dict(cfg, states, V, M, original_value_result)
        else:
            original_value_result = None
        area_ratio_check_normalised_results_csv_builder(results, test_name = cfg['filename'],  
                                                        intro_line = intro_line, normalised_by = cfg['normalise_results_by'],
                                                        original_value_result = original_value_result)
            
    else:
        area_ratio_check_normalised_results_csv_builder(results, test_name = cfg['filename'],  
                                       intro_line = intro_line, normalised_by = 'first value')        
    
    # and the pickle output
    area_ratio_check_pickle_data(cfg, results, test_name = cfg['filename'])    
           
    return cfg, states, V, M

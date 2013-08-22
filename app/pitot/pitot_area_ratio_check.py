#! /usr/bin/env python
"""
pitot_area_ratio_check.py: pitot area ratio check

This file includes functions that run through multiple specified
area ratios after a normal pitot run is completed.

Chris James (c.james4@uq.edu.au) - 20/08/13 

"""

import sys

from pitot_flow_functions import nozzle_expansion, conehead_calculation

def area_ratio_check(cfg, states, V, M):
    """Overarching area ratio check function."""
    
    # open a file to start saving our results
    area_ratio_output = open(cfg['filename']+'-area-ratio-check.csv',"w")  #csv_output file creation
    # print a line explaining the results
    intro_line_1 = "# Output of pitot area ratio checking program."
    area_ratio_output.write(intro_line_1 + '\n')
    if cfg['conehead']:
        intro_line_2 = "# area ratio, p8, T8, V8, M8, p10c, T10c, V10c"
    else:
        intro_line_2 = "# area ratio, p8, T8, V8, M8"
    area_ratio_output.write(intro_line_2 + '\n')
    
    # start by storing old area ratio so it can be retained later

    old_area_ratio = cfg['area_ratio']  
    
    print "Performing area ratio check by running through a series of different area ratios."
               
    for area_ratio in cfg['area_ratio_check_list']:
        # add the current area ratio
        cfg['area_ratio'] = area_ratio
        print 60*"-"
        print "Current area ratio = {0}.".format(area_ratio)
        # run the nozzle expansion
        cfg, states, V, M = nozzle_expansion(cfg, states, V, M)
        if cfg['conehead']: #do the conehead calculation if required
            cfg, states, V, M = conehead_calculation(cfg, states, V, M)
        # print some stuff to the screen
        print "V8 = {0} m/s, M8 = {1}.".format(V['s8'], M['s8'])
        print "State 8 (freestream at the nozzle exit):"
        states['s8'].write_state(sys.stdout)
        if cfg['conehead'] and cfg['conehead_completed']:
            print "V10c = {0} m/s.".format(V['s10c'])
            print "State 10c (surface of a 15 degree conehead):"
            states['s10c'].write_state(sys.stdout)
        elif cfg['conehead'] and not cfg['conehead_completed']:
            print "Conehead calculation failed so result is not being printed."
        
        #now add a new line to the output file
        #only prints the line to the csv if the conehead calc completed
        if cfg['conehead'] and cfg['conehead_completed']:        
            new_output_line = "{0},{1},{2},{3},{4},{5},{6},{7}"\
            .format(area_ratio, states['s8'].p, states['s8'].T,V['s8'],\
                    M['s8'], states['s10c'].p, states['s10c'].T,V['s10c'])
        elif not cfg['conehead']:
            new_output_line = "{0},{1},{2},{3},{4}"\
            .format(area_ratio, states['s8'].p, states['s8'].T,V['s8'], M['s8'])
        area_ratio_output.write(new_output_line + '\n')
    
    # close the output file
    area_ratio_output.close()        
    
    #return the original area ratio and values when we leave
    print 60*"-"
    print "Now returning original area ratio and values..."
    cfg['area_ratio'] = old_area_ratio
    # run the nozzle expansion
    cfg, states, V, M = nozzle_expansion(cfg, states, V, M)
    if cfg['conehead']: #do the conehead calculation if required
        cfg, states, V, M = conehead_calculation(cfg, states, V, M)
       
    return cfg, states, V, M
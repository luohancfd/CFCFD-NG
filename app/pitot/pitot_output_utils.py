#! /usr/bin/env python
"""
pitot_output_utils.py: pitot output utilities

This file collects all of the functions that the larger pitot program uses to
display and print its output.

Chris James (c.james4@uq.edu.au) - 07/05/13 

"""

from cfpylib.nm.zero_solvers import secant
# We base our calculation of gas properties upon calls to the NASA Glenn CEA code.
from cfpylib.gasdyn.cea2_gas import Gas, make_gas_from_name
from cfpylib.gasdyn.gas_flow import *
from cfpylib.gasdyn.ideal_gas_flow import p0_p, pitot_p

#to do some perfect gas stuff

import cfpylib.gasdyn.ideal_gas as pg

PRINT_STATUS = 1 #if print status is 1, some basic printouts are done
        
def txt_file_output(cfg, states, V, M):
    """Function that prints the txt output to screen and to a txt file.
    
    """

    txt_output = open(cfg['filename']+'.txt',"w")  #txt_output file creation
                
    if cfg['tunnel_mode'] == 'expansion-tube':
        version_printout = "Pitot Version: {0} doing an expansion tube calculation"\
        .format(cfg['VERSION_STRING'])
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel':
        version_printout = "Pitot Version: {0} doing a non-reflected shock tunnel calculation"\
        .format(cfg['VERSION_STRING'])
    elif cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        version_printout = "Pitot Version: {0} doing a reflected shock tunnel calculation"\
        .format(cfg['VERSION_STRING'])
    print version_printout
    txt_output.write(version_printout + '\n')
    
    if cfg['secondary']:
        description_sd = 'sd1 is secondary driver fill.'
        print description_sd
        txt_output.write(description_sd + '\n')   
    if cfg['tunnel_mode'] == 'expansion-tube':    
        description_1 = 'state 1 is shock tube fill. state 5 is acceleration tube fill.' 
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel' or cfg['tunnel_mode'] == 'reflected-shock-tunnel': 
        description_1 = 'state 1 is shock tube fill.' 
    print description_1
    txt_output.write(description_1 + '\n')
    if cfg['tunnel_mode'] == 'expansion-tube' and cfg['nozzle']:      
        description_2 = 'state 7 is expanded test gas entering the nozzle.'
    elif cfg['tunnel_mode'] == 'expansion-tube' and not cfg['nozzle']:
        description_2 = 'state 7 is expanded test gas entering the test section.'
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel' and cfg['nozzle']:
        description_2 = 'state 2 is shocked test gas entering the nozzle.'
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel' and not cfg['nozzle']:
        description_2 = 'state 2 is shocked test gas entering the test section.'    
    elif cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        description_2 = 'state 2 is the shock tube gas after it has been shocked once.'          
    print description_2
    txt_output.write(description_2 + '\n')
    
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        state5_description = 'state 5 is the stagnant shock tube gas after being processed by the reflected shock.'
        print state5_description
        txt_output.write(state5_description + '\n')
    
    if cfg['nozzle']:    
        description_3 = 'state 8 is test gas exiting the nozzle (using area ratio of {0}).'.format(cfg['area_ratio'])
        print description_3
        txt_output.write(description_3 + '\n')
    
    if cfg['shock_over_model']:    
        description_4 = 'state 10f is frozen shocked test gas flowing over the model.'
        print description_4
        txt_output.write(description_4 + '\n')  
    
        description_5 = 'state 10e is equilibrium shocked test gas flowing over the model.'
        print description_5
        txt_output.write(description_5 + '\n')
        
    if cfg['conehead']:
        description_6 = 'state 10c is conditions over {0} degree conehead in the test section.'.format(cfg['conehead_angle'])
        print description_6
        txt_output.write(description_6 + '\n')
        
    if cfg['wedge']:
        description_7 = 'state 10w is conditions over {0} degree wedge in the test section.'.format(cfg['wedge_angle'])
        print description_7
        txt_output.write(description_7 + '\n')
        
        
    if cfg['solver'] == 'eq':
        solver_printout = "Solver used is equilibrium."
    elif cfg['solver'] == 'pg':
        solver_printout = "Solver used is perfect gas."
    elif cfg['solver'] == 'pg-eq':
        solver_printout = "Solver used is pg-eq (state 1 is set as pg.)"
    print solver_printout
    txt_output.write(solver_printout + '\n')        
        
    facility_used = 'Facility is {0}.'.format(cfg['facility'])        
    print facility_used
    txt_output.write(facility_used + '\n')
    if cfg['solver'] == 'eq':
        test_gas_used = 'Test gas (state 1) is {0} (gamma = {1}, R = {2}, {3}).'.format(cfg['test_gas'],states['s1'].gam,states['s1'].R,states['s1'].reactants)
    elif cfg['solver'] == 'pg' or cfg['solver'] == 'pg-eq':
        test_gas_used = 'Test gas (state 1) is {0} (gamma = {1}, R = {2}).'.format(cfg['test_gas'],states['s1'].gam,states['s1'].R)
    print test_gas_used
    txt_output.write(test_gas_used + '\n')  
    if cfg['solver'] == 'eq':
        driver_gas_used = 'Driver gas is {0}.'.format(states['s4'].reactants)   
    else:
        if cfg['facility'] != 'custom':
            driver_gas_used = 'Driver gas is {0}.'.format(cfg['driver_gas'])
        else:
            driver_gas_used = 'Driver gas is {0}.'.format(cfg['driver_composition'])
    print driver_gas_used
    txt_output.write(driver_gas_used + '\n') 
            
    if cfg['shock_switch']:
        shock_warning1 = "NOTE: a reflected shock was done into the shock tube."
        print shock_warning1
        txt_output.write(shock_warning1 + '\n')
        
    if cfg['secondary'] and not cfg['shock_switch']:
        secondary_shockspeeds = "Vsd = {0:.2f} m/s, Msd1 = {1:.2f}".format(cfg['Vsd'],cfg['Msd1'])
    elif cfg['secondary'] and cfg['shock_switch']:
        secondary_shockspeeds = "Vsd = {0:.2f} m/s, Msd1 = {1:.2f}, Vr = {2:.2f}, Mr = {3:.2f}".\
        format(cfg['Vsd'],cfg['Msd1'], cfg['Vr'], cfg['Mr'])        
        print secondary_shockspeeds
        txt_output.write(secondary_shockspeeds + '\n')
    if cfg['tunnel_mode'] == 'expansion-tube':    
        shockspeeds = "Vs1 = {0:.2f} m/s, Ms1 = {1:.2f}, Vs2 = {2:.2f} m/s, Ms2 = {3:.2f}".\
        format(cfg['Vs1'],cfg['Ms1'],cfg['Vs2'],cfg['Ms2']) 
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel':
        shockspeeds = "Vs1 = {0:.2f} m/s, Ms1 = {1:.2f}".format(cfg['Vs1'],cfg['Ms1'])
    elif cfg['tunnel_mode'] == 'reflected-shock-tunnel':    
        shockspeeds = "Vs1 = {0:.2f} m/s, Ms1 = {1:.2f}, Vr = {2:.2f} m/s, Mr = {3:.2f}".\
        format(cfg['Vs1'],cfg['Ms1'],cfg['Vr'],cfg['Mr'])     
    print shockspeeds #prints above line in console
    txt_output.write(shockspeeds + '\n') #writes above line to txt_output file (input to write command must be a string)
                
    key = "{0:6}{1:11}{2:9}{3:6}{4:9}{5:6}{6:9}{7:8}{8:9}".format("state","P","T","a","V","M","rho","pitot","stgn")
    print key
    txt_output.write(key + '\n')
    
    units = "{0:6}{1:11}{2:9}{3:6}{4:9}{5:6}{6:9}{7:9}{8:9}".format("","Pa","K","m/s","m/s","","m^3/kg","kPa","MPa")
    print units
    txt_output.write(units + '\n')
    
    #new dictionaries here to add pitot and stagnation pressure calcs
    
    pitot = {} #pitot pressure dict
    p0 = {} #stagnation pressure dict
    
    def condition_printer(it_string):
        """Prints the values of a specified condition to the screen and to 
        the txt_output file. 
        
        I made a function of this so I didn't have to keep pasting the code in."""
        
        if states.has_key(it_string):
            
            if M[it_string] == 0:
                pitot[it_string] = states[it_string].p/1000.0
                p0[it_string] = states[it_string].p/1.0e6
            else:
                pitot[it_string] = pitot_p(states[it_string].p,M[it_string],states[it_string].gam)/1000.0
                #make total condition of relevant state for printing
                total_state = total_condition(states[it_string], V[it_string])
                p0[it_string] = total_state.p/1.0e6
            
            if states[it_string].p < 1.0e6: #change how the pressure is printed if it's too big, it keeps ruining the printouts!
                conditions = "{0:<6}{1:<11.7}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.5f}{7:<7.0f}{8:<9.1f}"\
                .format(it_string, states[it_string].p, states[it_string].T,
                        states[it_string].a,V[it_string],M[it_string],
                        states[it_string].rho, pitot[it_string], p0[it_string])
            else:
                conditions = "{0:<6}{1:<11.3e}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.5f}{7:<7.0f}{8:<9.1f}"\
                .format(it_string, states[it_string].p, states[it_string].T,
                        states[it_string].a,V[it_string],M[it_string],
                        states[it_string].rho, pitot[it_string], p0[it_string])
                    
            print conditions
            txt_output.write(conditions + '\n')

    #print the driver related stuff first
    
    condition_printer('s4')
    condition_printer('s3s')
    
    if cfg['secondary']: #need a separate printing thing for the secondary driver
        
        for i in range(1,4): #will do 1 - 3
        
            it_string = 'sd{0}'.format(i)
            condition_printer(it_string)
                    
    for i in range(1,4): #shock tube stuff
        
        it_string = 's{0}'.format(i)
        condition_printer(it_string)
    if cfg['tunnel_mode'] == 'expansion-tube':    
        for i in range(5,8): #acc tube extra states
            it_string = 's{0}'.format(i)
            try:
                condition_printer(it_string)
            except:
                #try again with ions off.
                states[it_string].with_ions = False
                condition_printer(it_string)
                states[it_string].with_ions = True
                
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
         condition_printer('s5')
            
    if cfg['nozzle']: #do nozzle calculations if asked to
        condition_printer('s8')
        
    #do the conditions over the model if asked
    if cfg['shock_over_model']:
        condition_printer('s10f')
        condition_printer('s10e')
            
    if cfg['conehead']:
        condition_printer('s10c')
        
    if cfg['wedge']:
        condition_printer('s10w')
                                           
    #some other useful calculations at the end
          
    states['test_section_total'] = total_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
    states['test_section_pitot'] = pitot_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
    
    cfg['stagnation_enthalpy'] = states['test_section_total'].h #J/kg
    if cfg['nozzle']:        
        stag_enth = 'The total enthalpy (Ht) leaving the nozzle is {0:<.5g} MJ/kg.'\
        .format(cfg['stagnation_enthalpy']/10**6)
    elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'expansion-tube':
        stag_enth = 'The total enthalpy (Ht) at the end of the acceleration tube (state 7) is {0:<.5g} MJ/kg.'\
        .format(cfg['stagnation_enthalpy']/10**6)
    elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'nr-shock-tunnel':
        stag_enth = 'The total enthalpy (Ht) at the end of the shock tube (state 2) is {0:<.5g} MJ/kg.'\
        .format(cfg['stagnation_enthalpy']/10**6)
    elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        stag_enth = 'The total enthalpy (Ht) in the stagnated region (state 5) is {0:<.5g} MJ/kg.'\
        .format(cfg['stagnation_enthalpy']/10**6)
    print stag_enth
    txt_output.write(stag_enth + '\n')
    
    #calculate flight equivalent velocity
    #for a description of why this is, refer to Bianca Capra's thesis page 104 - 105
    #Capra, B., Aerothermodynamic Simulation of Subscale Models of the FIRE II and
    #Titan Explorer Vehicles in Expansion Tubes, Ph.D. thesis, the University of Queens-
    #land, St. Lucia, Australia, 2006.
    cfg['u_eq'] = math.sqrt(2.0*cfg['stagnation_enthalpy']) 
    u_eq_print = 'The flight equivalent velocity (Ue) is {0:<.5g} m/s.'.format(cfg['u_eq'])
    print u_eq_print
    txt_output.write(u_eq_print + '\n')
    
    #if the test time calculation has been done, print it
    if cfg['calculate_test_time']: 
        basic_test_time_printout = 'Basic test time = {0:.2f} microseconds'.format(cfg['t_test_basic']*1.0e6)
        print  basic_test_time_printout
        txt_output.write(basic_test_time_printout + '\n')

    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        rst_1 = "At state 5 gam = {0}, R = {1}.".format(states['s5'].gam, states['s5'].R)
        print rst_1
        txt_output.write(rst_1 + '\n')
        
        rst_2 = "Species in state 5 at equilibrium:"
        print rst_2
        txt_output.write(rst_2 + '\n')
        
        rst_3 = '{0}'.format(states['s5'].species)
        print rst_3
        txt_output.write(rst_3 + '\n')        
    
    #added ability to get the species in the post-shock condition
    
    if cfg['shock_over_model']:
        species1 = 'species in the shock layer at equilibrium:'        
        print species1
        txt_output.write(species1 + '\n')
        
        species2 = '{0}'.format(states['s10e'].species)
        print species2
        txt_output.write(species2 + '\n')
        
    if cfg['mode'] == 'cea-printout' or cfg['mode'] == 'cea-txt-printout':
                        
        cea_printout_intro = "Printing gas state printouts for certain conditions..."
        print cea_printout_intro
        txt_output.write(cea_printout_intro + '\n')
        
        def gas_condition_printer(state, title, solver):
            """This is a function designed to mimic the gas printers
                from the gas object. I needed to be able to output the strings to 
                a text document as well as to the screen, so I just wrote my own
                function based on the other one (as it needed to do eq and pg objects).
                
                Chris James (c.james4@uq.edu.au) - 29-Jan-2013."""
                
            #make the strings we need
            
            if solver == 'eq' or solver == 'pg-eq':    
                line_one = '    p: {0:g} Pa, T: {1:g} K, rho: {2:g} kg/m**3, e: {3:g} J/kg, h: {4:g} J/kg, a: {5:g} m/s, s:{6:g} kJ/(kg.K)'\
                .format(state.p, state.T, state.rho, state.u, state.h, state.a, state.s)
                line_two = '    R: {0:g} J/(kg.K), gam: {1:g}, Cp: {2:g} J/(kg.K), mu: {3:g} Pa.s, k: {4:g} W/(m.K)'\
                .format(state.R, state.gam, state.cp, state.mu, state.k)
                line_three = '    species {0:s}: {1:s}'.format(state.outputUnits, str(state.species))
            elif solver == 'pg':
                line_one = '    p: {0:g} Pa, T: {1:g} K, rho: {2:g} kg/m**3, e: {3:g} J/kg, h: {4:g} J/kg, a: {5:g} m/s, s:{6:g} kJ/(kg.K)'\
                .format(state.p, state.T, state.rho, state.u, state.h, state.a, state.s)
                line_two = '    R: {0:g} J/(kg.K), gam: {1:g}, Cp: {2:g} J/(kg.K), mu: {3:g} Pa.s, k: {4:g} W/(m.K)'\
                .format(state.R, state.gam, state.C_p, state.mu, state.k)
                line_three = '    name: {0:s}'.format(state.name)
                
            #then start printing to screen /storing data
            
            intro_line = title
            print intro_line
            txt_output.write(intro_line + '\n')
            print line_one
            txt_output.write(line_one + '\n')
            print line_two
            txt_output.write(line_two + '\n')
            print line_three
            txt_output.write(line_three + '\n')
            
            return
                            
        #test section state
        gas_condition_printer(states[cfg['test_section_state']], 'Test section state ({0}):'\
        .format(cfg['test_section_state']), cfg['solver'])
        
        #test section total condition
        gas_condition_printer(states['test_section_total'],\
        'Test section state ({0}) total condition:'.format(cfg['test_section_state']), \
        cfg['solver'])
        
        #test section pitot condition
        gas_condition_printer(states['test_section_pitot'], \
        'Test section state ({0}) pitot condition:'.format(cfg['test_section_state']), \
        cfg['solver'])        
      
        if cfg['conehead']:
            #conehead condition in test section
            gas_condition_printer(states['s10c'], \
            "Conditions behind {0} degree conehead ('s10c'):".format(cfg['conehead_angle']), \
            cfg['solver'])
                           
        if cfg['shock_over_model']:
            #frozen normal shock in the test section
            gas_condition_printer(states['s10f'], \
            "Conditions behind frozen normal shock over test model ('s10f')", \
            cfg['solver'])
            #equilibrium normal shock in test condition
            gas_condition_printer(states['s10e'], \
            "Conditions behind equilibrium normal shock over test model ('s10e')", \
            cfg['solver'])
                   
    txt_output.close()
    
    return cfg, states, V, M
    
#----------------------------------------------------------------------------
    
def csv_file_output(cfg, states, V, M):
    """Function to do the csv prinouts for pitot."""
    
    csv_output = open(cfg['filename']+'.csv',"w")  #csv_output file creation
    
    if cfg['tunnel_mode'] == 'expansion-tube':
        csv_version_printout = "Pitot Version,{0},expansion-tube mode".format(cfg['VERSION_STRING'])
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel':
        csv_version_printout = "Pitot Version,{0},nr-shock-tunnel mode".format(cfg['VERSION_STRING'])
    elif cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        csv_version_printout = "Pitot Version,{0},reflected-shock-tunnel mode".format(cfg['VERSION_STRING']) 
    csv_output.write(csv_version_printout + '\n')
    
    if cfg['solver'] == 'eq':
        csv_solver_printout = "Solver,equilibrium."
    elif cfg['solver'] == 'pg':
        csv_solver_printout = "Solver,perfect gas"
    elif cfg['solver'] == 'pg-eq':
        csv_solver_printout = "Solver,pg eq"
    csv_output.write(csv_solver_printout + '\n')     
        
    csv_facility_used = 'Facility,{0}.'.format(cfg['facility'])        
    csv_output.write(csv_facility_used + '\n')
    
    csv_test_gas_used = 'Test gas (state 1),{0},gamma,{1},R,{2}'.format(cfg['test_gas'],states['s1'].gam,states['s1'].R)
    csv_output.write(csv_test_gas_used + '\n')  
    
    if cfg['solver'] == 'eq':
        csv_driver_gas_used = 'Driver gas,{0}.'.format(states['s4'].reactants)
    else:
        if cfg['facility'] != 'custom':
            csv_driver_gas_used = 'Driver gas,{0}.'.format(cfg['driver_gas'])
        else:
            csv_driver_gas_used = 'Driver gas,{0}.'.format(cfg['driver_composition'])
    csv_output.write(csv_driver_gas_used + '\n') 
            
    if cfg['secondary']:
        csv_secondary_shockspeeds = "Vsd,{0:.2f} m/s,Msd1,{1:.2f}".format(cfg['Vsd'],cfg['Msd1'])
        csv_output.write(csv_secondary_shockspeeds + '\n')
    if cfg['tunnel_mode'] == 'expansion-tube':
        csv_shockspeeds = "Vs1,{0:.2f} m/s,Ms1,{1:.2f},Vs2,{2:.2f} m/s,Ms2,{3:.2f}"\
        .format(cfg['Vs1'],cfg['Ms1'],cfg['Vs2'],cfg['Ms2']) 
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel':
        csv_shockspeeds = "Vs1,{0:.2f} m/s,Ms1,{1:.2f}".format(cfg['Vs1'],cfg['Ms1'])
    elif cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        csv_shockspeeds = "Vs1,{0:.2f} m/s,Ms1,{1:.2f},Vr,{2:.2f} m/s,Mr,{3:.2f}"\
        .format(cfg['Vs1'],cfg['Ms1'],cfg['Vr'],cfg['Mr'])         
    csv_output.write(csv_shockspeeds + '\n')
                
    csv_key = "{0:6},{1:11},{2:9},{3:6},{4:9},{5:6},{6:9},{7:8},{8:9}".format("state","P","T","a","V","M","rho","pitot","stgn")
    csv_output.write(csv_key + '\n')
    
    csv_units = "{0:6},{1:11},{2:9},{3:6},{4:9},{5:6},{6:9},{7:9},{8:9}".format("","Pa","K","m/s","m/s","","m^3/kg","kPa","MPa")
    csv_output.write(csv_units + '\n')
    
    #new dictionaries here to add pitot and stagnation pressure calcs
    
    pitot = {} #pitot pressure dict
    p0 = {} #stagnation pressure dict
           
    def csv_condition_printer(it_string):
        """Prints the values of a specified condition to the screen and to 
        the txt_output file. 
        
        I made a function of this so I didn't have to keep pasting the code in."""
        
        if states.has_key(it_string):
                
            if M[it_string] == 0:
                pitot[it_string] = states[it_string].p/1000.0
                p0[it_string] = states[it_string].p/1.0e6
            else:
                pitot[it_string] = pitot_p(states[it_string].p,M[it_string],states[it_string].gam)/1000.0
                p0[it_string] = p0_p(M[it_string], states[it_string].gam)*states[it_string].p/1.0e6
            
            csv_conditions = "{0:<6},{1:<11.7},{2:<9.1f},{3:<6.0f},{4:<9.1f},{5:<6.2f},{6:<9.4f},{7:<8.0f},{8:<9.1f}"\
            .format(it_string, states[it_string].p, states[it_string].T,
                    states[it_string].a,V[it_string],M[it_string],
                    states[it_string].rho, pitot[it_string], p0[it_string])

            csv_output.write(csv_conditions + '\n')

    #print the driver related stuff first
    csv_condition_printer('s4')
    csv_condition_printer('s3s')
    
    if cfg['secondary']: #need a separate printing thing for the secondary driver
        for i in range(1,4): #will do 1 - 3
            it_string = 'sd{0}'.format(i)
            csv_condition_printer(it_string)
                    
    for i in range(1,4): #shock tube stuff
        it_string = 's{0}'.format(i)
        csv_condition_printer(it_string)
    if cfg['tunnel_mode'] == 'expansion-tube':    
        for i in range(5,9): #acc tube and nozzle if it's there
            it_string = 's{0}'.format(i)
            csv_condition_printer(it_string)
            
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        csv_condition_printer('s5')
        
    #do the conditions over the model
    if cfg['shock_over_model']:
        csv_condition_printer('s10f')
        csv_condition_printer('s10e')
            
    if cfg['conehead']:
        csv_condition_printer('s10c')
        
    if cfg['wedge']:
        csv_condition_printer('s10w')
        
    if 'stagnation_enthalpy' not in cfg:
        #if stagnation enthalpy and u_eq not already calculated, calculate them here
        states['test_section_total'] = total_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
        states['test_section_pitot'] = pitot_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
        cfg['stagnation_enthalpy'] = states['test_section_total'].h #J/kg
        #calculate flight equivalent velocity
        #for a description of why this is, refer to Bianca Capra's thesis page 104 - 105
        #Capra, B., Aerothermodynamic Simulation of Subscale Models of the FIRE II and
        #Titan Explorer Vehicles in Expansion Tubes, Ph.D. thesis, the University of Queens-
        #land, St. Lucia, Australia, 2006.
        cfg['u_eq'] = math.sqrt(2.0*cfg['stagnation_enthalpy']) 
                                     
    csv_stag_enth = 'Ht,{0:<.5g} MJ/kg.'.format(cfg['stagnation_enthalpy']/10**6)
    csv_output.write(csv_stag_enth + '\n')
    
    csv_u_eq_print = 'Ue,{0:<.5g} m/s.'.format(cfg['u_eq'])
    csv_output.write(csv_u_eq_print + '\n')

    if cfg['calculate_test_time']: 
        csv_basic_test_time_printout = 'Basic test time,{0:.2f} microseconds'.format(cfg['t_test_basic']*1.0e6)
        csv_output.write(csv_basic_test_time_printout + '\n')   
    
    csv_output.close()
    
    return cfg, states, V, M
    
#----------------------------------------------------------------------------
    
def cleanup_function():
    """Function to clean up temp files if we want to."""
    
    import os
    
    if PRINT_STATUS: 
        print " "
        print "Removing temporary files and leaving the program."
    if os.path.isfile('thermo.inp'): os.remove('thermo.inp')
    if os.path.isfile('thermo.out'): os.remove('thermo.out')
    if os.path.isfile('thermo.lib'): os.remove('thermo.lib')
    if os.path.isfile('tmp.inp'): os.remove('tmp.inp')
    if os.path.isfile('tmp.out'): os.remove('tmp.out')
    if os.path.isfile('trans.inp'): os.remove('trans.inp')
    if os.path.isfile('trans.out'): os.remove('trans.out')
    if os.path.isfile('trans.lib'): os.remove('trans.lib')
    
    return

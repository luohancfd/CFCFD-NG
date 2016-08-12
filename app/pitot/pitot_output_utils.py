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
    
    if 's4i' in states:
        description_driver_fill = 'state 4i is the free piston driver fill condition.'
        print description_driver_fill
        txt_output.write(description_driver_fill + '\n')  
        
    description_driver = 'state 4 is the driver condition.'
    print description_driver
    txt_output.write(description_driver + '\n')  
     
    if cfg['secondary']:
        description_sd = 'state sd1 is secondary driver fill.'
        print description_sd
        txt_output.write(description_sd + '\n')   
    if cfg['tunnel_mode'] == 'expansion-tube':    
        description_1 = 'state 1 is shock tube fill. state 5 is acceleration tube fill.' 
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel' or cfg['tunnel_mode'] == 'reflected-shock-tunnel': 
        description_1 = 'state 1 is shock tube fill.' 
    print description_1
    txt_output.write(description_1 + '\n')   
    if cfg['tunnel_mode'] == 'expansion-tube':
        state2_description = "state 2 is the shocked test gas."
        print state2_description
        txt_output.write(state2_description + '\n')   
        
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
        
        state5_description = 'state 5_d is the stagnant driver gas after being processed by the reflected shock.'
        print state5_description
        txt_output.write(state5_description + '\n')
        
        state6_description = 'state 6 is the test gas condition at the throat (M = 1) of the de Lavel nozzle.'
        print state6_description
        txt_output.write(state6_description + '\n')
    
    if cfg['nozzle']:    
        description_3 = 'state 8 is test gas exiting the nozzle (using area ratio of {0}).'.format(cfg['area_ratio'])
        print description_3
        txt_output.write(description_3 + '\n')
        if cfg['frozen_nozzle_expansion']:
            frozen_nozzle_description = 'state 8 was found using a frozen nozzle expansion.'
            print frozen_nozzle_description
            txt_output.write(frozen_nozzle_description + '\n')            
    
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

    test_done = "Test is '{0}'".format(cfg['test'])
    print test_done
    txt_output.write(test_done + '\n')    
        
    facility_used = 'Facility is {0}.'.format(cfg['facility'])        
    print facility_used
    txt_output.write(facility_used + '\n')
    
    if cfg['custom_secondary_driver_gas']:
        if isinstance(states['sd1'], Gas):
            secondary_driver_gas_used = 'Custom secondary driver gas (state sd1) was used (gamma = {0}, R = {1}, {2} by {3}).'.format(states['sd1'].gam,states['sd1'].R,
                                                                                                                          states['sd1'].reactants, states['sd1'].outputUnits)
        elif isinstance(states['sd1'], pg.Gas):
            secondary_driver_gas_used = 'Custom secondary gas (state sd1) was used (gamma = {0}, R = {1}).'.format(states['sd1'].gam,states['sd1'].R)
    else:
        secondary_driver_gas_used = "Secondary gas (state sd1) is pure He."
        
    print secondary_driver_gas_used
    txt_output.write(secondary_driver_gas_used + '\n')    
    
    if cfg['solver'] == 'eq':
        test_gas_used = 'Test gas (state 1) is {0} (gamma = {1}, R = {2}, {3} by {4}).'.format(cfg['test_gas'],states['s1'].gam,states['s1'].R,
                                                                                        states['s1'].reactants, states['s1'].outputUnits)
    elif cfg['solver'] == 'pg' or cfg['solver'] == 'pg-eq':
        test_gas_used = 'Test gas (state 1) is {0} (gamma = {1}, R = {2}).'.format(cfg['test_gas'],states['s1'].gam,states['s1'].R)
    print test_gas_used
    txt_output.write(test_gas_used + '\n')  
    if cfg['custom_accelerator_gas']:
        if cfg['solver'] == 'eq':
            accelerator_gas_used = 'Custom Accelerator gas (state 5) was used (gamma = {0}, R = {1}, {2} by {3}).'.format(states['s5'].gam,states['s5'].R,
                                                                                                                          states['s5'].reactants, states['s5'].outputUnits)
        elif cfg['solver'] == 'pg' or cfg['solver'] == 'pg-eq':
            accelerator_gas_used = 'Custom Accelerator gas (state 5) was used (gamma = {0}, R = {1}).'.format(states['s5'].gam,states['s5'].R)
    else:
        accelerator_gas_used = "Accelerator gas (state 5) is Air."
        
    print accelerator_gas_used
    txt_output.write(accelerator_gas_used + '\n')
    
    if cfg['solver'] == 'eq':
        if cfg['facility'] != 'custom' and cfg['piston'] in ['Sangdi-1.8MPa', 'sangdi-1.8MPa','Sangdi-2.2MPa', 'sangdi-2.2MPa']:
            driver_gas_used = 'Driver gas is {0}.'.format(cfg['driver_gas'])  
        elif 'driver_pg_gam' not in cfg: # driver will be pg if it is!
            driver_gas_used = 'Driver gas is {0} (by {1}).'.format(states['s4'].reactants, states['s4'].outputUnits)  
        else:
            driver_gas_used = "Driver gas is a perfect gas with gam = {0} and R = {1}."\
            .format(cfg['driver_pg_gam'], cfg['driver_pg_R'])            
    else:
        if cfg['facility'] != 'custom':
            driver_gas_used = 'Driver gas is {0}.'.format(cfg['driver_gas'])
        else:
            if 'driver_composition' in cfg and cfg['solver'] in ['eq', 'pg-eq']:
                driver_gas_used = 'Driver gas is {0} (by {1}).'.format(cfg['driver_composition'], states['s4'].inputUnits)
            elif 'driver_composition' in cfg and cfg['solver'] in ['pg']:
                driver_gas_used = 'Driver gas is {0} (by {1}).'.format(cfg['driver_composition'], cfg['driver_inputUnits'])
            else:
                driver_gas_used = "Driver gas is a perfect gas with gam = {0} and R = {1}."\
                .format(cfg['driver_pg_gam'], cfg['driver_pg_R'])
    print driver_gas_used
    txt_output.write(driver_gas_used + '\n') 
            
    if cfg['shock_switch']:
        shock_warning_1 = "NOTE: a reflected shock was done into the shock tube."
        print shock_warning_1
        txt_output.write(shock_warning_1 + '\n')
    if cfg['rs_out_of_sd']:
        shock_warning_2 = "NOTE: a user specified reflected shock was done at the end of the secondary driver."
        print shock_warning_2
        txt_output.write(shock_warning_2 + '\n')  
    if cfg['rs_out_of_st']:
        shock_warning_3 = "NOTE: a user specified reflected shock was done at the end of the shock tube."
        print shock_warning_3
        txt_output.write(shock_warning_3 + '\n')   
        
    if cfg['secondary'] and not cfg['shock_switch']:
        secondary_shockspeeds = "Vsd = {0:.2f} m/s, Msd1 = {1:.2f}".format(cfg['Vsd'],cfg['Msd1'])
        print secondary_shockspeeds
        txt_output.write(secondary_shockspeeds + '\n')
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
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':    
        shockspeeds_2 = "Vr_d = {0:.2f} m/s, Mr_d = {1:.2f}".format(cfg['Vr_d'],cfg['Mr_d'])
        print shockspeeds_2
        txt_output.write(shockspeeds_2 + '\n')
    if cfg['rs_out_of_sd']:
        rs_out_of_sd = "Vr-sd = {0:.2f} m/s, Mr-sd = {1:.2f}".format(cfg['Vr-sd'],cfg['Mr-sd'])
        print rs_out_of_sd #prints above line in console
        txt_output.write(rs_out_of_sd + '\n') #writes above line to txt_output file (input to write command must be a string)            
    if cfg['rs_out_of_st']:
        rs_out_of_st = "Vr-st = {0:.2f} m/s, Mr-st = {1:.2f}".format(cfg['Vr-st'],cfg['Mr-st'])
        print rs_out_of_st #prints above line in console
        txt_output.write(rs_out_of_st + '\n') #writes above line to txt_output file (input to write command must be a string)        
                
    key = "{0:6}{1:11}{2:9}{3:6}{4:9}{5:6}{6:9}{7:8}{8:9}".format("state","P","T","a","V","M","rho","pitot","stgn")
    print key
    txt_output.write(key + '\n')
    
    units = "{0:6}{1:11}{2:9}{3:6}{4:9}{5:6}{6:9}{7:8}{8:9}".format("","Pa","K","m/s","m/s","","kg/m^3","kPa","MPa")
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
            
            # I decided to add a clone here to stop messing with the original object here.
            try:
                state = states[it_string].clone()
                clone_successful = True
            except Exception as e:
                print "Error {0}".format(str(e))   
                print "Failed to clone state {0} while trying to final total condition."
                if isinstance(states[it_string], Gas) and states[it_string].with_ions:
                    print "Trying again with ions turned off."
                    try:
                        states[it_string].with_ions = False
                        state = states[it_string].clone()
                        clone_successful = True
                    except Exception as e:
                         print "Error {0}".format(str(e))   
                         print "Failed to clone state {0} again while trying to final total condition." 
                         print "0.0 will be added as pitot and total pressures for this state."
                         pitot[it_string] = 0.0
                         p0[it_string] = 0.0
                         clone_successful = False
                else:
                     print "0.0 will be added as pitot and total pressures for this state."
                     pitot[it_string] = 0.0
                     p0[it_string] = 0.0
                     clone_successful = False
                                      
            if M[it_string] == 0:
                pitot[it_string] = state.p/1000.0
                p0[it_string] = state.p/1.0e6
            elif clone_successful:
                try:
                    pitot[it_string] = pitot_p(state.p,M[it_string],state.gam)/1000.0
                except:
                    try:
                        #try again with no ions turned on if it bails out
                        state.with_ions = False
                        pitot[it_string] = pitot_p(state.p,M[it_string],state.gam)/1000.0
                        state.with_ions = True
                    except:
                        # if this fails just give up and print 0.0
                        print "Failed to find pitot pressure for {0} will add 0.0 to print out.".format(it_string)
                        pitot[it_string] = 0.0
                try:
                    #make total condition of relevant state for printing
                    total_state = total_condition(state, V[it_string])
                    p0[it_string] = total_state.p/1.0e6
                except:
                    try:
                        #try again with no ions turned on if it bails out
                        state.with_ions = False
                        total_state = total_condition(state, V[it_string])
                        p0[it_string] = total_state.p/1.0e6
                        state.with_ions = True
                    except:
                        #failed again, we'll just leave it out
                        print "Failed to find total condition for {0} will add 0.0 to print out.".format(it_string)
                        p0[it_string] = 0.0
                    
            
            if state.p < 1.0e6: #change how the pressure is printed if it's too big, it keeps ruining the printouts!
                if state.rho < 0.01:
                    # also need something for when density is too small...
                    conditions = "{0:<6}{1:<11.7}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.2e}{7:<8.1f}{8:<9.3f}"\
                    .format(it_string, states[it_string].p, states[it_string].T, states[it_string].a,
                            V[it_string],M[it_string], states[it_string].rho, pitot[it_string], p0[it_string])
                else:
                    conditions = "{0:<6}{1:<11.7}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.5f}{7:<8.1f}{8:<9.3f}"\
                    .format(it_string, states[it_string].p, states[it_string].T, states[it_string].a,
                            V[it_string],M[it_string], states[it_string].rho, pitot[it_string], p0[it_string])
                    
            else:
                if state.rho < 0.01:
                    # also need something for when density is too small...
                    conditions = "{0:<6}{1:<11.3e}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.2e}{7:<8.1f}{8:<9.3f}"\
                    .format(it_string, states[it_string].p, states[it_string].T, states[it_string].a,
                            V[it_string],M[it_string], states[it_string].rho, pitot[it_string], p0[it_string])
                else:
                    conditions = "{0:<6}{1:<11.3e}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.5f}{7:<8.1f}{8:<9.3f}"\
                    .format(it_string, states[it_string].p, states[it_string].T, states[it_string].a,
                            V[it_string],M[it_string], states[it_string].rho, pitot[it_string], p0[it_string])                  
                    
            print conditions
            txt_output.write(conditions + '\n')
        
        return

    #print the driver related stuff first

    if 's4i' in states:
        condition_printer('s4i')
    
    condition_printer('s4')
    if M['s3s'] > 0.0:
        condition_printer('s3s')
    
    if cfg['secondary']: #need a separate printing thing for the secondary driver
        for i in range(1,4): #will do 1 - 3
            it_string = 'sd{0}'.format(i)
            condition_printer(it_string)
            if cfg['rs_out_of_sd'] and i == 2: 
                # need to add the reflected shock condition at state 2 if we did the normal shock here
                # this will print it at the right place
                it_string = 'sd{0}r'.format(i)
                condition_printer(it_string) 
            if cfg['sx_into_st'] and i == 2: 
                # need to add the reflected shock condition at state 2 if we did the normal shock here
                # this will print it at the right place
                it_string = 'sd{0}s'.format(i)
                condition_printer(it_string)            
                    
    for i in range(1,4): #shock tube stuff
        it_string = 's{0}'.format(i)
        condition_printer(it_string)
        if cfg['rs_out_of_st'] and i == 2: 
            # need to add the reflected shock condition at state 2 if we did the normal shock here
            # this will print it at the right place
            it_string = 's{0}r'.format(i)
            condition_printer(it_string)            
        
    if cfg['tunnel_mode'] == 'expansion-tube':    
        for i in range(5,8): #acc tube extra states
            if i == 6 and cfg['expand_to'] == 'p7':
                continue
            it_string = 's{0}'.format(i)
            condition_printer(it_string)
                                    
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
         condition_printer('s5')
         condition_printer('s5_d')
         if cfg['nozzle']:
             condition_printer('s6')             
            
    if cfg['nozzle']: #do nozzle calculations if asked to
        condition_printer('s8')
        
    #do the conditions over the model if asked
    if cfg['shock_over_model']:
        condition_printer('s10f')
        condition_printer('s10e')
            
    if cfg['conehead']:
        if 's10cf' in states:
            condition_printer('s10cf')
        if 's10c' in states:
            condition_printer('s10c')            
        
    if cfg['wedge']:
        if cfg['solver'] == 'pg':
            condition_printer('s10w')
        elif cfg['solver'] in ['eq', 'pg-eq']:
            condition_printer('s10wf')
            condition_printer('s10we')
        
    # added some extra code to calculate pre and post nozzle stagnation enthalpy 
    # if the user wants it
    if 'all_total' not in cfg:
        cfg['all_total'] = False
            
    if cfg['all_total'] and cfg['nozzle']:
        print "Calculating nozzle entry total condition as the user has asked for this."
        try:
            states['nozzle_entry_total'] = total_condition(states[cfg['nozzle_entry_state']], V[cfg['nozzle_entry_state']])
            states['nozzle_entry_pitot'] = pitot_condition(states[cfg['nozzle_entry_state']], V[cfg['nozzle_entry_state']])
            cfg['nozzle_entry_stagnation_enthalpy'] = states['nozzle_entry_total'].h - states['s1'].h #J/kg (take away the initial enthalpy in state 1 to get the change)
        except:
            try:
                # try again with ions turned off.
                states[cfg['nozzle_entry_state']].with_ions = False
                states['nozzle_entry_total'] = total_condition(states[cfg['nozzle_entry_state']], V[cfg['nozzle_entry_state']])
                states['nozzle_entry_pitot'] = pitot_condition(states[cfg['nozzle_entry_state']], V[cfg['nozzle_entry_state']])
                cfg['nozzle_entry_stagnation_enthalpy'] = states['nozzle_entry_total'].h - states['s1'].h #J/kg (take away the initial enthalpy in state 1 to get the change)
                states[cfg['nozzle_entry_state']].with_ions = True
            except:
                # just give up if it bails out again
                states['nozzle_entry_total'] = None
                states['nozzle_entry_pitot'] = None
                cfg['nozzle_entry_stagnation_enthalpy'] = None
        if cfg['nozzle_entry_stagnation_enthalpy']:        
            nozzle_entry_stag_enth = 'The total enthalpy (Ht) entering the nozzle is {0:<.5g} MJ/kg ({1} - s1).'\
            .format(cfg['nozzle_entry_stagnation_enthalpy']/10**6, cfg['nozzle_entry_state'])
            print nozzle_entry_stag_enth
            txt_output.write(nozzle_entry_stag_enth + '\n')
                
            
    #some other useful calculations at the end
    try:
        states['test_section_total'] = total_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
        states['test_section_pitot'] = pitot_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
        cfg['stagnation_enthalpy'] = states['test_section_total'].h - states['s1'].h #J/kg (take away the initial enthalpy in state 1 to get the change)
    except:
        try:
            # try again with ions turned off.
            states[cfg['test_section_state']].with_ions = False
            states['test_section_total'] = total_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
            states['test_section_pitot'] = pitot_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
            cfg['stagnation_enthalpy'] = states['test_section_total'].h - states['s1'].h #J/kg (take away the initial enthalpy in state 1 to get the change)
            states[cfg['test_section_state']].with_ions = True
        except:
            # just give up if it bails out again
            states['test_section_total'] = None
            states['test_section_pitot'] = None
            cfg['stagnation_enthalpy'] = None
    
    if cfg['stagnation_enthalpy']:        
        if cfg['nozzle']:        
            stag_enth = 'The total enthalpy (Ht) leaving the nozzle is {0:<.5g} MJ/kg (H8 - h1).'\
            .format(cfg['stagnation_enthalpy']/10**6)
        elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'expansion-tube':
            stag_enth = 'The total enthalpy (Ht) at the end of the acceleration tube (state 7) is {0:<.5g} MJ/kg (H7 - h1).'\
            .format(cfg['stagnation_enthalpy']/10**6)
        elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'nr-shock-tunnel':
            stag_enth = 'The total enthalpy (Ht) at the end of the shock tube (state 2) is {0:<.5g} MJ/kg (H2 - h1).'\
            .format(cfg['stagnation_enthalpy']/10**6)
        elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'reflected-shock-tunnel':
            stag_enth = 'The total enthalpy (Ht) in the stagnated region (state 5) is {0:<.5g} MJ/kg (H5 - h1).'\
            .format(cfg['stagnation_enthalpy']/10**6)
    else:
        stag_enth = "Was unable to find total condition. Therefore, unable to print stagnation enthalpy."
        
    print stag_enth  
    txt_output.write(stag_enth + '\n')
    
    if cfg['stagnation_enthalpy']:        
        if cfg['nozzle']:        
            stag_temp = 'The total temperature (Tt) leaving the nozzle is {0:<.5g} K.'.format(states['test_section_total'].T)
        elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'expansion-tube':
            stag_temp = 'The total temperature (Tt) at the end of the acceleration tube (state 7) is {0:<.5g} K.'.format(states['test_section_total'].T)
        elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'nr-shock-tunnel':
            stag_temp = 'The total temperature (Tt) at the end of the shock tube (state 2) is {0:<.5g} K.'.format(states['test_section_total'].T)
        elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'reflected-shock-tunnel':
            stag_temp = 'The total temperature (Tt) in the stagnated region (state 5) is {0:<.5g} K.'.format(states['test_section_total'].T)
        
    print stag_temp 
    txt_output.write(stag_temp + '\n')
    
    # I'm going to add freestream enthalpy (just the h component) to the output also
    # (take away the initial enthalpy in state 1 to get the change)
    cfg['freestream_enthalpy'] = states[cfg['test_section_state']].h - states['s1'].h #J/kg
    
    if cfg['nozzle']:        
        freestream_enth = 'The freestream enthalpy (h) leaving the nozzle is {0:<.5g} MJ/kg (h8 - h1).'\
        .format(cfg['freestream_enthalpy']/10**6)
    elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'expansion-tube':
        freestream_enth = 'The freestream enthalpy (h) at the end of the acceleration tube (state 7) is {0:<.5g} MJ/kg (h7 - h1).'\
        .format(cfg['freestream_enthalpy']/10**6)
    elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'nr-shock-tunnel':
        freestream_enth = 'The freestream enthalpy (h) at the end of the shock tube (state 2) is {0:<.5g} MJ/kg (h2 - h1).'\
        .format(cfg['freestream_enthalpy']/10**6)
    elif not cfg['nozzle'] and cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        freestream_enth = 'The freestream enthalpy (h) in the stagnated region (state 5) is {0:<.5g} MJ/kg (h5 - h1).'\
        .format(cfg['freestream_enthalpy']/10**6)
        
    print freestream_enth  
    txt_output.write(freestream_enth + '\n')
        
    if cfg['stagnation_enthalpy'] and cfg['stagnation_enthalpy'] >= 0.0:
        #calculate flight equivalent velocity
        #for a description of why this is, refer to Bianca Capra's thesis page 104 - 105
        #Capra, B., Aerothermodynamic Simulation of Subscale Models of the FIRE II and
        #Titan Explorer Vehicles in Expansion Tubes, Ph.D. thesis, the University of Queens-
        #land, St. Lucia, Australia, 2006.
        cfg['u_eq'] = math.sqrt(2.0*cfg['stagnation_enthalpy']) 
        u_eq_print = 'The flight equivalent velocity (Ue) is {0:<.5g} m/s.'.format(cfg['u_eq'])
    else:
        u_eq_print = "Unable to find equivalent velocity as stagnation enthalpy could not be found."
        
    print u_eq_print
    txt_output.write(u_eq_print + '\n')
    
    if cfg['wedge']:
        if 'beta_pg' in cfg:
            frozen_wedge = 'Frozen wedge beta angle is {0:.3f} degrees.'.format(math.degrees(cfg['beta_pg']))
        else:
            frozen_wedge = "Frozen wedge angle did not solve."
        print frozen_wedge
        txt_output.write(frozen_wedge + '\n')
        if cfg['solver'] in ['eq', 'pg-eq']:
            if 'beta_eq' in cfg:
                eq_wedge = 'Equilibrium wedge beta angle is {0:.3f} degrees.'.format(math.degrees(cfg['beta_eq']))
            else:
                eq_wedge = "Equilibrium wedge angle did not solve."
            print eq_wedge
            txt_output.write(eq_wedge + '\n')            

    #if the test time calculation has been done, print it
    if cfg['calculate_test_time']: 
        basic_test_time_printout = 'Basic test time = {0:.2f} microseconds'.format(cfg['t_test_basic']*1.0e6)
        print  basic_test_time_printout
        txt_output.write(basic_test_time_printout + '\n')

    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        rst_1 = "At state 5 gam = {0}, R = {1}.".format(states['s5'].gam, states['s5'].R)
        print rst_1
        txt_output.write(rst_1 + '\n')
        if cfg['solver'] in ['eq', 'pg-eq']:
            rst_2 = "Species in state 5 at equilibrium:"
            print rst_2
            txt_output.write(rst_2 + '\n')
            
            rst_3 = '{0}'.format(states['s5'].species)
            print rst_3
            txt_output.write(rst_3 + '\n')        
    
    #added ability to get the species in the post-shock condition
    
    if cfg['calculate_scaling_information']:
        scaling_intro = "Calculating scaling information for a represenative length scale of {0:.4f} m.".format(cfg['representative_length_scale'])
        print scaling_intro
        txt_output.write(scaling_intro + '\n')
        
        freestream_mu_string = 'Using a freestream ({0}) dynamic viscosity (mu) of {1:.4e} Pa.s.'.format(cfg['test_section_state'],
                                                                                                         states[cfg['test_section_state']].mu)
        print freestream_mu_string
        txt_output.write(freestream_mu_string + '\n')        

        rho_l_product_freestream_print = "Freestream ({0}) rhoL product is {1:.4e} kg/m**2.".format(cfg['test_section_state'], 
                                                                                                    cfg['rho_l_product_freestream'])
        print rho_l_product_freestream_print
        txt_output.write(rho_l_product_freestream_print + '\n')

        pressure_l_product_freestream_print = "Freestream ({0}) pL product is {1:.4f} Pa*m.".format(cfg['test_section_state'], 
                                                                                                    cfg['pressure_l_product_freestream'])
        print pressure_l_product_freestream_print
        txt_output.write(pressure_l_product_freestream_print + '\n')
        
        reynolds_number_freestream_print = "Freestream ({0}) Reynolds number is {1:.4f}.".format(cfg['test_section_state'], 
                                                                                                 cfg['reynolds_number_freestream'])                                                                                     
        
        print reynolds_number_freestream_print
        txt_output.write(reynolds_number_freestream_print + '\n')
    
        knudsen_number_freestream_print = "Freestream ({0}) Knudsen number is {1:.4e}.".format(cfg['test_section_state'], 
                                                                                               cfg['knudsen_number_freestream'])
                                                                                                 
        print knudsen_number_freestream_print
        txt_output.write(knudsen_number_freestream_print + '\n')            
        

        if cfg['shock_over_model']:
            state10e_mu_string = 'Using a test section post normal shock eq (s10e) dynamic viscosity (mu) of {0:.4e} Pa.s.'.format(states['s10e'].mu)
            print state10e_mu_string
            txt_output.write(state10e_mu_string + '\n')    
            
            rho_l_product_state10e_print = "Test section post normal shock eq (s10e) rhoL product is {0:.4e} kg/m**2.".format(cfg['rho_l_product_state10e'])
            print rho_l_product_state10e_print
            txt_output.write(rho_l_product_state10e_print + '\n')
    
            pressure_l_product_state10e_print = "Test section post normal shock eq (s10e) pL product is {0:.4f} Pa.m.".format(cfg['pressure_l_product_state10e'])
            print pressure_l_product_state10e_print
            txt_output.write(pressure_l_product_state10e_print + '\n')   

            reynolds_number_state10e_print = "Test section post normal shock eq (s10e) Reynolds number is {0:.4f}.".format(cfg['reynolds_number_state10e'])
                                                                                                 
            print reynolds_number_state10e_print
            txt_output.write(reynolds_number_state10e_print + '\n')  
     
            knudsen_number_state10e_print = "Test section post normal shock eq (s10e) Knudsen number is {0:.4e}.".format( cfg['knudsen_number_state10e'])
                                                                                                     
            print knudsen_number_state10e_print
            txt_output.write(knudsen_number_state10e_print + '\n')                                                                                                            
                                                                                                    
    
    if cfg['shock_over_model'] and 's10e' in states.keys():
        species1 = 'Species in the shock layer at equilibrium (s10e) (by {0}):'.format(states['s10e'].outputUnits)        
        print species1
        txt_output.write(species1 + '\n')
        
        species2 = '{0}'.format(states['s10e'].species)
        print species2
        txt_output.write(species2 + '\n')
        
    if cfg['test'] == 'fulltheory-pressure-ratios':
        pressure_ratio_line = "p2_p1 = {0}.".format(states['s2'].p / states['s1'].p)
        print pressure_ratio_line
        txt_output.write(pressure_ratio_line + '\n')

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
            
        print '-'*60
                   
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

    test_done = "Test,'{0}'".format(cfg['test'])
    csv_output.write(test_done + '\n')         
        
    csv_facility_used = 'Facility,{0}.'.format(cfg['facility'])        
    csv_output.write(csv_facility_used + '\n')
    
    csv_test_gas_used = 'Test gas (state 1),{0},gamma,{1},R,{2}'.format(cfg['test_gas'],states['s1'].gam,states['s1'].R)
    csv_output.write(csv_test_gas_used + '\n')  
    
    if cfg['custom_accelerator_gas']:
        csv_accelerator_gas_used = 'Custom Accelerator gas (state 5),gamma,{0},R,{1}.'.format(states['s5'].gam,states['s5'].R)
    else:
        csv_accelerator_gas_used = "Accelerator gas (state 5),Air."
    csv_output.write(csv_accelerator_gas_used + '\n')          
    
    if cfg['solver'] == 'eq':
        if cfg['facility'] != 'custom' and cfg['piston'] in ['Sangdi-1.8MPa', 'sangdi-1.8MPa','Sangdi-2.2MPa', 'sangdi-2.2MPa']:
            csv_driver_gas_used = 'Driver gas is {0}.'.format(cfg['driver_gas'])
        else:
            csv_driver_gas_used = 'Driver gas,{0}.'.format(states['s4'].reactants)
    else:
        if cfg['facility'] != 'custom':
            csv_driver_gas_used = 'Driver gas,{0}.'.format(cfg['driver_gas'])
        else:
            if 'driver_composition' in cfg:
                csv_driver_gas_used = 'Driver gas,{0}.'.format(cfg['driver_composition'])
            else:
                csv_driver_gas_used = 'Driver gas,custom pg, gam,{0},R,{1}'.format(cfg['driver_pg_gam'], cfg['driver_pg_R'])
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
    if cfg['rs_out_of_sd']:
        rs_out_of_sd = "Vr-sd,{0:.2f} m/s, Mr-sd,{1:.2f}".format(cfg['Vr-sd'],cfg['Mr-sd'])
        csv_output.write(rs_out_of_sd + '\n') #writes above line to txt_output file (input to write command must be a string)
    if cfg['rs_out_of_st']:
        rs_out_of_st = "Vr-st,{0:.2f} m/s, Mr-st,{1:.2f}".format(cfg['Vr-st'],cfg['Mr-st'])
        csv_output.write(rs_out_of_st + '\n') #writes above line to txt_output file (input to write command must be a string)
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':    
        shockspeeds_2 = "Vr_d,{0:.2f} m/s,Mr_d,{1:.2f}".format(cfg['Vr_d'],cfg['Mr_d'])
        csv_output.write(shockspeeds_2 + '\n')
                
    csv_key = "{0:6},{1:11},{2:9},{3:6},{4:9},{5:6},{6:9},{7:8},{8:9}".format("state","P","T","a","V","M","rho","pitot","stgn")
    csv_output.write(csv_key + '\n')
    
    csv_units = "{0:6},{1:11},{2:9},{3:6},{4:9},{5:6},{6:9},{7:9},{8:9}".format("","Pa","K","m/s","m/s","","kg/m^3","kPa","MPa")
    csv_output.write(csv_units + '\n')
    
    #new dictionaries here to add pitot and stagnation pressure calcs
    
    pitot = {} #pitot pressure dict
    p0 = {} #stagnation pressure dict
           
    def csv_condition_printer(it_string):
        """Prints the values of a specified condition to the screen and to 
        the txt_output file. 
        
        I made a function of this so I didn't have to keep pasting the code in."""
        
        if states.has_key(it_string):
            # I decided to add a clone here to stop messing with the original object here.
            try:
                state = states[it_string].clone()
                clone_successful = True
            except Exception as e:
                print "Error {0}".format(str(e))   
                print "Failed to clone state {0} while trying to final total condition."
                
                if isinstance(states[it_string], Gas) and states[it_string].with_ions:
                    print "Trying again with ions turned off."
                    try:
                        states[it_string].with_ions = False
                        state = states[it_string].clone()
                        clone_successful = True
                    except Exception as e:
                         print "Error {0}".format(str(e))   
                         print "Failed to clone state {0} again while trying to final total condition." 
                         print "0.0 will be added as pitot and total pressures for this state."
                         pitot[it_string] = 0.0
                         p0[it_string] = 0.0
                         clone_successful = False
                else:
                     print "0.0 will be added as pitot and total pressures for this state."
                     pitot[it_string] = 0.0
                     p0[it_string] = 0.0
                     clone_successful = False
                
            if M[it_string] == 0:
                pitot[it_string] = state.p/1000.0
                p0[it_string] = state.p/1.0e6
            elif clone_successful:
                try:
                    pitot[it_string] = pitot_p(state.p,M[it_string],state.gam)/1000.0
                except:
                    try:
                        #try again witcfg, states, V, Mh no ions turned on if it bails out
                        state.with_ions = False
                        pitot[it_string] = pitot_p(state.p,M[it_string],state.gam)/1000.0
                        state.with_ions = True
                    except:
                        # if this fails just give up and print 0.0
                        print "Failed to find pitot pressure for {0} will add 0.0 to csv print out.".format(it_string)
                        pitot[it_string] = 0.0
                try:
                    #make total condition of relevant state for printing
                    total_state = total_condition(state, V[it_string])
                    p0[it_string] = total_state.p/1.0e6
                except:
                    try:
                        #try again with no ions turned on if it bails out
                        state.with_ions = False
                        total_state = total_condition(state, V[it_string])
                        p0[it_string] = total_state.p/1.0e6
                        state.with_ions = True
                    except:
                        #failed again, we'll just leave it out
                        print "Failed to find total condition for {0} will add 0.0 to print out.".format(it_string)
                        p0[it_string] = 0.0
            
                
            if state.p < 1.0e6: #change how the pressure is printed if it's too big, it keeps ruining the printouts!  
                if state.rho < 0.01:
                    # also need something for when density is too small...
                    csv_conditions = "{0:<6},{1:<11.7},{2:<9.1f},{3:<6.0f},{4:<9.1f},{5:<6.2f},{6:<9.2e},{7:<8.1f},{8:<9.3f}"\
                    .format(it_string, states[it_string].p, states[it_string].T, states[it_string].a,
                            V[it_string],M[it_string], states[it_string].rho, pitot[it_string], p0[it_string])     
                else:
                    csv_conditions = "{0:<6},{1:<11.7},{2:<9.1f},{3:<6.0f},{4:<9.1f},{5:<6.2f},{6:<9.5f},{7:<8.1f},{8:<9.3f}"\
                    .format(it_string, states[it_string].p, states[it_string].T, states[it_string].a,
                            V[it_string],M[it_string], states[it_string].rho, pitot[it_string], p0[it_string])     
                    
            else:
                if state.rho < 0.01:
                    # also need something for when density is too small...
                    csv_conditions = "{0:<6},{1:<11.3e},{2:<9.1f},{3:<6.0f},{4:<9.1f},{5:<6.2f},{6:<9.2e},{7:<8.1f},{8:<9.3f}"\
                    .format(it_string, states[it_string].p, states[it_string].T, states[it_string].a,
                            V[it_string],M[it_string], states[it_string].rho, pitot[it_string], p0[it_string])     
                else:
                    csv_conditions = "{0:<6},{1:<11.3e},{2:<9.1f},{3:<6.0f},{4:<9.1f},{5:<6.2f},{6:<9.5f},{7:<8.1f},{8:<9.3f}"\
                    .format(it_string, states[it_string].p, states[it_string].T, states[it_string].a,
                            V[it_string],M[it_string], states[it_string].rho, pitot[it_string], p0[it_string])            
            
            csv_output.write(csv_conditions + '\n')
            
        return
        
    
    if 's4i' in states:
        csv_condition_printer('s4i')

    #print the driver related stuff first
    csv_condition_printer('s4')
    if M['s3s'] > 0.0:
        csv_condition_printer('s3s')
    
    if cfg['secondary']: #need a separate printing thing for the secondary driver
        for i in range(1,4): #will do 1 - 3
            it_string = 'sd{0}'.format(i)
            csv_condition_printer(it_string)
            if cfg['rs_out_of_sd'] and i == 2: 
                # need to add the reflected shock condition at state 2 if we did the normal shock here
                # this will print it at the right place
                it_string = 'sd{0}r'.format(i)
                csv_condition_printer(it_string)  
            if cfg['sx_into_st'] and i == 2: 
                # need to add the reflected shock condition at state 2 if we did the normal shock here
                # this will print it at the right place
                it_string = 'sd{0}s'.format(i)
                csv_condition_printer(it_string)    
                    
    for i in range(1,4): #shock tube stuff
        it_string = 's{0}'.format(i)
        csv_condition_printer(it_string)
        if cfg['rs_out_of_st'] and i == 2: 
            # need to add the reflected shock condition at state 2 if we did the normal shock here
            # this will print it at the right place
            it_string = 's{0}r'.format(i)
            csv_condition_printer(it_string)   
    if cfg['tunnel_mode'] == 'expansion-tube':    
        for i in range(5,8): #acc tube and nozzle if it's there
            if i == 6 and cfg['expand_to'] == 'p7':
                continue
            it_string = 's{0}'.format(i)
            csv_condition_printer(it_string)
            
    if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        csv_condition_printer('s5')
        csv_condition_printer('s5_d')
        if cfg['nozzle']:
            csv_condition_printer('s6') 
    
    if cfg['nozzle']:
        csv_condition_printer('s8')      
        
    #do the conditions over the model
    if cfg['shock_over_model']:
        csv_condition_printer('s10f')
        csv_condition_printer('s10e')
            
    if cfg['conehead']:
        if 's10cf' in states:
            csv_condition_printer('s10cf')
        if 's10c' in states:
            csv_condition_printer('s10c')
        
    if cfg['wedge']:
        if cfg['solver'] == 'pg':
            csv_condition_printer('s10w')
        elif cfg['solver'] in ['eq', 'pg-eq']:
            csv_condition_printer('s10wf')
            csv_condition_printer('s10we')
        
    if 'stagnation_enthalpy' not in cfg:
        #if stagnation enthalpy and u_eq not already calculated, calculate them here
        try:
            states['test_section_total'] = total_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
            states['test_section_pitot'] = pitot_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
            cfg['stagnation_enthalpy'] = states['test_section_total'].h - states['s1'].h #J/kg (take away the initial enthalpy in state 1 to get the change) #J/kg
        except:
            try:
                # try again with ions turned off.
                states[cfg['test_section_state']].with_ions = False
                states['test_section_total'] = total_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
                states['test_section_pitot'] = pitot_condition(states[cfg['test_section_state']], V[cfg['test_section_state']])
                cfg['stagnation_enthalpy'] = states['test_section_total'].h - states['s1'].h #J/kg (take away the initial enthalpy in state 1 to get the change) #J/kg
                states[cfg['test_section_state']].with_ions = True
            except:
                # just give up if it bails out again
                states['test_section_total'] = None
                states['test_section_pitot'] = None
                cfg['stagnation_enthalpy'] = None
        #calculate flight equivalent velocity
        #for a description of why this is, refer to Bianca Capra's thesis page 104 - 105
        #Capra, B., Aerothermodynamic Simulation of Subscale Models of the FIRE II and
        #Titan Explorer Vehicles in Expansion Tubes, Ph.D. thesis, the University of Queens-
        #land, St. Lucia, Australia, 2006.
        if cfg['stagnation_enthalpy']:
            cfg['u_eq'] = math.sqrt(2.0*cfg['stagnation_enthalpy']) 
    if cfg['stagnation_enthalpy']:                                 
        csv_stag_enth = 'Ht,{0:<.5g}, MJ/kg.'.format(cfg['stagnation_enthalpy']/10**6)
        csv_output.write(csv_stag_enth + '\n')
        
    csv_freestream_enth = 'h,{0:<.5g}, MJ/kg.'.format(cfg['freestream_enthalpy']/10**6)
    csv_output.write(csv_freestream_enth + '\n') 
    
    if cfg['stagnation_enthalpy']:   
        stag_temp = 'Tt, {0:<.5g}, K.'.format(states['test_section_total'].T)
        csv_output.write(stag_temp + '\n')
       
    if cfg['stagnation_enthalpy'] and cfg['stagnation_enthalpy'] >= 0.0 and cfg['u_eq']:           
        csv_u_eq_print = 'Ue,{0:<.5g}, m/s.'.format(cfg['u_eq'])
        csv_output.write(csv_u_eq_print + '\n')
        
    if cfg['wedge']:
        if 'beta_pg' in cfg:
            frozen_wedge = 'Frozen wedge beta angle, {0:.3f}, degrees.'.format(math.degrees(cfg['beta_pg']))
        else:
            frozen_wedge = "Frozen wedge is detached."
            csv_output.write(frozen_wedge + '\n')
        if cfg['solver'] in ['eq', 'pg-eq']:
            if 'beta_eq' in cfg:
                eq_wedge = 'Equilibrium wedge beta angle, {0:.3f}, degrees.'.format(math.degrees(cfg['beta_eq']))
            else:
                eq_wedge = "Equilibrium wedge is detached."
            csv_output.write(eq_wedge + '\n')    

    if cfg['calculate_test_time']: 
        csv_basic_test_time_printout = 'Basic test time,{0:.2f}, microseconds'.format(cfg['t_test_basic']*1.0e6)
        csv_output.write(csv_basic_test_time_printout + '\n')   
    
    csv_output.close()
    
    return cfg, states, V, M
    
def make_x_t_diagram(cfg, states, V, M, filename = 'x-t-diagram', show = False):
    """Function to make an x-t diagram if the user requests it.
       Currently only works for the X2 Expansion Tube but other facility values
       could be added in the future...
    """
    
    # X2 diaphragm and PCB locations are from 
    #Gildfind, D. E., Jacobs, P. A., and Morgan, R. G., "Vibration isolation in a free-
    #piston driven expansion tube facility," Shock Waves, Vol. 23, No. 5, 2013, pp. 431-438.
    
    transducer_list = ['sd1', 'sd2', 'sd3', 'st1', 'st2', 'st3', 
                       'at1', 'at2', 'at3', 'at4', 'at5', 'at6']
    # 0-m has been taken as the primary diaphragm location, locations are in m
    # tube exit is with the X2 nozzle, but it's good enough for now...
    loc_dict = {'sd1':2.577,'sd2':2.810, 'sd3':3.043,'st1':4.231,
                'st2':4.746, 'st3':5.260, 'at1':6.437,'at2':6.615,
                'at3':6.796,'at4':7.590, 'at5':7.846,'at6':8.096,
                'secondary diaphragm': 3.418, 'tertiary diaphragm': 5.976,
                'tube exit': 8.585}
                
    if cfg['facility'] != 'x2':
        print "Currently, this mode can only be used for the X2 Expansion Tube, will finish now..."
        return cfg
        
    print "Making x-t diagram of flow through the X2 Expansion Tube."
    
    shock_x_list = [0.0]
    shock_t_list = [0.0]
    
    transducer_times = {}
    
    with open(cfg['filename']+'-' + filename + '-summary.txt', 'w') as output_file:
    
        intro_line = "# x-t diagram summary made using Pitot Version {0}.".format(cfg['VERSION_STRING'])
        print intro_line
        output_file.write(intro_line + '\n')
    
        if cfg['secondary']:
                        
            # we have a secondary driver, so deal with it first...
            shock_x_list.append(loc_dict['secondary diaphragm'])
            shock_time = (shock_x_list[-1] - shock_x_list[-2]) / cfg['Vsd'] + shock_t_list[-1]
            shock_t_list.append(shock_time)
            
            print '-'*60            
            output_file.write('-'*60  + '\n')
            
            for transducer in ['sd1', 'sd2', 'sd3']:
                transducer_time = shock_t_list[-2] + (loc_dict[transducer] - shock_x_list[-2]) / cfg['Vsd']
                transducer_output =  "Shock will pass transducer {0} at t = {1} s ({2} microseconds)."\
                                     .format(transducer, transducer_time, transducer_time * 1.0e6)
                
                print transducer_output
                output_file.write(transducer_output + '\n')
                
                transducer_times[transducer] = transducer_time
                
            Vsd_line =  "Shock will reach end of the secondary driver tube at t = {0} s ({1} microseconds)."\
                        .format(shock_time, shock_time * 1.0e6)
            print Vsd_line
            output_file.write(Vsd_line + '\n')
            
            sec_drv_cs_x_list = [shock_x_list[-2], shock_x_list[-1]]
            cs_time = (sec_drv_cs_x_list[-1] - sec_drv_cs_x_list[-2]) / V['sd2'] + shock_t_list[-1]
            sec_drv_cs_t_list = [0.0, cs_time]
            
            # now do the shock tube...
            
            if cfg['tunnel_mode'] == 'expansion-tube': 
                shock_x_list.append(loc_dict['tertiary diaphragm'])
            else:
                # the shock tube runs right to the tube exit...
                shock_x_list.append(loc_dict['tube exit'])
            shock_time = (shock_x_list[-1] - shock_x_list[-2]) / cfg['Vs1'] + shock_t_list[-1]
            shock_t_list.append(shock_time)
            
            print '-'*60            
            output_file.write('-'*60  + '\n')
            
            if cfg['tunnel_mode'] == 'expansion-tube': 
                shk_tube_transducer_list = ['st1', 'st2', 'st3']
            else:
                shk_tube_transducer_list = ['st1', 'st2', 'st3', 'at1', 'at2', 'at3', 'at4', 'at5', 'at6']
        
            for transducer in shk_tube_transducer_list:
                transducer_time = shock_t_list[-2] + (loc_dict[transducer] - shock_x_list[-2]) / cfg['Vs1']
                transducer_output =  "Shock will pass transducer {0} at t = {1} s ({2} microseconds)."\
                                     .format(transducer, transducer_time, transducer_time * 1.0e6)
                
                print transducer_output
                output_file.write(transducer_output + '\n')
                
                transducer_times[transducer] = transducer_time
                
            Vs1_line =  "Shock will reach end of the shock tube at t = {0} s ({1} microseconds)."\
                        .format(shock_time, shock_time * 1.0e6)
            print Vs1_line
            output_file.write(Vs1_line + '\n')
            
        else:
            # do shock tube without secondary driver...
            if cfg['tunnel_mode'] == 'expansion-tube': 
                shock_x_list.append(loc_dict['secondary diaphragm'])
            else:
                # the shock tube runs right to the tube exit...
                shock_x_list.append(loc_dict['tube exit'])
            shock_time = (shock_x_list[-1] - shock_x_list[-2]) / cfg['Vs1'] + shock_t_list[-1]
            shock_t_list.append(shock_time)
            
            print '-'*60            
            output_file.write('-'*60  + '\n')
            
            if cfg['tunnel_mode'] == 'expansion-tube': 
                shk_tube_transducer_list = ['sd1', 'sd2', 'sd3']
            else:
                shk_tube_transducer_list = transducer_list
        
            for transducer in shk_tube_transducer_list:
                transducer_time = shock_t_list[-2] + (loc_dict[transducer] - shock_x_list[-2]) / cfg['Vs1']
                transducer_output =  "Shock will pass transducer {0} at t = {1} s ({2} microseconds)."\
                                    .format(transducer, transducer_time, transducer_time * 1.0e6)
            
                print transducer_output
                output_file.write(transducer_output + '\n')
                
                transducer_times[transducer] = transducer_time
            
            Vs1_line =  "Shock will reach end of the shock tube at t = {0} s ({1} microseconds)."\
                        .format(shock_time, shock_time * 1.0e6)
            print Vs1_line
            output_file.write(Vs1_line + '\n')
        
        # do the shock tube contact surface... 
        
        shk_tube_cs_x_list = [shock_x_list[-2], shock_x_list[-1]]
        cs_time = (shk_tube_cs_x_list[-1] - shk_tube_cs_x_list[-2]) / V['s2'] + shock_t_list[-2]
        shk_tube_cs_t_list = [shock_t_list[-2], cs_time]
        
        # now do the acc tube if we have an expansion tube...
        
        if cfg['tunnel_mode'] == 'expansion-tube': 
            shock_x_list.append(loc_dict['tube exit'])
            shock_time = (shock_x_list[-1] - shock_x_list[-2]) / cfg['Vs2'] + shock_t_list[-1]
            shock_t_list.append(shock_time) 
            
            print '-'*60            
            output_file.write('-'*60  + '\n')
            
            # do the shock tube contact surface...    
            acc_tube_cs_x_list = [shock_x_list[-2], shock_x_list[-1]]
            cs_time = (acc_tube_cs_x_list[-1] - acc_tube_cs_x_list[-2]) / V['s7'] + shock_t_list[-2]
            acc_tube_cs_t_list = [shock_t_list[-2], cs_time]
            
            if cfg['secondary']:
                acc_tube_transducer_list = ['at1', 'at2', 'at3', 'at4', 'at5', 'at6']
            else:
                acc_tube_transducer_list = ['st1', 'st2', 'st3', 'at1', 'at2', 'at3', 'at4', 'at5', 'at6']
            
            for transducer in acc_tube_transducer_list:
                transducer_time = shock_t_list[-2] + (loc_dict[transducer] - shock_x_list[-2]) /  cfg['Vs2']
                transducer_output =  "Shock will pass transducer {0} at t = {1} s ({2} microseconds)."\
                                    .format(transducer, transducer_time, transducer_time * 1.0e6)
                
                print transducer_output
                output_file.write(transducer_output + '\n')
                
                transducer_times[transducer] = transducer_time
                
            Vs2_line =  "Shock will reach end of the acceleration tube at t = {0} s ({1} microseconds)."\
                        .format(shock_time, shock_time * 1.0e6)
            print Vs2_line
            output_file.write(Vs2_line + '\n')
        
    # now let's do some plotting...
        
    import matplotlib.pyplot as mplt
    
    fig, ax1 = mplt.subplots()
    mplt.hold(True)
    
    # shocks
    ax1.plot(shock_x_list, shock_t_list, 'k-', label = 'shocks')
    # contact surfaces
    if cfg['secondary']:
        ax1.plot(sec_drv_cs_x_list, sec_drv_cs_t_list, 'k--')  
    ax1.plot(shk_tube_cs_x_list, shk_tube_cs_t_list, 'k--', label = 'contact surfaces')    
    if cfg['tunnel_mode'] == 'expansion tube': 
       ax1.plot(acc_tube_cs_x_list, acc_tube_cs_t_list, 'k--')
     
    y_min, y_max = ax1.get_ylim()
     
    if cfg['secondary'] and cfg['tunnel_mode'] == 'expansion-tube':
        mplt.vlines([loc_dict['secondary diaphragm'], loc_dict['tertiary diaphragm']], y_min, y_max, linestyle = 'dotted')
    elif cfg['secondary'] or cfg['tunnel_mode'] == 'expansion-tube':
        mplt.vlines([loc_dict['secondary diaphragm']], y_min, y_max, linestyle = 'dotted')
       
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('t (s)')
    ax1.set_title('x-t diagram for X2')
    ax1.legend(loc = 'best')
        
    print "Saving x-t diagram in pdf, and png formats under the filename {0}.*".format(filename)   
    mplt.savefig(filename + '.png',format='png', dpi=250)
    mplt.savefig(filename + '.pdf',format='pdf', dpi=250)
    
    if show:
        mplt.show()
    else:
        mplt.close(fig)
    
    print "x-t digram information has also been added to the cfg dictionary."
    
    cfg['shock_x_list'] = shock_x_list
    cfg['shock_t_list'] = shock_t_list  
    cfg['transducer_times'] = transducer_times
        
    return cfg 
    
#----------------------------------------------------------------------------
    
def cleanup_function():
    """Function to clean up temp files if we want to."""
    
    import os
    
    if PRINT_STATUS: 
        print "Removing temporary files and leaving the program."
        print "-"*60
    if os.path.isfile('thermo.inp'): os.remove('thermo.inp')
    if os.path.isfile('thermo.out'): os.remove('thermo.out')
    if os.path.isfile('thermo.lib'): os.remove('thermo.lib')
    if os.path.isfile('tmp.inp'): os.remove('tmp.inp')
    if os.path.isfile('tmp.out'): os.remove('tmp.out')
    if os.path.isfile('trans.inp'): os.remove('trans.inp')
    if os.path.isfile('trans.out'): os.remove('trans.out')
    if os.path.isfile('trans.lib'): os.remove('trans.lib')
    
    return

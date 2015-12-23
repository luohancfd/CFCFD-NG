#! /usr/bin/env python
"""
pitot_flow_functions.py: pitot flow functions

This file collects all of the functions that the larger pitot program uses to
process the gases as they travel through the tube.

Chris James (c.james4@uq.edu.au) - 07/05/13 

"""

from cfpylib.nm.zero_solvers import secant
# We base our calculation of gas properties upon calls to the NASA Glenn CEA code.
from cfpylib.gasdyn.cea2_gas import Gas, make_gas_from_name
from cfpylib.gasdyn.gas_flow import *
from cfpylib.gasdyn.ideal_gas_flow import p0_p, pitot_p

#to do some perfect gas stuff

import cfpylib.gasdyn.ideal_gas as pg

from pitot_input_utils import make_test_gas

PRINT_STATUS = 1 #if print status is 1, some basic printouts are done

def secondary_driver_calculation(cfg, states, V, M):
    """Function that does the expansion through the secondary driver 
    tube if required."""
    
    def error_in_velocity_s3s_to_sd3_expansion_pressure_iterator(psd1, state3s=states['s3s'], 
                                               V3sg=V['s3s'], statesd1=states['sd1'],
                                                statesd2=states['sd2'],Vsd=cfg['Vsd'],
                                                steps=cfg['secondary_driver_expansion_steps']):
        """Compute the velocity mismatch for a given pressure ratio across the 
        unsteady expansion in the driver into the secondary driver gas."""
        
        print "current guess for psd1 = {0} Pa".format(psd1)            
        
        statesd1.set_pT(psd1, 300.0) #set s1 at set pressure and ambient temp
        
        (Vsd2, Vsd2g) = normal_shock(statesd1, Vsd, statesd2)
        
        print "rho_sd2 = {0} m**3/kg.".format(statesd2.rho)
        
        #Across the contact surface, p3 == p2
        psd3 = statesd2.p
        
        # Across the expansion, we get a velocity, V5g.
        Vsd3g, statesd3 = finite_wave_dp('cplus', V3sg, state3s, psd3, steps=steps)
        
        print "rho_sd3 = {0} m**3/kg.".format(statesd3.rho)

        return (Vsd2g - Vsd3g)/Vsd2g
        
    def error_in_velocity_s3s_to_sd3_driver_expansion_shock_speed_iterator(Vsd,state3s=states['s3s'], 
                                               V3sg=V['s3s'], statesd1=states['sd1'],
                                                statesd2=states['sd2'],
                                                steps=cfg['secondary_driver_expansion_steps']):
        """Compute the velocity mismatch for a given shock speed with the driver
        unsteady expansion into the secondary driver gas behind it.
        
        Make sure you set the fill pressure in sd1 before you start!"""
        
        print '-'*60
        print "current guess for Vsd = {0} m/s".format(Vsd)            
        
        (Vsd2, Vsd2g) = normal_shock(statesd1, Vsd, statesd2)
               
        #Across the contact surface, p3 == p2
        psd3 = statesd2.p
        
        # Across the expansion, we get a velocity, V5g.
        Vsd3g, statesd3 = finite_wave_dp('cplus', V3sg, state3s, psd3, steps=steps)
        
        print "Current psd2 = {0} Pa, current psd3 = {1} Pa.".format(statesd2.p, statesd3.p)
        print "Current Vsd2g = {0} m/s, current Vsd3g = {1} m/s.".format(Vsd2g, Vsd3g)  

        return (Vsd2g - Vsd3g)/Vsd2g
        
    if PRINT_STATUS: print "Starting unsteady expansion of the primary driver gas into the secondary driver."
    
    if cfg['test'] == "fulltheory-shock": #get psd1 for our chosen shock speed
        cfg['psd1'] = secant(error_in_velocity_s3s_to_sd3_expansion_pressure_iterator, 5000.0, 6000.0, tol=1.0e-5,limits=[500.0,16000.0])
        if PRINT_STATUS: print "From secant solve: psd1 = {0} Pa".format(cfg['psd1'])
        #start using psd1 now, compute states sd1,sd2 and sd3 using the correct psd1
        if PRINT_STATUS: print "Now that psd1 is known, finding conditions at states sd2 and sd3."
        states['sd1'].set_pT(cfg['psd1'],cfg['T0'])
    
    elif cfg['test'] == "fulltheory-pressure" or cfg['test'] == 'fulltheory-pressure-ratios': #get Vsd for our chosen fill pressure
        if 'Vsd_guess_1' in cfg and 'Vsd_guess_2' in cfg:
            print "Using custom guesses for Vsd secant solver."
            print "('Vsd_guess_1' = {0} m/s and 'Vsd_guess_2' = {1} m/s)".\
                  format(cfg['Vsd_guess_1'], cfg['Vsd_guess_2'])
        elif 'Vsd_guess_1' not in cfg and 'Vsd_guess_2' not in cfg and cfg['tunnel_mode'] == 'expansion-tube':
            if states['s4'].T < 400.0:
                #lower guess for cold driver
                cfg['Vsd_guess_1'] = 1000.0; cfg['Vsd_guess_2'] = 1500.0
            else:
                cfg['Vsd_guess_1'] = 4000.0; cfg['Vsd_guess_2'] = 5000.0
        elif 'Vsd_guess_1' not in cfg and 'Vsd_guess_2' not in cfg and cfg['tunnel_mode'] == 'nr-shock-tunnel':
            if cfg['psd1'] < 300000.0:                
                cfg['Vsd_guess_1'] = 7000.0; cfg['Vsd_guess_2'] = 8000.0
            else:
                cfg['Vsd_guess_1'] = 3000.0; cfg['Vsd_guess_2'] = 4000.0
        else:
            cfg['Vsd_guess_1'] = 4000.0; cfg['Vsd_guess_2'] = 5000.0
            
        if 'Vsd_lower' not in cfg and 'Vsd_upper' not in cfg:
                cfg['Vsd_lower'] = 500.0; cfg['Vsd_upper'] = 17000.0                
        else:
            print "Using custom limits for Vsd secant solver."
            print "('Vsd_lower' = {0} m/s and 'Vsd_upper' = {1} m/s)".\
                  format(cfg['Vsd_lower'], cfg['Vsd_upper'])
        cfg['secondary_driver_secant_tol'] = 1.0e-5
        cfg['Vsd'] = secant(error_in_velocity_s3s_to_sd3_driver_expansion_shock_speed_iterator, 
                            cfg['Vsd_guess_1'], cfg['Vsd_guess_2'], 
                            tol=cfg['secondary_driver_secant_tol'],
                            limits=[cfg['Vsd_lower'], cfg['Vsd_upper']],
                            max_iterations = 100)
                            
        if cfg['Vsd'] == 'FAIL':
            print "Secant solver failed to settle on a Vs2 value after 100 iterations."
            print "Dropping tolerance from {0} to {1} and trying again..."\
                  .format(cfg['secondary_driver_secant_tol'], cfg['secondary_driver_secant_tol']*10.0)
            cfg['secondary_driver_secant_tol'] = cfg['secondary_driver_secant_tol'] * 10.0
            cfg['Vsd'] = secant(error_in_velocity_s3s_to_sd3_driver_expansion_shock_speed_iterator, 
                                cfg['Vsd_guess_1'], cfg['Vsd_guess_2'], 
                                tol=cfg['secondary_driver_secant_tol'],
                                limits=[cfg['Vsd_lower'], cfg['Vsd_upper']],
                                max_iterations = 100)
            if cfg['Vsd'] == 'FAIL':
                print "This still isn't working. Bailing out."
                raise Exception, "pitot_flow_functions.shock_tube_calculation() Run of pitot has failed in the secondary driver calculation."                              
                            
        if PRINT_STATUS: 
            print '-'*60
            print "From secant solve: Vsd = {0} m/s".format(cfg['Vsd'])
        #start using Vs1 now, compute states 1,2 and 3 using the correct Vs1
        if PRINT_STATUS: 
            print '-'*60
            print "Now that Vsd is known, finding conditions at states sd2 and sd3."
    
    (Vsd2, V['sd2']) = normal_shock(states['sd1'], cfg['Vsd'], states['sd2'])
    V['sd3'], states['sd3'] = finite_wave_dv('cplus', V['s3s'], states['s3s'], V['sd2'], steps=cfg['secondary_driver_expansion_steps'])

    if PRINT_STATUS: 
        print "state sd2: p = {0:.2f} Pa, T = {1:.2f} K.".format(states['sd2'].p, states['sd2'].T) 
        print "state sd3: p = {0:.2f} Pa, T = {1:.2f} K.".format(states['sd3'].p, states['sd3'].T) 

    #get mach numbers for the txt_output
    cfg['Msd1'] = cfg['Vsd']/states['sd1'].a
    M['sd2'] = V['sd2']/states['sd2'].a
    M['sd3']= V['sd3']/states['sd3'].a
    
    if PRINT_STATUS: print '-'*60
    
    return cfg, states, V, M
    
#----------------------------------------------------------------------------

def shock_tube_calculation(cfg, states, V, M):
    """Function that performs all of the shock tube calculations."""

    #--------------------- shock tube functions ---------------------------------
        
    def error_in_velocity_shock_tube_expansion_pressure_iterator(p1, 
                                               expanding_state=states[cfg['shock_tube_expansion']], 
                                               expansion_start_V=V[cfg['shock_tube_expansion']], 
                                                state1=states['s1'],
                                                state2=states['s2'],Vs1=cfg['Vs1'],
                                                solver=cfg['solver'], test_gas=cfg['test_gas'],
                                                steps=cfg['shock_tube_expansion_steps']):
        """Compute the velocity mismatch for a given pressure ratio across the 
        unsteady expansion in the shock tube."""
        
        print "current guess for p1 = {0} Pa".format(p1)            
        
        state1.set_pT(p1, 300.0) #set s1 at set pressure and ambient temp
        
        if solver == 'eq' or solver == 'pg':
            (V2, V2g) = normal_shock(state1, Vs1, state2)  
        elif solver == 'pg-eq': #if we're using the perfect gas eq solver, we want to make an eq state1, set it's T and p, clone it to state2_eq and then shock process it through CEA
            if gasName == 'mars' or gasName == 'co2' or gasName == 'venus':
                state1_eq, nothing1, nothing2, nothing3 = make_test_gas(test_gas) #make eq state1
            else: 
                state1_eq, nothing = make_test_gas(test_gas) #make eq state1
            state2 = state1_eq.clone() #clone state 1 to a new state2
            #set state 2 with state 1's properties as we're about to shock it:
            state2.p = state1.p; state2.T = state1.T
            state2.shock_process(Vs1) #shock it
            #now we're going to use continuity to set the gas velocity manually here:
            V1 = Vs1
            V2 = V1 * state1.rho / state2.rho
            V2g = V1 - V2
               
        #Across the contact surface, p3 == p2
        p3 = state2.p
        
        # Across the expansion, we get a velocity, V5g.
        V3g, state3 = finite_wave_dp('cplus', expansion_start_V, expanding_state, p3,steps=steps)

        return (V2g - V3g)/V2g
        
    def p1_reflected_iterator(p1, statesd2=states[cfg['shock_tube_expansion']], 
                              Vsd2g=V[cfg['shock_tube_expansion']], state1=states['s1'], 
                              state2=states['s2'], state3=states['s3'], Vs1=cfg['Vs1']):
        """Compute the velocity mismatch for a given pressure ratio across the reflected
            shock from the secondary driver into the shock tube. 
            
            Used when Vsd > Vs1"""
                    
        print "current guess for p1 = {0} Pa".format(p1)

        state1.set_pT(p1,300.0) #set s1 at set pressure and ambient temp
        
        V2, V2g = normal_shock(state1, Vs1, state2)
        
        def reflected_shock_speed_iterator(Vr,statesd2=statesd2,
                                           Vsd2g=Vsd2g, state3=state3,
                                           V2g=V2g):
            """iterates out the reflected shock speed that gives the right velocity behind
                a normal shock."""
    
            print "Current guess for Vr = {0} m/s".format(Vr)
            
            print Vsd2g-Vr
    
            V3, V3g = normal_shock(statesd2, Vsd2g-Vr, state3)
    
            return (V3g-(V2g+Vr))/V3g
        
        Vr = secant(reflected_shock_speed_iterator, 1000.0, 3000.0, tol=1.0e-4)
    
        V3,V3g = normal_shock(statesd2,Vsd2g-Vr,state3)
    
        return (state3.p-state2.p)/state2.p
        
    def error_in_velocity_shock_tube_expansion_shock_speed_iterator(Vs1,
                                                expanding_state=states[cfg['shock_tube_expansion']], 
                                                expansion_start_V=V[cfg['shock_tube_expansion']],                    
                                                state1=states['s1'], state2=states['s2'],
                                                solver=cfg['solver'], test_gas=cfg['test_gas'],
                                                steps=cfg['shock_tube_expansion_steps'],
                                                gas_guess = cfg['gas_guess']):
        """Compute the velocity mismatch for a given shock speed with the shock tube
        unsteady expansion behind it.
        
        Make sure you set the fill pressure in state 1 before you start!"""
        
        print '-'*60
        print "Current guess for Vs1 = {0} m/s".format(Vs1)            
        
        try:
            if solver == 'eq' or solver == 'pg':
                try:
                    if Vs1 > 8500 and solver == 'eq':
                        (V2, V2g) = normal_shock(state1, Vs1, state2, gas_guess)
                    else:
                        (V2, V2g) = normal_shock(state1, Vs1, state2)
                except:
                    print "Normal shock failed. Trying again with a gas guess."
                    (V2, V2g) = normal_shock(state1, Vs1, state2, gas_guess) 
                    
            elif solver == 'pg-eq': #if we're using the perfect gas eq solver, we want to make an eq state1, set it's T and p, clone it to state2_eq and then shock process it through CEA
                if test_gas == 'mars' or test_gas == 'co2' or test_gas == 'venus':
                    state1_eq, nothing1, nothing2, nothing3 = make_test_gas(test_gas) #make eq state1
                else: 
                    state1_eq, nothing = make_test_gas(test_gas) #make eq state1
                state2 = state1_eq.clone() #clone state 1 to a new state2
                #set state 2 with state 1's properties as we're about to shock it:
                state2.p = state1.p; state2.T = state1.T
                state2.shock_process(Vs1) #shock it
                #now we're going to use continuity to set the gas velocity manually here:
                V1 = Vs1
                V2 = V1 * state1.rho / state2.rho
                V2g = V1 - V2
        except Exception as e:
            print "Error {0}".format(str(e))
            raise Exception, "pitot_flow_functions.shock_tube_calculation() Normal shock calculation in the shock tube failed."   
        
        # checking to make sure normal shock worked properly and the pressure rose
        # as I had some issues with this with very low shock tube fill pressures (~1.0 Pa)
        if abs((state2.p - state1.p) / state2.p) < 0.10:
            print "For some reason p2 and p1 are too similar. Shock must have not occured properly."
            print "p1 = {0} Pa, p2 = {1} Pa.".format(state1.p, state2.p)
            raise Exception, "pitot_flow_functions.shock_tube_calculation() Normal shock calculation in the shock tube failed."
                     
        #Across the contact surface, p3 == p2
        p3 = state2.p
        
        # Across the expansion, we get a velocity, V5g.
        try:
            V3g, state3 = finite_wave_dp('cplus', expansion_start_V, expanding_state, p3, steps=steps)
        except Exception as e:
            print "Error {0}".format(str(e))
            print "Unsteady expansion into the shock tube failed. Trying again with ions turned off."
            expanding_state.no_ions = True
            try:
                V3g, state3 = finite_wave_dp('cplus', expansion_start_V, expanding_state, p3, steps=steps)
            except Exception as e:
                print "Error {0}".format(str(e))                
                raise Exception, "pitot_flow_functions.shock_tube_calculation() Unsteady expansion into the shock tube failed." 
                
        print "Current p2 = {0} Pa, current p3 = {1} Pa.".format(state2.p, state3.p)
        print "Current V2g = {0} m/s, current V3g = {1} m/s.".format(V2g, V3g)                
        
        return (V2g - V3g)/V2g

    def primary_shock_speed_reflected_iterator(Vs1,statesd2=states[cfg['shock_tube_expansion']],
                                               Vsd2g=V[cfg['shock_tube_expansion']],state1=states['s1'],
                                               state2=states['s2'],state3=states['s3']):
        """The other new iterator for pitot. Has another iterator inside it, will take
        AGES to run. Oh well, no other choice as far as I see."""
    
        print "Current guess for Vs1 = {0} m/s".format(Vs1)
        
        print "-"*60
        
        V2, V2g = normal_shock(state1, Vs1, state2)
    
        def reflected_shock_speed_iterator(Vr,statesd2=statesd2,
                                           Vsd2g=Vsd2g, state3=state3,
                                           V2g=V2g):
            """iterates out the reflected shock speed that gives the right velocity behind
                a normal shock."""
    
            print "Current guess for Vr = {0} m/s".format(Vr)
            
            print "-"*60            
            
            print "Vsd2g = {0} m/s.".format(Vsd2g)
            print "Mr = {0}.".format((Vsd2g-Vr)/statesd2.a)
            
            try:
                V3, V3g = normal_shock(statesd2, Vsd2g-Vr, state3)
            except Exception as e:
                print "Error {0}".format(str(e))
                raise Exception, "Normal shock calculation failed with the current Vr guess." 
            
            #print "V2g = {0}, V3g = {1}, V3 = {2}.".format(V2g, V3g, V3)
            
            velocity_difference = abs((V3g-V2g)/V3g)*100.0
            
            print "V2g = {0}, V3g = {1}, difference = {2} %.".format(V2g, V3g, velocity_difference)
            
            print "-"*60
    
            return (V3g-V2g)/V3g
    
        Vr = secant(reflected_shock_speed_iterator, -1000.0, -1200.0, tol=1.0e-4, limits=[-4000, 400.0])
    
        V3,V3g = normal_shock(statesd2,Vsd2g-Vr,state3)
        
        pressure_difference = abs((state3.p - state2.p) / state3.p) *100.0 
        
        print "p2 = {0} Pa, p3 = {1} Pa, difference = {2} %.".format(state2.p, state3.p, pressure_difference)
        
        print "-"*60
    
        return (state3.p-state2.p)/state3.p


    #----------------------- shock tube calculations --------------------------------                           
      
    if cfg['test'] == "fulltheory-shock": #get p1 for our chosen shock speed
        if cfg['secondary']: #if secondary, check which shock speed is greater, Vsd or Vs1
            if cfg['Vsd'] > cfg['Vs1']: #we get a shock into the shock tube, not an expansion
                cfg['shock_switch'] = True #flip over this switch
                if PRINT_STATUS: print "Vsd is a stronger shock than Vs1. Therefore a normal shock processes the secondary driver gas moving into the shock tube."
                cfg['p1'] = secant(p1_reflected_iterator, 1000.0, 100000.0, tol=1.0e-5,limits=[100.0,1000000.0])
            else: #same expansion as we're used to
                if PRINT_STATUS: print "Starting unsteady expansion of the secondary driver gas into the shock tube."
                cfg['p1'] = secant(error_in_velocity_shock_tube_expansion_pressure_iterator, 1000.0, 100000.0, tol=1.0e-4,limits=[100.0,1000000.0])
        else: #just do the expansion
            if PRINT_STATUS: print "Starting unsteady expansion of the driver gas into the shock tube."
            cfg['p1'] = secant(error_in_velocity_shock_tube_expansion_pressure_iterator, 1000.0, 100000.0, tol=1.0e-4,limits=[100.0,1500000.0])

        if PRINT_STATUS: print "From secant solve: p1 = {0} Pa".format(cfg['p1'])
        #start using p1 now, compute states 1,2 and 3 using the correct p1
        if PRINT_STATUS: print "Now that p1 is known, finding conditions at states 2 and 3."
        states['s1'].set_pT(cfg['p1'],cfg['T0'])
        
    elif cfg['test'] =="fulltheory-pressure" or cfg['test'] == 'fulltheory-pressure-ratios' \
        or cfg['test'] == 'theory-shock-tube-experiment-acc-tube': #get Vs1 for our chosen fill pressure
        if cfg['state2_no_ions']:
            # Need to turn ions off for state 2 here if it is required to make 
            # shock to state 2 work
            states['s2'].with_ions = False 
        if cfg['state3_no_ions']:
            # Need to turn ions off for state 2 here if it is required to make 
            # shock to state 2 work
            states['s3'].with_ions = False 
        if cfg['shock_switch']: #if we've been told to do a shock here instead of an expansion, do a shock instead of an expansion
            if PRINT_STATUS: print "The shock switch is turned on, therefore doing a shock here instead of the normal expansion... Turn this off if you didn't want it" 
            cfg['Vs1'] = secant(primary_shock_speed_reflected_iterator, 2000.0, 1500.0, tol=1.0e-6,limits=[500.0,10000.0])
        else: #just do the expansion
            if cfg['secondary']:
                if PRINT_STATUS: print "Starting unsteady expansion of the secondary driver gas into the shock tube."
            else:
                if PRINT_STATUS: print "Starting unsteady expansion of the driver gas into the shock tube."
            if 'Vs1_guess_1' not in cfg and 'Vs1_guess_2' not in cfg and cfg['tunnel_mode'] in ['expansion-tube', 'nr-shock-tunnel']:
                if cfg['secondary']:
                    cfg['Vs1_guess_1'] = cfg['Vsd'] + 2000.0; cfg['Vs1_guess_2'] = cfg['Vsd'] + 3000.0
                elif states['s4'].T < 400.0 and not cfg['secondary']:
                    # do a lower guess if we have a cold driver.
                    # Chris James (28/09/15)
                    cfg['Vs1_guess_1'] = 2000.0; cfg['Vs1_guess_2'] = 3000.0
                elif cfg['p1'] > 1000.0 and not cfg['secondary']:
                    cfg['Vs1_guess_1'] = 4000.0; cfg['Vs1_guess_2'] = 6000.0
                elif cfg['p1'] < 100.0 and not cfg['secondary']:
                    cfg['Vs1_guess_1'] = 10000.0; cfg['Vs1_guess_2'] = 12000.0
                else:
                    cfg['Vs1_guess_1'] = 6000.0; cfg['Vs1_guess_2'] = 8000.0
            elif 'Vs1_guess_1' not in cfg and 'Vs1_guess_2' not in cfg and cfg['tunnel_mode'] == 'reflected-shock-tunnel':   
                #start with a much lower speed guess in reflected-shock-tunnel mode
                if cfg['secondary']:
                    cfg['Vs1_guess_1'] = cfg['Vsd']+2000.0; cfg['Vs1_guess_2'] = cfg['Vsd']+3000.0
                else: 
                    cfg['Vs1_guess_1'] = 1000.0; cfg['Vs1_guess_2'] = 2000.0
            else:
                print "Using custom guesses for Vs1 secant solver."
                print "('Vs1_guess_1' = {0} m/s and 'Vs1_guess_2' = {1} m/s)".\
                      format(cfg['Vs1_guess_1'], cfg['Vs1_guess_2'])
                      
            if 'Vs1_lower' not in cfg and 'Vs1_upper' not in cfg:
                    cfg['Vs1_lower'] = 1000.0; cfg['Vs1_upper'] = 1000000.0                
            else:
                print "Using custom limits for Vs1 secant solver."
                print "('Vs1_lower' = {0} m/s and 'Vs1_upper' = {1} m/s)".\
                      format(cfg['Vs1_lower'], cfg['Vs1_upper'])
                      
            cfg['shock_tube_secant_tol'] = 1.0e-5
            cfg['Vs1'] = secant(error_in_velocity_shock_tube_expansion_shock_speed_iterator, 
                                cfg['Vs1_guess_1'], cfg['Vs1_guess_2'],
                                cfg['shock_tube_secant_tol'],
                                limits=[cfg['Vs1_lower'], cfg['Vs1_upper']],
                                max_iterations=100)
            if cfg['Vs1'] == 'FAIL':
                print "Secant solver failed to settle on a Vs1 value after 100 iterations."
                print "Dropping tolerance from {0} to {1} and trying again..."\
                      .format(cfg['shock_tube_secant_tol'], cfg['shock_tube_secant_tol']*10.0)
                cfg['shock_tube_secant_tol'] = cfg['shock_tube_secant_tol'] * 10.0
                cfg['Vs1'] = secant(error_in_velocity_shock_tube_expansion_shock_speed_iterator, 
                            cfg['Vs1_guess_1'], cfg['Vs1_guess_2'],
                            cfg['shock_tube_secant_tol'],
                            limits=[cfg['Vs1_lower'], cfg['Vs1_upper']],
                            max_iterations=100)
                if cfg['Vs1'] == 'FAIL':
                    print "Secant solver failed again with the lower tolerance."
                    print "Dropping tolerance from {0} to {1} and trying again..."\
                          .format(cfg['shock_tube_secant_tol'], cfg['shock_tube_secant_tol']*10.0)
                    cfg['shock_tube_secant_tol'] = cfg['shock_tube_secant_tol'] * 10.0
                    cfg['Vs1'] = secant(error_in_velocity_shock_tube_expansion_shock_speed_iterator, 
                                cfg['Vs1_guess_1'], cfg['Vs1_guess_2'],
                                cfg['shock_tube_secant_tol'],
                                limits=[cfg['Vs1_lower'], cfg['Vs1_upper']],
                                max_iterations=100)
                    if cfg['Vs1'] == 'FAIL':
                        print "Secant solver failed again with the lower tolerance."
                        print "Will try with the higher tolerance and 'state2_no_ions' turned on."
                        cfg['shock_tube_secant_tol'] = cfg['shock_tube_secant_tol'] / 10.0
                        cfg['state2_no_ions'] = True
                        states['s2'].with_ions = False
                        cfg['Vs1'] = secant(error_in_velocity_shock_tube_expansion_shock_speed_iterator, 
                                cfg['Vs1_guess_1'], cfg['Vs1_guess_2'],
                                cfg['shock_tube_secant_tol'],
                                limits=[cfg['Vs1_lower'], cfg['Vs1_upper']],
                                max_iterations=100)
                        if cfg['Vs1'] == 'FAIL':
                            print "This still isn't working. Bailing out."
                            raise Exception, "pitot_flow_functions.shock_tube_calculation() Run of pitot has failed in the shock tube calculation."  

        if PRINT_STATUS: 
            print '-' * 60
            print "From secant solve: Vs1 = {0} m/s".format(cfg['Vs1'])
        #start using Vs1 now, compute states 1,2 and 3 using the correct Vs1
        if PRINT_STATUS: 
            print '-' * 60
            print "Now that Vs1 is known, finding conditions at states 2 and 3."
    # first do the normal shock
    try:
        if cfg['Vs1'] > 8500 and 'solver' in ['eq', 'pg-eq']:
            (V2, V['s2']) = normal_shock(states['s1'], cfg['Vs1'], states['s2'], cfg['gas_guess'])
        else:
            (V2, V['s2']) = normal_shock(states['s1'], cfg['Vs1'], states['s2'])
    except:
        print "Normal shock failed. Trying again with a gas guess."
        (V2, V['s2']) = normal_shock(states['s1'], cfg['Vs1'], states['s2'], cfg['gas_guess'])
        
    if cfg['solver'] == 'pg-eq': #if we're using the pg-eq solver this is the point where we move from pg to eq gas objects
        if cfg['test_gas'] == 'mars' or cfg['test_gas'] == 'co2' or cfg['test_gas'] == 'venus':
            states['s1-eq'], nothing1, nothing2, nothing3 = make_test_gas(cfg['test_gas']) #make eq state1
        else: 
            states['s1-eq'], nothing = make_test_gas(cfg['test_gas']) #make eq state1
        states['s2'] = states['s1-eq'].clone() #get a new state2 ready   
        #set state 2 with state 1's properties as we're about to shock it:
        states['s2'].p = states['s1'].p; states['s2'].T = states['s1'].T
        states['s2'].shock_process(cfg['Vs1']) #shock it
        states['s7'] = states['s2'].clone() #need to turn s7 into an eq object too
        #now we're going to use continuity to set the gas velocity manually here:
        V1 = cfg['Vs1']
        V2 = V1 * states['s1'].rho / states['s2'].rho
        V['s2'] = V1 - V2
        
    if cfg['shock_switch'] and cfg['secondary']: #do a shock here if required
        V['s3'] = V['s2']
        #need to back out Vr again here, for now I was lazy and put the function back in:
        if PRINT_STATUS: print "Need to reiterate to find Vr again here..."   
        def reflected_shock_speed_iterator(Vr,statesd2=states['sd2'],
                                           Vsd2g=V[cfg['shock_tube_expansion']],
                                            state3=states['s3'], V2g=V['s2']):    
            """iterates out the reflected shock speed that gives the right velocity behind
                a normal shock."""
    
            print "Current guess for Vr = {0} m/s".format(Vr)
            
            print "-"*60            
            
            print "Vsd2g = {0} m/s.".format(Vsd2g)
            print "Mr = {0}.".format((Vsd2g-Vr)/statesd2.a)
            
            try:
                V3, V3g = normal_shock(statesd2, Vsd2g-Vr, state3)
            except Exception as e:
                print "Error {0}".format(str(e))
                raise Exception, "Normal shock calculation failed with the current Vr guess." 
            
            #print "V2g = {0}, V3g = {1}, V3 = {2}.".format(V2g, V3g, V3)
            
            velocity_difference = abs((V3g-V2g)/V3g)*100.0
            
            print "V2g = {0}, V3g = {1}, difference = {2} %.".format(V2g, V3g, velocity_difference)
            
            print "-"*60
    
            return (V3g-V2g)/V3g
        
        cfg['Vr'] = secant(reflected_shock_speed_iterator, -1000.0, -1200.0, tol=1.0e-4, limits=[-4000, 400.0])
        cfg['Mr'] = cfg['Vr']/states['sd2'].a
        
        (V3,Vjunk) = normal_shock(states['sd2'],V[cfg['shock_tube_expansion']] - cfg['Vr'],states['s3'])
    else:
        if cfg['tunnel_mode'] == 'reflected-shock-tunnel':
            #I seemed to have issues with using finite_wave_dv and the reflected shock tube mode
            # so I made it do finite_wave_dp instead
            V['s3'], states['s3'] = finite_wave_dp('cplus', V[cfg['shock_tube_expansion']], states[cfg['shock_tube_expansion']], states['s2'].p,cfg['shock_tube_expansion_steps'])
        else:
            V['s3'], states['s3'] = finite_wave_dv('cplus', V[cfg['shock_tube_expansion']], states[cfg['shock_tube_expansion']], V['s2'],cfg['shock_tube_expansion_steps'])
            # i'm adding the same stuff as the if statement above here incase there are issues with this expansion...
            if abs(V['s3'] - V['s2']) / V['s3'] > 0.10:
                print "For some reason V2 and V3 are too different. Expansion must have not occured properly."
                print "V2 = {0} Pa, V3 = {1} Pa.".format(V['s2'], V['s3'])
                print "Going to try expanding to a pressure now instead of a velocity."
                V['s3'], states['s3'] = finite_wave_dp('cplus', V[cfg['shock_tube_expansion']], states[cfg['shock_tube_expansion']], states['s2'].p,cfg['shock_tube_expansion_steps'])
    
    if cfg['tunnel_mode'] == 'expansion-tube' and cfg['V2_mirels_limit']:
        print "Setting V2 to Mirels limit of Vs1 value ({0} m/s)".format(cfg['Vs1'])
        V['s2'] = cfg['Vs1']
    
    if PRINT_STATUS: 
        print "state 2: p = {0:.2f} Pa, T = {1:.2f} K, V = {2:.2f} m/s.".format(states['s2'].p, states['s2'].T, V['s2']) 
        print "state 3: p = {0:.2f} Pa, T = {1:.2f} K, V = {2:.2f} m/s.".format(states['s3'].p, states['s3'].T, V['s3']) 
        if cfg['solver'] == 'eq' or cfg['solver'] == 'pg-eq':
            print 'species in state2 at equilibrium:'               
            print '{0}'.format(states['s2'].species)
            print 'state 2 gamma = {0}, state 2 R = {1}.'.format(states['s2'].gam,states['s2'].R)
        try:
            states['s2_total'] = total_condition(states['s2'], V['s2'])
            print 'The total enthalpy (Ht) at state 2 is {0:<.5g} MJ/kg (H2 - h1).'\
            .format((states['s2_total'].h - states['s1'].h)/1.0e6) #(take away the initial enthalpy in state 1 to get the change)
            cfg['Ht2'] = states['s2_total'].h - states['s1'].h
        except Exception as e:
            print e
            print "Failed to find total condition for state 2. Result will be set to None"
            cfg['Ht2'] = None
        
    if cfg['state2_no_ions']:
        # Turn with ions back on so it will be on for other states based on s7
        # if we turned it off to make the unsteady expansion work
        states['s2'].with_ions = True 
    if cfg['state3_no_ions']:
        # Turn with ions back on here if we turned it off
        states['s3'].with_ions = True 
    
    #get mach numbers for the txt_output
    cfg['Ms1'] = cfg['Vs1']/states['s1'].a
    M['s2'] = V['s2']/states['s2'].a
    M['s3']= V['s3']/states['s3'].a
    
    #erase the gas guess if we don't think we'll need it in the acc tube
    
    if ('p5' in cfg and cfg['solver'] == 'eq' or cfg['solver'] == 'pg-eq') \
    or (cfg['Vs2'] >= 9000.0 and cfg['solver'] == 'eq' or cfg['solver'] == 'pg-eq'): 
        # just use the gas guess if you don't know if you'll need it
        # and set an air guess too
        cfg['gas_guess_air'] = {'gam':1.35,'R':571.49}     
    else:
        cfg['gas_guess'] = None ; cfg['gas_guess_air'] = None
        
    if PRINT_STATUS: print '-'*60
        
    return cfg, states, V, M
    
#----------------------------------------------------------------------------
    
def shock_tube_rs_calculation(cfg, states, V, M):
    """This is a small function that performs a reflected shock on state 2
    when pitot is ask to simulate the effects of the reflected shock at the diaphragm.
    It will default to doing a full reflected shock or the user can mess around with the
    reflected shock Mach number to try 
    You add in this calc by adding 'perform_rs = True' to the cfg file.
    """
    
    if cfg['tunnel_mode'] != 'expansion-tube':
        print "Reflected shock calculation into the acceleration tube was asked for, but you are not running in expansion tube mode..."
        
        return cfg, states, V, M
        
    print "Doing reflected shock at the end of the shock tube that the user has asked for."
    
    #first build state2r as a clone of state 2
    states['s2r'] = states['s2'].clone()

    if cfg['Mr_st'] == "Maximum": # then perform the reflected shock
        cfg['Vr-st'] = reflected_shock(states['s2'], V['s2'], states['s2r'])
        cfg['Mr-st'] = (V['s2']+cfg['Vr-st'])/states['s2'].a #normally this would be V2 - V2r, but it's plus here as Vr has been left positive
        V['s2r'] = 0.0
        M['s2r']= V['s2r']/states['s2r'].a
    else: # perform it to a set strength
        cfg['Mr-st'] = cfg['Mr_st']
        cfg['Vr-st'] = cfg['Mr-st']*states['s2r'].a - V['s2']
        try:
            (V2r, V2rg) = normal_shock(states['s2'], cfg['Vr-st'] + V['s2'], states['s2r'])
            V['s2r'] = V['s2'] - V2rg
            M['s2r']= V['s2r']/states['s2r'].a           
        except Exception as e:
            print "Error {0}".format(str(e))
            raise Exception, "Reflected normal shock calculation at the end of the shocl tube failed."
           
    if PRINT_STATUS:
        print "Vr-st = {0} m/s, Mr-st = {1}.".format(cfg['Vr-st'], cfg['Mr-st'])
        print "state 2r: p = {0:.2f} Pa, T = {1:.2f} K.".format(states['s2r'].p, states['s2r'].T)
        print "Vs2r = {0} m/s, Ms2r = {1}.".format(V['s2r'], M['s2r'])
        if cfg['solver'] == 'eq' or cfg['solver'] == 'pg-eq':
            print 'species in state2r at equilibrium:'               
            print '{0}'.format(states['s2r'].species)
        try:
            states['s2r_total'] = total_condition(states['s2r'], V['s2r'])
            print 'The total enthalpy (Ht) at state 2r is {0:<.5g} MJ/kg (H2r - h1).'\
            .format((states['s2r_total'].h - states['s1'].h)/1.0e6) #(take away the initial enthalpy in state 1 to get the change)
            cfg['Ht2r'] = states['s2r_total'].h - states['s1'].h
        except Exception as e:
            print e
            print "Failed to find total condition for state 2r. Result will be set to None"
            cfg['Ht2r'] = None
            
    if PRINT_STATUS: print '-'*60
        
    return cfg, states, V, M  

#----------------------------------------------------------------------------
    
def morgan_diaphragm_analysis(cfg, states, V, M):
    """This is a thin diaphragm model from Morgan 1991
       "Double Diaphragm Driven Free Piston Expansion Tube.
       
       THIS IS CURRENTLY STILL BEING BUILT, DOES NOT WORK...
    """
    
    dt = 0.1
    diaphragm_mass = 0.191e3
    
    # model starts from a full reflection at the end of the shock tube, that then fades away
    states['s2r-mda'] = states['s2'].clone()
    
    cfg['Vr-mda'] = reflected_shock(states['s2'], V['s2'], states['s2r-mda'])
    cfg['Mr-mda'] = (V['s2']+cfg['Vr-mda'])/states['s2'].a #normally this would be V2 - V2r, but it's plus here as Vr has been left positive
    V['s2r-mda'] = 0.0
    M['s2r-mda']= V['s2r-mda']/states['s2r-mda'].a
    
    # now we start to step forward in time
    
    p2r_list = [states['s2r-mda'].p]
    Vr_list = [cfg['Vr-mda']]
    Mr_list = [cfg['Mr-mda']]
    Vd_list = [0.0]
    
    print "p2r = {0:.2f} Pa, Vr = {1:.2f} m/s, Mr = {2:.2f}, Vd = {3:.2f}."\
        .format(p2r_list[-1], Vr_list[-1], Mr_list[-1], Vd_list[-1])
       
    while Mr_list[-1] > 1:
        # first get diaphragm acceleration from pressure difference across diaphragm
        F = (p2r_list[-1] - states['s5'].p)*(math.pi*0.0425**2.0)
        a = F / diaphragm_mass #mass of a diaphragm on the bottom
        Vd = Vd_list[-1] + a*dt
        Vd_list.append(Vd)
        print "Vd = {0} m/s.".format(Vd)
        print '-'*60
        
        # now we need to find the Vr that makes V2 turn into Vd through the reflected shock
        # use a secant solver to make it happen
        def reflected_shock_speed_iterator(Vr, state2=states['s2'],
                                            state2r=states['s2r-mda'], V2g=V['s2'], Vd = Vd):    
            """iterates out the reflected shock speed that matches Vd"""
    
            print "Current guess for Vr = {0} m/s".format(Vr)
            
            print "-"*60            
            print "Mr = {0}.".format((V2g+Vr)/state2.a)
            
            try:
                V2r, V2rg = normal_shock(state2, V2g + Vr, state2r)
            except Exception as e:
                print "Error {0}".format(str(e))
                raise Exception, "Normal shock calculation failed with the current Vr guess." 
            
            #print "V2g = {0}, V3g = {1}, V3 = {2}.".format(V2g, V3g, V3)
            
            print "V2r = {0} m/s, V2rg = {1} m/s.".format(V2r, V2rg)
            
            V2r_lab = V2r - Vr
            
            velocity_difference = abs((V2r_lab-Vd)/V2r_lab)*100.0
            
            print "V2r_lab = {0}, Vd = {1}, difference = {2} %.".format(V2r_lab, Vd, velocity_difference)
            
            print "-"*60
    
            return (V2r_lab-Vd)/V2r_lab
        
        Vr = secant(reflected_shock_speed_iterator, Vr_list[-1], Vr_list[-1]-10.0, tol=1.0e-4,limits=[0.0,cfg['Vr-mda']])
        Vr_list.append(Vr)
        
        (V2r, V2rg) = normal_shock(states['s2'], Vr + V['s2'], states['s2r-mda'])
        
        Mr_list.append((V['s2'] + Vr)/states['s2'].a)
        
        p2r_list.append(states['s2r-mda'].p)
        
        print "p2r = {0:.2f} Pa, Vr = {1:.2f} m/s, Mr = {2:.2f}, Vd = {3:.2f}."\
        .format(p2r_list[-1], Vr_list[-1], Mr_list[-1], Vd_list[-1])
        print '-'*60
    
    return cfg, states, V, M
    
#----------------------------------------------------------------------------
    
def acceleration_tube_calculation(cfg, states, V, M):
    """Function that contains all of the acceleration tube calculations."""
            
    if PRINT_STATUS: print "Start unsteady expansion of the test gas into the acceleration tube."
    
    # need to know what our input state is to the acceleration tube, with a new
    # reflected shock at the end of the acc tube mode, it can change...
    
    if cfg['rs_out_of_st']: # state 2 after a reflected shock
        cfg['at_entry_state'] = 's2r'
    else: # this is just the normal one, state 2
        cfg['at_entry_state'] = 's2'

    #----------------------- acceleration tube functions -----------------------

    def error_in_velocity_s2_expansion_pressure_iterator(p5, state2=states[cfg['at_entry_state']], 
                                       V2g=V[cfg['at_entry_state']], state5=states['s5'],
                                        state6=states['s6'],Vs2=cfg['Vs2'],
                                        expansion_factor = cfg['expansion_factor'],
                                        ideal_gas_guess=cfg['gas_guess_air'],
                                        steps=cfg['acc_tube_expansion_steps']):
        """Compute the velocity mismatch for a given pressure ratio across the 
        unsteady expansion from state 2 to state 7."""
        
        print "current guess for p5 = {0} Pa".format(p5)
        
        state5.set_pT(p5,300.0) #set s1 at set pressure and ambient temp
        
        (V6, V6g) = normal_shock(state5, Vs2, state6,ideal_gas_guess)        
        
        #Across the contact surface, p3 == p2
        p7 = state6.p
        
        # Across the expansion, we get a velocity, V7g.
        V7g, state7 = finite_wave_dp('cplus', V2g, state2, p7, steps=steps)

        return (V6g - V7g)/V6g
        
    def error_in_pressure_s2_expansion_shock_speed_iterator(Vs2, state2=states[cfg['at_entry_state']], 
                                       V2g=V[cfg['at_entry_state']], state5=states['s5'],
                                        state6=states['s6'],
                                        high_temp_gas_guess=cfg['gas_guess_air'],
                                        steps=cfg['acc_tube_expansion_steps']):
        """Compute the velocity mismatch for a given shock speed in front of the 
        unsteady expansion from state 2 to state 7."""
        
        print '-'*60
        print "Current guess for Vs2 = {0} m/s".format(Vs2)
        
        # some extra code to try and get conditions above 19 km/s working with Pitot
        if Vs2 <= 19000.0 and cfg['solver'] in ['eq', 'pg-eq']:
            gas_guess = high_temp_gas_guess
        elif Vs2 > 19000.0 and cfg['solver'] in ['eq', 'pg-eq']:
            gas_guess = very_high_temp_gas_guess(Vs2)
        else:
            gas_guess = None
            
        momentum_change = 0.03375*1.0*100000.0
        
        (V6, V6g) = normal_shock(state5, Vs2, state6, gas_guess, momentum_change = momentum_change)
               
        #Across the contact surface, p3 == p2
        p7 = state6.p
                   
        # Across the expansion, we get a velocity, V7g.
        V7g, state7 = finite_wave_dp('cplus', V2g, state2, p7, steps=steps)
        
        print "Current p6 = {0} Pa, current p7 = {1} Pa.".format(state6.p, state7.p)
        print "Current V6g = {0} m/s, current V7g = {1} m/s.".format(V6g, V7g)      

        return (V6g - V7g)/V6g 
        
    def very_high_temp_gas_guess(Vs2):
        """A small function that will drop the gas guess further for the
           very high velocity gas giant entry conditions.
           
           Should only be used when Vs2 > 19,000 m/s.
        """
        
        print "Vs2 > 19 km/s a higher temperature gas guess gamma will be used."
        if Vs2 < 20500.0:
            gas_guess = {'gam':1.30, 'R':571.49}
        elif Vs2 < 21500.0:
            gas_guess = {'gam':1.26, 'R':571.49}
        elif Vs2 < 22500.0:
            gas_guess = {'gam':1.24, 'R':571.49}
        elif Vs2 < 23800.0:
            gas_guess = {'gam':1.22, 'R':571.49}
        elif Vs2 <= 24900.0:
            gas_guess = {'gam':1.20, 'R':571.49}
        elif Vs2 <= 26300.0:
            gas_guess = {'gam':1.18, 'R':571.49}
        elif Vs2 <= 27900.0:
            gas_guess = {'gam':1.16, 'R':571.49}
        elif Vs2 <= 29800.0:
            gas_guess = {'gam':1.14, 'R':571.49}
        elif Vs2 <= 32160.0:
            gas_guess = {'gam':1.12, 'R':571.49}
        elif Vs2 <= 34750.0:
            # this is the last value I could check to work, 
            # so anything above this I have made to code bail out on
            # Chris James (28/09/15)
            gas_guess = {'gam':1.10, 'R':571.49}
        elif Vs2 > 34750.0:
            print "No gas guess has been tested for this shock speed. Bailing out."
            raise Exception, "pitot_flow_functions.acceleration_tube_calculation() Run of pitot has failed in the acceleration tube calculation."
        
        return gas_guess
        
    #---------------------- acceleration tube calculations -----------------------
       
    if cfg['test'] == 'fulltheory-shock': #get p5 for our chosen shock speed
        #put two sets of limits here to try and make code that works will all conditions
        if cfg['Vs2'] > 4000.0 and 'p5_lower' not in cfg and 'p5_upper' not in cfg \
        and 'p5_guess_1' not in cfg and 'p5_guess_2' not in cfg:
            cfg['p5_lower'] = 1.0; cfg['p5_upper'] = 1000.0 #upper and lower limits (Pa)
            cfg['p5_guess_1'] = 10.0; cfg['p5_guess_2'] = 100.0 #first and second guesses for secant solver (Pa)
        elif cfg['Vs2'] <= 4000.0 and 'p5_lower' not in cfg and 'p5_upper' not in cfg \
        and 'p5_guess_1' not in cfg and 'p5_guess_2' not in cfg:
            cfg['p5_lower'] = 20.0; cfg['p5_upper'] = 1000.0 #upper and lower limits (Pa)
            cfg['p5_guess_1'] = 500.0; cfg['p5_guess_2'] = 1000.0 #first and second guesses for secant solver (Pa)
        cfg['p5'] = secant(error_in_velocity_s2_expansion_pressure_iterator, \
        cfg['p5_guess_1'], cfg['p5_guess_2'], tol=1.0e-5,limits=[cfg['p5_lower'],cfg['p5_upper']])
        if PRINT_STATUS: print "From secant solve: p5 = {0} Pa".format(cfg['p5'])
        
        #start using p5 now, compute states 5,6 and 7 using the correct p5
        if PRINT_STATUS: print "Now that p5 is known, finding conditions at state 5 and 6."
        states['s5'].set_pT(cfg['p5'],cfg['T0'])
        
    elif cfg['test'] == 'fulltheory-pressure' or cfg['test'] == 'fulltheory-pressure-ratios' \
    or cfg['test'] == 'experiment-shock-tube-theory-acc-tube': #compute the shock speed for the chosen fill pressure, uses Vs1 as starting guess
        #put two sets of limits here to try and make more stuff work
        if cfg['state7_no_ions']:
            # Need to turn ions off for state 2 here if it is required to make 
            # the unsteady expansion work (as state 2 is expanding into state 7)
            states['s2'].with_ions = False 
        if cfg['Vs1'] > 2000.0 and 'Vs2_lower' not in cfg and 'Vs2_upper' not in cfg:
            cfg['Vs2_lower'] = cfg['Vs1'] + 2000.0; cfg['Vs2_upper'] = 34750.0
        elif cfg['Vs1'] <= 2000.0 and 'Vs2_lower' not in cfg and 'Vs2_upper' not in cfg:
            cfg['Vs2_lower'] = cfg['Vs1'] + 1000.0; cfg['Vs2_upper'] = 34750.0
        if cfg['Vs1'] > 2000.0 and 'Vs2_guess_1' not in cfg and 'Vs2_guess_2' not in cfg:
            cfg['Vs2_guess_1'] = cfg['Vs1']+7000.0; cfg['Vs2_guess_2'] = cfg['Vs1']+8000.0
        elif cfg['Vs1'] <= 2000.0 and 'Vs2_guess_1' not in cfg and 'Vs2_guess_2' not in cfg:
            cfg['Vs2_guess_1'] = cfg['Vs1']+4000.0; cfg['Vs2_guess_2'] = cfg['Vs1']+5000.0
            
        cfg['acc_tube_secant_tol'] = 1.0e-5    
        cfg['Vs2'] = secant(error_in_pressure_s2_expansion_shock_speed_iterator, \
                            cfg['Vs2_guess_1'], cfg['Vs2_guess_2'], 
                            tol = cfg['acc_tube_secant_tol'],
                            limits=[cfg['Vs2_lower'],cfg['Vs2_upper']],
                                max_iterations = 100)
        if cfg['Vs2'] == 'FAIL':
            print "Acceleration tube secant solver did not converge after 100 iterations."
            raise Exception, "pitot_flow_functions.acceleration_tube_calculation() Run of pitot has failed in the acceleration tube calculation."                    
                  
        if PRINT_STATUS: 
            print '-'*60
            print "From secant solve: Vs2 = {0} m/s".format(cfg['Vs2'])     

        if PRINT_STATUS: 
            print '-'*60
            print "Now that Vs2 is known, finding conditions at state 6 and 7."
    elif cfg['test'] == 'experiment':
        if cfg['state7_no_ions']:
            # Need to turn ions off for state 2 here if it is required to make 
            # the unsteady expansion work (as state 2 is expanding into state 7)
            states['s2'].with_ions = False   
      
    # some extra code to try and get conditions above 19 km/s working with Pitot
    if cfg['Vs2'] <= 19000.0 and cfg['solver'] in ['eq', 'pg-eq']:
        gas_guess = cfg['gas_guess_air']
    elif cfg['Vs2'] > 19000.0 and cfg['solver'] in ['eq', 'pg-eq']:            
        gas_guess = very_high_temp_gas_guess(cfg['Vs2'])
    else:
        gas_guess = None
        
    if cfg['expand_to'] != 'p7':
        #if we're expanding to p7 we can't do this shock, so we skip it
        momentum_change = 0.03375*1.0*100000.0
        (V6, V['s6']) = normal_shock(states['s5'], cfg['Vs2'], states['s6'], gas_guess,  momentum_change = momentum_change)
        
    #do any modifications that were requested to the velocity behind the shock here 
    # new if statement here as we now have the ability to expand to a pressure if required - CMJ (16/09/15)
    if cfg['expand_to'] == 'flow-behind-shock' or cfg['expand_to'] == 'shock-speed':
        if cfg['expand_to'] == 'flow-behind-shock':
            print "State 7 is being expanded to V6 ({0}) multiplied by an expansion factor of {1}.".format(V['s6'], cfg['expansion_factor'])
            acc_tube_expand_to_V = V['s6']*cfg['expansion_factor']
        elif cfg['expand_to'] == 'shock-speed':
            acc_tube_expand_to_V = cfg['Vs2']*cfg['expansion_factor'] 
            print "State 7 is being expanded to the shock speed of Vs2 ({0} m/s) multiplied by an expansion factor of {1}."\
            .format(cfg['Vs2'], cfg['expansion_factor'])
        try:
            V['s7'], states['s7'] = finite_wave_dv('cplus', V[cfg['at_entry_state']], states[cfg['at_entry_state']], acc_tube_expand_to_V, steps=cfg['acc_tube_expansion_steps'])
        except Exception as e:
            print "Finding state7 failed. Trying again with 'state7_no_ions' turned on."
            cfg['state7_no_ions'] = True
            states['s2'].with_ions = False
            V['s7'], states['s7'] = finite_wave_dv('cplus', V[cfg['at_entry_state']], states[cfg['at_entry_state']], acc_tube_expand_to_V, steps=cfg['acc_tube_expansion_steps'])
    elif cfg['expand_to'] == 'p7':
        print "State 7 is being expanded to a specified p7 value of {0} Pa.".format(cfg['p7'])
        try:
            V['s7'], states['s7'] = finite_wave_dp('cplus', V[cfg['at_entry_state']], states[cfg['at_entry_state']], cfg['p7'], steps=cfg['acc_tube_expansion_steps'])
        except Exception as e:
            print "Finding state7 failed. Trying again with 'state7_no_ions' turned on."
            cfg['state7_no_ions'] = True
            states['s2'].with_ions = False
            V['s7'], states['s7'] = finite_wave_dp('cplus', V[cfg['at_entry_state']], states[cfg['at_entry_state']], cfg['p7'], steps=cfg['acc_tube_expansion_steps'])   
        
        cfg['Vs2'] = V['s7']
    if cfg['state7_no_ions']:
        # Turn with ions back on so it will be on for other states based on s7
        # if we turned it off to make the unsteady expansion work
        states['s7'].with_ions = True 
    
    #get mach numbers for the txt_output
    cfg['Ms2'] = cfg['Vs2']/states['s5'].a
    if not cfg['expand_to'] == 'p7': #no Vs2 or state 6 if we expand to a pressure to find state 7
        M['s6'] = V['s6']/states['s6'].a
    M['s7']= V['s7']/states['s7'].a
    
    if PRINT_STATUS:
        print '-'*60
        if not cfg['expand_to'] == 'p7': #no Vs2 or state 6 if we expand to a pressure to find state 7
            print "state 6: p = {0:.2f} Pa, T = {1:.2f} K, V = {2:.2f} m/s.".format(states['s6'].p, states['s6'].T,  V['s6']) 
        print "state 7: p = {0:.2f} Pa, T = {1:.2f} K. V = {2:.2f} m/s.".format(states['s7'].p, states['s7'].T,  V['s7'])
        if cfg['solver'] == 'eq' or cfg['solver'] == 'pg-eq':
            print 'species in state7 at equilibrium:'               
            print '{0}'.format(states['s7'].species)
    
    if PRINT_STATUS: print '-'*60
    
    return cfg, states, V, M
    
#----------------------------------------------------------------------------

def rs_calculation(cfg, states, V, M):
    """This is a small function that performs a reflected shock on state 2
    when pitot is used to simulate a nr-shock-tunnel.
    I basically wrote it as an easy way to simulate a rst in pitot without
    having to change too much of the code.
    You add in this calc by adding 'perform_rs = True' to the cfg file.
    """
    
    if cfg['tunnel_mode'] == 'expansion-tube':
        print "Reflected shock calculation was asked for, but nothing was done as tube is being simulated in expansion tube mode."
        
        return cfg, states, V, M
    
    #first build state 5 as a clone of state 2
    states['s5'] = states['s2'].clone()
    # then perform the reflected shock
    cfg['Vr'] = reflected_shock(states['s2'], V['s2'], states['s5'])
    cfg['Mr'] = (V['s2']+cfg['Vr'])/states['s2'].a #normally this would be V2 - Vr, but it's plus here as Vr has been left positive
    V['s5'] = 0.0
    M['s5']= V['s5']/states['s5'].a
    
    # now we do the same to state 3 so we can check the tailoring
    
    states['s5_d'] = states['s3'].clone()
    # then perform the reflected shock
    cfg['Vr_d'] = reflected_shock(states['s3'], V['s3'], states['s5_d'])
    cfg['Mr_d'] = (V['s3']+cfg['Vr_d'])/states['s3'].a #normally this would be V2 - Vr, but it's plus here as Vr has been left positive
    V['s5_d'] = 0.0
    M['s5_d']= V['s5_d']/states['s5_d'].a
    
    if PRINT_STATUS: 
        print "state 5: p = {0:.2f} Pa, T = {1:.2f} K.".format(states['s5'].p, states['s5'].T)
        if cfg['solver'] == 'eq' or cfg['solver'] == 'pg-eq':
            print '{0}'.format(states['s5'].species)
        print "state 5 driver ('5_d'): p = {0:.2f} Pa, T = {1:.2f} K.".format(states['s5_d'].p, states['s5_d'].T)
        
        
    return cfg, states, V, M
    
#----------------------------------------------------------------------------    
    
def test_section_setup(cfg, states, V, M):
    """Function to setup what is the test section state. It will also do the
    steady expansion through the nozzle if it detects that it needs to."""
                      
    if cfg['tunnel_mode'] == 'expansion-tube' and cfg['nozzle']:
        cfg['nozzle_entry_state'] = 's7'
        cfg, states, V, M == nozzle_expansion(cfg, states, V, M)
        cfg['test_section_state'] = 's8'
    elif cfg['tunnel_mode'] == 'expansion-tube' and not cfg['nozzle']:
        cfg['test_section_state'] = 's7' 
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel' and cfg['nozzle']:
        if PRINT_STATUS: print "Starting steady expansion through the nozzle."
        cfg['nozzle_entry_state'] = 's2'
        cfg, states, V, M == nozzle_expansion(cfg, states, V, M)
        cfg['test_section_state'] = 's8'
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel' and not cfg['nozzle']:            
        cfg['test_section_state'] = 's2'
    elif cfg['tunnel_mode'] == 'reflected-shock-tunnel' and cfg['nozzle']:
        if PRINT_STATUS: print "Starting steady expansion through the nozzle."
        cfg['nozzle_entry_state'] = 's5'
        cfg, states, V, M == nozzle_expansion(cfg, states, V, M)
        cfg['test_section_state'] = 's8'
    elif cfg['tunnel_mode'] == 'reflected-shock-tunnel' and not cfg['nozzle']:            
        cfg['test_section_state'] = 's5'        
            
    if PRINT_STATUS and cfg['nozzle']: print '-'*60
    
    return cfg, states, V, M

#----------------------------------------------------------------------------

def nozzle_expansion(cfg, states, V, M):
    """Function that just does the nozzle steady expansion.
    
    I built it to be used with the changing area ratio code.
    
    """

    if PRINT_STATUS: print "Starting steady expansion through the nozzle."
    try:
        (V['s8'], states['s8']) = steady_flow_with_area_change(states[cfg['nozzle_entry_state']], V[cfg['nozzle_entry_state']],cfg['area_ratio'])
        M['s8']= V['s8']/states['s8'].a
    except Exception as e:
        print "Nozzle expansion failed. Going to try again with no ions."
        try:
            states['s7'].with_ions = False
            (V['s8'], states['s8']) = steady_flow_with_area_change(states[cfg['nozzle_entry_state']], V[cfg['nozzle_entry_state']], cfg['area_ratio'])
            M['s8']= V['s8']/states['s8'].a
            states['s8'].with_ions = True
            print "Nozzle expansion sucessful with ions turned off."
        except Exception as e:
            print "Error {0}".format(str(e))
            raise Exception, "pitot_flow_functions.nozzle_expansion(): Run of pitot failed in the nozzle expansion calculation."
    
    print "state 8: p = {0:.2f} Pa, T = {1:.2f} K.".format(states['s8'].p, states['s8'].T)        
    if cfg['solver'] == 'eq' or cfg['solver'] == 'pg-eq':
        print 'species in state8 at equilibrium:'               
        print '{0}'.format(states['s8'].species)
    
    return cfg, states, V, M     
    
#----------------------------------------------------------------------------
    
def shock_over_model_calculation(cfg, states, V, M):
    """Function that takes the cfg, states, V and M dictionaries
    and does the shock over model calculations if the user requests it.
    The changed dictionaries are then returned.
    
    """
    
    if PRINT_STATUS: print "Starting frozen normal shock calculation over the test model."
    
    try:
        states['s10f'] = states[cfg['test_section_state']].clone()
    except Exception as e:
        print "Error {0}".format(str(e))
        print "Failed to clone test section gas state."
        print "Trying again with ions turned off."
        states[cfg['test_section_state']].with_ions = False
        try:
            states['s10f'] = states[cfg['test_section_state']].clone() 
            states[cfg['test_section_state']].with_ions = True
            print "Managed to clone test section state with ions turned off."
        except Exception as e:
            print "Error {0}".format(str(e))
            print "Failed to clone test section gas state."
            states[cfg['test_section_state']].with_ions = True
            if 's10f' in states.keys():
                del states['s10f']
            print "Frozen normal shock calculation over the test model failed."
            print "Result will not be printed."                
    
    if 's10f' in states.keys():
        try:
            (V10, V['s10f']) = shock_ideal(states[cfg['test_section_state']], V[cfg['test_section_state']], states['s10f'])
            M['s10f']= V['s10f']/states['s10f'].a
        except Exception as e:
            print "Error {0}".format(str(e))
            print "Frozen normal shock calculation over the test model failed."
            print "Result will not be printed."
            if 's10f' in states.keys():
                del states['s10f']
                
    if PRINT_STATUS: print '-'*60                
    if PRINT_STATUS: print "Starting equilibrium normal shock calculation over the test model."  
    
    if cfg['solver'] == 'eq' or cfg['solver'] == 'pg-eq': 
        try:
            states['s10e'] = states[cfg['test_section_state']].clone()
            states['s10e'].with_ions = True
        except Exception as e:
            print "Error {0}".format(str(e))
            print "Failed to clone test section gas state."
            print "Trying again with ions turned off."
            states[cfg['test_section_state']].with_ions = False
            try:
                states['s10e'] = states[cfg['test_section_state']].clone() 
                states[cfg['test_section_state']].with_ions = True
                states['s10e'].with_ions = True
                print "Managed to clone test section state with ions turned off."
            except Exception as e:
                print "Error {0}".format(str(e))
                print "Failed to clone test section gas state."
                states[cfg['test_section_state']].with_ions = True
                if 's10e' in states.keys():
                    del states['s10e']
                print "Equilibrium normal shock calculation over the test model failed."
                print "Result will not be printed." 
        if 's10e' in states.keys():
            try:
                (V10, V['s10e']) = normal_shock(states[cfg['test_section_state']], V[cfg['test_section_state']], states['s10e'])
                M['s10e']= V['s10e']/states['s10e'].a
                print states['s10e'].species
                if abs((states['s10e'].p - states[cfg['test_section_state']].p) / states['s10e'].p) < 0.10:
                    print "For some reason p10e and p{0} are too similar. Shock must have not occured properly.".format(cfg['test_section_state'][1])
                    print "p{0} = {1} Pa, p10e = {2} Pa."\
                    .format(cfg['test_section_state'][1], states[cfg['test_section_state']].p, states['s10e'].p)
                    raise Exception, "pitot_flow_functions.shock_over_model_calculation() Eq shock over model calculation failed."
            except Exception as e:
                print "Error {0}".format(str(e))
                print "Equilibrium normal shock calculation over the test model failed."
                if cfg['gas_guess']: 
                    print "Will try again with a high temperature gas guess."
                    try:
                        (V10, V['s10e']) = normal_shock(states[cfg['test_section_state']], 
                                                        V[cfg['test_section_state']], 
                                                        states['s10e'], cfg['gas_guess'])
                        M['s10e']= V['s10e']/states['s10e'].a
                        if abs((states['s10e'].p - states[cfg['test_section_state']].p) / states['s10e'].p) < 0.10:
                            print "For some reason p10e and p{0} are too similar. Shock must have not occured properly.".format(cfg['test_section_state'][1])
                            print "p{0} = {1} Pa, p10e = {2} Pa."\
                            .format(cfg['test_section_state'][1], states[cfg['test_section_state']].p, states['s10e'].p)
                            raise Exception, "pitot_flow_functions.shock_over_model_calculation() Eq shock over model calculation failed."
                    except Exception as e:
                        print "Error {0}".format(str(e))
                        print "Result will not be printed."
                        if 's10e' in states.keys():
                            del states['s10e']
                else:
                    print "Result will not be printed."
                    if 's10e' in states.keys():
                        del states['s10e']
                
    elif cfg['solver'] == 'pg': #we need to make a cea2 gas object to do this equilibrium calculaiton if every other gas object is pg
        if cfg['test_gas'] != 'custom':
            states[cfg['test_section_state']+'eq'] = make_test_gas(cfg['test_gas'])[0]
        else: # need to do something slightly different if we have a custom test gas
            states[cfg['test_section_state']+'eq'] = Gas(cfg['test_gas_composition'],
                                                     inputUnits=cfg['test_gas_inputUnits'], 
                                                     outputUnits='moles', with_ions=True)                                       
        states[cfg['test_section_state']+'eq'].set_pT(states[cfg['test_section_state']].p,states[cfg['test_section_state']].T)
        states['s10e'] = states[cfg['test_section_state']+'eq'].clone()
        try:
            (V10, V['s10e']) = normal_shock(states[cfg['test_section_state']+'eq'], V[cfg['test_section_state']], states['s10e'])
            M['s10e']= V['s10e']/states['s10e'].a
        except Exception as e:
            print "Error {0}".format(str(e))
            print "Equilibrium normal shock calculation over the test model failed."
            print "Result will not be printed."
            if 's10e' in states.keys():
                del states['s10e']
        
    if PRINT_STATUS: print '-'*60
                
    return cfg, states, V, M
    
#----------------------------------------------------------------------------
    
def wedge_calculation(cfg, states, V, M):
    """Function that takes the cfg, states, V and M dictionaries
    and does a wedge calculation if the user requests it.
    It takes the wedge angle from the config dictionary and builds a new
    state, state 10w.
    
    """
    
    if PRINT_STATUS: print "Starting calculation of conditions behind a {0} degree wedge in the test section.".format(cfg['wedge_angle']) 
    
    #--------------------------------------------------------------------------
    print "Starting frozen wedge calculation."
    
    try:
        # we take the MM and gamma of the test section state and use that to specify
        # a perfect gas version of the test section state to work with
        pg_test_section_state = pg.Gas(Mmass=states[cfg['test_section_state']].Mmass,
                                gamma=states[cfg['test_section_state']].gam, name='pg_test_section_state')
        pg_test_section_state.set_pT(states[cfg['test_section_state']].p, 
                                     states[cfg['test_section_state']].T)
        states['s10w'] = pg_test_section_state.clone()
    except Exception as e:
        print "Error {0}".format(str(e))
        print "Failed to make perfect gas test section state."
        if 's10w' in states.keys():
            del states['s10w']
            print "Frozen wedge calculation failed."
            print "Result will not be printed."
            
    if 's10w' in states.keys(): 
        # start by getting the beta angle over the wedge
        cfg['wedge_angle_radians'] = math.radians(cfg['wedge_angle'])
        try:
            cfg['beta_pg'] = beta_oblique(pg_test_section_state, V[cfg['test_section_state']], cfg['wedge_angle_radians'])
            print "Beta = {0} degrees.".format(math.degrees(cfg['beta_pg']))
            # now get the wedge surface conditions
            cfg['wedge_angle_calculated'], V['s10w'], states['s10w'] = theta_oblique(pg_test_section_state, V[cfg['test_section_state']], cfg['beta_pg'])
            wedge_angle_calculated_degrees = math.degrees(cfg['wedge_angle_calculated'])
            wedge_angle_error = ((cfg['wedge_angle']-wedge_angle_calculated_degrees)/cfg['wedge_angle'])*100.0
            # need to check the calculated theta
            print "Selected wedge angle is {0} degrees, calculated wedge angle is {1} degrees, error is {2} %"\
                  .format(cfg['wedge_angle'], wedge_angle_calculated_degrees, wedge_angle_error)
            if wedge_angle_error > 1.0:
                print "Wedge angle error is too large. Going to throw this result out."
                if 's10w' in states.keys():
                    del states['s10w']
            if 's10w' in states.keys():        
                M['s10w']= V['s10w']/states['s10w'].a
        except Exception as e:
            print "Error {0}".format(str(e))
            print "Frozen wedge calculation failed."
            print "Result will not be printed."
            if 's10w' in states.keys():
                del states['s10w']
            
    if PRINT_STATUS: print '-'*60
        
    #--------------------------------------------------------------------------
    # now do the eq wedge calc if we're in eq mode
        
    if cfg['solver'] == 'eq' or cfg['solver'] == 'pg-eq': 
        print "Starting equilibrium wedge calculation."
        if 's10w' in states.keys():
            # first give the frozen calc frozen notation
            states['s10wf'] = states['s10w']
            V['s10wf'] = V['s10w']
            M['s10wf'] = M['s10w']
            
        try:
            states['s10we'] = states[cfg['test_section_state']].clone()
        except Exception as e:
            print "Error {0}".format(str(e))
            print "Failed to clone test section gas state."
            print "Trying again with ions turned off."
            states[cfg['test_section_state']].with_ions = False
            try:
                states['s10we'] = states[cfg['test_section_state']].clone() 
                states[cfg['test_section_state']].with_ions = True
                states['s10we'].with_ions = True
                print "Managed to clone test section state with ions turned off."
            except Exception as e:
                print "Error {0}".format(str(e))
                print "Failed to clone test section gas state."
                states[cfg['test_section_state']].with_ions = True
                if 's10we' in states.keys():
                    del states['s10we']
                print "Equilibrium wedge calculation failed."
                print "Result will not be printed." 
                if 's10we' in states.keys():   
                    del states['s10we'] 
            
        if 's10we' in states.keys():
            # start by getting the beta angle over the wedge
            cfg['wedge_angle_radians'] = math.radians(cfg['wedge_angle'])
            try:
                cfg['beta_eq'] = beta_oblique(states[cfg['test_section_state']], V[cfg['test_section_state']], cfg['wedge_angle_radians'])
                print "Beta = {0} degrees.".format(math.degrees(cfg['beta_eq']))
                # now get the wedge surface conditions
                cfg['wedge_angle_calculated'], V['s10we'], states['s10we'] = theta_oblique(states[cfg['test_section_state']], V[cfg['test_section_state']], cfg['beta_eq'])
                wedge_angle_calculated_degrees = math.degrees(cfg['wedge_angle_calculated'])
                wedge_angle_error = ((cfg['wedge_angle']-wedge_angle_calculated_degrees)/cfg['wedge_angle'])*100.0
                # need to check the calculated theta
                print "Selected wedge angle is {0} degrees, calculated wedge angle is {1} degrees, error is {2} %"\
                .format(cfg['wedge_angle'], wedge_angle_calculated_degrees, wedge_angle_error)
                if wedge_angle_error > 1.0:
                    print "Wedge angle error is too large. Going to throw this result out."
                    if 's10we' in states.keys():
                        del states['s10we']
                if 's10we' in states.keys(): 
                    M['s10we']= V['s10we']/states['s10we'].a
            except Exception as e:
                print "Error {0}".format(str(e))
                print "Equilibrium wedge calculation failed."
                print "Result will not be printed."
                if 's10we' in states.keys():
                    del states['s10we']
                
    return cfg, states, V, M
    
#----------------------------------------------------------------------------
    
def conehead_calculation(cfg, states, V, M):
    """Function that takes the cfg, states, V and M dictionaries, does a 
    calculation for a conehead in the test section at a specified angle
    and then returns the changes cfg, states, V, M dictionaries."""
     
    cfg['conehead_completed'] = True #variable we'll set to false if the calculation fails
    if 'conehead_no_ions' not in cfg: #add this variable and set it to false if the user has not used it
        cfg['conehead_no_ions'] = False
    cfg['conehead_no_ions_required'] = False #variable we'll set to True if we need to turn ions off in the conehead calc 
    
    if PRINT_STATUS: print 'Starting taylor maccoll conehead calculation on {0} degree conehead.'.format(cfg['conehead_angle'])
    
    if cfg['conehead_no_ions']: states[cfg['test_section_state']].with_ions = False    
    
    try:
        shock_angle = beta_cone(states[cfg['test_section_state']], V[cfg['test_section_state']], math.radians(cfg['conehead_angle']))
    except Exception as e:
        print "Error {0}".format(str(e))
        print "beta_cone function bailed out while trying to find a shock angle."
        print "will try again with ions turned off in the calculation."
        cfg['conehead_no_ions_required'] = True
        cfg['conehead_no_ions'] = True        
        states[cfg['test_section_state']].with_ions = False
        try:
            shock_angle = beta_cone(states[cfg['test_section_state']], V[cfg['test_section_state']], math.radians(cfg['conehead_angle']))
        except Exception as e:
            print "Error {0}".format(str(e))    
            print "beta_cone function bailed out while trying to find a shock angle."
            print "Stopping here. Try another nozzle area ratio."
            print "Result will not show state 10c."
            cfg['conehead_completed'] = False
            return cfg, states, V, M  
            
            
    if PRINT_STATUS: print "Shock angle over cone:", math.degrees(shock_angle)
    # Reverse the process to get the flow state behind the shock and check the surface angle is correct
    try:    
        delta_s, V['s10c'], states['s10c'] = theta_cone(states[cfg['test_section_state']], V[cfg['test_section_state']], shock_angle)
    except Exception as e:
        print "Error {0}".format(str(e))
        print "theta_cone bailed out while trying to find cone surface conditions."
        print "Stopping here. Try another nozzle area ratio."
        print "Result will not show state 10c."
        cfg['conehead_completed'] = False
        return cfg, states, V, M
        
    M['s10c'] = V['s10c']/states['s10c'].a
    if PRINT_STATUS: print "Surface angle should be the same.....: 15deg = ", math.degrees(delta_s), "deg"
    #if PRINT_STATUS: print "\nConehead surface conditions:"
    #if PRINT_STATUS: states['s10c'].write_state(sys.stdout)
    # Need to check whether the pressure are the same
    if PRINT_STATUS: print "Computed conehead pressure is {0} Pa".format(states['s10c'].p)
    
    # turn the ions back on at the end
    if cfg['conehead_no_ions']: states[cfg['test_section_state']].with_ions = True
    if cfg['conehead_no_ions_required']: cfg['conehead_no_ions'] = False
    
    if PRINT_STATUS: print '-'*60
    
    return cfg, states, V, M
    
#----------------------------------------------------------------------------
        
def test_time_calculator(cfg, states, V):
    """Function that takes the cfg, states and V dictionaries and will return
    the cfg dictionary with basic test time calculated if it is able to."""
    
    #This was based off code written by David Gildfind that was based off a paper
    #by Allan Paull and Ray Stalker. I've only considered the basic test time case here
    # A. Paull & R. J. Stalker, "Test Flow Disturbances in an Expansion Tube", J. Fluid Mech. (1992), vol. 245, pp. 493-521 (p499).
    
    cfg['calculate_test_time'] = False #defaults to false unless we have enough info to calculate it
    
    if cfg['facility'] == 'x2' and cfg['tunnel_mode'] == 'expansion-tube':
        #a few different tunnel scenarios give different lengths
        #all the distances are taken from my L1D run file, written by David Gildfind
        #0m is the primary diaphragm burst location
        cfg['calculate_test_time'] = True
        if not cfg['secondary'] and not cfg['nozzle']:
            distances = [3.418, 8.979] #secondary diaphragm, and then end of acc tube (m)
        elif cfg['secondary'] and not cfg['nozzle']:
            distances = [3.418, 5.976, 8.979] #secondary diaphragm, tertiary diaphragm, end of acc tube (m)
        elif not cfg['secondary'] and cfg['nozzle']:
            distances = [3.418, 8.585] #secondary diaphragm, entrance to nozzle (m)
        elif cfg['secondary'] and cfg['nozzle']:
            distances = [3.418, 5.976, 8.585] #secondary diaphragm, tertiary, entrance to nozzle (m)
    
        t_start = 0.0 #time of primary diaphragm burst
    
        if cfg['secondary']: #if secondary we have a third tube to consider
            #calculate some lengths
            L_sec_drv = distances[0]
            L_shk_tube = distances[1] - distances[0]
            L_acc_tube = distances[2] - distances[1]
            
            #shocks
            t_inc_shock_sd = t_start + L_sec_drv/cfg['Vsd']
            t_inc_shock_st = t_inc_shock_sd + L_shk_tube/cfg['Vs1']
            t_inc_shock_at = t_inc_shock_st + L_acc_tube/cfg['Vs2']
            
            #contact surfaces
            t_cs_sd = t_start + L_sec_drv/V['sd3']
            t_cs_st = t_inc_shock_sd + L_shk_tube/V['s3']
            t_cs_at = t_inc_shock_st + L_acc_tube/V['s7']

        
        else: #so we're going straight from primary driver to shock tube now
            #calculate some lengths
            L_shk_tube = distances[0]
            L_acc_tube = distances[1] - distances[0]
            
            #shocks
            t_inc_shock_st = t_start + L_shk_tube/cfg['Vs1']
            t_inc_shock_at = t_inc_shock_st +  L_acc_tube/cfg['Vs2']
            
            #contact surfaces
            t_cs_sd = t_start + L_shk_tube/V['s3']
            t_cs_at = t_inc_shock_st +  L_acc_tube/V['s7']
            
        
        #now we can actually calculate the basic test time...
        # (borrowed from the procedure David Gildfind used in his PhD)
        
        #now we just need to calculate the time taken for the unsteady 
        #expansion to hit the end of the acc tube
        
        t_final_usx = t_inc_shock_st + L_acc_tube/(V['s7']-states['s7'].a)
        
        cfg['t_test_basic'] = t_final_usx - t_cs_at
        
    return cfg

#! /usr/bin/env python 
"""
pitot.py: Equilibrium expansion tube condition builder

ADD INTRO TO THE PROGRAM

Getting the program set up
--------------------------
pitot.py is not a stand-alone file.
It comes as part of the cfcfd3 compressible-flow collection and
depends upon functions from the cfpylib library to do the specific 
calculations.

Someday I'll make a proper makefile for this guy, but for now,
just copy pitot.py into your e3bin.

You may then call upon pitot.py so long as you have suitable
enviroment variables set, as per the installation instructions
for Eilmer3.


Some History
------------
Since the dawn of the time, Richard Morgan used the program pitot to
do basic calculations and create new expansion tube conditions.
The code was originally written in GWBasic, and as an exercise in 
programming (and to hopefully make the life of my colleagues and I
easier) I learnt some basic GWBasic syntax and ported the program to
Python at the start of 2012. It has since become an ongoing project
and this latest version has completely rebuilt the code from the 
ground up as an new program that makes use of the classes and functions
available to me as part of cfpylib inside the cfcfd code collection.

.. Author: Chris James (c.james4@uq.edu.au)
   The Centre for Hypersonics
   The University of QLD, St Lucia.

.. Versions:
   when dinosaurs were still around: RGM wrote the original code in GWBasic
   Summer 2011 - 2012: ported to Python as a direct copy by myself
   first half of 2012: added new parts to the code to make it more friendly
   October 2012: started full scale rebuild of code as an exercise in expansion
       tubes for myself at the start of my postgrad, and to build a program
       that is more useful for myself and others
"""

#--------------------- intro stuff --------------------------------------

#this import pulled from estcj with some additions

import sys, os, math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from cfpylib.nm.zero_solvers import secant
# We base our calculation of gas properties upon calls to the NASA Glenn CEA code.
from cfpylib.gasdyn.cea2_gas import Gas, make_gas_from_name
from cfpylib.gasdyn.gas_flow import *
from cfpylib.gasdyn.ideal_gas_flow import *

VERSION_STRING = "19-Oct-2012"

DEBUG_PITOT = False

PRINT_STATUS = 1 #if print status is 1, some basic printouts are done

#define any functions needed
      
def is_valid(command,valid_commands):
    """Prompts the user for a valid command, and only returns a valid command.

    is_valid_command(float,tuple<floats>) -> float

    Preconditions: input must be a single number."""
    
    check = False
    
    while check == False:
        for i in valid_commands:
            if i == command:
                check = True
    
        if check == False:
             print 'That is an invalid command. Valid commands are {0}'.format(valid_commands)
             command = str(raw_input('Try again? '))             
    return command
    
def make_test_gas(gasName, outputUnits='massf'):
    """
    Manufacture a Gas object for the test gas from a small library of options.
    
    Also has workable gamma's and R's for the test gases at high temperature stored
    in a dictionary of the form {'gam' : gam, 'R' : R}.
    
    I stole this from the make_gas_from_name function in cea2_gas.py and added
    my own gases - Chris James

    :param gasName: one of the names for the special cases set out below
    """
    if gasName.lower() == 'air':
        return Gas({'Air':1.0,}, outputUnits=outputUnits, trace=1.0e-4), {'gam':1.35,'R':571.49}
    elif gasName.lower() == 'air5species':
        return Gas(reactants={'N2':0.79, 'O2':0.21, 'N':0.0, 'O':0.0, 'NO':0.0}, 
                   inputUnits='moles', onlyList=['N2','O2','N','O','NO'],
                   outputUnits=outputUnits), {'gam':1.35,'R':571.49}
    elif gasName.lower() == 'n2':
        return Gas(reactants={'N2':1.0, 'N':0.0}, onlyList=['N2', 'N'],
                   outputUnits=outputUnits),{'gam': 1.36,'R':593.56}
    elif gasName.lower() == 'titan':
        return Gas(reactants={'N2':0.95, 'CH4':0.05}, inputUnits='moles',
                   outputUnits=outputUnits), {'gam':1.35,'R':652.03}      
    elif gasName.lower() == 'mars':
        return Gas(reactants={'CO2':0.96, 'N2':0.04}, inputUnits='moles',
                   outputUnits=outputUnits)    
    elif gasName.lower() == 'co2':
        return Gas(reactants={'CO2':1.0}, outputUnits=outputUnits)
    elif gasName.lower() == 'gasgiant_h2ne':
        return Gas(reactants={'H2':0.85, 'Ne':0.15}, inputUnits='moles',
                   outputUnits=outputUnits),{'gam':1.15,'R':3245.1}
    elif gasName.lower() == 'gasgiant_h2he':
        return Gas(reactants={'H2':0.85, 'He':0.15}, inputUnits='moles',
                   outputUnits=outputUnits), None
    else:
        raise Exception, 'make_gas_from_name(): unknown gasName: %s' % gasName               
    
#--------------------------- property dictionaries -------------------------

#primary driver conditions, sorted into a dictionary with a gas object at the
#right conditions, and then mach number at the change over from steady to
#unsteady expansion, this was based on calcs done by RGM

primary_driver = dict([('He:1.0', [Gas({'He':1.0},inputUnits='moles'),2.15]),
                       ('He:0.80,Ar:0.20',[Gas({'He':0.8,'Ar':0.2},inputUnits='moles'),1]),
                        ('He:0.90,Ar:0.10',[Gas({'He':0.9,'Ar':0.1},inputUnits='moles'),1.59])])
  
def main():

    import optparse
    
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--config', dest='config', default='x2-sec-nozzle',
                  choices=['x2','x2-sec-nozzle','x2-sec','x2-nozzle'],
                  help=("tunnel configuration to use. "
                        "currently only set up for X2. "
                        "x2 = x2 with no secondary driver, no nozzle; "
                        "x2-sec = x2 with secondary driver, no nozzle; "
                        "x2-nozzle = x2 with no secondary driver, nozzle; "
                        "x2-sec-nozzle = x2 with secondary driver, nozzle "))
    op.add_option('--driver_gas', dest='driver_gas', default='He:1.0',
                  choices=['He:0.80,Ar:0.20', 'He:0.90,Ar:0.10','He:1.0'],
                  help=("driver gas configuration: "
                        "'He:0.80,Ar:0.20'; " "'He:0.90,Ar:0.10'; " "'He:1.0' "
                        "default is currently 100% He"))
    op.add_option('--test_gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'n2', 'titan', 
                           'gasgiant_h2ne', 'gasgiant_h2he'],
                  help=("name of test gas: "
                        "air; " "air5species; " "n2; " "titan; " "gasgiant_h2ne; "
                        "gasgiant_h2he" "default is air"))
    op.add_option('--Vs1', dest='Vs1', type='float', default=None,
                  help=("first shock speed, in m/s"))
    op.add_option('--Vs2', dest='Vs2', type='float', default=None,
                  help=("second shock speed, in m/s"))
    op.add_option('--Vs3', dest='Vs3', type='float', default=None,
                  help=("third shock speed, in m/s "
                        "not needed if secondary driver isn't used"))   
    op.add_option('--filename',dest='filename',type='string',default='x2run.txt',
                  help=("filename the result will be saved to "
                          "defaults to x2run.txt"))
    
    opt, args = op.parse_args()
    
    config = opt.config
    driver_gas = opt.driver_gas
    gasName = opt.gasName
    filename =opt.filename
    
    Vs1 = opt.Vs1
    Vs2 = opt.Vs2
    Vs3 = opt.Vs3
    
    bad_input = False
    
    if Vs1 is None:
        print "Need to supply a float value for Vs1."
        bad_input = True
        
    if Vs2 is None:
        print "Need to supply a float value for Vs2."
        bad_input = True
    
    if config == 'x2-sec-nozzle':
        if Vs3 is None:
            print "Need to supply a float value for Vs3."
            bad_input = True
        secondary = True
        nozzle = True    
    elif config == 'x2-sec':
        if Vs3 is None:
            print "Need to supply a float value for Vs3."
            bad_input = True
        secondary = True
        nozzle = False    
    elif config == 'x2-nozzle':
        secondary = False
        nozzle = True
    else:
        secondary = False
        nozzle = False
    
    if bad_input: #bail out here if you end up having issues with your input
        return -2
        
    if PRINT_STATUS: print "Let's get started, shall we."
    if PRINT_STATUS: 
        print 'selected Vs1 = {0} m/s'.format(Vs1)
        print 'selected Vs2 = {0} m/s'.format(Vs2)
        if secondary:
            print 'selected Vs3 = {0} m/s'.format(Vs3)
    
    
    #CURRENTLY RUNNING WITHOUT ANY OF THE REFLECTED SHOCK TYPE STUFF 
            
    states = {} #states dictionary that we'll fill up later
    V = {} #same for velocity
    M = {} #same for Mach number
    
    #a couple of switches
    
    ref = 0 #ref = 1 allows shock reflection at D2
    stand = 1 #if stand = 1, calculate pressure for standing shock in region 2
    REV = 0
    
    #initial condition I want to test will be Umar's condition
    #(100%He driver, no secondary driver)
    
    #state 4, piston burst pressure
    states['s4']=primary_driver[driver_gas][0]
    states['s4'].set_pT(2.79e7,2700.0)
    V['s4']=0.0
    M['s4']=0.0
    
    M['s11']=primary_driver['He:0.90,Ar:0.10'][1]

    #state11, piston after steady expansion at throat of driver
    
    if PRINT_STATUS: print "Start steady expansion of driver gas."     
    
    (states['s11'], V['s11']) = expand_from_stagnation(1.0/(p0_p(M['s11'],states['s4'].gam)),states['s4'])
    
    #build the gas objects for all the other sections based on knowledge of what is what
    
    T0 = 300.0 #atmospheric temperature (K), used for any section starting at ambient
    p0 = 101300.0 #atmospheric p (Pa)
    
    #start with shock tube and acc tube, start with atmospheric p and T
    
    if secondary: #state 1 is pure He secondary driver
        states['s1'] =  Gas({'He':1.0,})
        states['s1'].set_pT(p0,T0)
        V['s1']=0.0
        M['s1']=0.0
    else: #state 1 is shock tube
        states['s1'], gas_guess = make_test_gas(gasName)
        states['s1'].set_pT(p0,T0)
        V['s1']=0.0
        M['s1']=0.0
        
    if secondary: #state 5 is shock tube
        states['s5'], gas_guess = make_test_gas(gasName)
        states['s5'].set_pT(p0,T0)
        V['s5']=0.0
        M['s5']=0.0
    else: #state 5 is acceleration tube
        states['s5'] = Gas({'Air':1.0,})
        states['s5'].set_pT(p0,T0)
        V['s5']=0.0
        M['s5']=0.0
        
    if secondary: #state 8 is acceleration tube
        states['s8'] = Gas({'Air':1.0,})
        states['s8'].set_pT(p0,T0)
        V['s8']=0.0
        M['s8']=0.0
    
    #now let's clone for the states derived from these
    
    states['s3'] = states['s11'].clone() #3 will be 11 after unsteady expansion
    states['s2'] = states['s1'].clone() #2 is 1 shock heated
    states['s6'] = states['s5'].clone() #6 is 5 shock heated
    states['s7'] = states['s2'].clone() #7 is 2 after unsteady expansion
    
    if secondary:
        states['s9'] = states['s8'].clone() #9 is 8 shock heated
        states['s10'] = states['s6'].clone() #10 is 6 after unsteady expansion
      
    arbitrary = True #dummy variable
    
    while arbitrary:
        
        #stole some code from PJ and am modifying it for my needs
        
        print "Start unsteady expansion of the driver gas."
        # For the unsteady expansion of the test driver into the tube, regulation of the amount
        # of expansion is determined by the shock-processed gas in the next section.
        # Across the contact surface between these gases, the pressure and velocity
        # have to match so we set up some trials of various pressures and check 
        # that velocities match.
        
        def error_in_velocity_s11_driver_expansion(p1, state11=states['s11'], 
                                                   V11g=V['s11'], state1=states['s1'],
                                                    state2=states['s2'],Vs1=Vs1):
            """Compute the velocity mismatch for a given pressure ratio across the 
            unsteady expansion in the driver."""
            
            print "current guess for p1 = {0} Pa".format(p1)            
            
            state1.set_pT(p1,300.0) #set s1 at set pressure and ambient temp
            
            (V2, V2g) = normal_shock(state1, Vs1, state2)
            
            #Across the contact surface, p3 == p2
            p3 = state2.p
            
            # Across the expansion, we get a velocity, V5g.
            V3g, state3 = finite_wave_dp('cplus', V11g, state11, p3)

            return (V2g - V3g)/V2g        
        
        #get p1 for our chosen shock speed
        p1 = secant(error_in_velocity_s11_driver_expansion, 1000.0, 100000.0, tol=1.0e-3,limits=[100.0,1000000.0])
        if PRINT_STATUS: print "From secant solve: p1 = {0} Pa".format(p1)
        
        #start using p1 now, compute states 1,2 and 3 using the correct p1
        if PRINT_STATUS: print "Once p1 is known, find conditions at state 1 and 2."
        states['s1'].set_pT(p1,T0)
        (V2, V['s2']) = normal_shock(states['s1'], Vs1, states['s2'])
        V['s3'], states['s3'] = finite_wave_dv('cplus', V['s11'], states['s11'], V['s2'])
        
        #get mach numbers for the output
        Ms1 = Vs1/states['s1'].son
        M['s2'] = V['s2']/states['s2'].son
        M['s3']= V['s3']/states['s3'].son
        
        #do ideal gas guess for gas constants if needed
        
        if Vs2 >= 9000.0:
            ideal_gas_guess = gas_guess
        else:
            ideal_gas_guess = None
            
        if secondary: 
            if PRINT_STATUS: print "Start unsteady expansion of the secondary driver gas."
        else: 
            if PRINT_STATUS: print "Start unsteady expansion of the test gas."

        def error_in_velocity_s2_expansion(p5, state2=states['s2'], 
                                           V2g=V['s2'], state5=states['s5'],
                                            state6=states['s6'],Vs2=Vs2,
                                            ideal_gas_guess=ideal_gas_guess):
            """Compute the velocity mismatch for a given pressure ratio across the 
            unsteady expansion from state 2 to state 7."""
            
            print "current guess for p5 = {0} Pa".format(p5)
            
            state5.set_pT(p5,300.0) #set s1 at set pressure and ambient temp
            
            (V6, V6g) = normal_shock(state5, Vs2, state6,ideal_gas_guess)
            
            #Across the contact surface, p3 == p2
            p7 = state6.p
            
            # Across the expansion, we get a velocity, V7g.
            V7g, state7 = finite_wave_dp('cplus', V2g, state2, p7)

            return (V6g - V7g)/V6g        
        
        #get p5 for our chosen shock speed
        if secondary: #do a wider guess here        
            p5 = secant(error_in_velocity_s2_expansion, 50.0, 100000.0, tol=1.0e-3,limits=[50.0,100000.0])
        else: #don't bother
            p5 = secant(error_in_velocity_s2_expansion, 10.0, 1000.0, tol=1.0e-3,limits=[5.0,1000.0])
        if PRINT_STATUS: print "From secant solve: p5 = {0} Pa".format(p5)
        
        #start using p5 now, compute states 5,6 and 7 using the correct p5
        if PRINT_STATUS: print "Once p5 is known, find conditions at state 5 and 6."
        states['s5'].set_pT(p5,T0)
        (V6, V['s6']) = normal_shock(states['s5'], Vs2, states['s6'],ideal_gas_guess)
        V['s7'], states['s7'] = finite_wave_dv('cplus', V['s2'], states['s2'], V['s6'])
        
        #get mach numbers for the output
        Ms2 = Vs2/states['s5'].son
        M['s6'] = V['s6']/states['s6'].son
        M['s7']= V['s7']/states['s7'].son
        
        #do ideal gas guess for gas constants if needed
        
        if secondary:
            if Vs3 >= 10000.0:
                ideal_gas_guess = gas_guess
            else:
                ideal_gas_guess = None
                
            if PRINT_STATUS: print "Start unsteady expansion of the test gas."
    
            def error_in_velocity_s6_expansion(p8, state6=states['s6'], 
                                               V6g=V['s6'], state8=states['s8'],
                                                state9=states['s9'],Vs3=Vs3,
                                                ideal_gas_guess=ideal_gas_guess):
                """Compute the velocity mismatch for a given pressure ratio across the 
                unsteady expansion from state 6 to state 10."""
                
                print "current guess for p8 = {0} Pa".format(p8)
                
                state8.set_pT(p8,300.0) #set s1 at set pressure and ambient temp
                
                (V9, V9g) = normal_shock(state8, Vs3, state9,ideal_gas_guess)
                
                #Across the contact surface, p10 = p9
                p10 = state9.p
                
                # Across the expansion, we get a velocity, V7g.
                V10g, state10 = finite_wave_dp('cplus', V6g, state6, p10)
    
                return (V9g - V10g)/V9g        
            
            #get p8 for our chosen shock speed
            p8 = secant(error_in_velocity_s6_expansion, 10.0, 100.0, tol=1.0e-3,limits=[5.0,100.0])
            if PRINT_STATUS: print "From secant solve: p8 = {0} Pa".format(p8)
            
            #start using p8 now, compute states 8,9 and 10 using the correct p8
            if PRINT_STATUS: print "Once p8 is known, find conditions at state 8 and 9."
            states['s8'].set_pT(p8,T0)
            (V9, V['s9']) = normal_shock(states['s8'], Vs3, states['s9'],ideal_gas_guess)
            V['s10'], states['s10'] = finite_wave_dv('cplus', V['s6'], states['s6'], V['s9'])
            
            #get mach numbers for the output
            Ms2 = Vs3/states['s8'].son
            M['s9'] = V['s9']/states['s9'].son
            M['s10']= V['s10']/states['s10'].son        
        
        #do the nozzle calc up to s15 now
        if nozzle:
            if PRINT_STATUS: print "Start steady expansion through the nozzle."
            if secondary: 
                nozzle_entry_state = 's10'
            else:
                nozzle_entry_state = 's7'
            (V['s15'], states['s15']) = steady_flow_with_area_change(states[nozzle_entry_state], V[nozzle_entry_state], 2.5)
            M['s15']= V['s15']/states['s15'].son
        
        #do normal shock over model
        if PRINT_STATUS: print "Start normal shock calculation over the test model."  
        states['s16'] = states['s15'].clone()
        (V16, V['s16']) = normal_shock(states['s15'], V['s15'], states['s16'])
        M['s16']= V['s16']/states['s16'].son
        
        #--------------------------- output --------------------------------

        output = open(filename,"w")  #output file creation
        
        #added a line that tells the user easily what the important lines of result are
        
        if secondary:
            important_lines1 = 'Line 1 is secondary driver fill. Line 4 is reservoir. Line 5 is shock tube fill.'
            print important_lines1
            output.write(important_lines1 + '\n')   
            
            important_lines2 = 'Line 10 is acceleration tube fill. Line 8 is test gas entering nozzle.'
            print important_lines2
            output.write(important_lines2 + '\n')
            
            important_lines3 = 'Line 15 is test gas exiting nozzle (using area ratio of 2.5).'
            print important_lines3
            output.write(important_lines3 + '\n')
            
        else:
            important_lines1 = 'State 1 is shock tube fill. State 4 is reservoir at diaphragm burst.'
            print important_lines1
            output.write(important_lines1 + '\n')   
            
            important_lines2 = 'State 5 is acceleration tube fill. State 7 is test gas entering nozzle.'
            print important_lines2
            output.write(important_lines2 + '\n')
            
            important_lines3 = 'State 15 is test gas exiting nozzle (using area ratio of 2.5).'
            print important_lines3
            output.write(important_lines3 + '\n')
            
            important_lines4 = 'State 16 is shocked test gas flowing over the model.'
            print important_lines4
            output.write(important_lines4 + '\n')
        
        test_gas_used = 'Test gas is {0}'.format(states['s1'].reactants)        
        print test_gas_used
        output.write(test_gas_used + '\n')        
        
        shockspeeds = "Vs1= {0} m/s,Ms1= {1:.2f} ,Vs2= {2} m/s,Ms2= {3:.2f}".format(Vs1,Ms1,Vs2,Ms2) 
        print shockspeeds #prints above line in console
        output.write(shockspeeds + '\n') #writes above line to output file (input to write command must be a string)
        
        gas = "Driver g/R = {0:.2f}, {1}, test gas {2:.2f}, {3}, accel {4:.2f}, {5}"\
        .format(states['s4'].gam, states['s4'].R,
                states['s1'].gam, states['s1'].R,
                states['s5'].gam, states['s5'].R)
        print gas
        output.write(gas + '\n')
        
        key = "{0:7}{1:10}{2:9}{3:9}{4:9}{5:9}{6:9}{7:9}{8:9}".format("Region","P","T","a","V","M","rho","pitot","stgn")
        print key
        output.write(key + '\n')
        
        units = "{0:7}{1:10}{2:9}{3:9}{4:9}{5:9}{6:9}{7:9}{8:9}".format("","Pa","K","m/s","m/s","","m^3/kg","kPa","MPa")
        print units
        output.write(units + '\n')
        
        #new dictionaries here to add pitot and stagnation pressure calcs
        
        pitot = {} #pitot pressure dict
        p0 = {} #stagnation pressure dict
        
        for i in range(1,17): #will cover 1-16
            """if M[I]== 0:
                XPP=0;PS=0
            else:
                #little bit here to cover where the nozzle code was fucking up
                if not isnan(P[I]):
                    XPP=int(P[I]*pitot(M[I],G[I])/1000)
                    PS=P[I]/(IP(M[I],G[I]))/10**6 #isentropic stgn pressure (MPa)
                else:
                    XPP=0;PS=0"""
            
            it_string = 's{0}'.format(i)
                    
            if states.has_key(it_string):
                
                if M[it_string]==0:
                    pitot[it_string] = 0
                    p0[it_string] = 0
                else:
                    pitot[it_string] = pitot_p(states[it_string].p,M[it_string],states[it_string].gam)/1000.0
                    p0[it_string] = p0_p(M[it_string], states[it_string].gam)*states[it_string].p/1.0e6
                
                conditions = "{0:<7}{1:<10.7g}{2:<9.1f}{3:<9.0f}{4:<9.0f}{5:<9.2f}{6:<9.4f}{7:<9.0f}{8:<9.1f}"\
                .format(i, states[it_string].p, states[it_string].T,
                        states[it_string].son,V[it_string],M[it_string],
                        states[it_string].rho, pitot[it_string], p0[it_string])
                        
                print conditions
                output.write(conditions + '\n')
            
        #added some other useful calculations to the end
               
        #Cp of test gas (J/KgK) (gamma*R/(gamma-1))
        
        """Cp_test_gas = (gases[test_gas][0]*gases[test_gas][1])/(gases[test_gas][0]-1)

        if secondary == 'y':
            stagnation_enthalpy_accn = (((0.5)*U[8]**2.0)+(Cp_test_gas*T[8]))/10**6 #MJ/kg
        else:
            stagnation_enthalpy_accn = (((0.5)*U[7]**2.0)+(Cp_test_gas*T[7]))/10**6 #MJ/kg
        
        stagnation_enthalpy = (((0.5)*U[15]**2.0)+(Cp_test_gas*T[15]))/10**6 #MJ/kg

        stag_enth_accn = 'The stagnation enthalpy (Ht) entering the nozzle is {0:<.5g} MJ/kg.'.format(stagnation_enthalpy_accn)
                
        print stag_enth_accn
        output.write(stag_enth_accn + '\n')     

        # calculate flight equivalent speed
        u_eq = math.sqrt(2.0*stagnation_enthalpy_accn)

        u_eq_print = 'The flight equivalent speed is {0:<.5g} km/s.'.format(u_eq)
        print u_eq_print
        output.write(u_eq_print + '\n')

        stag_enth = 'The stagnation enthalpy (Ht) leaving the nozzle is {0:<.5g} MJ/kg.'.format(stagnation_enthalpy)
                
        print stag_enth
        output.write(stag_enth + '\n')
        
        if stand == 1: #changed this if statement to make it print standing shock stuff if stand =1
            P[15]=P[2]*press(M[2],G[2])
            standing = "Pressure behind standing shock in region 2= {0} kPa".format(int(P[15]/1000))
            print standing
            output.write(standing + '\n')
        
        reflected = "Reflected shock speed= {0:.2f} m/s, ref = {1}".format(UR,ref)
        print reflected
        output.write(reflected + '\n')
        
        if REV == 1:
            REV1 = "Reverse shock needed to match velocity behind incident shock"
            print REV1
            output.write(REV1 + '\n')
            REV2 = "reverse shock M= {0}, absolute downstream velocity= {1} m/s".format(int(MR1*100)/100, int(U[11]-MR1*a[11]))
            print REV2
            output.write(REV2 + '\n')"""
            
        output.close()
        
        change_tuple = ('Vs1','Vs2','Vs3', 'quit')
        
        change = is_valid(raw_input('What do you want to change? {0}'.format(change_tuple)), change_tuple)
            
        if change == 'quit':
            arbitrary == False
            break
            #quit()
                    
        elif change == 'Vs1':
            difference = int(raw_input('How much do you want to change it by? '))
            Vs1 += difference
            print 'Vs1 = {0} m/s'.format(Vs1)
            
        elif change == 'Vs2':    
            difference = int(raw_input('How much do you want to change it by? '))
            Vs2 += difference
            print 'Vs2 = {0} m/s'.format(Vs2)
        
        else:
            difference = int(raw_input('How much do you want to change it by? '))
            Vs3 += difference
            print 'Vs3 = {0} m/s'.format(Vs3)
                           
if __name__ == '__main__':
    main()

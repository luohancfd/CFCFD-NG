#! /usr/bin/env python

"""
pitot.py: Equilibrium expansion tube simulator

This program can be used to estimate the flow conditions for a shock-processed 
flow for an expansion tube. It is currently setup to work primarily with the 
X2 and X3 expansion tube's at the University of QLD, but I daresay it could be
expanded (excuse the pun) to use other facilities without too much effort.

The gas is assumed to remain in thermochemical equilibrium and the flow 
processing is done in decoupled quasi-one-dimensional wave processes such as 
shock waves and expansion fans.

The aim of the program is for it to be versatile, and it can do a number of
different calculations and configurations using different command line arguments:

* expansion tube with and without a pure He secondary driver tube and/or nozzle
* different driver conditions currently (100%He,80%He:20%Ar,90%He:10%Ar for X2, 60%He:40%Ar for X3)
* choosing shock speeds to give certain fill pressures in the tube
* choosing fill pressures to give certain shock speeds
* inputing experimental data of fill pressures and shock speeds to 'fill in the gaps' so to speak
* the code has been setup to do high total pressure conditions where a shock
  propagates into the shock tube from the secondary driver (instead of the usual
  expansion). I don't think this feature is working very well at the moment though.. 
 
When run as an application, this program takes its input as
command line arguments, performs the requested calculations and prints a table
with all of the results to the screen. Any changes you want to make can then 
be made to certain parameters, or you can quit the program. 
These results are also conveniently printed to a text file.

To see what specific inputs are required, start the program as::

$ pitot.py --help

Which particular input parameters you need to supply depends on the
chosen task, however, a fully theoretical basic x2 condition can be ran by using::

$ pitot.py --facility x2 --test fulltheory-pressure --config nozzle --driver_gas He:1.0 --test_gas air --p1 3000.0 --p5 10.0

due to the some of the default values this can be shortened down to this once 
you get the feel of it everything:

$ pitot.py --driver_gas He:1.0 --test_gas air --p1 3000.0 --p5 10.0
    
The full output is a bit too much to include here, but you should see that the
stagnation enthalpy leaving the nozzle is 62.87 MJ/kg, and the velocity exiting
the nozzle (V in state 8) should be 10639.5 m/s.

A brief description of the states can be found below:
    
State 4: primary driver at the burst of the steel diaphragm (27.9MPa,2700K), 
    taken to be basically stagnated (M~0)
State 3s: primary driver after steady expansion to the throat mach number 
    (M=1 for 80%He, M=1.59 for 90%He,M=2.15 for 100%He, based on area ratio stuff done by Richard)

If secondary driver is used:
State sd1 : secondary driver fill condition (pure He)
State sd2: shocked secondary driver gas
State sd3: state 11 unsteadily expanding into the secondary driver tube

State 1: shock tube fill condition (test gas)
State 2: shocked test gas
State 3: state 11 or state sd2 unsteadily expanding into the shock tube

State 5: acceleration tube fill condition (air)
State 6: shocked acceleration tube gas
State 7: test gas unsteadily expanding into the acceleration tube

State 8: nozzle exit condition after steady expansion through the nozzle
State 10: post shock condition at the stagnation point of the model


Getting the program set up
--------------------------
pitot.py is not a stand-alone file.
It comes as part of the cfcfd3 compressible-flow collection and
depends upon functions from the cfpylib library to do the specific 
calculations.

Someday I'll make a proper makefile for this guy, but for now,
just make sure cfpylib is in your Python $PATH.

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
    24-Oct-2012: changed notation of each state to fall in line with the notation
        used by David Gildfind in his work, I thought it made more sense as
        it means the shock tube and acc tubes always have the same notation
        (regardless of whether or not secondary driver is in place). Added one of
        my own conditions as an example too. Also added ability to use pressures
        and shock speeds from testing with pitot.
    30-Oct-2012: modified some of the notation to add x3 as a facility. Tested that all
        of the x2 simulations still worked, and then uploaded the new version.
    8-Nov-2012: added the final stuff for the scramjet conditions, I don't know
        how well it's working now, but for better or for worse, it's in there.
    29-Nov-2012: Did some little changes to the printout at the end of the test
        to make it clearer and some of the basic syntax.
    17-Dec-2012: Added two more test gases to the code.
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
from cfpylib.gasdyn.ideal_gas_flow import p0_p, pitot_p

VERSION_STRING = "19-Dec-2012"

DEBUG_PITOT = False

PRINT_STATUS = 1 #if print status is 1, some basic printouts are done

#------------------- define any functions needed ---------------------------
      
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
    
def make_test_gas(gasName, outputUnits='moles'):
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
                   inputUnits='moles', onlyList=['N2','O2','N','O','NO'],with_ions=False,
                   outputUnits=outputUnits), {'gam':1.35,'R':571.49}
    elif gasName.lower() == 'n2':
        return Gas(reactants={'N2':1.0, 'N':0.0}, onlyList=['N2', 'N'], with_ions=True,
                   outputUnits=outputUnits), {'gam': 1.36,'R':593.56}
    elif gasName.lower() == 'titan':
        return Gas(reactants={'N2':0.95, 'CH4':0.05}, inputUnits='moles', with_ions=True,
                   outputUnits=outputUnits), {'gam':1.35,'R':652.03}      
    elif gasName.lower() == 'mars':
        return Gas(reactants={'CO2':0.96, 'N2':0.04}, inputUnits='moles', with_ions=True,
                   outputUnits=outputUnits)    
    elif gasName.lower() == 'co2':
        return Gas(reactants={'CO2':1.0}, outputUnits=outputUnits)
    elif gasName.lower() == 'gasgiant_h215ne':
        return Gas(reactants={'H2':0.85, 'Ne':0.15}, inputUnits='moles', with_ions=True,
                   outputUnits=outputUnits),{'gam':1.15,'R':3245.1}
    elif gasName.lower() == 'gasgiant_h240ne':
        return Gas(reactants={'H2':0.6, 'Ne':0.4}, inputUnits='moles', with_ions=True,
                   outputUnits=outputUnits),{'gam':1.5,'R':1443.2}
    elif gasName.lower() == 'gasgiant_h285ne':
        return Gas(reactants={'H2':0.15, 'Ne':0.85}, inputUnits='moles', with_ions=True,
                   outputUnits=outputUnits),{'gam':1.5,'R':547.8}
    elif gasName.lower() == 'gasgiant_h215he':
        return Gas(reactants={'H2':0.85, 'He':0.15}, inputUnits='moles', with_ions=True,
                   outputUnits=outputUnits), {'gam':1.2,'R':6303.2}
    elif gasName.lower() == 'gasgiant_h210he':
        return Gas(reactants={'H2':0.9, 'He':0.10}, inputUnits='moles', with_ions=True,
                   outputUnits=outputUnits), None    
    elif gasName.lower() == 'gasgiant_h210ne':
        return Gas(reactants={'H2':0.9, 'Ne':0.10}, inputUnits='moles', with_ions=True,
                   outputUnits=outputUnits), None                  
    else:
        raise Exception, 'make_gas_from_name(): unknown gasName: %s' % gasName               
    
#--------------------------- property dictionaries -------------------------

#primary driver conditions, sorted into a dictionary with a gas object at the
#right conditions, and then mach number at the change over from steady to
#unsteady expansion, this was based on calcs done by RGM

primary_driver_x2 = dict([('He:1.0', [Gas({'He':1.0},inputUnits='moles'),2.15]),
                       ('He:0.80,Ar:0.20',[Gas({'He':0.8,'Ar':0.2},inputUnits='moles'),1]),
                        ('He:0.90,Ar:0.10',[Gas({'He':0.9,'Ar':0.1},inputUnits='moles'),1.59])])
                        
primary_driver_x3 = dict([('He:0.60,Ar:0.40',[Gas({'He':0.6,'Ar':0.4},inputUnits='moles'),2.23])])
  
def main():

#---------------------- getting the inputs set up --------------------------

    import optparse
    
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--facility', dest='facility', default='x2',
                  choices=['x2','x3'],
                  help=("facility to use; "
                        "x2 = the x2 expansion tube; "
                        "x3 = the x3 expansion tube; "
                        "defaults to x2 "))
    op.add_option('--test', dest='test', default='fulltheory-pressure',
                  choices=['fulltheory-shock','fulltheory-pressure','test'],
                  help=("type of test to run. "
                        "fulltheory-shock = fully theoretical run where fill pressures are found from set shock speeds "
                        "fulltheory-pressure = fully theoretical run where shock speeds are found from set fill pressures "
                        "test = partially theoretical run where both shock speeds and fill pressures are specified based on tunnel data "))
                        
    op.add_option('--config', dest='config', default='nozzle',
                  choices=['basic','sec-nozzle','sec','nozzle'],
                  help=("tunnel configuration to use. "
                        "basic = no secondary driver, no nozzle; "
                        "sec = with secondary driver, no nozzle; "
                        "nozzle = with no secondary driver, nozzle; "
                        "sec-nozzle = with secondary driver, nozzle "))
    op.add_option('--driver_gas', dest='driver_gas', default='He:1.0',
                  choices=['He:0.80,Ar:0.20', 'He:0.90,Ar:0.10','He:1.0',
                           'He:0.60,Ar:0.40'],
                  help=("driver gas configuration: "
                        "'He:0.80,Ar:0.20'; " 
                        "'He:0.90,Ar:0.10'; " 
                        "'He:1.0' "
                        "'He:0.60,Ar:0.40'"
                        "all by moles "
                        "default is currently 100% He"))
    op.add_option('--test_gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'n2', 'titan', 
                           'gasgiant_h215ne', 'gasgiant_h215he',
                           'gasgiant_h240ne','gasgiant_h285ne', 
                           'gasgiant_h210he', 'gasgiant_h210ne'],
                  help=("name of test gas: "
                        "air; " "air5species; " "n2; " "titan; " "gasgiant_h215ne; "
                        "gasgiant_h215he; " "gasgiant_h240ne; " "gasgiant_h285ne; " 
                        "gasgiant_h210he; "  "gasgiant_h210ne; " 
                        "default is air"))
    op.add_option('--Vs1', dest='Vs1', type='float', default=None,
                  help=("first shock speed, in m/s"))
    op.add_option('--Vs2', dest='Vs2', type='float', default=None,
                  help=("second shock speed, in m/s"))
    op.add_option('--Vsd', dest='Vsd', type='float', default=None,
                  help=("third shock speed, in m/s "
                        "not needed if secondary driver isn't used"))   
    op.add_option('--p1', dest='p1', type='float', default=None,
                  help=("shock tube fill pressure, in Pa"))
    op.add_option('--p5', dest='p5', type='float', default=None,
                  help=("acceleration tube fill pressure, in Pa"))
    op.add_option('--psd1', dest='psd1', type='float', default=None,
                  help=("secondary driver fill pressure, in Pa "
                        "not needed if secondary driver isn't used"))   
    op.add_option('--area_ratio', dest='ar', type='float', default=None,
                  help=("nozzle area ratio"
                        "default value is 2.5"))
    op.add_option('--conehead', dest='conehead', default=None,
                  help=("switch to calculate conehead pressure over a 15 degree conehead "
                        "default is don't do it, type yes or anything else to make it happen "))                    
    op.add_option('--filename',dest='filename',type='string',default='x2run.txt',
                  help=("filename the result will be saved to "
                          "defaults to x2run.txt"))
    op.add_option('--expansion',dest='expansion',type='float',default='1.0',
                  help=("Use this to expand the test gas further into the acceleration tube "
                        "Basically just multiplies V7 by this percentage before the final expansion is done "
                        "Defaults to 1.0 "))
    op.add_option('--shock_switch',dest='shock_switch', default=None,
                  help=("Used to trigger a shock instead of an expansion from the secondary driver gas moving into the shock tube"))
    
    opt, args = op.parse_args()
    
    facility = opt.facility
    test = opt.test
    config = opt.config
    driver_gas = opt.driver_gas
    gasName = opt.gasName
    filename =opt.filename
    
    Vs1 = opt.Vs1
    Vs2 = opt.Vs2
    Vsd = opt.Vsd
    
    p1 = opt.p1
    p5 = opt.p5
    psd1 = opt.psd1
    
    area_ratio = opt.ar
    conehead = opt.conehead
    shock_switch = opt.shock_switch
    expansion = opt.expansion

    bad_input = False

    if not area_ratio: #if they haven't specified their own area ratio

        if facility == 'x2':
            area_ratio = 2.5 #effective area ratio x2's nozzle normally has
        elif facility == 'x3':
            area_ratio = 5.8 #geometric area ratio of x3's nozzle
    
    if test == 'fulltheory-shock': #shows we're going to solve fill pressures from shock speeds
        solver = 'shock_speeds'
        if p1: #if they specify both pressures and shock speefs, bail out
            print "You need to choose either pressures or shock speeds to solve for. Bailing out here."
            bad_input = True
    
        if Vs1 is None:
            print "Need to supply a float value for Vs1."
            bad_input = True
        
        if Vs2 is None:
            print "Need to supply a float value for Vs2."
            bad_input = True
    
        if config == 'sec-nozzle':
            if Vsd is None:
                print "Need to supply a float value for Vsd."
                bad_input = True
            secondary = True
            nozzle = True    
        elif config == 'sec':
            if Vsd is None:
                print "Need to supply a float value for Vsd."
                bad_input = True
            secondary = True
            nozzle = False    
        elif config == 'nozzle':
            secondary = False
            nozzle = True
        else:
            secondary = False
            nozzle = False
        
        if bad_input: #bail out here if you end up having issues with your input
            return -2
            
        if PRINT_STATUS: print "Let's get started, shall we."
        if PRINT_STATUS: 
            if secondary:
                print 'selected Vsd = {0} m/s'.format(Vsd)
            print 'selected Vs1 = {0} m/s'.format(Vs1)
            print 'selected Vs2 = {0} m/s'.format(Vs2)

        
    elif test == 'fulltheory-pressure':  #shows we're going to solve shock speeds from fill pressures
        solver = 'fill_pressures'
        if Vs1: #if they specify both pressures and shock speeds, bail out
            print "You need to choose either pressures or shock speeds to solve for. Bailing out here."
            bad_input = True
    
        if p1 is None:
            print "Need to supply a float value for p1."
            bad_input = True
        
        if p5 is None:
            print "Need to supply a float value for p5."
            bad_input = True
    
        if config == 'sec-nozzle':
            if psd1 is None:
                print "Need to supply a float value for psd1."
                bad_input = True
            secondary = True
            nozzle = True    
        elif config == 'sec':
            if psd1 is None:
                print "Need to supply a float value for psd1."
                bad_input = True
            secondary = True
            nozzle = False    
        elif config == 'nozzle':
            secondary = False
            nozzle = True
        else:
            secondary = False
            nozzle = False

        if bad_input: #bail out here if you end up having issues with your input
            return -2

        if PRINT_STATUS: print "Let's get started, shall we."
        if PRINT_STATUS: 
            if secondary:
                print 'selected secondary driver fill pressure (psd1) = {0} Pa.'.format(psd1)

            print 'selected shock tube fill pressure (p1) = {0} Pa.'.format(p1)
            print 'selected acceleration tube fill pressure (p5) = {0} Pa.'.format(p5)
            
    elif test == 'test': #test based on real experimental data
        
        if Vs1 is None:
            print "Need to supply a float value for Vs1."
            bad_input = True
        
        if Vs2 is None:
            print "Need to supply a float value for Vs2."
            bad_input = True
            
        if p1 is None:
            print "Need to supply a float value for p1."
            bad_input = True
        
        if p5 is None:
            print "Need to supply a float value for p5."
            bad_input = True 
    
        if config == 'sec-nozzle':
            if Vsd is None:
                print "Need to supply a float value for Vsd."
                bad_input = True
            if psd1 is None:
                print "Need to supply a float value for psd1."
                bad_input = True
            if Vsd > Vs1:
                shock_switch = True
            secondary = True
            nozzle = True    
        elif config == 'sec':
            if Vsd is None:
                print "Need to supply a float value for Vsd."
                bad_input = True
            if psd1 is None:
                print "Need to supply a float value for psd1."
                bad_input = True
            if Vsd > Vs1:
                shock_switch = True
            secondary = True
            nozzle = False    
        elif config == 'nozzle':
            secondary = False
            nozzle = True
        else:
            secondary = False
            nozzle = False
        
        if bad_input: #bail out here if you end up having issues with your input
            return -2
            
        if PRINT_STATUS: print "Let's get started, shall we."
        if PRINT_STATUS: 
            if secondary:
                print 'selected Vsd = {0} m/s'.format(Vsd)
            print 'selected Vs1 = {0} m/s'.format(Vs1)
            print 'selected Vs2 = {0} m/s'.format(Vs2)
            if secondary:
                print 'selected secondary driver fill pressure (psd1) = {0} Pa.'.format(psd1)
            print 'selected shock tube fill pressure (p1) = {0} Pa.'.format(p1)
            print 'selected acceleration tube fill pressure (p5) = {0} Pa.'.format(p5)
            
    else: #test input must be bad
        "You haven't specified a relevant test condition. Bailing out..."
        bad_input = True
    
        if bad_input: #bail out here if you end up having issues with your input
            return -2


      
    arbitrary = True #dummy variable
    
    while arbitrary:
        
     #---------------- building initial states ----------------------------------    
               
        states = {} #states dictionary that we'll fill up later
        V = {} #same for velocity
        M = {} #same for Mach number
        
        
        #state 4, piston burst pressure
        
        if facility == 'x2':
            states['s4']=primary_driver_x2[driver_gas][0].clone()
            states['s4'].set_pT(2.79e7,2700.0)
            V['s4']=0.0
            M['s4']=0.0
        
            M['s3s']=primary_driver_x2[driver_gas][1]
            
        elif facility == 'x3':
            states['s4']=primary_driver_x3[driver_gas][0],clone()
            states['s4'].set_pT(2.79e7,2700.0)
            V['s4']=0.0
            M['s4']=0.0
        
            M['s3s']=primary_driver_x3[driver_gas][1]
    
        #state3s, piston after steady expansion at throat of driver
        
        if PRINT_STATUS: print "Start steady expansion of driver gas."     
        
        (states['s3s'], V['s3s']) = expand_from_stagnation(1.0/(p0_p(M['s3s'],states['s4'].gam)),states['s4'])
        
        #build the gas objects for all the other sections based on knowledge of what is what
        
        T0 = 300.0 #atmospheric temperature (K), used for any section starting at ambient
        p0 = 101300.0 #atmospheric p (Pa)
        
        #start with shock tube and acc tube, start with atmospheric p and T
        
        if secondary: #state sd1 is pure He secondary driver
            states['sd1'] =  Gas({'He':1.0,},outputUnits='moles')
            if psd1:
                states['sd1'].set_pT(psd1,T0)
            else:
                states['sd1'].set_pT(p0,T0)
            V['sd1']=0.0
            M['sd1']=0.0
    
        #state 1 is shock tube
        states['s1'], gas_guess = make_test_gas(gasName)
        if p1:
            states['s1'].set_pT(p1,T0)
        else:
            states['s1'].set_pT(p0,T0)
        V['s1']=0.0
        M['s1']=0.0
            
    
        #state 5 is acceleration tube
        states['s5'] = Gas({'Air':1.0,},outputUnits='moles')
        if p5:
            states['s5'].set_pT(p5,T0)
        else:
            states['s5'].set_pT(p0,T0)
        V['s5']=0.0
        M['s5']=0.0
               
        #now let's clone for the states derived from these
        
        if secondary:
            states['sd2'] = states['sd1'].clone() #sd2 is sd1 shock heated
            states['sd3'] = states['s3s'].clone() #sd3 is s11 after unsteady
            states['s3'] = states['sd2'].clone() #3 will be sd2 after unsteady expansion
        else:
            states['s3'] = states['s3s'].clone() #3 will be 11 after unsteady expansion
            
        states['s2'] = states['s1'].clone() #2 is 1 shock heated
        states['s6'] = states['s5'].clone() #6 is 5 shock heated
        states['s7'] = states['s2'].clone() #7 is 2 after unsteady expansion   
        
    #--------- unsteady expansion of driver gas-----------------------------
        
        print "Start unsteady expansion of the driver gas."
        # For the unsteady expansion of the test driver into the tube, regulation of the amount
        # of expansion is determined by the shock-processed gas in the next section.
        # Across the contact surface between these gases, the pressure and velocity
        # have to match so we set up some trials of various pressures and check 
        # that velocities match.
        
        if secondary:
            
            def error_in_velocity_s3s_to_sd3_expansion_pressure_iterator(psd1, state3s=states['s3s'], 
                                                       V3sg=V['s3s'], statesd1=states['sd1'],
                                                        statesd2=states['sd2'],Vsd=Vsd):
                """Compute the velocity mismatch for a given pressure ratio across the 
                unsteady expansion in the driver into the secondary driver gas."""
                
                print "current guess for psd1 = {0} Pa".format(psd1)            
                
                statesd1.set_pT(psd1,300.0) #set s1 at set pressure and ambient temp
                
                (Vsd2, Vsd2g) = normal_shock(statesd1, Vsd, statesd2)
                
                #Across the contact surface, p3 == p2
                psd3 = statesd2.p
                
                # Across the expansion, we get a velocity, V5g.
                Vsd3g, statesd3 = finite_wave_dp('cplus', V3sg, state3s, psd3)
    
                return (Vsd2g - Vsd3g)/Vsd2g
                
            def error_in_velocity_s3s_to_sd3_driver_expansion_shock_speed_iterator(Vsd,state3s=states['s3s'], 
                                                       V3sg=V['s3s'], statesd1=states['sd1'],
                                                        statesd2=states['sd2']):
                """Compute the velocity mismatch for a given shock speed with the driver
                unsteady expansion into the secondary driver gas behind it.
                
                Make sure you set the fill pressure in sd1 before you start!"""
                
                print "current guess for Vsd = {0} m/s".format(Vsd)            
                
                (Vsd2, Vsd2g) = normal_shock(statesd1, Vsd, statesd2)
                
                #Across the contact surface, p3 == p2
                psd3 = statesd2.p
                
                # Across the expansion, we get a velocity, V5g.
                Vsd3g, statesd3 = finite_wave_dp('cplus', V3sg, state3s, psd3)
    
                return (Vsd2g - Vsd3g)/Vsd2g  
            
            if test == "fulltheory-shock": #get psd1 for our chosen shock speed
                psd1 = secant(error_in_velocity_s3s_to_sd3_expansion_pressure_iterator, 5000.0, 6000.0, tol=1.0e-4,limits=[1000.0,1000000.0])
                if PRINT_STATUS: print "From secant solve: psd1 = {0} Pa".format(psd1)
                #start using psd1 now, compute states sd1,sd2 and sd3 using the correct psd1
                if PRINT_STATUS: print "Once p1 is known, find conditions at state 1 and 2."
                states['sd1'].set_pT(psd1,T0)
            
            elif test == "fulltheory-pressure": #get Vsd for our chosen fill pressure
                Vsd = secant(error_in_velocity_s3s_to_sd3_driver_expansion_shock_speed_iterator, 4000.0, 5000.0, tol=1.0e-4,limits=[1000.0,1000000.0])
                if PRINT_STATUS: print "From secant solve: Vsd = {0} m/s".format(Vsd)
                #start using Vs1 now, compute states 1,2 and 3 using the correct Vs1
                if PRINT_STATUS: print "Once Vsd is known, find conditions at state sd1 and sd2."
            
            (Vsd2, V['sd2']) = normal_shock(states['sd1'], Vsd, states['sd2'])
            V['sd3'], states['sd3'] = finite_wave_dv('cplus', V['s3s'], states['s3s'], V['sd2'])
        
            #get mach numbers for the output
            Msd1 = Vsd/states['sd1'].son
            M['sd2'] = V['sd2']/states['sd2'].son
            M['sd3']= V['sd3']/states['sd3'].son
        
        if secondary:
            shock_tube_expansion = 'sd2' #state sd2 expands into shock tube
        else:
            shock_tube_expansion = 's3s' #otherwise state 3s is expanding into shock tube
            
#--------------------- shock tube-------------------------------------------

#--------------------- shock tube functions ---------------------------------
            
        def error_in_velocity_shock_tube_expansion_pressure_iterator(p1, 
                                                   expanding_state=states[shock_tube_expansion], 
                                                   expansion_start_V=V[shock_tube_expansion], 
                                                    state1=states['s1'],
                                                    state2=states['s2'],Vs1=Vs1):
            """Compute the velocity mismatch for a given pressure ratio across the 
            unsteady expansion in the shock tube."""
            
            print "current guess for p1 = {0} Pa".format(p1)            
            
            state1.set_pT(p1,300.0) #set s1 at set pressure and ambient temp
            
            (V2, V2g) = normal_shock(state1, Vs1, state2)
            
            #Across the contact surface, p3 == p2
            p3 = state2.p
            
            # Across the expansion, we get a velocity, V5g.
            V3g, state3 = finite_wave_dp('cplus', expansion_start_V, expanding_state, p3)

            return (V2g - V3g)/V2g
            
        def p1_reflected_iterator(p1, statesd2=states[shock_tube_expansion], 
                                  Vsd2g=V[shock_tube_expansion], state1=states['s1'], 
                                  state2=states['s2'], state3=states['s3'], Vs1=Vs1):
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
        
                V3, V3g = normal_shock(statesd2, Vsd2g+Vr, state3)
        
                return (V3-(V2g+Vr))/V3
            
            Vr = secant(reflected_shock_speed_iterator, 100.0, 200.0, tol=1.0e-4)
        
            V3,V3g = normal_shock(statesd2,Vsd2g+Vr,state3)
        
            return (state3.p-state2.p)/state2.p
                
            
            #Across the contact surface, p3 == p2
            p3 = state2.p
            
            # therefore we know the required pressure ratio across our shock
            
            print 'shocking state psd2 = {0} Pa, p3 = {1} Pa'.format(shocking_state.p,p3)
            
            (Vs1found, V3, V3g, state3) = normal_shock_p2p1(shocking_state, p3/shocking_state.p)
            
            print 'Vs1 found is {0} m/s'.format(Vs1found)
            
            print 'V2g = {0} m/s, V3g = {1} m/s'.format(V2g,V3g)
            
            print 'V3 = {0} m/s'.format(V3)
            
            return (V2g - V3g)/V2g        
            
        def error_in_velocity_shock_tube_expansion_shock_speed_iterator(Vs1,
                                                    expanding_state=states[shock_tube_expansion], 
                                                    expansion_start_V=V[shock_tube_expansion],                    
                                                    state1=states['s1'], state2=states['s2']):
            """Compute the velocity mismatch for a given shock speed with the shock tube
            unsteady expansion behind it.
            
            Make sure you set the fill pressure in state 1 before you start!"""
            
            print "current guess for Vs1 = {0} m/s".format(Vs1)            
            
            (V2, V2g) = normal_shock(state1, Vs1, state2)
            
            #Across the contact surface, p3 == p2
            p3 = state2.p
            
            # Across the expansion, we get a velocity, V5g.
            V3g, state3 = finite_wave_dp('cplus', expansion_start_V, expanding_state, p3)
        
            return (V2g - V3g)/V2g

        def primary_shock_speed_reflected_iterator(Vs1,statesd2=states[shock_tube_expansion],
                                                   Vsd2g=V[shock_tube_expansion],state1=states['s1'],
                                                   state2=states['s2'],state3=states['s3']):
            """The other new iterator for pitot. Has another iterator inside it, will take
            AGES to run. Oh well, no other choice as far as I see."""
        
            print "Current guess for Vs1 = {0} m/s".format(Vs1)
            
            V2, V2g = normal_shock(state1, Vs1, state2)
        
            def reflected_shock_speed_iterator(Vr,statesd2=statesd2,
                                               Vsd2g=Vsd2g, state3=state3,
                                               V2g=V2g):
                """iterates out the reflected shock speed that gives the right velocity behind
                    a normal shock."""
        
                print "Current guess for Vr = {0} m/s".format(Vr)
        
                V3, V3g = normal_shock(statesd2, Vsd2g+Vr, state3)
        
                return (V3-(V2g+Vr))/V3
        
            Vr = secant(reflected_shock_speed_iterator, 100.0, 200.0, tol=1.0e-4)
        
            V3,V3g = normal_shock(statesd2,Vsd2g+Vr,state3)
        
            return (state3.p-state2.p)/state2.p


#----------------------- shock tube calculations --------------------------------                           
            
        if test == "fulltheory-shock": #get p1 for our chosen shock speed
            if secondary: #if secondary, check which shock speed is greater, Vsd or Vs1
                if Vsd > Vs1: #we get a shock into the shock tube, not an expansion
                    shock_switch = True #flip over this switch
                    if PRINT_STATUS: print "Vsd is a stronger shock than Vs1. Therefore a normal shock processes the secondary driver gas moving into the shock tube."
                    p1 = secant(p1_reflected_iterator, 1000.0, 100000.0, tol=1.0e-3,limits=[100.0,1000000.0])
                else: #same expansion as we're used to
                    if PRINT_STATUS: print "Start unsteady expansion of the secondary driver gas into the shock tube."
                    p1 = secant(error_in_velocity_shock_tube_expansion_pressure_iterator, 1000.0, 100000.0, tol=1.0e-3,limits=[100.0,1000000.0])
            else: #just do the expansion
                if PRINT_STATUS: print "Start unsteady expansion of the driver gas into the shock tube."
                p1 = secant(error_in_velocity_shock_tube_expansion_pressure_iterator, 1000.0, 100000.0, tol=1.0e-3,limits=[100.0,1000000.0])

            if PRINT_STATUS: print "From secant solve: p1 = {0} Pa".format(p1)
            #start using p1 now, compute states 1,2 and 3 using the correct p1
            if PRINT_STATUS: print "Once p1 is known, find conditions at state 1 and 2."
            states['s1'].set_pT(p1,T0)
            
        elif test =="fulltheory-pressure": #get Vs1 for our chosen fill pressure
            if shock_switch: #if we've been told to do a shock here instead of an expansion, do a shock instead of an expansion
                if PRINT_STATUS: print "The shock switch is turned on, therefore doing a shock here instead of the normal expansion... Turn this off if you didn't want it"
                if PRINT_STATUS: print "Current settings with the shock switch turned on here HALF the selected p1 and p5 to match experiment."
                if PRINT_STATUS: print "the fill pressure's will be returned back to normal when the results are printed at the end..."
                states['s1'].p = states['s1'].p/2.0
                states['s5'].p = states['s5'].p/2.0
                Vs1 = secant(primary_shock_speed_reflected_iterator, 2000.0, 4000.0, tol=1.0e-3,limits=[500.0,1000000.0])
            else: #just do the expansion
                Vs1 = secant(error_in_velocity_shock_tube_expansion_shock_speed_iterator, 2000.0, 4000.0, tol=1.0e-3,limits=[1000.0,1000000.0])
            if PRINT_STATUS: print "From secant solve: Vs1 = {0} m/s".format(Vs1)
            #start using Vs1 now, compute states 1,2 and 3 using the correct Vs1
            if PRINT_STATUS: print "Once Vs1 is known, find conditions at state 1 and 2."
        
        (V2, V['s2']) = normal_shock(states['s1'], Vs1, states['s2'])
        if shock_switch: #do a shock here if required
            V['s3'] = V['s2']
            #need to back out Vr again here, for now I was lazy and put the function back in:
            if PRINT_STATUS: print "Need to reiterate to find Vr again here..."   
            def reflected_shock_speed_iterator(Vr,statesd2=states['sd2'],
                                               Vsd2g=V[shock_tube_expansion],
                                                state3=states['s3'], V2g=V['s2']):
                """iterates out the reflected shock speed that gives the right velocity behind
                    a normal shock."""
        
                print "Current guess for Vr = {0} m/s".format(Vr)
        
                V3, V3g = normal_shock(statesd2, Vsd2g+Vr, state3)
        
                return (V3-(V2g+Vr))/V3
            
            Vr = secant(reflected_shock_speed_iterator, 100.0, 200.0, tol=1.0e-4)
            (V3,Vjunk) = normal_shock(states['sd2'],V[shock_tube_expansion]+Vr,states['s3'])
        else:
            V['s3'], states['s3'] = finite_wave_dv('cplus', V[shock_tube_expansion], states[shock_tube_expansion], V['s2'])
        
        #get mach numbers for the output
        Ms1 = Vs1/states['s1'].son
        M['s2'] = V['s2']/states['s2'].son
        M['s3']= V['s3']/states['s3'].son
        
        #do ideal gas guess for gas constants if needed
        
        if p5: #just use the gas guess if you don't know
            ideal_gas_guess = gas_guess
        elif Vs2 >= 9000.0:
            ideal_gas_guess = gas_guess
        else:
            ideal_gas_guess = None

#--------------------- acceleration tube -----------------------------------
                    
        if PRINT_STATUS: print "Start unsteady expansion of the test gas into the acceleration tube."

#----------------------- acceleration tube functions -----------------------

        def error_in_velocity_s2_expansion_pressure_iterator(p5, state2=states['s2'], 
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
            
        def error_in_pressure_s2_expansion_shock_speed_iterator(Vs2, state2=states['s2'], 
                                           V2g=V['s2'], state5=states['s5'],
                                            state6=states['s6'],
                                            ideal_gas_guess=ideal_gas_guess):
            """Compute the velocity mismatch for a given shock speed in front of the 
            unsteady expansion from state 2 to state 7."""
            
            print "current guess for Vs2 = {0} m/s".format(Vs2)
            
            (V6, V6g) = normal_shock(state5, Vs2, state6, ideal_gas_guess)
            
            #Across the contact surface, p3 == p2
            p7 = state6.p
                       
            # Across the expansion, we get a velocity, V7g.
            V7g, state7 = finite_wave_dp('cplus', V2g, state2, p7)

            return (V6g - V7g)/V6g    
            
#---------------------- acceleration tube calculations -----------------------
        
        if test == 'fulltheory-shock': #get p5 for our chosen shock speed
            #put two sets of starting guesses here to try and make code that works will all conditions
            if Vs2 > 4000.0:
                p5 = secant(error_in_velocity_s2_expansion_pressure_iterator, 10.0, 30.0, tol=1.0e-4,limits=[1.0,1000.0])
            else:
                p5 = secant(error_in_velocity_s2_expansion_pressure_iterator, 400.0, 500.0, tol=1.0e-4,limits=[20.0,1000.0])
            if PRINT_STATUS: print "From secant solve: p5 = {0} Pa".format(p5)
            
            #start using p5 now, compute states 5,6 and 7 using the correct p5
            if PRINT_STATUS: print "Once p5 is known, find conditions at state 5 and 6."
            states['s5'].set_pT(p5,T0)
            
        elif test == 'fulltheory-pressure': #compute the shock speed for the chosen fill pressure, uses Vs1 as starting guess
            #put two sets of starting guesses here to try and make more stuff work
            if Vs1 < 2000.0:
                Vs2 = secant(error_in_pressure_s2_expansion_shock_speed_iterator, Vs1+3000.0, 12000.0, tol=1.0e-3,limits=[Vs1,100000.0])
            else:
                Vs2 = secant(error_in_pressure_s2_expansion_shock_speed_iterator, Vs1+3000.0, 12000.0, tol=1.0e-3,limits=[Vs1,100000.0])
            if PRINT_STATUS: print "From secant solve: Vs2 = {0} m/s".format(Vs2)
            #start using Vs1 now, compute states 1,2 and 3 using the correct Vs1
            if PRINT_STATUS: print "Once Vs2 is known, find conditions at state 5 and 6."            
            
        (V6, V['s6']) = normal_shock(states['s5'], Vs2, states['s6'],ideal_gas_guess)
        V['s7'], states['s7'] = finite_wave_dv('cplus', V['s2'], states['s2'], (V['s6']*expansion))
        
        #get mach numbers for the output
        Ms2 = Vs2/states['s5'].son
        M['s6'] = V['s6']/states['s6'].son
        M['s7']= V['s7']/states['s7'].son
        
#------------------- finishing off any other needed calculations -------------
                          
        #do the nozzle calc up to s8 now
        if nozzle:
            if PRINT_STATUS: print "Start steady expansion through the nozzle."
            (V['s8'], states['s8']) = steady_flow_with_area_change(states['s7'], V['s7'], area_ratio)
            M['s8']= V['s8']/states['s8'].son
            state_over_model = 's8'
        else:
            state_over_model = 's7'             
        
        #do normal shock over model
        
        if PRINT_STATUS: print "Start frozen normal shock calculation over the test model."  
        states['s10f'] = states[state_over_model].clone()
        (V10, V['s10f']) = shock_ideal(states[state_over_model], V[state_over_model], states['s10f'])
        M['s10f']= V['s10f']/states['s10f'].son        
        
        if PRINT_STATUS: print "Start equilibrium normal shock calculation over the test model."  
        states['s10e'] = states[state_over_model].clone()
        (V10, V['s10e']) = normal_shock(states[state_over_model], V[state_over_model], states['s10e'])
        M['s10e']= V['s10e']/states['s10e'].son
        
        if conehead:
            if PRINT_STATUS: print 'Doing taylor maccoll conehead calculation on 15 degree conehead.'
            shock_angle = beta_cone(states[state_over_model], V[state_over_model], math.radians(15))
            if PRINT_STATUS: print "\nShock angle over cone:", math.degrees(shock_angle)
            # Reverse the process to get the flow state behind the shock and check the surface angle is correct
            delta_s, V['s10c'], states['s10c'] = theta_cone(states[state_over_model], V[state_over_model], shock_angle)
            M['s10c'] = V['s10c']/states['s10c'].son
            if PRINT_STATUS: print "Surface angle should be the same.....: 15deg = ", math.degrees(delta_s), "deg"
            if PRINT_STATUS: print "\nConehead surface conditions:"
            if PRINT_STATUS: states['s10c'].write_state(sys.stdout)
            # Need to check whether the pressure are the same
            if PRINT_STATUS: print "Computed conehead pressure is {0} Pa".format(states['s10c'].p)

        if test == "fulltheory-pressure" and shock_switch: #restore the fill pressure's back to normal
            states['s1'].p = states['s1'].p*2.0
            states['s5'].p = states['s5'].p*2.0
            
        if test == "fulltheory-shock" and shock_switch: #cut fill pressure's in half at the end
            states['s1'].p = states['s1'].p*2.0
            states['s5'].p = states['s5'].p*2.0
                
        
      
#--------------------------- output --------------------------------

        output = open(filename,"w")  #output file creation
                    
        print " "
        
        version_printout = "Pitot Version: {0}".format(VERSION_STRING)
        print version_printout
        output.write(version_printout + '\n')
        
        if secondary:
            description_sd = 'sd1 is secondary driver fill.'
            print description_sd
            output.write(description_sd + '\n')   
            
        description_1 = 'state 1 is shock tube fill. state 5 is acceleration tube fill.'    
        print description_1
        output.write(description_1 + '\n')
            
        description_2 = 'state 7 is expanded test gas entering the nozzle.'
        print description_2
        output.write(description_2 + '\n')
            
        description_3 = 'state 8 is test gas exiting the nozzle (using area ratio of {0}).'.format(area_ratio)
        print description_3
        output.write(description_3 + '\n')
            
        description_4 = 'state 10f is frozen shocked test gas flowing over the model.'
        print description_4
        output.write(description_4 + '\n')  
        
        description_5 = 'state 10e is equilibrium shocked test gas flowing over the model.'
        print description_5
        output.write(description_5 + '\n')        
            
        facility_used = 'Facility is {0}.'.format(facility)        
        print facility_used
        output.write(facility_used + '\n')
        
        test_gas_used = 'Test gas is {0} (gamma = {1}, R = {2}, {3}).'.format(gasName,states['s1'].gam,states['s1'].R,states['s1'].reactants)
        print test_gas_used
        output.write(test_gas_used + '\n')  
        
        driver_gas_used = 'Driver gas is {0}.'.format(driver_gas)       
        print driver_gas_used
        output.write(driver_gas_used + '\n') 
                
        if shock_switch:
            shock_warning1 = "NOTE: a reflected shock was done into the shock tube and as such,"
            print shock_warning1
            output.write(shock_warning1 + '\n')
            
            shock_warning2 = "fill pressure's have been artifically modified by the code to match with xpt."
            print shock_warning2
            output.write(shock_warning2 + '\n')

        if secondary:
            secondary_shockspeeds = "Vsd = {0:.2f} m/s, Msd1 = {1:.2f}".format(Vsd,Msd1)
            print secondary_shockspeeds
            output.write(secondary_shockspeeds + '\n')
        
        shockspeeds = "Vs1 = {0:.2f} m/s, Ms1 = {1:.2f} ,Vs2 = {2:.2f} m/s, Ms2 = {3:.2f}".format(Vs1,Ms1,Vs2,Ms2) 
        print shockspeeds #prints above line in console
        output.write(shockspeeds + '\n') #writes above line to output file (input to write command must be a string)
                
        key = "{0:6}{1:11}{2:9}{3:6}{4:9}{5:6}{6:9}{7:8}{8:9}".format("state","P","T","a","V","M","rho","pitot","stgn")
        print key
        output.write(key + '\n')
        
        units = "{0:6}{1:11}{2:9}{3:6}{4:9}{5:6}{6:9}{7:9}{8:9}".format("","Pa","K","m/s","m/s","","m^3/kg","kPa","MPa")
        print units
        output.write(units + '\n')
        
        #new dictionaries here to add pitot and stagnation pressure calcs
        
        pitot = {} #pitot pressure dict
        p0 = {} #stagnation pressure dict
        
        def condition_printer(it_string):
            """Prints the values of a specified condition to the screen and to 
            the output file. 
            
            I made a function of this so I didn't have to keep pasting the code in."""
            
            if states.has_key(it_string):
                    
                if M[it_string] == 0:
                    pitot[it_string] = 0
                    p0[it_string] = 0
                else:
                    pitot[it_string] = pitot_p(states[it_string].p,M[it_string],states[it_string].gam)/1000.0
                    p0[it_string] = p0_p(M[it_string], states[it_string].gam)*states[it_string].p/1.0e6
                
                conditions = "{0:<6}{1:<11.7}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.4f}{7:<8.0f}{8:<9.1f}"\
                .format(it_string, states[it_string].p, states[it_string].T,
                        states[it_string].son,V[it_string],M[it_string],
                        states[it_string].rho, pitot[it_string], p0[it_string])
                        
                print conditions
                output.write(conditions + '\n')   
                
        #print the driver related stuff first
        
        condition_printer('s4')
        condition_printer('s3s')
        
        if secondary: #need a separate printing thing for the secondary driver
            
            for i in range(1,4): #will do 1 - 3
            
                it_string = 'sd{0}'.format(i)
                condition_printer(it_string)
                        
        for i in range(1,4): #shock tube stuff
            
            it_string = 's{0}'.format(i)
            condition_printer(it_string)
            
        for i in range(5,9): #acc tube and nozzle if it's there
        
            it_string = 's{0}'.format(i)
            condition_printer(it_string)
            
        #do the conditions over the model
        condition_printer('s10f')
        condition_printer('s10e')
                
        if conehead:
            condition_printer('s10c')
                                       
        #some other useful calculations at the end
              
        total = total_condition(states[state_over_model], V[state_over_model])
        
        stagnation_enthalpy = total.h #J/kg
        if nozzle:        
            stag_enth = 'The stagnation enthalpy (Ht) leaving the nozzle is {0:<.5g} MJ/kg.'.format(stagnation_enthalpy/10**6)
        else:
            stag_enth = 'The stagnation enthalpy (Ht) at the end of the acceleration tube {0:<.5g} MJ/kg.'.format(stagnation_enthalpy/10**6)
        print stag_enth
        output.write(stag_enth + '\n')
        
        #calculate flight equivalent speed
        u_eq = math.sqrt(2.0*stagnation_enthalpy) #supposedly why this is is discussed in Bianca's thesis
        u_eq_print = 'The flight equivalent speed is {0:<.5g} m/s.'.format(u_eq)
        print u_eq_print
        output.write(u_eq_print + '\n')
        
        #added ability to get the species in the post-shock condition
        
        show_species = True
        
        if show_species:
            species1 = 'species in the shock layer at equilibrium:'        
            print species1
            output.write(species1 + '\n')
            
            species2 = '{0}'.format(states['s10e'].species)
            print species2
            output.write(species2 + '\n')
        
        #added taylor maccoll conehead calcs if people want to use them, appropriated from Fabs (thankyou)
                    
        output.close()
        
#--------------------------- user commands -----------------------------------
        
        if test == 'fulltheory-shock':
            if secondary:
                change_tuple = ('Vsd','Vs1','Vs2', 'ar', 'quit')
            else:
                change_tuple = ('Vs1','Vs2','ar','quit')
                
        elif test == 'fulltheory-pressure':
            if secondary:
                change_tuple = ('psd1','p1','p5','ar', 'quit')
            else:
                change_tuple = ('p1','p5','ar', 'quit')
                
        else:
            change_tuple = ('ar', 'quit')
            
        change = is_valid(raw_input('What do you want to change? {0}'.format(change_tuple)), change_tuple)
                
        if change == 'quit':
            print "Removing temporary files and leaving the program."
            if os.path.isfile('thermo.inp'): os.remove('thermo.inp')
            if os.path.isfile('thermo.out'): os.remove('thermo.out')
            if os.path.isfile('thermo.lib'): os.remove('thermo.lib')
            if os.path.isfile('tmp.inp'): os.remove('tmp.inp')
            if os.path.isfile('tmp.out'): os.remove('tmp.out')
            if os.path.isfile('trans.inp'): os.remove('trans.inp')
            if os.path.isfile('trans.out'): os.remove('trans.out')
            if os.path.isfile('trans.lib'): os.remove('trans.lib')
            arbitrary == False
            break
            #quit()
                    
        elif change == 'Vs1':
            difference = float(raw_input('How much do you want to change it by? '))
            Vs1 += difference
            print 'Vs1 = {0} m/s'.format(Vs1)
            
        elif change == 'Vs2':    
            difference = float(raw_input('How much do you want to change it by? '))
            Vs2 += difference
            print 'Vs2 = {0} m/s'.format(Vs2)
            
        elif change == 'ar':
            area_ratio = float(raw_input('What do you want the new nozzle area ratio to be? '))
            print 'New nozzle area ratio is {0}'.format(area_ratio)
            
        elif change == 'Vsd':
            difference = float(raw_input('How much do you want to change it by? '))
            Vsd += difference
            print 'Vsd = {0} m/s'.format(Vsd)
    
        elif change == 'p1':
            difference = float(raw_input('How much do you want to change it by? '))
            p1 += difference
            print 'p1 = {0} Pa'.format(p1)
            
        elif change == 'p5':    
            difference = float(raw_input('How much do you want to change it by? '))
            p5 += difference
            print 'p5 = {0} Pa'.format(p5)
                        
        elif change == 'psd1':
            difference = float(raw_input('How much do you want to change it by? '))
            psd1 += difference
            print 'psd1 = {0} m/s'.format(psd1)


#------------------------ running stuff----------------------------------------
                           
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "pitot - Equilibrium expansion tube simulator"
        print "start with --help for help with inputs"
        print "There are a few demos here that you can run to see how the program works:"
        print "demo_s: demo of Umar's high speed air condition where shock speeds are guessed to find the fill pressure's he used."
        print "demo_p: demo of Umar's high speed air condition where the fill pressure's are specified."
        print "hadas85-fulltheory: fully theoretical demo of Hadas' 8.5km/s titan condition where fill pressure's are specified."
        print "hadas85-tunnel: run through of Hadas' 8.5km/s titan condition with both fill pressure's and shock speed's from xpt specified."
        print "chrishe-fulltheory: full theoretical demo of one of my He conditions with secondary driver."
        print "dave-scramjet-p: a fully theoretical test of one of David Gildfind's scramjet conditions that iterates through fill pressures."
        print "dave-scramjet-s: a fully theoretical test of one of David Gildfind's scramjet conditions that iterates through shock speeds."
        print "dave-scramjet-tunnel: a recreation of one of Dave's scramjet conditions using test data."
        print "x3: basic x3 scramjet condition example."
        demo = raw_input("Type one of the demo commands for a demo run or anything else to quit ")
        if demo == 'demo_s':
            print "This is a demo of pitot recreating Umar Sheikh's high speed air condition where shock speeds are guessed to find the fill pressure's he used"
            print "Fill pressure's we are aiming for are p1 = 3000 Pa, p5 = 10 Pa"
            print " "
            sys.argv = ['pitot.py', '--test','fulltheory-shock', '--Vs1','5645.0',
                        '--Vs2','11600.0','--filename','demo_s.txt']
            main()
            
        elif demo == "demo_p":
            print "This is demo of pitot recreating Umar Sheikh's high speed air condition where fill pressure's are specified."
            print " "            
            sys.argv = ['pitot.py', '--p1','3000.0','--p5','10.0','--filename','demo_p.txt']
            main()
            
        elif demo == 'hadas85-fulltheory':
            print "This is the fully theoretical demo of pitot recreating Hadas' 8.5 km/s titan condition."
            print " "
            sys.argv = ['pitot.py', '--driver_gas','He:0.80,Ar:0.20','--test_gas','titan', 
                        '--p1','3200.0','--p5','10.0', '--ar','3.0', '--conehead','yes'
                        '--filename','hadas85-fulltheory']
            main()  
            
        elif demo == 'hadas85-tunnel':
            print "This is the demo of pitot recreating Hadas' 8.5 km/s titan condition with experimental shock speeds specified."
            print " "
            sys.argv = ['pitot.py', '--test','test', '--driver_gas','He:0.80,Ar:0.20',
                        '--test_gas','titan', '--p1','3200.0','--p5','10.0', 
                        '--Vs1','4100.0','--Vs2','8620.0','--ar','3.0', 
                        '--filename','hadas85-tunnel.txt']
            main()
            
        elif demo == 'chrishe-fulltheory':
            print "This is the demo of pitot recreating my 14 km/s 85%H2:15%He condition fully theoretically."
            print " "
            sys.argv = ['pitot.py', '--test','fulltheory-pressure', '--config',
                        'sec-nozzle','--driver_gas','He:1.0', '--test_gas',
                        'gasgiant_h215he','--psd1','17500','--p1','4700',
                        '--p5','4.0','--filename','chrishe-fulltheory.txt']
            main()
            
        elif demo == 'dave-scramjet-p':
            print "This is the demo of pitot recreating one of Dave Gildfind's scramjet conditions that iterates through fill pressures."
            print " "
            sys.argv = ['pitot.py', '--test','fulltheory-pressure', '--config',
                        'sec-nozzle','--driver_gas','He:0.80,Ar:0.20', '--test_gas',
                        'air','--psd1','100000','--p1','486000',
                        '--p5','1500.0','--filename','dave-scramjet-p.txt','--shock_switch','True']
            main()
            
                    
        elif demo == 'dave-scramjet-s':
            print "This is the demo of pitot recreating one of Dave Gildfind's scramjet conditions that iterates through shock speeds."
            print " "
            sys.argv = ['pitot.py', '--test','fulltheory-shock', '--config',
                        'sec-nozzle','--driver_gas','He:0.80,Ar:0.20', '--test_gas',
                        'air','--Vsd','4290.0','--Vs1','1588.0',
                        '--Vs2','3424.0','--filename','dave-scramjet-s.txt']
            main()
        
        elif demo == 'dave-scramjet-tunnel':
            print "This is the demo of pitot recreating one of Dave Gildfind's scramjet condition from test data."
            print " "
            sys.argv = ['pitot.py', '--test','test', '--config',
                        'sec','--driver_gas','He:0.80,Ar:0.20', '--test_gas',
                        'air','--Vsd','4178.0','--Vs1','1417.0',
                        '--Vs2','3264.0','--psd1','100000','--p1','690800','--p5','288.2',
                        '--filename','dave-scramjet-tunnel.txt']
            main()
            
        elif demo == 'x3':
            print "This is the demo of pitot recreating one of the basic x3 conditions."
            print " "
            sys.argv = ['pitot.py', '--facility', 'x3', '--test','fulltheory-pressure', '--config',
                        'sec','--driver_gas','He:0.60,Ar:0.40', '--test_gas',
                        'air','--psd1','133000','--p1','73000',
                        '--p5','210.0','--filename','x3.txt']
            main()
    
    else:
        main()

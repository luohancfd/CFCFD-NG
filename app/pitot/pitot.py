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
    
The full txt_output is a bit too much to include here, but you should see that the
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

You may then call upon pitot.py as long as you have suitable
enviroment variables set, as per the installation instructions
for Eilmer3.


Some History
------------
Since the dawn of time, Richard Morgan used the program pitot to perform
basic calculations and create new expansion tube conditions.
The code was originally written in GWBasic, and as an exercise in 
programming (and to hopefully make the life of my colleagues and I
easier) I learnt some basic GWBasic syntax and ported the program to
Python at the start of 2012. It has since become an ongoing project
and this latest version has completely rebuilt the code from the 
ground up as a new program that makes use of the classes and functions
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
    15-Jan-2013: added basic test time calculation from the work of Paull and Stalker (1992) to the code
    16-Jan-2013: added a second csv output to the code that should be handy for pulling the data into spreadsheets etc.
    17-Jan-2013: started implementing changes suggested by Peter Jacobs and Fabian Zander:
        Changed the variable 'state_over_model' to 'test_section_state' so it's less ambiguous.
        Did a general tidy up of my print statements and comments.
        Cut out the section at the end of the code where variables can be changed, the user can do this externally.
    22-Jan-2013: added perfect gas solver to the code so it can now do equilibrium and perfect gas calculations
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

#to do some perfect gas stuff

import cfpylib.gasdyn.ideal_gas as pg

VERSION_STRING = "22-Jan-2013"

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
    op.add_option('--mode', dest='mode', default='printout',
                 choices=['printout','return'],
                 help=("software mode; "
                        "printout = normal run, prints out a summary to the screen, a txt file and a csv file at the end of the program, then quits; "
                        "return = simpler run, useful if pitot is to be used inside a bigger program, returns a set of values at the end of the run, then quits; "
                        "defaults to printout "))
    op.add_option('--solver', dest='solver', default='eq',
                 choices=['eq','pg'],
                 help=("solver to use; "
                        "eq = equilibrium calculations using CEA code; "
                        "pg = perfect gas solver; "
                        "defaults to eq "))
    op.add_option('--facility', dest='facility', default='x2',
                  choices=['x2','x3'],
                  help=("facility to use; "
                        "x2 = the x2 expansion tube; "
                        "x3 = the x3 expansion tube; "
                        "defaults to x2 "))
    op.add_option('--test', dest='test', default='fulltheory-pressure',
                  choices=['fulltheory-shock','fulltheory-pressure','experiment'],
                  help=("type of test to run. "
                        "fulltheory-shock = fully theoretical run where fill pressures are found from set shock speeds "
                        "fulltheory-pressure = fully theoretical run where shock speeds are found from set fill pressures "
                        "experiment = partially theoretical run where both shock speeds and fill pressures are specified based on experimental data "))                   
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
    op.add_option('--shock_over_model', dest='shock_over_model', default=None,
                  help=("switch to calculate conditions along the stagnation streamline " 
                        "behind the normal shock over the test model "
                        "default is don't do it, type yes or anything else to make it happen "))   
    op.add_option('--filename',dest='filename',type='string',default='x2run',
                  help=("filename the result will be saved to "
                          "defaults to x2run"
                          "no need to specify .txt or .csv units as this is added by the code"))
    op.add_option('--expansion',dest='expansion',type='float',default='1.0',
                  help=("Use this to expand the test gas further into the acceleration tube "
                        "Basically just multiplies V7 by this percentage before the final expansion is done "
                        "Defaults to 1.0 "))
    op.add_option('--shock_switch',dest='shock_switch', default=None,
                  help=("Used to trigger a shock instead of an expansion from the secondary driver gas moving into the shock tube"))
    
    opt, args = op.parse_args()
    
    mode = opt.mode
    solver = opt.solver
    facility = opt.facility
    test = opt.test
    config = opt.config
    driver_gas = opt.driver_gas
    gasName = opt.gasName
    filename =opt.filename
    
    Vs1 = opt.Vs1
    Vs2 = opt.Vs2
    Vsd = opt.Vsd
    shock_over_model = opt.shock_over_model
    
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
            
        if PRINT_STATUS: print "Let's get started, shall we:"
        if PRINT_STATUS: 
            if secondary:
                print 'Selected Vsd = {0} m/s'.format(Vsd)
            print 'Selected Vs1 = {0} m/s'.format(Vs1)
            print 'Selected Vs2 = {0} m/s'.format(Vs2)

        
    elif test == 'fulltheory-pressure':  #shows we're going to solve shock speeds from fill pressures
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

        if PRINT_STATUS: print "Let's get started, shall we:"
        if PRINT_STATUS: 
            if secondary:
                print 'Selected secondary driver fill pressure (psd1) = {0} Pa.'.format(psd1)

            print 'Selected shock tube fill pressure (p1) = {0} Pa.'.format(p1)
            print 'Selected acceleration tube fill pressure (p5) = {0} Pa.'.format(p5)
            
    elif test == 'experiment': #test based on real experimental data
        
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
            
        if PRINT_STATUS: print "Let's get started, shall we:"
        if PRINT_STATUS: 
            if secondary:
                print 'Selected Vsd = {0} m/s'.format(Vsd)
            print 'Selected Vs1 = {0} m/s'.format(Vs1)
            print 'Selected Vs2 = {0} m/s'.format(Vs2)
            if secondary:
                print 'Selected secondary driver fill pressure (psd1) = {0} Pa.'.format(psd1)
            print 'Selected shock tube fill pressure (p1) = {0} Pa.'.format(p1)
            print 'Selected acceleration tube fill pressure (p5) = {0} Pa.'.format(p5)
            
    else: #test input must be bad
        print "You haven't specified a relevant test condition. Bailing out..."
        bad_input = True
    
        if bad_input: #bail out here if you end up having issues with your input
            return -2

  
#---------------- building initial states ---------------------------------- 
    #here the states are made as CEA2 gas states, and then the ideal gam and MW are pulled out if required   
    #and the gas objects are redefined to be perfect gas
       
    states = {} #states dictionary that we'll fill up later
    V = {} #same for velocity
    M = {} #same for Mach number
    
    
    #state 4 is diaphragm burst state (taken to basically be a total condition)
    
    if facility == 'x2':
        states['s4']=primary_driver_x2[driver_gas][0].clone()
        p4 = 2.79e7; T4 = 2700.0 #Pa, K
        states['s4'].set_pT(p4,T4)
        V['s4']=0.0
        M['s4']=0.0
        
        if solver == 'pg': #make perfect gas object, and then re-set the state
            states['s4']=pg.Gas(Mmass=states['s4'].Mmass,
                                        gamma=states['s4'].gam, name='s4')
            states['s4'].set_pT(p4,T4)
    
        M['s3s']=primary_driver_x2[driver_gas][1]
        
    elif facility == 'x3':
        states['s4']=primary_driver_x3[driver_gas][0].clone()
        p4 = 2.79e7; T4 = 2700.0 #Pa, K
        states['s4'].set_pT(p4,T4)
        V['s4']=0.0
        M['s4']=0.0
        
        if solver == 'pg': #make perfect gas object if asked to do so, and then re-set the state
            states['s4']=pg.Gas(Mmass=states['s4'].Mmass,
                                        gamma=states['s4'].gam, name='s4')
            states['s4'].set_pT(p4, T4)
    
        M['s3s']=primary_driver_x3[driver_gas][1]

    #state3s is driver gas after steady expansion at the throat between 
    #the primary driver and the next section
        
    (states['s3s'], V['s3s']) = expand_from_stagnation(1.0/(p0_p(M['s3s'],states['s4'].gam)),states['s4'])
    
    #build the gas objects for all the other sections based on knowledge of what is what
    
    T0 = 300.0 #atmospheric temperature (K), used for any section starting at ambient
    p0 = 101300.0 #atmospheric pressure (Pa)
    
    #start with shock tube and acc tube, start with atmospheric p and T
    
    if secondary: #state sd1 is pure He secondary driver
        states['sd1'] =  Gas({'He':1.0,},outputUnits='moles')
        if not psd1: #set atmospheric state if a pressure was not specified
            psd1 = p0
        states['sd1'].set_pT(psd1,T0)
        if solver == 'pg': #make perfect gas object if asked to do so, then re-set the gas state
            states['sd1']=pg.Gas(Mmass=states['sd1'].Mmass,
                                 gamma=states['sd1'].gam, name='sd1')
            states['sd1'].set_pT(psd1,T0)
        V['sd1']=0.0
        M['sd1']=0.0

    #state 1 is shock tube
    states['s1'], gas_guess = make_test_gas(gasName)
    if not p1: #set atmospheric state if a pressure was not specified
        p1 = p0
    states['s1'].set_pT(p1,T0)
    if solver == 'pg': #make perfect gas object if asked to do so, then re-set the gas state
        states['s1']=pg.Gas(Mmass=states['s1'].Mmass,
                            gamma=states['s1'].gam, name='s1')
        states['s1'].set_pT(p1,T0)
    V['s1']=0.0
    M['s1']=0.0
        

    #state 5 is acceleration tube
    states['s5'] = Gas({'Air':1.0,},outputUnits='moles')
    if not p5: #set atmospheric state if a pressure was not specified
        p5 = p0
    states['s5'].set_pT(p5,T0)
    if solver == 'pg': #make perfect gas object if asked to do so, then re-set the gas state
        states['s5']=pg.Gas(Mmass=states['s5'].Mmass,
                            gamma=states['s5'].gam, name='s5')
        states['s5'].set_pT(p5,T0)
    V['s5']=0.0
    M['s5']=0.0
    #now let's clone the states we just defined to get the states derved from these
   
    if secondary:
        states['sd2'] = states['sd1'].clone() #sd2 is sd1 shock heated
        states['sd3'] = states['s3s'].clone() #sd3 is s3s after unsteady expansion
        states['s3'] = states['sd2'].clone() #s3 will be sd2 after unsteady expansion
    else:
        states['s3'] = states['s3s'].clone() #s3 will be s3s after unsteady expansion
        
    states['s2'] = states['s1'].clone() #s2 is s1 shock heated
    states['s6'] = states['s5'].clone() #s6 is s5 shock heated
    states['s7'] = states['s2'].clone() #7s is s2 after unsteady expansion   
       
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
    
        #get mach numbers for the txt_output
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
                p1 = secant(p1_reflected_iterator, 1000.0, 100000.0, tol=1.0e-4,limits=[100.0,1000000.0])
            else: #same expansion as we're used to
                if PRINT_STATUS: print "Start unsteady expansion of the secondary driver gas into the shock tube."
                p1 = secant(error_in_velocity_shock_tube_expansion_pressure_iterator, 1000.0, 100000.0, tol=1.0e-4,limits=[100.0,1000000.0])
        else: #just do the expansion
            if PRINT_STATUS: print "Start unsteady expansion of the driver gas into the shock tube."
            p1 = secant(error_in_velocity_shock_tube_expansion_pressure_iterator, 1000.0, 100000.0, tol=1.0e-4,limits=[100.0,1000000.0])

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
            Vs1 = secant(primary_shock_speed_reflected_iterator, 2000.0, 4000.0, tol=1.0e-4,limits=[500.0,1000000.0])
        else: #just do the expansion
            Vs1 = secant(error_in_velocity_shock_tube_expansion_shock_speed_iterator, 4000.0, 6000.0, tol=1.0e-3,limits=[1000.0,1000000.0])
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
    
    #get mach numbers for the txt_output
    Ms1 = Vs1/states['s1'].son
    M['s2'] = V['s2']/states['s2'].son
    M['s3']= V['s3']/states['s3'].son
    
    #do ideal gas guess for gas constants if needed
    
    if p5: #just use the gas guess if you don't know
        ideal_gas_guess = gas_guess
        ideal_gas_guess_air = {'gam':1.35,'R':571.49}
    elif Vs2 >= 9000.0:
        ideal_gas_guess = gas_guess
        ideal_gas_guess_air = {'gam':1.35,'R':571.49}
    else:
        ideal_gas_guess = None
        ideal_gas_guess_air = None
    
    if solver == 'pg': #if using the pg solver the ideal gas guess causes problems, so just get rid of it!
        ideal_gas_guess = None
        ideal_gas_guess_air = None

#--------------------- acceleration tube -----------------------------------
                
    if PRINT_STATUS: print "Start unsteady expansion of the test gas into the acceleration tube."

#----------------------- acceleration tube functions -----------------------

    def error_in_velocity_s2_expansion_pressure_iterator(p5, state2=states['s2'], 
                                       V2g=V['s2'], state5=states['s5'],
                                        state6=states['s6'],Vs2=Vs2,
                                        ideal_gas_guess=ideal_gas_guess_air):
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
                                        ideal_gas_guess=ideal_gas_guess_air):
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
            Vs2 = secant(error_in_pressure_s2_expansion_shock_speed_iterator, Vs1+7000.0, 15100.0, tol=1.0e-5,limits=[Vs1+1000.0,100000.0])
        else:
            Vs2 = secant(error_in_pressure_s2_expansion_shock_speed_iterator, Vs1+7000.0, 15100.0, tol=1.0e-5,limits=[Vs1+1000.0,100000.0])
        if PRINT_STATUS: print "From secant solve: Vs2 = {0} m/s".format(Vs2)
        #start using Vs1 now, compute states 1,2 and 3 using the correct Vs1
        if PRINT_STATUS: print "Once Vs2 is known, find conditions at state 5 and 6."            
        
    (V6, V['s6']) = normal_shock(states['s5'], Vs2, states['s6'],ideal_gas_guess)
    V['s7'], states['s7'] = finite_wave_dv('cplus', V['s2'], states['s2'], (V['s6']*expansion))
    
    #get mach numbers for the txt_output
    Ms2 = Vs2/states['s5'].son
    M['s6'] = V['s6']/states['s6'].son
    M['s7']= V['s7']/states['s7'].son
    
#------------------- finishing off any other needed calculations -------------
                      
    #-------------- do the nozzle calc up to s8 now --------------------
    if nozzle:
        if PRINT_STATUS: print "Start steady expansion through the nozzle."
        (V['s8'], states['s8']) = steady_flow_with_area_change(states['s7'], V['s7'], area_ratio)
        M['s8']= V['s8']/states['s8'].son
        test_section_state = 's8'
    else:
        test_section_state = 's7'             
    
    #-------------- do normal shock over model if asked to --------------------------
    
    if shock_over_model:
        if PRINT_STATUS: print "Start frozen normal shock calculation over the test model."  
        states['s10f'] = states[test_section_state].clone()
        (V10, V['s10f']) = shock_ideal(states[test_section_state], V[test_section_state], states['s10f'])
        M['s10f']= V['s10f']/states['s10f'].son

        if PRINT_STATUS: print "Start equilibrium normal shock calculation over the test model."  
        if solver == 'eq': 
            states['s10e'] = states[test_section_state].clone()
            (V10, V['s10e']) = normal_shock(states[test_section_state], V[test_section_state], states['s10e'])
            M['s10e']= V['s10e']/states['s10e'].son
        elif solver == 'pg': #we need to make a cea2 gas object to do this equilibrium calculaiton if every other gas object is pg
            states[test_section_state+'eq'] = make_test_gas(gasName)[0]
            states[test_section_state+'eq'].set_pT(states[test_section_state].p,states[test_section_state].T)
            states['s10e'] = states[test_section_state+'eq'].clone()
            (V10, V['s10e']) = normal_shock(states[test_section_state+'eq'], V[test_section_state], states['s10e'])
            M['s10e']= V['s10e']/states['s10e'].son
    
    #------------- do conehead calculations if asked to do so -----------
    
    if conehead:
        if PRINT_STATUS: print 'Doing taylor maccoll conehead calculation on 15 degree conehead.'
        shock_angle = beta_cone(states[test_section_state], V[test_section_state], math.radians(15))
        if PRINT_STATUS: print "\nShock angle over cone:", math.degrees(shock_angle)
        # Reverse the process to get the flow state behind the shock and check the surface angle is correct
        delta_s, V['s10c'], states['s10c'] = theta_cone(states[test_section_state], V[test_section_state], shock_angle)
        M['s10c'] = V['s10c']/states['s10c'].son
        if PRINT_STATUS: print "Surface angle should be the same.....: 15deg = ", math.degrees(delta_s), "deg"
        if PRINT_STATUS: print "\nConehead surface conditions:"
        if PRINT_STATUS: states['s10c'].write_state(sys.stdout)
        # Need to check whether the pressure are the same
        if PRINT_STATUS: print "Computed conehead pressure is {0} Pa".format(states['s10c'].p)
        
    #-------------- if the normal shock thing was done, fix it up before print out ----------

    if test == "fulltheory-pressure" and shock_switch: #restore the fill pressure's back to normal
        states['s1'].p = states['s1'].p*2.0
        states['s5'].p = states['s5'].p*2.0
        
    if test == "fulltheory-shock" and shock_switch: #cut fill pressure's in half at the end
        states['s1'].p = states['s1'].p*2.0
        states['s5'].p = states['s5'].p*2.0
        
    #----------- test time calculations -------------------------------
    
    #This was based off code written by David Gildfind that was based off a paper
    #by Allan Paull and Ray Stalker. I've only considered the basic test time case here
    # A. Paull & R. J. Stalker, "Test Flow Disturbances in an Expansion Tube", J. Fluid Mech. (1992), vol. 245, pp. 493-521 (p499).
    
    
    if facility == 'x2':
        #a few different tunnel scenarios give different lengths
        #all the distances are taken from my L1D run file, written by David Gildfind
        #0m is the primary diaphragm burst location
        if config == 'basic':
            distances = [3.418, 8.979] #secondary diaphragm, and then end of acc tube (m)
        elif config == 'sec':
            distances = [3.418, 5.976, 8.979] #secondary diaphragm, tertiary diaphragm, end of acc tube (m)
        elif config == 'nozzle':
            distances = [3.418, 8.585] #secondary diaphragm, entrance to nozzle (m)
        elif config == 'sec-nozzle':
            distances = [3.418, 5.976, 8.585] #secondary diaphragm, tertiary, entrance to nozzle (m)
    
        t_start = 0.0 #time of primary diaphragm burst
    
        if secondary: #if secondary we have a third tube to consider
            #calculate some lengths
            L_sec_drv = distances[0]
            L_shk_tube = distances[1] - distances[0]
            L_acc_tube = distances[2] - distances[1]
            
            #shocks
            t_inc_shock_sd = t_start + L_sec_drv/Vsd
            t_inc_shock_st = t_inc_shock_sd + L_shk_tube/Vs1
            t_inc_shock_at = t_inc_shock_st + L_acc_tube/Vs2
            
            #contact surfaces
            t_cs_sd = t_start + L_sec_drv/V['sd3']
            t_cs_st = t_inc_shock_sd + L_shk_tube/V['s3']
            t_cs_at = t_inc_shock_st + L_acc_tube/V['s7']

        
        else: #so we're going straight from primary driver to shock tube now
            #calculate some lengths
            L_shk_tube = distances[0]
            L_acc_tube = distances[1] - distances[0]
            
            #shocks
            t_inc_shock_st = t_start + L_shk_tube/Vs1
            t_inc_shock_at = t_inc_shock_st +  L_acc_tube/Vs2
            
            #contact surfaces
            t_cs_sd = t_start + L_shk_tube/V['s3']
            t_cs_at = t_inc_shock_st +  L_acc_tube/V['s7']
            
        
        #now we can actually calculate the basic test time...
        # (borrowed from the procedure David Gildfind used in his PhD)
        
        #now we just need to calculate the time taken for the unsteady 
        #expansion to hit the end of the acc tube
        
        t_final_usx = t_inc_shock_st + L_acc_tube/(V['s7']-states['s7'].son)
        
        t_test_basic = t_final_usx - t_cs_at
    
    if mode == 'printout':
     
    #--------------------------- txt_output --------------------------------
    
        txt_output = open(filename+'.txt',"w")  #txt_output file creation
                    
        print " "
        
        version_printout = "Pitot Version: {0}".format(VERSION_STRING)
        print version_printout
        txt_output.write(version_printout + '\n')
        
        if secondary:
            description_sd = 'sd1 is secondary driver fill.'
            print description_sd
            txt_output.write(description_sd + '\n')   
            
        description_1 = 'state 1 is shock tube fill. state 5 is acceleration tube fill.'    
        print description_1
        txt_output.write(description_1 + '\n')
            
        description_2 = 'state 7 is expanded test gas entering the nozzle.'
        print description_2
        txt_output.write(description_2 + '\n')
            
        description_3 = 'state 8 is test gas exiting the nozzle (using area ratio of {0}).'.format(area_ratio)
        print description_3
        txt_output.write(description_3 + '\n')
        
        if shock_over_model:    
            description_4 = 'state 10f is frozen shocked test gas flowing over the model.'
            print description_4
            txt_output.write(description_4 + '\n')  
        
            description_5 = 'state 10e is equilibrium shocked test gas flowing over the model.'
            print description_5
            txt_output.write(description_5 + '\n')
            
        if solver == 'eq':
            solver_printout = "Solver used is equilibrium."
        elif solver == 'pg':
            solver_printout = "Solver used is perfect gas."
        print solver_printout
        txt_output.write(solver_printout + '\n')        
            
        facility_used = 'Facility is {0}.'.format(facility)        
        print facility_used
        txt_output.write(facility_used + '\n')
        if solver == 'eq':
            test_gas_used = 'Test gas is {0} (gamma = {1}, R = {2}, {3}).'.format(gasName,states['s1'].gam,states['s1'].R,states['s1'].reactants)
        elif solver == 'pg':
            test_gas_used = 'Test gas is {0} (gamma = {1}, R = {2}).'.format(gasName,states['s1'].gam,states['s1'].R)
        print test_gas_used
        txt_output.write(test_gas_used + '\n')  
        
        driver_gas_used = 'Driver gas is {0}.'.format(driver_gas)       
        print driver_gas_used
        txt_output.write(driver_gas_used + '\n') 
                
        if shock_switch:
            shock_warning1 = "NOTE: a reflected shock was done into the shock tube and as such,"
            print shock_warning1
            txt_output.write(shock_warning1 + '\n')
            
            shock_warning2 = "fill pressure's have been artifically modified by the code to match with xpt."
            print shock_warning2
            txt_output.write(shock_warning2 + '\n')
    
        if secondary:
            secondary_shockspeeds = "Vsd = {0:.2f} m/s, Msd1 = {1:.2f}".format(Vsd,Msd1)
            print secondary_shockspeeds
            txt_output.write(secondary_shockspeeds + '\n')
        
        shockspeeds = "Vs1 = {0:.2f} m/s, Ms1 = {1:.2f} ,Vs2 = {2:.2f} m/s, Ms2 = {3:.2f}".format(Vs1,Ms1,Vs2,Ms2) 
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
                    pitot[it_string] = 0
                    p0[it_string] = 0
                else:
                    pitot[it_string] = pitot_p(states[it_string].p,M[it_string],states[it_string].gam)/1000.0
                    #make total condition of relevant state for printing
                    total_state = total_condition(states[it_string], V[it_string])
                    p0[it_string] = total_state.p/1.0e6
                
                if states[it_string].p < 1.0e6: #change how the pressure is printed if it's too big, it keeps ruining the printouts!
                    conditions = "{0:<6}{1:<11.7}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.5f}{7:<7.0f}{8:<9.1f}"\
                    .format(it_string, states[it_string].p, states[it_string].T,
                            states[it_string].son,V[it_string],M[it_string],
                            states[it_string].rho, pitot[it_string], p0[it_string])
                else:
                    conditions = "{0:<6}{1:<11.3e}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.5f}{7:<7.0f}{8:<9.1f}"\
                    .format(it_string, states[it_string].p, states[it_string].T,
                            states[it_string].son,V[it_string],M[it_string],
                            states[it_string].rho, pitot[it_string], p0[it_string])
                        
                print conditions
                txt_output.write(conditions + '\n')
    
                
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
            
        #do the conditions over the model if asked
        if shock_over_model:
            condition_printer('s10f')
            condition_printer('s10e')
                
        if conehead:
            condition_printer('s10c')
                                       
        #some other useful calculations at the end
              
        total = total_condition(states[test_section_state], V[test_section_state])
        
        stagnation_enthalpy = total.h #J/kg
        if nozzle:        
            stag_enth = 'The total enthalpy (Ht) leaving the nozzle is {0:<.5g} MJ/kg.'.format(stagnation_enthalpy/10**6)
        else:
            stag_enth = 'The total enthalpy (Ht) at the end of the acceleration tube {0:<.5g} MJ/kg.'.format(stagnation_enthalpy/10**6)
        print stag_enth
        txt_output.write(stag_enth + '\n')
        
        #calculate flight equivalent velocity
        #for a description of why this is, refer to Bianca Capra's thesis page 104 - 105
        #Capra, B., Aerothermodynamic Simulation of Subscale Models of the FIRE II and
        #Titan Explorer Vehicles in Expansion Tubes, Ph.D. thesis, the University of Queens-
        #land, St. Lucia, Australia, 2006.
        u_eq = math.sqrt(2.0*stagnation_enthalpy) 
        u_eq_print = 'The flight equivalent velocity (Ue) is {0:<.5g} m/s.'.format(u_eq)
        print u_eq_print
        txt_output.write(u_eq_print + '\n')
        
        #if the test time calculation has been done, print it
        if t_test_basic: 
            basic_test_time_printout = 'Basic test time = {0:.2f} microseconds'.format(t_test_basic*1.0e6)
            print  basic_test_time_printout
            txt_output.write(basic_test_time_printout + '\n')
        
        #added ability to get the species in the post-shock condition
        
        if shock_over_model:
            species1 = 'species in the shock layer at equilibrium:'        
            print species1
            txt_output.write(species1 + '\n')
            
            species2 = '{0}'.format(states['s10e'].species)
            print species2
            txt_output.write(species2 + '\n')
        
        #added taylor maccoll conehead calcs if people want to use them, appropriated from Fabs (thankyou)
                    
        txt_output.close()
        
    #------------------------- cut down csv output -----------------------------
        
        csv_output = open(filename+'.csv',"w")  #csv_output file creation
              
        csv_version_printout = "Pitot Version,{0}".format(VERSION_STRING)
        csv_output.write(csv_version_printout + '\n')
        
        if solver == 'eq':
            csv_solver_printout = "Solver,equilibrium."
        elif solver == 'pg':
            csv_solver_printout = "Solver,perfect gas"
        csv_output.write(csv_solver_printout + '\n')     
            
        csv_facility_used = 'Facility,{0}.'.format(facility)        
        csv_output.write(csv_facility_used + '\n')
        
        csv_test_gas_used = 'Test gas,{0},gamma,{1},R,{2}'.format(gasName,states['s1'].gam,states['s1'].R)
        csv_output.write(csv_test_gas_used + '\n')  
        
        csv_driver_gas_used = 'Driver gas,{0}.'.format(driver_gas)       
        csv_output.write(csv_driver_gas_used + '\n') 
                
        if secondary:
            csv_secondary_shockspeeds = "Vsd,{0:.2f} m/s,Msd1,{1:.2f}".format(Vsd,Msd1)
            csv_output.write(csv_secondary_shockspeeds + '\n')
        
        csv_shockspeeds = "Vs1,{0:.2f} m/s,Ms1,{1:.2f},Vs2,{2:.2f} m/s,Ms2,{3:.2f}".format(Vs1,Ms1,Vs2,Ms2) 
        csv_output.write(csv_shockspeeds + '\n')
                
        csv_key = "{0:6},{1:11},{2:9},{3:6},{4:9},{5:6},{6:9},{7:8},{8:9}".format("state","P","T","a","V","M","rho","pitot","stgn")
        csv_output.write(csv_key + '\n')
        
        csv_units = "{0:6},{1:11},{2:9},{3:6},{4:9},{5:6},{6:9},{7:9},{8:9}".format("","Pa","K","m/s","m/s","","m^3/kg","kPa","MPa")
        csv_output.write(csv_units + '\n')
               
        def csv_condition_printer(it_string):
            """Prints the values of a specified condition to the screen and to 
            the txt_output file. 
            
            I made a function of this so I didn't have to keep pasting the code in."""
            
            if states.has_key(it_string):
                    
                if M[it_string] == 0:
                    pitot[it_string] = 0
                    p0[it_string] = 0
                else:
                    pitot[it_string] = pitot_p(states[it_string].p,M[it_string],states[it_string].gam)/1000.0
                    p0[it_string] = p0_p(M[it_string], states[it_string].gam)*states[it_string].p/1.0e6
                
                csv_conditions = "{0:<6},{1:<11.7},{2:<9.1f},{3:<6.0f},{4:<9.1f},{5:<6.2f},{6:<9.4f},{7:<8.0f},{8:<9.1f}"\
                .format(it_string, states[it_string].p, states[it_string].T,
                        states[it_string].son,V[it_string],M[it_string],
                        states[it_string].rho, pitot[it_string], p0[it_string])
    
                csv_output.write(csv_conditions + '\n')
    
                
        #print the driver related stuff first
        
        csv_condition_printer('s4')
        csv_condition_printer('s3s')
        
        if secondary: #need a separate printing thing for the secondary driver
            
            for i in range(1,4): #will do 1 - 3
            
                it_string = 'sd{0}'.format(i)
                csv_condition_printer(it_string)
                        
        for i in range(1,4): #shock tube stuff
            
            it_string = 's{0}'.format(i)
            csv_condition_printer(it_string)
            
        for i in range(5,9): #acc tube and nozzle if it's there
        
            it_string = 's{0}'.format(i)
            csv_condition_printer(it_string)
            
        #do the conditions over the model
        if shock_over_model:
            csv_condition_printer('s10f')
            csv_condition_printer('s10e')
                
        if conehead:
            csv_condition_printer('s10c')
                                       
    
        csv_stag_enth = 'Ht,{0:<.5g} MJ/kg.'.format(stagnation_enthalpy/10**6)
        csv_output.write(csv_stag_enth + '\n')
        
        csv_u_eq_print = 'Ue,{0:<.5g} m/s.'.format(u_eq)
        csv_output.write(csv_u_eq_print + '\n')
    
        if t_test_basic: 
            csv_basic_test_time_printout = 'Basic test time,{0:.2f} microseconds'.format(t_test_basic*1.0e6)
            csv_output.write(csv_basic_test_time_printout + '\n')   
        
        csv_output.close()
         
        if solver == 'eq': 
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
        
    elif mode == 'return': #return a few values and then quit
        if solver == 'eq':
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
        
        if secondary:
            return states, V, M, Vsd, Vs1, Vs2
        else:
            return states, V, M, Vs1, Vs2

                    

#------------------------ running stuff----------------------------------------
                           
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "pitot - Equilibrium expansion tube simulator"
        print "start with --help for help with inputs"
        print "There are a few demos here that you can run to see how the program works:"
        print "demo-s-eq: demo of Umar's high speed air condition where shock speeds are guessed to find the fill pressure's he used."
        print "demo-p-eq: equilibrium demo of Umar's high speed air condition where the fill pressure's are specified."
        print "demo-p-pg: perfect gas demo of Umar's high speed air condition where the fill pressure's are specified."
        print "hadas85-full-theory-eq: equilibrium fully theoretical demo of Hadas' 8.5km/s titan condition where fill pressure's are specified."
        print "hadas85-full-theory-pg: perfect gas fully theoretical demo of Hadas' 8.5km/s titan condition where fill pressure's are specified."
        print "hadas85-experiment-eq: equilibrium run through of Hadas' 8.5km/s titan condition with both fill pressure's and shock speed's from xpt specified."
        print "hadas85-experiment-pg: perfect gas run through of Hadas' 8.5km/s titan condition with both fill pressure's and shock speed's from xpt specified."
        print "chrishe-full-theory-eq: fully theoretical equilibrium demo of one of my He conditions with secondary driver."
        print "chrishe-full-theory-pg: fully theoretical perfect gas demo of one of my He conditions with secondary driver."
        print "dave-scramjet-p: a fully theoretical test of one of David Gildfind's scramjet conditions that iterates through fill pressures."
        print "dave-scramjet-s: a fully theoretical test of one of David Gildfind's scramjet conditions that iterates through shock speeds."
        print "dave-scramjet-tunnel: a recreation of one of Dave's scramjet conditions using test data."
        print "x3: basic x3 scramjet condition example."
        demo = raw_input("Type one of the demo commands for a demo run or anything else to quit ")
        
        if demo == 'demo-s-eq':
            print "This is an equilibrium demo of pitot recreating Umar Sheikh's high speed air condition where shock speeds are guessed to find the fill pressure's he used"
            print "Fill pressure's we are aiming for are p1 = 3000 Pa, p5 = 10 Pa"
            print " "
            sys.argv = ['pitot.py', '--solver','eq','--test','fulltheory-shock', '--Vs1','5645.0',
                        '--Vs2','11600.0','--filename',demo]
            main()
            
        elif demo == "demo-p-eq":
            print "This is equilibrium demo of pitot recreating Umar Sheikh's high speed air condition where fill pressure's are specified."
            print " "            
            sys.argv = ['pitot.py', '--p1','3000.0','--p5','10.0','--filename',demo]
            main()
            
        elif demo == "demo-p-pg":
            print "This is perfect gas demo of pitot recreating Umar Sheikh's high speed air condition where fill pressure's are specified."
            print " "            
            sys.argv = ['pitot.py', '--solver','pg', '--p1','3000.0','--p5','10.0','--filename',demo]
            main()
            
        elif demo == 'hadas85-full-theory-eq':
            print "This is the equilibrium fully theoretical demo of pitot recreating Hadas' 8.5 km/s titan condition."
            print " "
            sys.argv = ['pitot.py', '--driver_gas','He:0.80,Ar:0.20','--test_gas','titan', 
                        '--p1','3200.0','--p5','10.0', '--ar','3.0', '--conehead','yes'
                        '--filename',demo]
            main() 
            
        elif demo == 'hadas85-full-theory-pg':
            print "This is the perfect gas fully theoretical demo of pitot recreating Hadas' 8.5 km/s titan condition."
            print " "
            sys.argv = ['pitot.py', '--driver_gas','He:0.80,Ar:0.20','--test_gas','titan', 
                        '--p1','3200.0','--p5','10.0', '--ar','3.0', '--conehead','yes'
                        '--filename',demo, '--solver','pg']
            main()  
            
        elif demo == 'hadas85-experiment-eq':
            print "This is the equilibrium demo of pitot recreating Hadas' 8.5 km/s titan condition with experimental shock speeds specified."
            print " "
            sys.argv = ['pitot.py', '--test','experiment', '--driver_gas','He:0.80,Ar:0.20',
                        '--test_gas','titan', '--p1','3200.0','--p5','10.0', 
                        '--Vs1','4100.0','--Vs2','8620.0','--ar','3.0', 
                        '--filename',demo,'--solver','eq']
            main()
            
        elif demo == 'hadas85-experiment-pg':
            print "This is the perfect gas demo of pitot recreating Hadas' 8.5 km/s titan condition with experimental shock speeds specified."
            print " "
            sys.argv = ['pitot.py', '--test','experiment', '--driver_gas','He:0.80,Ar:0.20',
                        '--test_gas','titan', '--p1','3200.0','--p5','10.0', 
                        '--Vs1','4100.0','--Vs2','8620.0','--ar','3.0', 
                        '--filename',demo,'--solver','pg']
            main()
            
        elif demo == 'chrishe-full-theory-eq':
            print "This is the equilibrium demo of pitot recreating my 14 km/s 85%H2:15%He condition fully theoretically."
            print " "
            sys.argv = ['pitot.py', '--test','fulltheory-pressure', '--config',
                        'sec-nozzle','--driver_gas','He:1.0', '--test_gas',
                        'gasgiant_h215he','--psd1','17500','--p1','4700',
                        '--p5','6.37','--filename',demo, '--shock_over_model','True']
            main()
            
        elif demo == 'chrishe-full-theory-pg':
            print "This is the perfect gas demo of pitot recreating my 14 km/s 85%H2:15%He condition fully theoretically."
            print " "
            sys.argv = ['pitot.py', '--solver','pg', '--test','fulltheory-pressure', '--config',
                        'sec-nozzle','--driver_gas','He:1.0', '--test_gas',
                        'gasgiant_h215he','--psd1','17500','--p1','4700',
                        '--p5','6.37','--filename',demo, '--shock_over_model','True']
            main()
            
        elif demo == 'dave-scramjet-p':
            print "This is the demo of pitot recreating one of Dave Gildfind's scramjet conditions that iterates through fill pressures."
            print " "
            sys.argv = ['pitot.py', '--test','fulltheory-pressure', '--config',
                        'sec-nozzle','--driver_gas','He:0.80,Ar:0.20', '--test_gas',
                        'air','--psd1','100000','--p1','486000',
                        '--p5','1500.0','--filename',demo,'--shock_switch','True']
            main()
            
                    
        elif demo == 'dave-scramjet-s':
            print "This is the demo of pitot recreating one of Dave Gildfind's scramjet conditions that iterates through shock speeds."
            print " "
            sys.argv = ['pitot.py', '--test','fulltheory-shock', '--config',
                        'sec-nozzle','--driver_gas','He:0.80,Ar:0.20', '--test_gas',
                        'air','--Vsd','4290.0','--Vs1','1588.0',
                        '--Vs2','3424.0','--filename',demo]
            main()
        
        elif demo == 'dave-scramjet-tunnel':
            print "This is the demo of pitot recreating one of Dave Gildfind's scramjet condition from test data."
            print " "
            sys.argv = ['pitot.py', '--test','test', '--config',
                        'sec','--driver_gas','He:0.80,Ar:0.20', '--test_gas',
                        'air','--Vsd','4178.0','--Vs1','1417.0',
                        '--Vs2','3264.0','--psd1','100000','--p1','690800','--p5','288.2',
                        '--filename',demo]
            main()
            
        elif demo == 'x3':
            print "This is the demo of pitot recreating one of the basic x3 conditions."
            print " "
            sys.argv = ['pitot.py', '--facility', 'x3', '--test','fulltheory-pressure', '--config',
                        'sec','--driver_gas','He:0.60,Ar:0.40', '--test_gas',
                        'air','--psd1','133000','--p1','73000',
                        '--p5','210.0','--filename',demo]
            main()
    
    else:
        main()

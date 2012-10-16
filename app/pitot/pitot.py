from __future__ import division #makes everything fp division

"""
import of RGM's Pitot program from the ancient GWBasic programming language to Python by Chris James (41429341, uqcjame4)

programmer contact: Chris James (c.james4@uq.edu.au)

thanks to Fabs for the idea, and PJ for helping me get to the source code

original code completed 9/1/12. Nozzle code written by Chris James (based on
NACA compressible flow tables) added the same day.

V2 with better printing created by Chris James start of April 2012.

V3 (the interactive version), created by Chris James late April 2012. V3 based on
knowledge I gained while completing the CSSE1001 course.

V4 I fixed up some little parts of the nozzle flow code

V5 got the 10%He test condition working, added stagnation enthalpy out of nozzle calculation

V6 added a much better mach number finder for the nozzle code, should be much more accurate
and not crash sometimes now. Serves me right for reusing crap code I wrote a year ago
also added a couple of little calculations done by Elise Fahy (e.fahy@uq.edu.au) Chris James 28/08/12

V7 I chose to try and start trying to make pitot fall in line with some of the
stuff in cfcfd's cfpylib. Should make this program a lot lighter and easier
to follow for my colleagues and I. October 2012.

1 preshot conditions in primary driver
2 shock heated secondary driver gas
3 expanded primary driver after unsteady expansion 
4 pre shot conditions in reservoir
5 pre shot conditions in shock tube
6 shock heated shock tube gas
7 expanded secondary driver gas
8 test gas after shock heating and unsteady expansion
9 shock heated accelerator gas
10 pre shot condition in acceleration tube
11 start of unsteady expansion (used for different configurations)
12 
13
14
15 test flow at the end of the nozzle (using area ratio of 2.5) (new results line written by Chris James)

NEW GASES CAN BE ADDED BY ADDING THEIR GAMMA AND R VALUES TO THE DICTIONARY 'gases'

"""

#this import pulled from estcj with some additions

import sys, os, math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from cfpylib.nm.zero_solvers import secant
# We base our calculation of gas properties upon calls to the NASA Glenn CEA code.
from cfpylib.gasdyn.cea2_gas import Gas, make_gas_from_name
from cfpylib.gasdyn.gas_flow import *
from cfpylib.gasdyn.ideal_gas_flow import *

#define any functions needed
      
def is_valid(command,valid_commands):
    """Prompts the user for a valid command, and only returns a valid command.

    is_valid_command(float,tuple<floats>) -> float

    Precondtions: input must be a single number."""
    
    check = False
    
    while check == False:
    
        for i in valid_commands:
            
        #print i
            
            if i == command:
                
                check = True
            
        if check == False:
             print 'That is an invalid command. Valid commands are {0}'.format(valid_commands)
             command = str(raw_input('Try again? '))
                
    return command

    
    
#------- property dictionaries ----------

#(I used the 4th edition of Cengel and Boles' Thermodynamics textbook
#to get the majority of these gas values. I calculated mixture values using
#my own gas mixture program that uses elementary physics to get gamma and R
#of gas mixtures. I probably included that file with this program.)

#primary driver conditions, sorted into a dictionary with gamma, R, and then mach number at throat

primary_driver = dict([('100', [Gas({'He':1.0},inputUnits='moles'),2.15]),
                       ('80',[Gas({'He':0.8,'Ar':0.2},inputUnits='moles'),1]),
                        ('90',[Gas({'He':0.9,'Ar':0.1},inputUnits='moles'),1.59])])
  
def main():

    #CURRENTLY RUNNING WITHOUT ANY OF THE REFLECTED SHOCK TYPE STUFF 
    
    secondary = 'n'
    
    #decided to cut out all the automation for this rebuilt so I can work with
    #easier for now
    
    #I've managed to cut down the storage of data a lot now, I just have a
    #dictionary for the gas state, u and M in each section
    
    states = {} #states dictionary that we'll fill up later
    V = {}
    M = {}
    
    #a couple of switches
    
    ref = 0 #ref = 1 allows shock reflection at D2
    stand = 1 #if stand = 1, calculate pressure for standing shock in region 2
    REV = 0
    
    #initial condition I want to test will be Umar's condition
    #(100%He driver, no secondary driver)
    
    #state 4, piston burst pressure
    states['s4']=primary_driver['100'][0]
    states['s4'].set_pT(2.79e7,2700.0)
    V['s4']=0.0
    M['s4']=0.0
    
    M['s11']=primary_driver['100'][1]

    #state11, piston after steady expansion at throat of driver
    
    (states['s11'], V['s11']) = expand_from_stagnation(1.0/(p0_p(M['s11'],states['s4'].gam)),states['s4'])
    
    #build the gas objects for all the other sections based on knowing what is in each
    
    T0 = 300.0 #atmospheric temperature (K), used for any section starting at ambient
    p0 = 101300.0 #atmospheric p (Pa)
    
    #start with shock tube and acc tube, start with atmospheric p and T
    
    states['s1'] = Gas({'Air':1.0,})
    states['s1'].set_pT(p0,T0)
    V['s1']=0.0
    M['s1']=0.0
    
    states['s5'] = Gas({'Air':1.0,})
    states['s5'].set_pT(p0,T0)
    V['s5']=0.0
    M['s5']=0.0
    
    #now let's do some clone for the states derived from these
    
    states['s3'] = states['s11'].clone() #3 will be 11 after unsteady expansion
    
    states['s2'] = states['s1'].clone() #2 is 1 shock heated
    
    states['s6'] = states['s5'].clone() #6 is 5 shock heated
    
    states['s7'] = states['s2'].clone() #7 is 2 after unsteady expansion
    
    #states['s15'] = states['s7'].clone() #15 is 7 after steady expansion through nozzle
    
    #states['s16'] = states['s15'].clone() #16 is 15 after being processed by normal shock
    
    print 'Enter shock speeds below to find solution.'
    
    vs1 = raw_input('Vs1? '); vs2 = raw_input('Vs2? '); #vs3 = raw_input('Vs3? ')
    
    Vs1 = int(vs1);Vs2 = int(vs2);#Us3 = int(us3) #shock speeds (m/s) (mess around with these three speeds to iteratively find a suitable condition. a couple of short tips: you don't want your acceleration tube fill pressure much lower than 10Pa, takes too long to fill. You probably never want more than an atmosphere in your secondary driver or shock tube either.)
      
    arbitrary = True #dummy variable
    
    while arbitrary:
        
        #stole some code from PJ and am modifying it for my needs
        
        print "\nNow do unsteady expansion..."
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
            
            print "current guess for p1 = {0}".format(p1)            
            
            state1.set_pT(p1,300.0) #set s1 at set pressure and ambient temp
            
            (V2, V2g) = normal_shock(state1, Vs1, state2)
            
            #Across the contact surface, p3 == p2
            p3 = state2.p
            
            # Across the expansion, we get a velocity, V5g.
            V3g, state3 = finite_wave_dp('cplus', V11g, state11, p3)

            return (V2g - V3g)/V2g        
        
        #get p1 for our chosen shock speed
        p1 = secant(error_in_velocity_s11_driver_expansion, 1000.0, 100000.0, tol=1.0e-3)
        print "From secant solve: p1 = {0}".format(p1)
        
        #start using p1 now, compute states 1,2 and 3 using the correct p1
        states['s1'].set_pT(p1,T0)
        (V2, V['s2']) = normal_shock(states['s1'], Vs1, states['s2'])
        V['s3'], states['s3'] = finite_wave_dv('cplus', V['s11'], states['s11'], V['s2'])
        
        #states['s1'].write_state(sys.stdout)            
        #states['s2'].write_state(sys.stdout)    
        #states['s3'].write_state(sys.stdout)    
        
        #get mach numbers for the output
        Ms1 = Vs1/states['s1'].son
        M['s2'] = V['s2']/states['s2'].son
        M['s3']= V['s3']/states['s3'].son

        def error_in_velocity_s2_expansion(p5, state2=states['s2'], 
                                                   V2g=V['s2'], state5=states['s5'],
                                                    state6=states['s6'],Vs2=Vs2):
            """Compute the velocity mismatch for a given pressure ratio across the 
            unsteady expansion from state 2 to state 7."""
            
            print "current guess for p5 = {0}".format(p5)
            
            state5.set_pT(p5,300.0) #set s1 at set pressure and ambient temp
            
            (V6, V6g) = normal_shock(state5, Vs2, state6)
            
            #Across the contact surface, p3 == p2
            p6 = state6.p
            
            # Across the expansion, we get a velocity, V7g.
            V7g, state7 = finite_wave_dp('cplus', V2g, state2, p6)

            return (V6g - V7g)/V6g        
        
        #get p1 for our chosen shock speed
        p5 = secant(error_in_velocity_s2_expansion, 10.0, 100.0, tol=1.0e-3)
        print "From secant solve: p5 = {0}".format(p5)
        
        #start using p5 now, compute states 5,6 and 7 using the correct p5
        states['s5'].set_pT(p5,T0)
        (V6, V['s6']) = normal_shock(states['s5'], Vs2, states['s6'])
        V['s7'], states['s7'] = finite_wave_dv('cplus', V['s2'], states['s2'], V['s6'])
        
        #get mach numbers for the output
        Ms2 = Vs2/states['s2'].son
        M['s6'] = V['s6']/states['s6'].son
        M['s7']= V['s7']/states['s7'].son
        
        #do the nozzle calc up to s15 now
        
        (V['s15'], states['s15']) = steady_flow_with_area_change(states['s7'], V['s7'], 2.5)
        M['s15']= V['s15']/states['s15'].son
               
                     
        #----------------- output -------------------

        output = open('x2run.txt',"w")  #output file creation
        
        #added a line that tells the user easily what the important lines of result are
        
        if secondary == 'y':
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
            important_lines1 = 'Line 1 is shock tube fill. Line 4 is reservoir.'
            print important_lines1
            output.write(important_lines1 + '\n')   
            
            important_lines2 = 'Line 5 is acceleration tube fill. Line 7 is test gas entering nozzle.'
            print important_lines2
            output.write(important_lines2 + '\n')
            
            important_lines3 = 'Line 15 is test gas exiting nozzle (using area ratio of 2.5).'
            print important_lines3
            output.write(important_lines3 + '\n')
        
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
        
        for i in range(1,16): #will cover 1-16
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
                
                conditions = "{0:<7}{1:<10.7g}{2:<9.0f}{3:<9.0f}{4:<9.0f}{5:<9.2f}{6:<9.4f}{7:<9.0f}{8:<9.1f}"\
                .format(i, states[it_string].p, states[it_string].T,
                        states[it_string].son,V[it_string],M[it_string],
                        states[it_string].rho, 0,0)
                        
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
        
        change_tuple = ('VS1','VS2','VS3', 'quit')
        
        change = is_valid(raw_input('What do you want to change? {0}'.format(change_tuple)), change_tuple)
            
        if change == 'quit':
            arbitrary == False
            break
            #quit()
                    
        elif change == 'VS1':
            difference = int(raw_input('How much do you want to change it by? '))
            Vs1 += difference
            print 'VS1 = {0}m/s'.format(Vs1)
            
        elif change == 'VS2':    
            difference = int(raw_input('How much do you want to change it by? '))
            Us2 += difference
            print 'VS2 = {0}m/s'.format(Vs2)
        
        else:
            difference = int(raw_input('How much do you want to change it by? '))
            Us3 += difference
            print 'VS3 = {0}m/s'.format(Vs3)
        
                   
main()

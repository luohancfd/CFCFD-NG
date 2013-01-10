"""
Pitot classic: cleanup and conversion of the original GWBasic code by Richard Morgan

Not alot to say here, I decided to try and clean up the original version of Pitot
written by Richard so it could still exist as an archival piece, but in a nicer
version for anyone that wants to mess with it. Mainly done for Richard's benefit/
some tradition. I've tried to use functions from the code library where possible, 
and I have updated the tunnel syntax so it matches the new version of pitot, and the
syntax that was agreed on by Richard, PJ and I last a couple of months ago when I
was developing the last version of Pitot.

29-Nov-12: I decided to add the code from the new pitot that does the normal shock
calcuation at the end of the tunnel, on the stagnation streamline over the test piece,
to this code as well, as it's very handy for my own testing. I also wanted to add a 
version string to the results printout so people know what version of the code they
used to create certain results.

Chris James (c.james4@uq.edu.au) 29 Nov 2012

"""

VERSION_STRING = "29-Nov-2012"

#start with some useful imports
import os, sys, math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from cfpylib.gasdyn.ideal_gas_flow import *
from cfpylib.gasdyn.cea2_gas import Gas, make_gas_from_name
from cfpylib.gasdyn.gas_flow import shock_ideal, normal_shock

#define our functions

def press (M,G):
    """function to find static pressure ratio across the shock
    inputs are mach number (M) and gamma (G)"""
    
    Pratio = (2*G*M**2-(G-1))/(G+1)
    return Pratio
    
def dens (M,G):
    """function to find density ratio across the shock
    inputs are mach number (M) and gamma (G)"""
    
    rhoratio = (G+1)*M**2/((G-1)*M**2+2)
    return rhoratio

def isen(a,G):
    """function to find isentropic pressure ratio from sound speed
    inputs are sound speed ratio (a) and gamma (G)"""
    
    IPratio = a**(2*G/(G-1))
    return IPratio
    
def ad(G,au,Uu,Ud):
    """function to find the downstream sound speed across U-a wave
    inputs are gamma (G), upstream sound speed (au), upstream velocity (Uu), and downstream velocity (Ud)"""
    
    ad = au+(G-1.0)*(Uu-Ud)/2.0
    return ad
    
def IP(M,G):
    """function to find the static to total pressure ratio
    inputs are mach number (M), and gamma (G)"""
    
    ratio = (1.0+(G-1.0)*M**2/2.0)**(G/(1.0-G))
    return ratio

def ar (P,G):
    """function to find the sound speed ratio from the isentropic pressure ratio
    inputs are pressure ratio (p), and gamma (g)"""
    
    ratio = P**((G-1)/2*G)
    return ratio

def pitot(M,G):
    """function to calculate pitot pressure RATIO using the rayleigh pitot formula
    inputs are mach number (M) and gamma (G)"""
    
    rayleigh = M**2.0*((G+1.0)/2.0)**((G+1.0)/(G-1.0))*(G-(G-1.0)/2.0*M**-2.0)**(1.0/(1.0-G))
    return rayleigh
       
def nozzleMfinder (gamma, Aratio, Min, Moutguess):
    """Function that will find M exiting a supersonic nozzle when gamma and the area
    ratio (A/A*) are known at the exit of the nozzle.
    
    Uses the bisection method, does a max of 100 iterations.
    
    An initial guess for M is needed. (I'd recommend Min +/- 5 depending on
    whether you expect M to go down or up).
    
    Output is the iteratively solved Mach number for that area ratio.
    
    Reference: A/A* equation from NACA report 1135
    
    Mfinder (float, float, float, float) -> float
    
    Chris James (c.james4@uq.edu.au) 28/08/12"""
    
    #lambda function below to get the Aratio (A/A*)
    
    AonAstar = lambda gamma, M: 1.0/((((gamma+1.0)/2.0)**((gamma+1.0)/(2.0*(gamma-1.0))))\
    *(M*(1.0+((gamma-1.0)/2.0)*M**2.0)**-((gamma+1.0)/(2.0*(gamma-1.0))))) #A/A*
    
    
    AonAstarin = AonAstar(gamma,Min) #get value coming in
    AonAstarout = AonAstarin*Aratio #get value coming out based on area ratio specified
    
    #print AonAstarin, AonAstarout    
    
    f = lambda gamma, M: AonAstar(gamma,M) - AonAstarout    
    
    left = Min #left guess is the input Mach number to the nozzle
    right = Moutguess #right guess is the Mach number out guess
    tol = 0.00000001 #tolerance
    counter = 0 #iteration counter
    
    #print f(gamma,left),f(gamma,right)
    
    if f(gamma,left)*f(gamma,right) < 0.0: #check the root is in the interval
        
        while True:
            counter += 1
            mid = 0.5*(left+right)
            #print mid
            if abs(left-right)<tol or counter > 10000:
                break
            if f(gamma,mid)*f(gamma,left)<0.0:
                right = mid
            else:
                left = mid
        
        #print 'the Mach number is {0}, A/A* is {1}'.format(mid,AonAstar(gamma,mid))
        return mid #the mach number for that area ratio
    else:
        #print 'no solution in the interval'
        return

def nozzleTratio(Min,Mout,gamma):
    """Function to find the ratio of static temperature (T) through a nozzle when
    the Mach number in (Min) and the Mach number out of the nozzle (Mout) and 
    gamma are known. Uses the isentropic stagnation temperature ratio in and out
    of the nozzle to get the temperature ratio across the whole nozzle, knowing
    that stagnation temperature is constant through the nozzle.
    
    returns Tout/Tin
    
    reference: NACA report 1135
    
    nozzleTratio (float, float, float) -> float
    
    Chris James (c.james4@uq.edu.au)"""

    stagnationTratioin = (1.0+((gamma-1.0)/2.0)*Min**2.0)**(-1.0) #Tin/Ttin
    
    stagnationTratioout = (1.0+((gamma-1.0)/2.0)*Mout**2.0)**(-1.0) #Tout/Ttout
    
    Tratio = stagnationTratioout/stagnationTratioin

    return Tratio

def nozzlepratio(Min,Mout,gamma):
    """Function to find the ratio of static pressure (p) through a nozzle when
    the Mach number in (Min) and the Mach number out of the nozzle (Mout) and 
    gamma are known. Uses the isentropic stagnation pressure ratio in and out
    of the nozzle to get the pressure ratio across the whole nozzle, assuming
    that stagnation pressure is constant through the nozzle.
    
    returns pout/pin
    
    reference: NACA report 1135
    
    nozzlepratio (float, float,float) -> float
    
    Chris James (c.james4@uq.edu.au)"""

    stagnationpratioin = (1.0+((gamma-1.0)/2.0)*Min**2.0)**(-gamma/(gamma-1.0)) #pin/ptin
    
    stagnationpratioout = (1.0+((gamma-1.0)/2.0)*Mout**2.0)**(-gamma/(gamma-1.0)) #pout/ptout
    
    pratio = stagnationpratioout/stagnationpratioin

    return pratio

def nozzlerhoratio(Min,Mout,gamma):   
    """Function to find the ratio of static density (rho) through a nozzle when
    the Mach number in (Min) and the Mach number out of the nozzle (Mout) and 
    gamma are known. Uses the isentropic stagnation density ratio in and out
    of the nozzle to get the density ratio across the whole nozzle, assuming
    that stagnation density is constant through the nozzle.
    
    returns rhoout/rhoin
    
    reference: NACA report 1135
    
    nozzlerhoratio (float, float,float) -> float
    
    Chris James (c.james4@uq.edu.au)"""

    stagnationrhoratioin = (1.0+((gamma-1.0)/2.0)*Min**2.0)**(-1.0/(gamma-1.0)) #rhoin/rhotin
    
    stagnationrhoratioout = (1.0+((gamma-1.0)/2.0)*Mout**2.0)**(-1.0/(gamma-1.0)) #rhoout/rhotout
    
    rhoratio = stagnationrhoratioout/stagnationrhoratioin

    return rhoratio
    
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

#piston configuration conditions, sorted into a dicionary with reservoir temp and then pressure
#at the burst of the primary steel diaphragm

pistons = dict(x2_lwp = [2700.0, 2.79e7])

#primary driver conditions, sorted into a dictionary with gamma, R, and then mach number at throat

primary_driver = dict([('100', [1.667,2077,2.15]),('80',[1.667,742.9,1.0]),('90',[1.67,1095,1.59])])

#different test gases, sorted into a dictionary with gamma and then R
#the third value is now a proper gas object for the normal shock calc

gases = {}
gases['He'] = [1.667,2077.0, Gas(reactants={'He':1.0}, inputUnits='moles', 
                with_ions=True, outputUnits='moles')]
gases['air'] = [1.4,286.0, Gas(reactants={'air':1.0}, inputUnits='moles', 
                with_ions=True, outputUnits='moles')]
gases['90H210He'] = [1.42,3754.2, Gas(reactants={'H2':0.9, 'He':0.1}, 
                inputUnits='moles', with_ions=True, outputUnits='moles')]
gases['90H210Ne'] = [1.42,2169.3, Gas(reactants={'H2':0.9, 'Ne':0.1}, 
                inputUnits='moles', with_ions=True, outputUnits='moles')]
gases['85H215He'] = [1.443,3593.0, Gas(reactants={'H2':0.85, 'He':0.15}, 
                inputUnits='moles', with_ions=True, outputUnits='moles')]
gases['85H215Ne'] = [1.443,1753.7, Gas(reactants={'H2':0.85, 'Ne':0.15}, 
                inputUnits='moles', with_ions=True, outputUnits='moles')]
gases['15H285Ne'] = [1.6277,476.25, Gas(reactants={'H2':0.15, 'Ne':0.85}, 
                inputUnits='moles', with_ions=True, outputUnits='moles')]
gases['60H240Ne'] = [1.5098,895.68, Gas(reactants={'H2':0.6, 'Ne':0.4}, 
                inputUnits='moles', with_ions=True, outputUnits='moles')]
gases['venus'] = [1.2929,191.36, Gas(reactants={'CO2':0.965, 'N2':0.035}, 
                inputUnits='moles', with_ions=True, outputUnits='moles')]
gases['H2'] = [1.405,4124.0, Gas(reactants={'H2':1.0}, inputUnits='moles', 
                with_ions=True, outputUnits='moles')]
gases['N2'] = [1.40,296.8, Gas(reactants={'N2':1.0}, inputUnits='moles', 
                with_ions=True, outputUnits='moles')]
gases['mars'] = [1.293,191.71, Gas(reactants={'CO2':0.96, 'N2':0.04}, 
                inputUnits='moles', with_ions=True, outputUnits='moles')]
gases['titan'] = [1.395,303.29, Gas(reactants={'N2':0.95, 'CH4':0.05}, 
                inputUnits='moles', with_ions=True, outputUnits='moles')]

  
def main():
    """top level function!"""
    
    #dictionaries to store useful stuff
    
    # pressure, mach number, temperature, sound speed, velocity, density, R, gamma
    
    P={};M={};T={};a={};V={};rho={};RU2={};G={}
    
    #a couple of switches
    
    ref = 0 #ref = 1 allows shock reflection at D2
    stand = 1 #if stand = 1, calculate pressure for standing shock in region 2
    REV = 0
    
    #inputs (asks user for some of these inputs)
         
    print 'Driver condition'
    
    piston = is_valid(raw_input('Piston configuration? {0}'.format(pistons.keys())),pistons.keys())
    
    primary = is_valid(str(raw_input('Percentage of Helium in the primary driver? {0} '.format(primary_driver.keys()))),primary_driver.keys())
    
    GP = primary_driver[primary][0]; RP = primary_driver[primary][1]
       
    secondary_input = is_valid(str(raw_input('Is the secondary driver being used (y/n)? ')),['y','n'])
    if secondary_input == 'y': 
        secondary = True   
    else:
        secondary = None
       
    test_gas = is_valid(raw_input('Test gas? {0} '.format(gases.keys())), gases.keys())
    
    if secondary:
        
        GS = gases['He'][0]; RS = gases['He'][1]
        
    GT = gases[test_gas][0]; RT = gases[test_gas][1]
        
    #air as accelorator gas
    GA = gases['air'][0]; RA = gases['air'][1]
    
    filename = raw_input('filename? ')
    
    if filename == '':
        filename = 'x2run.txt'
               
    print 'Enter shock speeds (in m/s) below to find solution.'
    
    if secondary: Vsd = int(raw_input('Vsd? '))
    
    Vs1 = int(raw_input('Vs1? '))
    Vs2 = int(raw_input('Vs2? ')) 
    
    T0 = 300.0;T['s4'] = pistons[piston][0]; P['s4'] = pistons[piston][1] #ambient + reservoir temp (K), reservoir pressure (Pa) (this is the current lightweight piston driver condition used for all shots in X2)
    M['s3s'] = primary_driver[str(primary)][2] #Mach numbers terminating steady expansions (used for orifice plate) set M[11] to 0 for constant area shock tunnel (not normally used), M[11] =1 is no orifice plate, sonic throat. normal configuration, theoretical optimum. set M[11] to 2.15 to use pure He driver orifice plate 
    
    #put gammas into dictionaries
    if secondary:
        G['sd1'] = GS; G['sd2'] = GS; G['sd3'] = GP; G['s3'] = GS
    else:
        G['s3'] = GP
    G['s4'] = GP; G['s3s'] = GP; G['s1'] = GT; G['s2'] = GT
    G['s7'] = GT; G['s5'] = GA; G['s6'] = GA; G['s13'] = GA; G['s8'] = GT
    
    #add ambient temperature to arrays where required
    T['sd1' ]= T0; T['s1'] = T0; T['s5'] = T0
    
    M['s4'] = 0; V['s4']= 0; rho['s4'] = 0
    
    arbitrary = True #dummy variable
    
    while arbitrary:
        
        #output file creation

        output = open(filename,"w")
    
        #----------------- starting off ----------------------
        
        #calculate the sound speeds known
        if secondary:
            a['sd1'] = math.sqrt(T['sd1']*GS*RS)
        a['s4'] = math.sqrt(T['s4']*GP*RP) 
        a['s1'] = math.sqrt(T['s1']*GT*RT)
        a['s5'] = math.sqrt(T['s5']*GA*RA)
        a['s3s'] = a['s4']*IP(M['s3s'],GP)**((GP-1.0)/2.0/GP) 
        P['s3s'] = P['s4']*IP(M['s3s'],GP)
        V['s3s'] = a['s3s']*M['s3s'] 
        T['s3s'] = T['s4']*(a['s3s']/a['s4'])**2.0
        rho['s3s'] = P['s3s']/RP/T['s3s']
        RU2['s3s'] = rho['s3s']*V['s3s']**2.0
        
        #---------------------- secondary driver ------------------------------
        
        if secondary: #if secondary driver is in use
            #conditions behind shock         
            Msd = Vsd/a['sd1']
            V['sd2'] = Vsd*(1.0-u2_u1(Msd,GS))
            T['sd2'] = T['sd1']*T2_T1(Msd,GS)
            a['sd2'] = math.sqrt(T['sd2']*GS*RS)
            M['sd2'] = V['sd2']/a['sd2']
            
            #unsteady expansion using Vsd2 = Vsd3
            V['sd3'] = V['sd2']
            
            a['sd3'] = ad(GP, a['s3s'], V['s3s'], V['sd3'])
            P['sd3'] = P['s3s']*isen((a['sd3']/a['s3s']),GP) 
            T['sd3'] = T['s3s']*(a['sd3']/a['s3s'])**2.0
            M['sd3'] = V['sd3']/a['sd3'] 
            rho['sd3'] = P['sd3']/RP/T['sd3'] 
            RU2['sd3']=rho['sd3']*V['sd3']**2.0
            
            #can now get conditions at state sd2 now we have Psd3
            P['sd2'] = P['sd3']
            rho['sd2'] = P['sd2']/RS/T['sd2']
            RU2['sd2'] = rho['sd2']*V['sd2']**2.0
            #and then finally the fill pressure...
            P['sd1'] = P['sd2']/p2_p1(Msd,GS)
            M['sd1']=0; V['sd1']=0; rho['sd1'] = P['sd1']/RS/T['sd1']; RU2['sd1'] = 0

        #-------------------- shock tube ------------------------------
        
        #conditions behind shock         
        Ms1 = Vs1/a['s1']
        V['s2'] = Vs1*(1.0 - u2_u1(Ms1,GT))
        T['s2'] = T['s1']*T2_T1(Ms1,GT)
        a['s2'] = math.sqrt(T['s2']*GT*RT)
        M['s2'] = V['s2']/a['s2']
        
        #unsteady expansion using Vs2 = Vs3
        V['s3'] = V['s2']
        
        if secondary: #expansion of secondary driver gas into shock tube
            a['s3'] = ad(GS, a['sd2'], V['sd2'], V['s3'])
            P['s3'] = P['sd2']*isen((a['s3']/a['sd2']),GS) 
            T['s3'] = T['sd2']*(a['s3']/a['sd2'])**2.0
            rho['s3'] = P['s3']/RS/T['s3'] 
        else: #expansion of primary driver gas into shock tube
            a['s3'] = ad(GP, a['s3s'], V['s3s'], V['s3'])
            P['s3'] = P['s3s']*isen((a['s3']/a['s3s']),GP) 
            T['s3'] = T['s3s']*(a['s3']/a['s3s'])**2.0
            rho['s3'] = P['s3']/RP/T['s3'] 
        M['s3'] = V['s3']/a['s3'] 
        RU2['s3']=rho['s3']*V['s3']**2.0
        
        #can now get conditions at state s2 now we have Ps3
        P['s2'] = P['s3']
        rho['s2'] = P['s2']/RT/T['s2']
        RU2['s2'] = rho['s2']*V['s2']**2.0
        #and then finally the fill pressure...
        P['s1'] = P['s2']/p2_p1(Ms1,GT)
        M['s1']=0; V['s1']=0; rho['s1'] = P['s1']/RT/T['s1']; RU2['s1'] = 0
        
        #------------- acceleration tube -------------------------------
        
        #conditions behind shock         
        Ms2 = Vs2/a['s5']
        V['s6'] = Vs2*(1.0 - u2_u1(Ms2,GA))
        T['s6'] = T['s5']*T2_T1(Ms2,GA)
        a['s6'] = math.sqrt(T['s6']*GA*RA)
        M['s6'] = V['s6']/a['s6']
        
        #unsteady expansion using Vs6 = Vs7
        V['s7'] = V['s6']
        
        a['s7'] = ad(GT, a['s2'], V['s2'], V['s7'])
        P['s7'] = P['s2']*isen((a['s7']/a['s2']),GT) 
        T['s7'] = T['s2']*(a['s7']/a['s2'])**2.0
        rho['s7'] = P['s7']/RA/T['s7'] 
        M['s7'] = V['s7']/a['s7'] 
        RU2['s7']=rho['s7']*V['s7']**2.0
        
        #can now get conditions at state s6 now we have Ps7
        P['s6'] = P['s7']
        rho['s6'] = P['s6']/RA/T['s6']
        RU2['s6'] = rho['s6']*V['s6']**2.0
        #and then finally the fill pressure...
        P['s5'] = P['s6']/p2_p1(Ms2,GA)
        M['s5']=0; V['s5']=0; rho['s5'] = P['s5']/RA/T['s5']; RU2['s5'] = 0
               
        #--------------- nozzle calculations--------------------
        
        #currently using area ratio of 2.5!
        
        M['s8'] = nozzleMfinder (GT, 2.5  , M['s7'], M['s7']+5.0)
        
        T['s8'] = T['s7']*nozzleTratio(M['s7'],M['s8'],GT)
        P['s8'] = P['s7']*nozzlepratio(M['s7'],M['s8'],GT)
        V['s8'] = V['s7']*(M['s8']/M['s7'])*(T['s8']/T['s7'])**(1.0/2.0)
        rho['s8'] = rho['s7']*nozzlerhoratio(M['s7'],M['s8'],GT)
        a['s8'] = V['s8']/M['s8']
        
        #------------ normal shock calculation over the model -----------
        
        #NOTE: this part of the code actually does use chemistry for the equilibrium part...
        #it also uses the states dictionary syntax from pitot, so I can use the pitot printing function        
        
        #build states dictionary:
        states = {}
        
        #first build a gas object at the nozzle conditions
        
        states['s8'] = gases[test_gas][2].clone()
        states['s8'].set_pT(P['s8'],T['s8'])
        
        #and then the two different state 10 conditions that we will then shock
        
        states['s10f'] = states['s8'].clone()
        states['s10e'] = states['s8'].clone()
        
        #do a frozen and then an equilibrium normal shock on these conditions:
        
        V10f, V10fg = shock_ideal(states['s8'],V['s8'],states['s10f'])
        V['s10f'] = V10fg; M['s10f'] = V['s10f']/states['s10f'].son
        
        V10e, V10eg = normal_shock(states['s8'],V['s8'],states['s10e'])
        V['s10e'] = V10eg; M['s10e'] = V['s10e']/states['s10e'].son
                     
        #------------------------ output -------------------------------------
        
        #pillaged from new pitot!
        
        print " "
        
        version_printout = "Pitot Classic Version: {0}".format(VERSION_STRING)
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
            
        description_3 = 'state 8 is test gas exiting the nozzle (using area ratio of {0}).'.format(2.5)
        print description_3
        output.write(description_3 + '\n')
        
        description_4 = 'state 10f is frozen shocked test gas flowing over the model.'
        print description_4
        output.write(description_4 + '\n')  
        
        description_5 = 'state 10e is equilibrium shocked test gas flowing over the model.'
        print description_5
        output.write(description_5 + '\n')     
        
        test_gas_used = 'Test gas is {0} (gamma = {1}, R = {2}, {3}).'.format(test_gas,GT,RT,gases[test_gas][2].reactants)
        print test_gas_used
        output.write(test_gas_used + '\n')  
        
        driver_gas_used = 'Driver gas has {0}%He (gamma = {1}, R = {2})'.format(primary,GP,RP)        
        print driver_gas_used
        output.write(driver_gas_used + '\n')  

        if secondary:
            secondary_shockspeeds = "Vsd = {0:.2f} m/s, Msd = {1:.2f}".format(Vsd,Msd)
            print secondary_shockspeeds
            output.write(secondary_shockspeeds + '\n')
        
        shockspeeds = "Vs1= {0:.2f} m/s, Ms1= {1:.2f} ,Vs2= {2:.2f} m/s, Ms2= {3:.2f}".format(Vs1,Ms1,Vs2,Ms2) 
        print shockspeeds #prints above line in console
        output.write(shockspeeds + '\n') #writes above line to output file (input to write command must be a string)
        
        key = "{0:6}{1:11}{2:9}{3:6}{4:9}{5:6}{6:9}{7:8}{8:9}".format("state","P","T","a","V","M","rho","pitot","stgn")
        print key
        output.write(key + '\n')
        
        units = "{0:6}{1:11}{2:9}{3:6}{4:9}{5:6}{6:9}{7:9}{8:9}".format("","Pa","K","m/s","m/s","","m^3/kg","kPa","MPa")
        print units
        output.write(units + '\n')
        
        #new dictionaries here to add pitot and stagnation pressure calcs
        
        pitot_pressure = {} #pitot pressure dict
        p0 = {} #stagnation pressure dict
        
        def condition_printer(it_string):
            """Prints the values of a specified condition to the screen and to 
            the output file. 
            
            I made a function of this so I didn't have to keep pasting the code in."""
            
            if M.has_key(it_string):
                    
                if M[it_string] == 0:
                    pitot_pressure[it_string] = 0
                    p0[it_string] = 0
                else:
                    pitot_pressure[it_string] = P[it_string]*pitot(M[it_string],G[it_string])/1000.0
                    p0[it_string] = p0_p(M[it_string], G[it_string])*P[it_string]/1.0e6
                
                conditions = "{0:<6}{1:<11.7g}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.4f}{7:<8.0f}{8:<9.1f}"\
                .format(it_string, P[it_string], T[it_string], a[it_string],
                        V[it_string],M[it_string],rho[it_string], 
                        pitot_pressure[it_string], p0[it_string])
                        
                print conditions
                output.write(conditions + '\n')
                
        def condition_printer_pitot(it_string):
            """The states based condition printer from pitot, used for the
            normal shock stuff at the end of the tunnel!"""
            
            if states.has_key(it_string):
                    
                if M[it_string] == 0:
                    pitot_pressure[it_string] = 0
                    p0[it_string] = 0
                else:
                    pitot_pressure[it_string] = pitot_p(states[it_string].p,M[it_string],states[it_string].gam)/1000.0
                    p0[it_string] = p0_p(M[it_string], states[it_string].gam)*states[it_string].p/1.0e6
                
                conditions = "{0:<6}{1:<11.7}{2:<9.1f}{3:<6.0f}{4:<9.1f}{5:<6.2f}{6:<9.4f}{7:<8.0f}{8:<9.1f}"\
                .format(it_string, states[it_string].p, states[it_string].T,
                        states[it_string].son,V[it_string],M[it_string],
                        states[it_string].rho, pitot_pressure[it_string], p0[it_string])
                        
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
        condition_printer_pitot('s10f')
        condition_printer_pitot('s10e')
            
                
        Cp_test_gas = (gases[test_gas][0]*gases[test_gas][1])/(gases[test_gas][0]-1.0)
        
        stagnation_enthalpy = (((0.5)*V['s7']**2.0)+(Cp_test_gas*T['s7'])) #MJ/kg

        stag_enth = 'The stagnation enthalpy (Ht) at the end of the acceleration tube {0:<.5g} MJ/kg.'.format(stagnation_enthalpy/1.0e6)
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
    
        output.close()
               
        #------------------------------ user commands---------------------------
        
        if secondary:
            change_tuple = ('Vsd','Vs1','Vs2', 'quit')
        else:
            change_tuple = ('Vs1','Vs2', 'quit')
        
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
                    
        elif change == 'Vsd':
            difference = int(raw_input('How much do you want to change it by? '))
            Vsd += difference
            print 'Vsd = {0}m/s'.format(Vsd)
            
        elif change == 'Vs1':    
            difference = int(raw_input('How much do you want to change it by? '))
            Vs1 += difference
            print 'Vs1 = {0}m/s'.format(Vs1)
        
        else:
            difference = int(raw_input('How much do you want to change it by? '))
            Vs2 += difference
            print 'Vs2 = {0}m/s'.format(Vs2)

        
                   
main()

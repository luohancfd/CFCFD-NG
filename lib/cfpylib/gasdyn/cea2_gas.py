"""
cea2_gas.py: Thermodynamic properties of a gas mixture in chemical equilibrium.

It interfaces to the CEA code by writing a small input file,
running the CEA code as a child process and then reading the results
from the CEA plot file.

See the report::

    Bonnie J. McBride and Sanford Gordon
    Computer Program for Calculation of Complex Chemical Equilibrium
    Compositions and Applications II. Users Manual and Program
    Description. NASA Reference Publication 1311, June 1996.

for details of the input and output file formats.

.. Author: 
   PA Jacobs
   Institute of Aerodynamics and Flow Technology
   The German Aerospace Center, Goettingen.

.. Versions:
   24-Dec-02: First code.
   10-May-04: Updated for a mix of species.
   06-Feb-05: renamed to cea_gas.py
   28-Feb-08: Added a get_eq_massf() access function.
   28-Fed-08: Major changes to allow proper calculation at high temps.
   11-Dec-08: Addition of basic incident Shock function
   RJG and DFP have made these later changes which are significant.
"""

import sys, string, math, os, subprocess, re
import copy

# -------------------------------------------------------------------
# First, global data.

DEBUG_GAS = 0
R_universal = 8314.0;  # J/kgmole.K
CEA_COMMAND_NAME = 'cea2'

# -------------------------------------------------------------------
# Second, utility functions.

def locate_executable_file(name):
    """
    Locates an executable file, if available somewhere on the PATH.

    :param name: may be a simple file name or fully-qualified path.
    :returns: the full program name, if is is found and is executable,
        else None.
    """
    def is_exe(path):
        return os.path.exists(path) and os.access(path, os.X_OK)

    head, tail = os.path.split(name)
    if head:
        # If there is a head component, we may have been given
        # full path to the exe_file.
        if is_exe(name): return name
    else:
        # We've been given the name of the program
        # without the fully-qualified path in front,
        # now search the PATH for the program.
        for path in os.environ["PATH"].split(os.pathsep):
            fullName = os.path.join(path, name)
            if is_exe(fullName): return fullName
    return None

def run_cea_program(jobName):
    """
    Runs the CEA program on the specified job.
    """
    inpFile = jobName + '.inp'
    if os.path.exists(inpFile):
        if DEBUG_GAS >= 2:
            print 'Start cea program on job %s...' % jobName
        p = subprocess.Popen(CEA_COMMAND_NAME, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate(jobName + '\n')
        return_code = p.wait()
        if DEBUG_GAS >= 2:
            print ('%s finished job %s.' % (CEA_COMMAND_NAME,jobName))
    else:
        raise Exception, 'The file %s is not present.' % inpFile

def get_cea2_float( value_str ):
    """
    Clean up the CEA2 short-hand notation for exponential format.
    """
    if value_str.find("-")>0:
    	value_str = value_str.replace("-","e-")
    if value_str.find("+")>0:
    	value_str = value_str.replace("+","e+")
    return float(value_str)
   
# ----------------------------------------------------------------

class Gas(object):
    """
    Provides the equation of state for the gas.
    """
    def __init__(self, gasName, speciesList=[], massfList=[], use_out_file=False):
        """
        Set up a new obects, from either a name of species list.
        """
	if locate_executable_file(CEA_COMMAND_NAME) is None:
            print "Could not find the executable program %s" % CEA_COMMAND_NAME
            print "The chemical equilibrium-analysis program is external"
            print "to the cfcfd3 code collection and needs to be obtained from NASA Glenn."
            print "Quitting the current program because we really can't do anything further."
            sys.exit()
	# ----------------------------------------------------------------
        self.gasName = gasName           # see method EOS for possible names
        self.species = speciesList       # species names as per CEA database
        self.nsp     = len(speciesList)  # number of species
        self.massf   = massfList         # user-specified mass fractions
        self.eq_massf = copy.copy(self.massf) # will overwrite this list later
        self.eq_molef = copy.copy(self.massf) # also to be written over
	# print self.massf
        # self.set_from_pAndT(1.0e5,1000.0)
	# check for electron to set with ions flag
	self.with_ions = False
	for sp in speciesList:
	    if sp=="e-": self.with_ions = True
        # self.with_ions = True # force it on for test PJ 18-Jan-2011
        self.p = 1.0e5 # need a value in case EOS gets called
        self.T = 300.0 # likewise
	self.EOS(speciesOut=1, thermoProps=1, transProps=1, 
                 problemType='pT', use_out_file=use_out_file)
        return
	
    def get_eq_massf(self, isp):
	"""
        Returns equil mass-fraction for species at index isp.
        """
	return self.eq_massf[isp]

    def get_eq_molef(self, isp):
	"""
        Returns equil mole-fraction for species at index isp.
        """
	return self.eq_molef[isp]  

    def set_from_pAndT(self, p, T, speciesOut=0, thermoProps=1,
                       transProps=1, use_out_file=False):
        """
        Fills out gas state from given p and T.
        """
        self.p = p;        # pressure, Pa
        self.T = T;        # temperature, K
        flag = self.EOS(speciesOut, thermoProps, transProps, use_out_file=use_out_file)
        return flag

    def write_state(self,strm):
        """
        Writes the gas state data to the specified stream.
        """
        strm.write('    p: %g Pa, T: %g K, rho: %g kg/m**3, e: %g J/kg, a: %g m/s\n'
                   % (self.p, self.T, self.rho, self.u, self.son) )
        strm.write('    R: %g J/(kg.K), gam: %g, Cp: %g J/(kg.K), mu: %g Pa.s, k: %g W/(m.K)\n'
                   % (self.R, self.gam, self.cp, self.mu, self.k) )
        for i in range(self.nsp):
            strm.write('species[%d] is %s: specified-massf %e, eq-massf %e, eq-molef %e\n' %
                       (i, self.species[i], self.massf[i], self.eq_massf[i], self.eq_molef[i]))

    def EOS(self, speciesOut=0, thermoProps=1, transProps=1, problemType='pT', use_out_file=False):
        """
        Computes the gas state, taking into account the high-temperature effects.

        It does this by writing a suitable input file for the CEA code,
        calling that code and then pulling the relevant results out of the CEA plt file.

        :param self: the gas state to be filled in
	:param speciesOut: a flag where (-1) do not plot species,
                           (0) requests massfs,
                           (1) requests molefs
	:param thermoProps: a flag where (0) does nothing,
                            (1) request Cv, Cp, entropy, a, gamma etc
        :param problemType: a string specifying input parameters: 'pT', 'rhoT', 'ps'
            
        :returns: None, but does update the contents of the gas state as a side-effect.
        """
        if DEBUG_GAS >= 2:
            print 'EOS: Write temporary input file.'
        inp_file_name = 'tmp.inp'
        fp = open(inp_file_name, 'w')
        fp.write('# %s generated by cea2_gas.py\n' % inp_file_name)
        if problemType == 'rhoT':
            if self.with_ions:
                fp.write('problem case=estcj tv ions\n')
            else:
                fp.write('problem case=estcj tv\n')
            fp.write('   rho,kg/m**3 %e\n' % self.rho)
            fp.write('   t(k)        %e\n' % self.T)
            if DEBUG_GAS >= 2:
                print 'EOS: input to CEA2 rho: %g, T: %g' % (self.rho, self.T)
        elif problemType == 'pT':
            if self.with_ions:
                fp.write('problem case=estcj tp ions\n')
            else:
                fp.write('problem case=estcj tp\n')
            fp.write('   p(bar)      %e\n' % (self.p / 1.0e5) )
            fp.write('   t(k)        %e\n' % self.T)
            if DEBUG_GAS >= 2:
                print 'EOS: input to CEA2 p: %g, T: %g' % (self.p, self.T)
        elif problemType == 'ps':
            if self.with_ions:
                fp.write('problem case=estcj ps ions\n')
            else:
                fp.write('problem case=estcj ps\n')
            fp.write('   p(bar)      %e\n' % (self.p / 1.0e5) )
            fp.write('   s/r         %e\n' % (self.s / R_universal) )
            if DEBUG_GAS >= 2:
                print 'EOS: input to CEA2 p: %g, s: %g' % (self.p, self.s)
        else:
            print 'Invalid problemType: %s' % problemType
            sys.exit( 0 );
        # Select the gas components.
        # Look for valid gas names here...
        fp.write('reac\n')
        if self.gasName == 'n2':
            fp.write('   name= N2  wtf= 1.\n')
            fp.write('   name= N   wtf= 0.\n')
            fp.write('only N2 N')
        elif self.gasName == 'co2':
            fp.write('   name= CO2  wtf= 1.\n')
        elif self.gasName == 'air':
            fp.write('   name= Air  wtf= 1.\n')
        elif self.gasName == 'h2ne':
            fp.write('   name= H2  moles= 0.15\n')
            fp.write('   name= Ne  moles= 0.85\n')
        elif self.gasName == 'air5species':
            # Assume air is 79% N2, 21% O2 by mole fraction
            fp.write('   name= N2  wtf= 0.767082\n')
            fp.write('   name= O2  wtf= 0.232918\n')
            fp.write('only N2 O2 N O NO')
        elif self.gasName == 'mix':
            for i in range(self.nsp):
                if self.massf[i] > 0.0:
                    fp.write('   name= %s  wtf= %g \n'%(self.species[i], self.massf[i]))
	    fp.write('only')
	    for i in range(self.nsp):
                fp.write(' %s' % self.species[i])
        fp.write('\n')
        if speciesOut == 0 :                    
            fp.write('output massf trans\n')    # gives us mass fractions (can use trace=1.e-10 here)
        else :
            fp.write('output trans\n')          # gives us mole fractions (can use trace=1.e-10 here)
        fp.write('plot')			       # only plotting mass fractions
        if speciesOut > -1:
            for i in range(self.nsp):
            	fp.write(' %s' % self.species[i])
        if thermoProps:
            fp.write(' h u s cp gam son rho p t')
        if transProps:
            fp.write(' vis cond')
        fp.write('\n')
        fp.write('end\n') 
        fp.close()
        #
        # Make sure that the database input files are in the working dir
        #
        if not os.path.exists('thermo.inp'):
            print 'Copying thermo.inp to the current working directory'
            os.system("cp %s/e3bin/thermo.inp ." % ( os.getenv("HOME") ) )
            print 'Copying trans.inp to the current working directory'
            os.system("cp %s/e3bin/trans.inp ." % ( os.getenv("HOME") ) )
        #
        # Make sure that binary versions of the database files exist.
        #
        if not os.path.exists('thermo.lib'):
            print 'Make the binary database for thermodynamic properties'
            run_cea_program('thermo')
            print 'Make the binary database for transport properties'
            run_cea_program('trans')
        #
        # Now, run the cea program in anger.
        run_cea_program('tmp')
        #
        # Pick out the interesting bits from the plotfile.
        # These data will be on the first non-comment line.
        if use_out_file==False:
            fp = open('tmp.plt', 'r')
            L = string.split(fp.readline())
            while L[0] == "#":
                # Skip over comment lines
                L = string.split(fp.readline())
            fp.close()
            #
            # Fill out equilibrium mass fractions...
	    # print "len(L) = %d, len(species) = %d" % ( len(L), len(self.species) )
            expected_entries = len(self.species)
            if thermoProps:
                expected_entries += 9
            if transProps:
                expected_entries += 2
            if len(L)!=expected_entries:
                print "cea .plt file returned %d values, while %d values were requested." % ( len(L), expected_entries )
                print "This probably means too many values were requested."
                print "Try calling the EOS function with 'use_out_file=True'."
                sys.exit()
            if speciesOut == 0 : 
                for i in range(self.nsp):
                    self.eq_massf[i] = float(L[i])
                next = self.nsp
            elif speciesOut == 1:
                for i in range(self.nsp):
                    self.eq_molef[i] = float(L[i])
                next = self.nsp
            else:
                next = 0
            #
            # Fill out thermo properties if requested
            if thermoProps:
                self.h = float(L[self.nsp]) * 1.0e3;    # enthalpy, J/kg
                self.u = float(L[self.nsp+1]) * 1.0e3;  # internal energy, J/kg
                self.e = self.u
                self.s = float(L[self.nsp+2]) * 1.0e3;  # entropy, J/kg.K
                self.cp = float(L[self.nsp+3]) * 1.0e3; # specific heat, const p,
                self.C_p = self.cp
                self.gam = float(L[self.nsp+4])         # ratio of specific-heats
                self.son = float(L[self.nsp+5])         # sound speed, m/s
                self.a = self.son
                self.rho = float(L[self.nsp+6])         # density, kg/m^3
                self.p = float(L[self.nsp+7]) * 1.0e5   # pressure, Pa
                self.T = float(L[self.nsp+8])           # temperature, Kelvin
                self.R = self.p / (self.rho * self.T);  # gas constant, J/kg.K
                self.C_v = self.C_p - self.R            # specific heat, const volume
                next += 9
            if transProps:
	            self.mu =  float(L[next]) * 1.0e-4     # dynamic viscosity, Pa.s
    	            self.k  =  float(L[next+1]) * 1.0e-1   # W/(m.K)
        else:
	    # use the .out file as this allows more species to be included
            fp = open('tmp.out', 'r')
            lines = fp.readlines()
            fp.close()
            species_data = dict()
            pT_thermo_props = False
            conductivity_found = False
            for s in self.species:
                species_data.setdefault(s,0.0)
            for line in lines:
               if line=="\n": continue
               if line.find("PRODUCTS WHICH WERE CONSIDERED BUT WHOSE")>=0:
                   break
               tks = line.split()
               if line.find("THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED")>=0:
                   pT_thermo_props = True
               elif pT_thermo_props==True:
                   #
                   # Fill out thermo properties if requested
                   if thermoProps:
                       if line.find("H, KJ/KG")>=0:
                           self.h = get_cea2_float(tks[2]) * 1.0e3
                       elif line.find("U, KJ/KG")>=0:
                           self.u = get_cea2_float(tks[2]) * 1.0e3
                           self.e = self.u
                       elif line.find("S, KJ/(KG)(K)")>=0:
                           self.s = get_cea2_float(tks[2]) * 1.0e3
                       elif line.find("Cp, KJ/(KG)(K)")>=0:
                           self.cp = get_cea2_float(tks[2]) * 1.0e3
                           self.C_p = self.cp
                       elif line.find("GAMMAs")>=0:
                           self.gam = get_cea2_float(tks[1])
                       elif line.find("SON VEL,M/SEC")>=0:
                           self.son = get_cea2_float(tks[2])
                           self.a = self.son
                       elif line.find("P, BAR")>=0:
                           self.p = get_cea2_float(tks[2]) * 1.0e5
                           # print "p = ", self.p
                       elif line.find("T, K")>=0:
                           self.T = get_cea2_float(tks[2])
                           # print "T = ", self.T
                       elif line.find("RHO, KG/CU M")>=0:
                           self.rho = get_cea2_float(tks[3])
                           # print "rho = ", self.rho
                   #
                   # get species data
                   for s in self.species:
                       if tks[0]==s or tks[0]==("*"+s):
                           species_data[s] = get_cea2_float(tks[1])
                           # print "%s = %e" % ( s, species_data[s] )
                   #
                   # Fill out trans properties if requested
                   if transProps:
                       if line.find("VISC,MILLIPOISE")>=0:
                           self.mu = get_cea2_float(tks[1]) * 1.0e-4
                           # print "mu = ", self.mu
                       elif conductivity_found==False and line.find("CONDUCTIVITY")>=0 and len(tks)==2:
                           self.k = get_cea2_float(tks[1]) * 1.0e-1
                           # print "k = ", self.k
                           # want to use the first conductivity value (for equilibrium reaction)
                           conductivity_found = True
            #
            # Fill out species properties in requested form
            for i,val in enumerate(species_data.values()):
                if speciesOut == 0:
                    self.eq_massf[i] = val
                else: 
                    self.eq_molef[i] = val
            #
            # Calculation remaining thermo properties if requested
            if thermoProps:
                self.R = self.p / (self.rho * self.T);  # gas constant, J/kg.K
                self.C_v = self.C_p - self.R            # specific heat, const volume
        return
	
    def Shock(self, Us, speciesOut=0, use_out_file=False):
        """
        Computes the equilibrium post-shock gas state using CEA.

        :param self: the gas state to be filled in
	:param speciesOut: a flag where (0) requests massfs, (1) requests molefs

        :returns: the equilibrium post-shock (eps) gas-state as a new class
        """
        if DEBUG_GAS >= 2:
            print 'EOS: Write temporary input file.'
        inp_file_name = 'shock.inp'
        fp = open(inp_file_name, 'w')
        fp.write('# %s generated by cea_gas.py\n' % inp_file_name)
	if self.with_ions:
	    fp.write('problem shock inc eq ions\n')
	else:
	    fp.write('problem shock inc eq\n')
        fp.write('   p(bar)      %e\n' % (self.p / 1.0e5) )
        fp.write('   t(k)        %e\n' % self.T)
	fp.write('   u1          %e\n' % Us)
        # Select the gas components.
        # Look for valid gas names here...
        fp.write('reac\n')
        if self.gasName == 'n2':
            fp.write('   name= N2  wtf= 1.\n')
            fp.write('   name= N   wtf= 0.\n')
            fp.write('only N2 N')
        elif self.gasName == 'co2':
            fp.write('   name= CO2  wtf= 1.\n')
        elif self.gasName == 'air':
            fp.write('   name= Air  wtf= 1.\n')
        elif self.gasName == 'h2ne':
            fp.write('   name= H2  moles= 0.15\n')
            fp.write('   name= Ne  moles= 0.85\n')
        elif self.gasName == 'air5species':
            # Assume air is 79% N2, 21% O2 by mole fraction
            fp.write('   name= N2  wtf= 0.767082\n')
            fp.write('   name= O2  wtf= 0.232918\n')
            fp.write('only N2 O2 N O NO')
        elif self.gasName == 'mix':
            for i in range(len(self.species)):
                if self.massf[i] > 0.0:
                    fp.write('    name= %s  wtf= %g \n'%(self.species[i], self.massf[i]))
	    fp.write('only ')
	    for i in range(len(self.species)):
                fp.write(' %s' % self.species[i])
        fp.write('\n')
        if speciesOut == 0 :                    
            fp.write('output massf trans\n')    # gives us mass fractions (can use trace=1.e-10 here)
        else :
            fp.write('output trans\n')          # gives us mole fractions (can use trace=1.e-10 here)
        fp.write('plot')			       # plotting T, p and mass fractions
	fp.write(' t p rho')
        for i in range(len(self.species)):
            fp.write(' %s' % self.species[i])
        fp.write(' vis cond\n')
        fp.write('end\n') 
        fp.close()
        #
        # Make sure that the database input files are in the working dir
        #
        if not os.path.exists('thermo.inp'):
            print 'Copying thermo.inp to the current working directory'
            os.system("cp %s/e3bin/thermo.inp ." % ( os.getenv("HOME") ) )
            print 'Copying trans.inp to the current working directory'
            os.system("cp %s/e3bin/trans.inp ." % ( os.getenv("HOME") ) )
        #
        # Make sure that binary versions of the database files exist.
        #
        if not os.path.exists('thermo.lib'):
            print 'Make the binary database for thermodynamic properties'
            run_cea_program('thermo')
            print 'Make the binary database for transport properties'
            run_cea_program('trans')
        #
        # Now, run the cea program in anger.
        run_cea_program('shock')
        #
        # Create a instance of Gas to store the equilibrium post-shock gas-state
        eps = Gas( "mix", self.species, self.massf, use_out_file )
        #
        # Pick out the interesting bits from the plotfile.
        # These data will be on the first non-comment line.
        if use_out_file==False:
            # use the .plt file as per usual
            fp = open('shock.plt', 'r')
            L = string.split(fp.readline())
            while L[0] == "#":
                # Skip over comment lines
                L = string.split(fp.readline())
            fp.close()
            expected_entries = len(self.species) + 3 + 2
            if len(L)!=expected_entries:
                print "cea returned %d values, while %d values were requested." % ( len(L), expected_entries )
                print "This probably means too many values were requested."
                print "Try calling the Shock function with 'use_out_file=True'."
                sys.exit()
            # Fill out the equilibrium post-shock gas-state
	    eps.T = float(L[0])
	    eps.p = float(L[1])*1.0e5
            eps.rho = float(L[2])
	    # print "T = %f K, p = %f bar " % ( eps.T, eps.p*1.0e-5 ) 
            if speciesOut == 0 : 
                for i in range(3,3+len(self.species)):
	            # print "eq_massf[%d] = %f" % ( i-3, float(L[i]) )
                    eps.eq_massf[i-3] = float(L[i])
            else :
                for i in range(3,3+len(self.species)):
                    eps.eq_molef[i-3] = float(L[i])
	    eps.mu =  float(L[3+len(self.species)]) * 1.0e-4
	    eps.k  =  float(L[3+len(self.species)+1]) * 1.0e-1
        else:
	    # use the .out file as this allows more species to be included
            fp = open('shock.out', 'r')
            lines = fp.readlines()
            fp.close()
            incident_shock_data = False
            species_data = dict()
            conductivity_found = False
            for s in self.species:
                species_data.setdefault(s,0.0)
            for line in lines:
               if line=="\n": continue
               if line.find("PRODUCTS WHICH WERE CONSIDERED BUT WHOSE")>=0:
                   break
               tks = line.split()
               if line.find("SHOCKED GAS (2)--INCIDENT--EQUILIBRIUM")>=0:
                   incident_shock_data = True
               elif incident_shock_data==True:
                   # thermo data
                   if line.find("P, BAR")>=0:
                       eps.p = get_cea2_float(tks[2]) * 1.0e5
                       # print "p = ", eps.p
                   elif line.find("T, K")>=0:
                       eps.T = get_cea2_float(tks[2])
                       # print "T = ", eps.T
                   elif line.find("RHO, KG/CU M")>=0:
                       eps.rho = get_cea2_float(tks[3])
                       # print "rho = ", eps.rho
                   # species data
                   for s in self.species:
                       if tks[0]==s or tks[0]==("*"+s):
                           species_data[s] = get_cea2_float(tks[1])
                           # print "%s = %e" % ( s, species_data[s] )
                   # transport data
                   if line.find("VISC,MILLIPOISE")>=0:
                       eps.mu = get_cea2_float(tks[1]) * 1.0e-4
                       # print "mu = ", eps.mu
                   elif conductivity_found==False and line.find("CONDUCTIVITY")>=0 and len(tks)==2:
                       eps.k = get_cea2_float(tks[1]) * 1.0e-1
                       # print "k = ", eps.k
                       # want to use the first conductivity value (for equilibrium reaction)
                       conductivity_found = True
            #
            # Fill out species properties in requested form
            for i,val in enumerate(species_data.values()):
                if speciesOut == 0:
                    eps.eq_massf[i] = val
                else: 
                    eps.eq_molef[i] = val
	return eps

# --------------------------------------------------------------

if __name__ == '__main__':
    print 'Test the Gas class...'
    print 'Default constructor with air as the test gas...'
    a = Gas('air')
    a.write_state(sys.stdout)
    print 'and the same air at a higher temperature'
    a.set_from_pAndT(100.0e3,4000.0)
    a.write_state(sys.stdout)
    print 'Try a user-specified mix of Helium, N2 and N'
    b = Gas('mix',['N2','N','He'], [1.0, 0.0, 0.0])
    b.write_state(sys.stdout)
    print 'and the same initial mix at a higher temperature'
    b.set_from_pAndT(263.79e3,5719.0)
    b.write_state(sys.stdout)
    print 'End of test.'


#! /usr/bin/env python
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
   19-Feb-12: some refactoring and general clean-up by PJ
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
    outFile = jobName + '.out'
    pltFile = jobName + '.plt'
    if os.path.exists(inpFile):
        if DEBUG_GAS >= 2:
            print 'cea2_gas: Start cea program on job %s...' % jobName
        # We should remove the results files from previous runs. 
        if os.path.exists(pltFile): os.remove(pltFile)
        if os.path.exists(outFile): os.remove(outFile)
        p = subprocess.Popen(CEA_COMMAND_NAME, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate(jobName + '\n')
        return_code = p.wait()
        if DEBUG_GAS >= 2:
            print('cea2_gas: %s finished job %s.' % (CEA_COMMAND_NAME,jobName))
        if return_code != 0:
            print('cea2_gas: return-code from cea2 program is nonzero.')
            raise Exception,  'cea2-return-code = %d' % return_code
        fp = open(outFile, 'r')
        outFileText = fp.read()
        outFileIsBad = False
        # Look for the summary table header.
        if outFileText.find('THERMODYNAMIC PROPERTIES') == -1: outFileIsBad = True
        if outFileIsBad:
            print('cea2_gas: the output file seems incomplete; you should go check.')
            raise Exception, 'cea2_gas: detected badness in cea2 output file.'
    else:
        raise Exception, 'cea2_gas: The file %s is not present.' % inpFile

def get_cea2_float(token_list):
    """
    Clean up the CEA2 short-hand notation for exponential format.

    CEA2 seems to write exponential-format numbers in a number of ways:

    | 1.023-2
    | 1.023+2
    | 1.023 2
    """
    if len(token_list) == 0:
        value_str = '0.0'
    elif len(token_list) == 1:
        value_str = token_list[0]
        if value_str.find("-")>0:
            value_str = value_str.replace("-","e-")
        if value_str.find("+")>0:
            value_str = value_str.replace("+","e+")
    elif len(token_list) == 2:
        value_str = token_list[0] + 'e+' + token_list[1]
    else:
        print "get_cea2_float(): too many tokens:", token_list
        value_str = '0.0'
    return float(value_str)
   
# ----------------------------------------------------------------

class Gas(object):
    """
    Provides the equation of state for the gas.
    """
    def __init__(self, gasName, speciesList=[], massfList=[], use_out_file=False):
        """
        Set up a new obects, from either a name of species list.

        :param gasName: The generic name of the gas:
            'n2', 'co2', 'air', 'h2ne', 'air5species', 'mix'
        :param speciesList: if gasName is 'mix', this list of strings is used
            to specify the species in the mix
        :param massfList: if gasName is 'mix', this list of floats is used
            to specify the mass-fractions of the species
        :param use_out_file: flag to indicate which cea2 output file to scan:
            False==plt file (with limited species); True==text output file. 
            Note that the text output file allows more species to be handled
            then the .plt file but that there may be small differences in the 
            least significant digits.  The text file seems to have more digits
            written for a floating-point number.

        .. todo: Need to think more carefully about the mass-fraction, mole-fraction data.
           Presently, if mass-fractions are requested, mole-fractions do not get computed.
        """
	if locate_executable_file(CEA_COMMAND_NAME) is None:
            print "Could not find the executable program %s" % CEA_COMMAND_NAME
            print "The chemical equilibrium-analysis program is external"
            print "to the cfcfd3 code collection and needs to be obtained from NASA Glenn."
            print "Quitting the current program because we really can't do anything further."
            sys.exit()
	# ----------------------------------------------------------------
        self.gasName = gasName                # see method EOS for possible names
        self.species = copy.copy(speciesList) # species names as per CEA database
        self.nsp     = len(speciesList)       # number of species
        self.massf   = copy.copy(massfList)   # user-specified mass fractions
        self.eq_massf = copy.copy(self.massf) # will overwrite this list later
        self.eq_molef = copy.copy(self.massf) # also to be written over
        self.use_out_file = use_out_file      # True == use text output to enable more species.
                                              # False == use plt file
	self.with_ions = "e-" in speciesList
        self.p = 1.0e5 # need a value in case EOS gets called
        self.T = 300.0 # likewise
        self.Us = 0.0 # m/s
	self.EOS(speciesOut=1, thermoProps=1, transProps=1, problemType='pT')
        return

    def clone(self):
        """
        Clone the current Gas object to make another, just the same.

        :returns: the new Gas object.
        """
        other = Gas(self.gasName, self.species, self.massf, self.use_out_file)
        other.p = self.p
        other.T = self.T
        other.Us = self.Us
	other.EOS(speciesOut=1, thermoProps=1, transProps=1, problemType='pT')
        return other

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

    def set_from_pAndT(self, p, T, speciesOut=0, thermoProps=1, transProps=1):
        """
        Fills out gas state from given p and T.
        """
        self.p = p;        # pressure, Pa
        self.T = T;        # temperature, K
        flag = self.EOS(speciesOut, thermoProps, transProps)
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

    def write_cea2_input_file(self, speciesOut, thermoProps, transProps, problemType):
        """
        Set up a problem-description file for CEA2.

        :param speciesOut: format for species fractions: 
            0=mass fractions, 1=mole fractions, -1=do not report species fractions 
	:param thermoProps: a flag:
            0=don't request thermoprops, 1=request Cv, Cp, entropy, a, gamma etc
	:param transProps: a flag:
            0=don't request transport props, 1=request viscosity and thermal-conductivity
        :param problemType: a string specifying type of CEA analysis: 
            'pT', 'rhoT', 'ps', 'shock'
        :returns: None
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
            assert self.rho > 0.0
            assert self.T > 0.0
            fp.write('   rho,kg/m**3 %e\n' % self.rho)
            fp.write('   t(k)        %e\n' % self.T)
            if DEBUG_GAS >= 2:
                print 'EOS: input to CEA2 rho: %g, T: %g' % (self.rho, self.T)
        elif problemType == 'pT':
            if self.with_ions:
                fp.write('problem case=estcj tp ions\n')
            else:
                fp.write('problem case=estcj tp\n')
            assert self.p > 0.0, self.T > 0.0
            fp.write('   p(bar)      %e\n' % (self.p / 1.0e5) )
            fp.write('   t(k)        %e\n' % self.T)
            if DEBUG_GAS >= 2:
                print 'EOS: input to CEA2 p: %g, T: %g' % (self.p, self.T)
        elif problemType == 'ps':
            if self.with_ions:
                fp.write('problem case=estcj ps ions\n')
            else:
                fp.write('problem case=estcj ps\n')
            assert self.p > 0.0
            fp.write('   p(bar)      %e\n' % (self.p / 1.0e5) )
            fp.write('   s/r         %e\n' % (self.s / R_universal) )
            if DEBUG_GAS >= 2:
                print 'EOS: input to CEA2 p: %g, s: %g' % (self.p, self.s)
        elif problemType == 'shock':
            if self.with_ions:
                fp.write('problem shock inc eq ions\n')
            else:
                fp.write('problem shock inc eq\n')
            assert self.p > 0.0, self.T > 0.0
            fp.write('   p(bar)      %e\n' % (self.p / 1.0e5) )
            fp.write('   t(k)        %e\n' % self.T)
            fp.write('   u1          %e\n' % self.Us)
        else:
            raise Exception, 'cea2_gas: Invalid problemType: %s' % problemType
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
                    assert self.massf[i] >= 0.0
                    fp.write('   name= %s  wtf= %g \n'%(self.species[i], self.massf[i]))
	    fp.write('only')
	    for i in range(self.nsp):
                fp.write(' %s' % self.species[i])
        fp.write('\n')
        if speciesOut == 0 :                    
            fp.write('output massf trans\n') # gives us mass fractions (can use trace=1.e-10 here)
        else :
            fp.write('output trans\n') # gives us mole fractions (can use trace=1.e-10 here)
        fp.write('plot') # only plotting mass fractions
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
        return

    def scan_cea2_dot_plt_file(self, speciesOut, thermoProps, transProps):
        """
        Scan the plot file generated by CEA2 and extract our gas-properties data.

        :param speciesOut: format for species fractions: 
            0=mass fractions, 1=mole fractions, -1=do not report species fractions 
	:param thermoProps: a flag:
            0=don't request thermoprops, 1=request Cv, Cp, entropy, a, gamma etc
	:param transProps: a flag:
            0=don't request transport props, 1=request viscosity and thermal-conductivity
        :returns: None, but does update the contents of the gas state as a side-effect.
        """
        # Pick out the interesting bits from the plotfile.
        # These data will be on the first non-comment line.
        fp = open('tmp.plt', 'r')
        L = string.split(fp.readline())
        while L[0] == "#":
            # Skip over comment lines
            L = string.split(fp.readline())
        # We should now have the first, non-comment line.
        fp.close()
        expected_entries = len(self.species)
        if thermoProps:
            expected_entries += 9
        if transProps:
            expected_entries += 2
        if len(L) != expected_entries:
            print "cea .plt file returned %d values, while %d values were requested." \
                % (len(L), expected_entries)
            print "This probably means too many values were requested."
            print "Try calling the EOS function with 'use_out_file=True'."
            raise Exception, 'cea2_gas: wrong number of data elements in .plt file'
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
            self.mu = float(L[next]) * 1.0e-4       # dynamic viscosity, Pa.s
            self.k  = float(L[next+1]) * 1.0e-1     # W/(m.K)
        return

    def scan_cea2_dot_out_file(self, speciesOut, thermoProps, transProps):
        """
        Scan the output text file generated by CEA2 and extract our gas-properties data.

        :param speciesOut: format for species fractions: 
            0=mass fractions, 1=mole fractions, -1=do not report species fractions 
	:param thermoProps: a flag:
            0=don't request thermoprops, 1=request Cv, Cp, entropy, a, gamma etc
	:param transProps: a flag:
            0=don't request transport props, 1=request viscosity and thermal-conductivity
        :returns: None, but does update the contents of the gas state as a side-effect.
        """
        # use the .out file as this allows more species to be included
        fp = open('tmp.out', 'r')
        lines = fp.readlines()
        fp.close()
        species_data = dict()
        thermo_props_found = False
        conductivity_found = False
        incident_shock_data = False
        for s in self.species:
            species_data.setdefault(s,0.0)
        for line in lines:
           if line=="\n": continue
           if line.find("PRODUCTS WHICH WERE CONSIDERED BUT WHOSE")>=0:
               break
           tks = line.split()
           if line.find("THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED")>=0:
               thermo_props_found = True
           elif line.find("SHOCKED GAS (2)--INCIDENT--EQUILIBRIUM")>=0:
               incident_shock_data = True
           elif thermo_props_found or incident_shock_data:
               # Fill out thermo properties if requested
               if thermoProps:
                   if line.find("H, KJ/KG")>=0:
                       self.h = get_cea2_float(tks[2:]) * 1.0e3
                   elif line.find("U, KJ/KG")>=0:
                       self.u = get_cea2_float(tks[2:]) * 1.0e3
                       self.e = self.u
                   elif line.find("S, KJ/(KG)(K)")>=0:
                       self.s = get_cea2_float(tks[2:]) * 1.0e3
                   elif line.find("Cp, KJ/(KG)(K)")>=0:
                       self.cp = get_cea2_float(tks[2:]) * 1.0e3
                       self.C_p = self.cp
                   elif line.find("GAMMAs")>=0:
                       self.gam = get_cea2_float(tks[1:])
                   elif line.find("SON VEL,M/SEC")>=0:
                       self.son = get_cea2_float(tks[2:])
                       self.a = self.son
                   elif line.find("P, BAR")>=0:
                       self.p = get_cea2_float(tks[2:]) * 1.0e5
                       # print "p = ", self.p
                   elif line.find("T, K")>=0:
                       self.T = get_cea2_float(tks[2:])
                       # print "T = ", self.T
                   elif line.find("RHO, KG/CU M")>=0:
                       self.rho = get_cea2_float(tks[3:])
                       # print "rho = ", self.rho
               # Scan tokens for species data
               for s in self.species:
                   if tks[0]==s or tks[0]==("*"+s):
                       species_data[s] = get_cea2_float(tks[1:])
                       # print "%s = %e" % ( s, species_data[s] )
               # Fill out transport properties if requested
               if transProps:
                   if line.find("VISC,MILLIPOISE")>=0:
                       self.mu = get_cea2_float(tks[1:]) * 1.0e-4
                       # print "mu = ", self.mu
                   elif conductivity_found==False and line.find("CONDUCTIVITY")>=0 and len(tks)==2:
                       self.k = get_cea2_float(tks[1:]) * 1.0e-1
                       # print "k = ", self.k
                       # want to use the first conductivity value (for equilibrium reaction)
                       conductivity_found = True
        # Fill out species properties in requested form
        for i,val in enumerate(species_data.values()):
            if speciesOut == 0:
                self.eq_massf[i] = val
            else: 
                self.eq_molef[i] = val
        # Calculate remaining thermo properties, if requested
        if thermoProps:
            self.R = self.p / (self.rho * self.T);  # gas constant, J/kg.K
            self.C_v = self.C_p - self.R            # specific heat, const volume
        return

    def EOS(self, speciesOut=0, thermoProps=1, transProps=1, problemType='pT'):
        """
        Computes the gas state, taking into account the high-temperature effects.

        It does this by writing a suitable input file for the CEA code,
        calling that code and then extracting the relevant results from 
        the CEA output or plot file.

        :param self: the gas state to be filled in
        :param speciesOut: format for species fractions: 
            0=mass fractions, 1=mole fractions, -1=do not report species fractions 
	:param thermoProps: a flag:
            0=don't request thermoprops, 1=request Cv, Cp, entropy, a, gamma etc
	:param transProps: a flag:
            0=don't request transport props, 1=request viscosity and thermal-conductivity
        :param problemType: a string specifying the type of CEA analysis:
            'pT', 'rhoT', 'ps', shock
        :returns: None, but does update the contents of the gas state as a side-effect.
        """
        self.write_cea2_input_file(speciesOut, thermoProps, transProps, problemType)
        # Make sure that the database input files are in the working dir
        if not os.path.exists('thermo.inp'):
            print 'Copying thermo.inp to the current working directory'
            os.system("cp %s/e3bin/thermo.inp ." % ( os.getenv("HOME") ) )
            print 'Copying trans.inp to the current working directory'
            os.system("cp %s/e3bin/trans.inp ." % ( os.getenv("HOME") ) )
        # Make sure that binary versions of the database files exist.
        if not os.path.exists('thermo.lib'):
            print 'Make the binary database for thermodynamic properties'
            run_cea_program('thermo')
            print 'Make the binary database for transport properties'
            run_cea_program('trans')
        # Now, run the cea program on the actual job.
        run_cea_program('tmp')
        # ... and scan the output to extract the gas-properties data.
        if self.use_out_file:
            self.scan_cea2_dot_out_file(speciesOut, thermoProps, transProps)
        else:
            self.scan_cea2_dot_plt_file(speciesOut, thermoProps, transProps)
        return

    def shock_process(self, Us, speciesOut=0):
        """
        Compute the gas state after being processed by an incident shock.

        :param Us: shock speed into quiescent gas, m/s
        :param speciesOut: format for species fractions: 
            0=mass fractions, 1=mole fractions, -1=do not report species fractions 
        :returns: a reference to the post-shock gas state (self)

        .. This recovers (approximately) Dan's original Shock function.
        """
        self.Us = Us
        self.EOS(speciesOut=speciesOut, thermoProps=1, transProps=1, problemType='shock')
        return self

# --------------------------------------------------------------

if __name__ == '__main__':
    print 'Test the Gas class...'
    #
    print '\nDefault constructor with air as the test gas...'
    a = Gas('air')
    a.write_state(sys.stdout)
    print 'and the same air at a higher temperature'
    a.set_from_pAndT(100.0e3,4000.0)
    a.write_state(sys.stdout)
    #
    print '\nTry a user-specified mix of Helium, N2 and N'
    b = Gas('mix', ['N2','N','He'], [1.0, 0.0, 0.0])
    b.write_state(sys.stdout)
    print 'and the same initial mix at a higher temperature'
    b.set_from_pAndT(263.79e3,5719.0)
    b.write_state(sys.stdout)
    #
    print '\nStart again with low-T air as the test gas'
    a = Gas('air', use_out_file=True)
    c = a.clone()
    c.write_state(sys.stdout)
    print 'and shock process it'
    c.shock_process(4000.0)
    c.write_state(sys.stdout)
    #
    print 'End of test.'


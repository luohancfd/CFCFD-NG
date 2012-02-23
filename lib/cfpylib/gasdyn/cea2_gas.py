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
from copy import copy

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

def run_cea_program(jobName,checkTableHeader=True):
    """
    Runs the CEA program on the specified job.

    :param jobName: string that is used to construct input and output file names
    :param checkTableHeader: boolean flag to activate checking of output file 
        table header.  We use this as a test to see if the cea2 program has run
        the job succesfully.
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
        if checkTableHeader:
            # Look for the summary table header 
            if outFileText.find('THERMODYNAMIC PROPERTIES') == -1:
            	outFileIsBad = True
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
        print "get_cea2_float(): too many tokens (expected one or two, only):", token_list
        value_str = '0.0'
    return float(value_str)
   
# ----------------------------------------------------------------

class Gas(object):
    """
    Provides the equation of state for the gas.
    """
    def __init__(self, reactants=[], fractions=[], onlyList=[],
                 inputUnits='massf', outputUnits='massf', 
                 with_ions=False, trace=1.0e-10):
        """
        Set up a new obects, from either a name of species list.

        :param reactants: list of strings is used to specify the reactants in the mix.
            The names are as per the CEA database.
            Note that other chemical species may be added to the mix by cea2.
        :param fractions: list of floats is used to specify the mass- 
            or mole-fractions of the reactants
        :param onlyList: list of strings limiting the species in the mix.
        :param inputUnits: string 'moles' or 'massf'
        :param outputUnits: string 'moles' or 'massf'
        :param with_ions: boolean flag indicating whether electrons and ions 
            should be included in the mix
        :param trace: fraction below which a species will be neglected in the output
        """
	if locate_executable_file(CEA_COMMAND_NAME) is None:
            print "Could not find the executable program %s" % CEA_COMMAND_NAME
            print "The chemical equilibrium-analysis program is external"
            print "to the cfcfd3 code collection and needs to be obtained from NASA Glenn."
            print "Quitting the current program because we really can't do anything further."
            sys.exit()
	# ----------------------------------------------------------------
        self.reactants = copy(reactants)
        self.fractions = copy(fractions)
        self.inputUnits = inputUnits
        self.outputUnits = outputUnits
        self.onlyList = copy(onlyList)
        # Make a full list of all species that we expect to have in the mix.
        # FIX_ME this won't work as cea2 will add whatever it likes...
        self.species = copy(reactants)
        self.speciesFraction = copy(fractions)
        for s in self.onlyList:
            if s not in self.species: 
                self.species.append(s)
                self.speciesFraction.append(0.0)
        self.nsp = len(self.species)
	self.with_ions = with_ions or ("e-" in self.species)
        self.trace = trace
        self.Us = 0.0 # m/s
        self.have_run_cea = False
        return

    def clone(self):
        """
        Clone the current Gas object to make another, just the same.

        :returns: the new Gas object.
        """
        other = Gas(self.reactants, self.fractions, self.onlyList, 
                    self.inputUnits, self.outputUnits, self.with_ions)
        if self.have_run_cea:
            other.p = self.p
            other.T = self.T
            other.Us = self.Us
            other.EOS(problemType='pT', transProps=True)
        return other

    def set_pT(self, p, T, transProps=True):
        """
        Fills out gas state from given p and T.
        """
        self.p = p;        # pressure, Pa
        self.T = T;        # temperature, K
        flag = self.EOS(problemType='pT', transProps=transProps)
        return flag

    def write_state(self, strm):
        """
        Writes the gas state data to the specified stream.
        """
        strm.write('    p: %g Pa, T: %g K, rho: %g kg/m**3, e: %g J/kg, a: %g m/s\n'
                   % (self.p, self.T, self.rho, self.u, self.son) )
        strm.write('    R: %g J/(kg.K), gam: %g, Cp: %g J/(kg.K), mu: %g Pa.s, k: %g W/(m.K)\n'
                   % (self.R, self.gam, self.cp, self.mu, self.k) )
        for i in range(self.nsp):
            strm.write('%s fraction for species[%d] %s %e\n' %
                       (self.outputUnits, i, self.species[i], self.speciesFraction[i]))

    def write_cea2_input_file(self, problemType, transProps):
        """
        Set up a problem-description file for CEA2.

        :param problemType: a string specifying type of CEA analysis that is requested: 
            'pT', 'rhoT', 'ps', 'shock'
	:param transProps: a boolean flag:
            False=don't request transport props, True=request viscosity and thermal-conductivity
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
        fp.write('reac\n')
        for i in range(len(self.reactants)):
            if self.fractions[i] > 0.0:
                fp.write('   name= %s  %s=%g\n' % (self.species[i], self.inputUnits, self.fractions[i]))
        #
        if len(self.onlyList) > 0:
            fp.write('only %s\n' % (' '.join(self.onlyList)))
        #
        fp.write('output')
        if self.outputUnits == 'massf': fp.write(' MASSF')
        fp.write(' trace=%e' % self.trace)
        if transProps: fp.write(' trans')
        fp.write('\n')
        #
        fp.write('end\n') 
        fp.close()
        return

    def scan_cea2_dot_out_file(self, transProps):
        """
        Scan the output text file generated by CEA2 and extract our gas-properties data.

	:param transProps: a boolean flag:
            False=don't request transport props, True=request viscosity and thermal-conductivity
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
           if line.find("PRODUCTS WHICH WERE CONSIDERED BUT WHOSE")>=0: break
           if line.find("THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED")>=0:
               thermo_props_found = True
           elif line.find("SHOCKED GAS (2)--INCIDENT--EQUILIBRIUM")>=0:
               incident_shock_data = True
           elif thermo_props_found or incident_shock_data:
               tks = line.split()
               # Fill out thermo properties
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
        # Calculate remaining thermo properties
        self.R = self.p / (self.rho * self.T);  # gas constant, J/kg.K
        self.C_v = self.C_p - self.R            # specific heat, const volume
        return

    def EOS(self, problemType='pT', transProps=True):
        """
        Computes the gas state, taking into account the high-temperature effects.

        It does this by writing a suitable input file for the CEA code,
        calling that code and then extracting the relevant results from 
        the CEA output or plot file.

        :param self: the gas state to be filled in
        :param problemType: a string specifying the type of CEA analysis:
            'pT', 'rhoT', 'ps', shock
	:param transProps: a boolean flag:
            False=don't request transport props, True=request viscosity and thermal-conductivity
        :returns: None, but does update the contents of the gas state as a side-effect.
        """
        # Make sure that the database input files are in the working dir
        if not os.path.exists('thermo.inp'):
            print 'Copying thermo.inp to the current working directory'
            os.system("cp %s/e3bin/thermo.inp ." % ( os.getenv("HOME") ) )
            print 'Copying trans.inp to the current working directory'
            os.system("cp %s/e3bin/trans.inp ." % ( os.getenv("HOME") ) )
        # Make sure that binary versions of the database files exist.
        if not os.path.exists('thermo.lib'):
            print 'Make the binary database for thermodynamic properties'
            run_cea_program('thermo',checkTableHeader=False)
            print 'Make the binary database for transport properties'
            run_cea_program('trans',checkTableHeader=False)
        # Now, run the cea program on the actual job.
        self.write_cea2_input_file(problemType, transProps)
        run_cea_program('tmp')
        self.scan_cea2_dot_out_file(transProps)
        self.have_run_cea = True
        return

    def shock_process(self, Us):
        """
        Compute the gas state after being processed by an incident shock.

        :param Us: shock speed into quiescent gas, m/s
        :returns: a reference to the post-shock gas state (self)

        .. This recovers (approximately) Dan's original Shock function.
        """
        self.Us = Us
        self.EOS(problemType='shock', transProps=True)
        return self

# --------------------------------------------------------------

if __name__ == '__main__':
    print 'Test the Gas class...'
    #
    print '\nDefault constructor with air as the test gas...'
    a = Gas(['Air'], [1.0,], outputUnits='moles')
    a.set_pT(100.0e3, 300.0)
    a.write_state(sys.stdout)
    print 'and the same air at a higher temperature'
    a.set_pT(100.0e3, 4000.0)
    a.write_state(sys.stdout)

    sys.exit()
    #
    print '\nAir-5-species for Luke: 79% N2, 21% O2 by mole fraction.'
    a = Gas(reactants=['N2','O2','N','O','NO'], 
            fractions=[0.79,0.21,0.0,0.0,0.0],
            inputUnits='moles', outputUnits='massf',
            onlyList=['N2','O2','N','O','NO'])
    a.set_pT(100.0e3, 300.0)
    a.write_state(sys.stdout)
    print 'and the same air at a higher temperature'
    a.set_pT(100.0e3, 4000.0)
    a.write_state(sys.stdout)
    #
    print '\nTry a user-specified mix of Helium, N2 and N'
    b = Gas(['N2','N','He'], [1.0, 0.0, 0.0])
    b.set_pT(100.0e3, 300.0)
    b.write_state(sys.stdout)
    print 'and the same initial mix at a higher temperature'
    b.set_pT(263.79e3, 5719.0)
    b.write_state(sys.stdout)
    #
    print '\nStart again with low-T air as the test gas'
    a = Gas(['Air'], [1.0,]); a.set_pT(100.0e3, 300.0)
    c = a.clone()
    c.write_state(sys.stdout)
    print 'and shock process it'
    c.shock_process(4000.0)
    c.write_state(sys.stdout)
    #
    print 'End of test.'


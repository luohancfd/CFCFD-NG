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
   PA Jacobs RJ Gollan and DF Potter
   Institute of Aerodynamics and Flow Technology
   The German Aerospace Center, Goettingen.
   and
   School of Mechanical Engineering
   The University of Queensland

.. Versions:
   24-Dec-02: First code.
   10-May-04: Updated for a mix of species.
   06-Feb-05: renamed to cea_gas.py
   28-Feb-08: Added a get_eq_massf() access function.
   28-Fed-08: Major changes to allow proper calculation at high temps.
   11-Dec-08: Addition of basic incident Shock function
   19-Feb-12: some refactoring, simplification and general clean-up
"""

import sys, string, math, os, subprocess, re
from copy import copy

# -------------------------------------------------------------------
# First, global data.

DEBUG_GAS = 0
R_universal = 8314.0;  # J/kgmole.K

# Set name for cea executable. If we are not on a Windows
# machine then assume we are on a Linux-like machine
if sys.platform.startswith('win'):
    CEA_COMMAND_NAME = 'fcea2.exe'
else:
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
        #
        # At the highest level of estcj we have added
        # e3bin and local estcj path to sys.path. Searching 
        # over sys.path  ensures that estcj/cea2_gas will 
        # work on Windows machines. Luke D. 24-May-12
        for path in sys.path:
            fullName = os.path.join(path, name)
            if is_exe(fullName): return fullName
        # Note that sys.path is initialized from PYTHONPATH,
        # at least on linux machines,
        # so we might need to search the PATH as well.  PJ 25-Jul-12
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
        if value_str.find("****") >= 0:
            # We have one of the dodgy strings such as ******e-3
            # CEA2 seems to write such for values like 0.0099998
            # We should be able to tolerate one of these, at most,
            # because we should be able to back out the intended
            # value from knowledge of the rest of the list.
            return None
        if value_str.find("-") > 0:
            value_str = value_str.replace("-","e-")
        if value_str.find("+") > 0:
            value_str = value_str.replace("+","e+")
    elif len(token_list) == 2:
        value_str = token_list[0] + 'e+' + token_list[1]
    else:
        print "get_cea2_float(): too many tokens (expected one or two, only):", token_list
        value_str = '0.0'
    try:
        value = float(value_str)
    except:
        print "Cannot make a float from this string: ", value_str
        sys.exit(-1)
    return value
   
# ----------------------------------------------------------------

class Gas(object):
    """
    Provides the equation of state for the gas.
    """
    def __init__(self, reactants={}, onlyList=[],
                 inputUnits='massf', outputUnits='massf', 
                 with_ions=False, trace=1.0e-6):
        """
        Set up a new object, from either a name of species list.

        :param reactants: dictionary of reactants and their mixture fractions
            The keys used to specify the reactants in the mix 
            and the (float) values are their mass- or mole-fractions.
            The names are as per the CEA database.
            Note that other chemical species may be added to the mix by cea2.
        :param onlyList: list of strings limiting the species in the mix.
        :param inputUnits: string 'moles' or 'massf'
        :param outputUnits: string 'moles' or 'massf'
        :param with_ions: boolean flag indicating whether electrons and ions 
            should be included in the mix
        :param trace: fraction below which a species will be neglected in CEA
        """
	if locate_executable_file(CEA_COMMAND_NAME) is None:
            print "Could not find the executable program %s" % CEA_COMMAND_NAME
            print "The chemical equilibrium-analysis program is external"
            print "to the cfcfd3 code collection and needs to be obtained from NASA Glenn."
            print "Quitting the current program because we really can't do anything further."
            sys.exit()
	# ----------------------------------------------------------------
        assert inputUnits == 'moles' or inputUnits == 'massf'
        assert outputUnits == 'moles' or outputUnits == 'massf'
        self.reactants = copy(reactants)
        self.inputUnits = inputUnits
        self.outputUnits = outputUnits
        self.onlyList = copy(onlyList)
        self.species = {} # will be read from CEA2 output
	self.with_ions = with_ions or ('e-' in self.reactants.keys()) or ('e-' in self.onlyList)
        self.trace = trace
        self.Us = 0.0 # m/s
        self.have_run_cea = False
        return

    def clone(self, newOutputUnits=None):
        """
        Clone the current Gas object to make another, just the same.

        :returns: the new Gas object.
        """
        if newOutputUnits == None: newOutputUnits = self.outputUnits
        other = Gas(self.reactants, self.onlyList, self.inputUnits,
                    newOutputUnits, self.with_ions, self.trace)
        if self.have_run_cea:
            other.p = self.p
            other.T = self.T
            other.Us = self.Us
            other.trace = self.trace
            other.EOS(problemType='pT', transProps=True)
        return other

    def set_pT(self, p, T, transProps=True):
        """
        Fills out gas state from given pressure and temperature.

        :param p: pressure, Pa
        :param T: temperature, K
        """
        self.p = p; self.T = T
        return self.EOS(problemType='pT', transProps=transProps)

    def set_rhoT(self, rho, T, transProps=True):
        """
        Fills out gas state from given density and temperature.

        :param rho: density, kg/m**3
        :param T: temperature, K
        """
        self.rho = rho; self.T = T
        return self.EOS(problemType='rhoT', transProps=transProps)

    def set_rhoe(self, rho, e, transProps=True):
        """
        Fills out gas state from given density and internal energy.

        :param rho: density, kg/m**3
        :param e: internal energy of mixture, J/kg
        """
        self.rho = rho; self.e = e
        return self.EOS(problemType='rhoe', transProps=transProps)

    def set_ps(self, p, s, transProps=True):
        """
        Fills out gas state from given pressure and specific entropy.

        :param p: pressure, Pa
        :param s: entropy, J/(kg.K)
        """
        self.p = p; self.s = s
        return self.EOS(problemType='ps', transProps=transProps)

    def set_ph(self, p, h, transProps=True):
        """
        Fills out gas state from given pressure and enthalpy.

        :param p: pressure, Pa
        :param h: enthalpy, J/kg
        """
        self.p = p; self.h = h
        return self.EOS(problemType='ph', transProps=transProps)

    def write_state(self, strm):
        """
        Writes the gas state data to the specified stream.
        """
        strm.write('    p: %g Pa, T: %g K, rho: %g kg/m**3, e: %g J/kg, h: %g J/kg, a: %g m/s, s: %g J/(kg.K)\n'
                   % (self.p, self.T, self.rho, self.e, self.h, self.a, self.s) )
        strm.write('    R: %g J/(kg.K), gam: %g, Cp: %g J/(kg.K), mu: %g Pa.s, k: %g W/(m.K)\n'
                   % (self.R, self.gam, self.cp, self.mu, self.k) )
        strm.write('    species %s: %s\n' % (self.outputUnits, str(self.species)) )
        return

    def get_fractions(self, speciesList):
        """
        Gets a list of mole- or mass-fractions for the specified species.

        :param speciesList: the species names for which we want a list of fractions.
        :returns: list of floats representing the fractions of each species in the mix
            Note that the mass-fractions or mole-fractions are returned, based on
            the value of outputUnits in the Gas object.
        """
        fractionList = []
        for s in speciesList:
            if s in self.species.keys(): 
                fractionList.append(self.species[s])
            else:
                fractionList.append(0.0)
        return fractionList

    def write_cea2_input_file(self, problemType, transProps):
        """
        Set up a problem-description file for CEA2.

        :param problemType: a string specifying type of CEA analysis that is requested: 
            'pT', 'rhoT', 'rhoe', 'ps', 'shock'
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
        elif problemType == 'rhoe':
            if self.with_ions:
                fp.write('problem case=estcj vu ions\n')
            else:
                fp.write('problem case=estcj vu\n')
            assert self.rho > 0.0
            fp.write('   rho,kg/m**3 %e\n' % self.rho)
            fp.write('   u/r         %e\n' % (self.e / R_universal) )
            if DEBUG_GAS >= 2:
                print 'EOS: input to CEA2 rho: %g, e: %g' % (self.rho, self.e)
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
                print 'EOS: input to CEA2 p: %g, s/r: %g' % (self.p, self.s)
        elif problemType == 'ph':
            if self.with_ions:
                fp.write('problem case=estcj ph ions\n')
            else:
                fp.write('problem case=estcj ph\n')
            assert self.p > 0.0
            fp.write('   p(bar)      %e\n' % (self.p / 1.0e5) )
            fp.write('   h/r         %e\n' % (self.h / R_universal) )
            if DEBUG_GAS >= 2:
                print 'EOS: input to CEA2 p: %g, h/r: %g' % (self.p, self.s)
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
        for s in self.reactants.keys():
            f = self.reactants[s]
            if f > 0.0:
                if self.inputUnits == 'moles':
                    fp.write('   name= %s  moles=%g' % (s, f))
                else:
                    fp.write('   name= %s  wtf=%g' % (s, f))
                if problemType in ['ph', 'rhoe']: fp.write(' t=300')
                fp.write('\n')
        #
        if len(self.onlyList) > 0:
            fp.write('only %s\n' % (' '.join(self.onlyList)))
        #
        fp.write('output')
        if self.outputUnits == 'massf': fp.write(' massf')
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
        thermo_props_found = False
        conductivity_found = False
        incident_shock_data = False
        for line in lines:
            if line=="\n": continue
            if line.find("PRODUCTS WHICH WERE CONSIDERED BUT WHOSE")>=0: break
            if (line.find("THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED")>=0 or
                line.find("THERMODYNAMIC EQUILIBRIUM COMBUSTION PROPERTIES AT ASSIGNED")>=0):
                thermo_props_found = True
            elif line.find("SHOCKED GAS (2)--INCIDENT--EQUILIBRIUM")>=0:
                incident_shock_data = True
            elif thermo_props_found or incident_shock_data:
                tokens = line.split()
                # Fill out thermo properties
                if line.find("H, KJ/KG")>=0:
                    self.h = get_cea2_float(tokens[2:]) * 1.0e3
                elif line.find("U, KJ/KG")>=0:
                    self.e = get_cea2_float(tokens[2:]) * 1.0e3
                elif line.find("S, KJ/(KG)(K)")>=0:
                    self.s = get_cea2_float(tokens[2:]) * 1.0e3
                elif line.find("Cp, KJ/(KG)(K)")>=0:
                    self.cp = get_cea2_float(tokens[2:]) * 1.0e3
                    self.C_p = self.cp
                elif line.find("GAMMAs")>=0:
                    self.gam = get_cea2_float(tokens[1:])
                elif line.find("M, (1/n)")>=0:
                    self.Mmass = get_cea2_float(tokens[2:])
                elif line.find("SON VEL,M/SEC")>=0:
                    self.a = get_cea2_float(tokens[2:])
                elif line.find("P, BAR")>=0:
                    self.p = get_cea2_float(tokens[2:]) * 1.0e5
                    # print "p = ", self.p
                elif line.find("T, K")>=0:
                    self.T = get_cea2_float(tokens[2:])
                    # print "T = ", self.T
                elif line.find("RHO, KG/CU M")>=0:
                    self.rho = get_cea2_float(tokens[3:])
                    # print "rho = ", self.rho
                # Fill out transport properties if requested
                if transProps:
                    if line.find("VISC,MILLIPOISE")>=0:
                        self.mu = get_cea2_float(tokens[1:]) * 1.0e-4
                        # print "mu = ", self.mu
                    elif conductivity_found==False and line.find("CONDUCTIVITY")>=0 and len(tokens)==2:
                        self.k = get_cea2_float(tokens[1:]) * 1.0e-1
                        # print "k = ", self.k
                        # want to use the first conductivity value (for equilibrium reaction)
                        conductivity_found = True
                else:
                    self.mu = 0.0
                    self.k = 0.0
                # Get the shock specific parameters if appropriate
                if incident_shock_data:
                    if line.find("U2, M/SEC")>=0:
                        self.u2 = get_cea2_float(tokens[2:])
        # Calculate remaining thermo properties
        self.R = R_universal / self.Mmass  # gas constant, J/kg.K
        self.C_v = self.C_p - self.R       # specific heat, const volume
        # Check for small or zero pressure value printed by CEA2; 
        # it may have underflowed when printed in bars.
        if self.p < 1000.0:
            self.p = self.rho * self.R * self.T
        #
        # Scan lines again, this time looking for species fractions.
        species_fractions_found = False
        # Re-initialise the species list/fractions so that we ensure that there is
        # no 'left-over' information from last time
        self.species = {}    
        for line in lines:
            line = line.strip()
            if len(line) == 0: continue
            if line.find('MOLE FRACTIONS') >= 0:
                species_fractions_found = True
                continue
            if line.find('MASS FRACTIONS') >= 0:
                species_fractions_found = True
                continue
            if line.find('* THERMODYNAMIC PROPERTIES FITTED') >= 0: break
            if species_fractions_found:
                tokens = line.split()
                s = tokens[0].replace('*', '')
                self.species[s] = get_cea2_float(tokens[1:])
                # print "%s = %e" % (s, self.species[s])
        # Now check for any None values, where CEA2 wrote a dodgy format float.
        dodgyCount = 0
        sumFractions = 0.0
        for s in self.species.keys():
            if self.species[s] == None:
                dodgyCount += 1
                dodgySpecies = s
            else:
                sumFractions += self.species[s]
        if dodgyCount > 1:
            print "Cannot evaluate species fractions"
            print "because there are too many dodgy values"
            sys.exit(-1)
        # but we can recover one missing value.
        if dodgyCount == 1:
            self.species[dodgySpecies] = 1.0 - sumFractions
        return

    def EOS(self, problemType='pT', transProps=True):
        """
        Computes the gas state, taking into account the high-temperature effects.

        It does this by writing a suitable input file for the CEA code,
        calling that code and then extracting the relevant results from 
        the CEA output or plot file.

        :param self: the gas state to be filled in
        :param problemType: a string specifying the type of CEA analysis:
            'pT', 'rhoT', 'rhoe', 'ps', shock
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

def make_gas_from_name(gasName, outputUnits='massf'):
    """
    Manufacture a Gas object from a small library of options.

    :param gasName: one of the names for the special cases set out below
    :returns: as Gas object
    """
    if gasName.lower() == 'air':
        return Gas({'Air':1.0,}, outputUnits=outputUnits, trace=1.0e-4)
    elif gasName.lower() == 'air-ions':
        return Gas({'Air':1.0,}, outputUnits=outputUnits, trace=1.0e-4,
                   with_ions=True)
    elif gasName.lower() == 'air5species':
        return Gas(reactants={'N2':0.79, 'O2':0.21}, inputUnits='moles',
                   onlyList=['N2','O2','N','O','NO'],
                   outputUnits=outputUnits)
    elif gasName.lower() == 'air7species':
        return Gas(reactants={'N2':0.79, 'O2':0.21}, inputUnits='moles', 
                   onlyList=['N2','O2','N','O','NO','NO+','e-'],
                   outputUnits=outputUnits, with_ions=True)
    elif gasName.lower() == 'air11species':
        return Gas(reactants={'N2':0.79, 'O2':0.21}, inputUnits='moles', 
                   onlyList=['N2','O2','N','O','NO','N+','O+','N2+','O2+','NO+','e-'],
                   outputUnits=outputUnits, with_ions=True, trace=1.0e-30)
    elif gasName.lower() == 'air13species':
        return Gas(reactants={'N2':0.7811, 'O2':0.2095, 'Ar':0.0093}, inputUnits='moles', 
                   onlyList=['N2','O2','Ar','N','O','NO','Ar+','N+','O+','N2+','O2+','NO+','e-'],
                   outputUnits=outputUnits, with_ions=True, trace=1.0e-30)
    elif gasName.lower() == 'n2':
        return Gas(reactants={'N2':1.0, 'N':0.0}, onlyList=['N2','N'],
                   outputUnits=outputUnits)
    elif gasName.lower() == 'n2-ions':
        return Gas(reactants={'N2':1.0, 'N':0.0}, 
                   onlyList=['N2','N','N2+','N+','e-'],
                   outputUnits=outputUnits, with_ions=True)
    elif gasName.lower() == 'co2':
        return Gas(reactants={'CO2':1.0}, 
                   onlyList=['CO2','C2','C','CO','O2','O'],
                   outputUnits=outputUnits)
    elif gasName.lower() == 'co2-ions':
        return Gas(reactants={'CO2':1.0}, 
                   onlyList=['CO2','C2','C','CO','O2','O','C+','CO+','O2+','O+','e-'],
                   outputUnits=outputUnits, with_ions=True)
    elif gasName.lower() == 'mars-basic':
        return Gas(reactants={'CO2':0.97,'N2':0.03}, inputUnits='massf', 
                   onlyList=['C','C2','CN','CO','CO2','N','N2','NO','O','O2'],
                   outputUnits=outputUnits)
    elif gasName.lower() == 'mars-trace':
        return Gas(reactants={'CO2':0.9668,'N2':0.0174,'O2':0.0011,'Ar':0.0147}, inputUnits='massf', 
                   onlyList=['C','C2','CN','CO','CO2','N','N2','NO','O','O2','Ar'],
                   outputUnits=outputUnits)
    elif gasName.lower() == 'mars-trace-ions':
        return Gas(reactants={'CO2':0.9668,'N2':0.0174,'O2':0.0011,'Ar':0.0147}, inputUnits='massf', 
                   onlyList=['C','C2','CN','CO','CO2','N','N2','NO','O','O2','Ar',
                             'C+','CO+','NO+','O+','O2+','e-'],
                   outputUnits=outputUnits, with_ions=True)
    elif gasName.lower() == 'h2ne':
        return Gas(reactants={'H2':0.85, 'Ne':0.15}, inputUnits='moles',
                   onlyList=['H2','H','Ne'],
                   outputUnits=outputUnits)
    elif gasName.lower() == 'h2ne-ions':
        return Gas(reactants={'H2':0.85, 'Ne':0.15}, inputUnits='moles',
                   onlyList=['H2','H','Ne','H+','e-'],
                   outputUnits=outputUnits, with_ions=True)
    elif gasName.lower() == 'jupiter-like':
        return Gas(reactants={'H2':0.15, 'Ne':0.85}, inputUnits='moles',
                   onlyList=['H2','H','Ne','H+','e-'],
                   outputUnits=outputUnits, with_ions=True)
    elif gasName.lower() == 'titan-like':
        return Gas(reactants={'N2':0.95,'CH4':0.05}, inputUnits='moles',
                   onlyList=['N2','CH4','CH3','CH2','CH','C2','H2','CN','NH','HCN','N','C','H'],
                   outputUnits=outputUnits, with_ions=False)
    elif gasName.lower() == 'titan-like-ions':
        return Gas(reactants={'N2':0.95,'CH4':0.05}, inputUnits='moles',
                   onlyList=['N2','CH4','CH3','CH2','CH','C2','H2','CN','NH','HCN','N','C','H',
                             'N2+','CN+','N+','C+','H+','e-'],
                   outputUnits=outputUnits, with_ions=True)
    elif gasName.lower() == 'ar':
        return Gas(reactants={'Ar':1.0, 'Ar+':0.0, 'e_minus':0.0},
                   inputUnits='moles', outputUnits=outputUnits, 
                   with_ions=True, trace=1.0e-16)
    elif gasName.lower() == 'kr':
        return Gas(reactants={'Kr':1.0, 'Kr+':0.0, 'e_minus':0.0},
                   inputUnits='moles', outputUnits=outputUnits, 
                   with_ions=True, trace=1.0e-16)
    else:
        raise Exception, 'make_gas_from_name(): unknown gasName: %s' % gasName

def list_gas_names():
    """
    :returns: the list of gases available in make_gas_from_name()
    """
    return ['air', 'air-ions', 'air5species', 'air7species', 'air11species',
            'air13species', 'n2', 'n2-ions', 'co2', 'co2-ions', 'mars-trace', 'mars-basic',
            'h2ne', 'h2ne-ions', 'jupiter-like', 'titan-like', 'titan-like-ions', 'ar', 'kr']
    
def make_reactants_dictionary( species_list ):
    """
    Creates the CEA reactants dictionary from a list of species
    in the lib/gas format
    :param species_list: lib/gas species list
    """
    nsp = len(species_list)
    reactants = dict()
    for sp in species_list:
        # replace names containing '_plus' with '+' 
        sp = sp.replace("_plus","+")
        # replace names containing '_minus' with '-' 
        sp = sp.replace("_minus","-")
	reactants.setdefault(sp,0.0)
    return reactants

def get_species_composition( sp, species_data ):
    """
    Creates a list of mass or mole fractions for a species
    in lib/gas form from the CEA species_data dictionary
    :param sp: a single lib/gas species
    :param species_data: the CEA species_data dictionary
    """
    # replace names containing '_plus' with '+' 
    if ( sp.find("_plus")>=0 ): sp = sp[0:sp.find("_plus")] + "+"
    # replace names containing '_minus' with '-' 
    if ( sp.find("_minus")>=0 ): sp = sp[0:sp.find("_minus")] + "-"
    if sp in species_data.keys():
	    return species_data[sp]
    else:
	    return 0.0
	    
def get_with_ions_flag( species_list ):
    """
    Determines the 'with_ions' flag from a list of species
    in the lib/gas format
    :param species_list: lib/gas species list
    """
    for sp in species_list:
        if sp.find("_plus")>=0: return True
        if sp.find("_minus")>=0: return True
    return False

# --------------------------------------------------------------

if __name__ == '__main__':
    print 'Test/demonstrate the Gas class...'
    #
    print '\nDefault constructor with Air as the test gas.'
    a = Gas({'Air':1.0,}, outputUnits='moles')
    a.set_pT(100.0e3, 300.0)
    a.write_state(sys.stdout)
    print 'and the same Air at a higher temperature'
    a.set_pT(100.0e3, 4000.0)
    a.write_state(sys.stdout)
    #
    print '\nCheck enthalpy specification'
    b = make_gas_from_name('air', outputUnits='moles')
    b.set_ph(a.p, a.h)
    b.write_state(sys.stdout)
    #
    print '\nCheck internal-energy specification'
    b = make_gas_from_name('air', outputUnits='moles')
    b.set_rhoe(a.rho, a.e)
    b.write_state(sys.stdout)
    #
    print '\nAir-5-species for nenzfr: 79% N2, 21% O2 by mole fraction.'
    a = Gas(reactants={'N2':0.79, 'O2':0.21, 'N':0.0, 'O':0.0, 'NO':0.0}, 
            inputUnits='moles', outputUnits='massf',
            onlyList=['N2','O2','N','O','NO'])
    a.set_pT(100.0e3, 300.0)
    a.write_state(sys.stdout)
    print 'and isentropically compress to a higher pressure'
    a.set_ps(10.0e6, a.s)
    a.write_state(sys.stdout)
    #
    print '\nTry an odd mix of Helium, N2 and N'
    b = Gas({'N2':1.0, 'N':0.0, 'He':0.0})
    b.set_pT(100.0e3, 300.0)
    b.write_state(sys.stdout)
    print 'and the same initial mix and volume at a higher temperature'
    b.set_rhoT(b.rho, 5000.0)
    b.write_state(sys.stdout)
    #
    print '\nStart again with low-T air as the test gas'
    a = Gas({'Air':1.0,}); a.set_pT(100.0e3, 300.0)
    a.write_state(sys.stdout)
    print 'clone it, changing species-fraction units'
    c = a.clone(newOutputUnits='moles')
    c.write_state(sys.stdout)
    print 'and shock process it'
    c.shock_process(4000.0)
    c.write_state(sys.stdout)
    #
    print 'End of test.'


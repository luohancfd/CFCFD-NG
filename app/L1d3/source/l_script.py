#! /usr/bin/env python
"""
Python program to write the input parameter file for L1d3.

It is intended for the user to define their particular facility and
flow in terms of the data objects defined in this module.  As part of
its initialization, this program will execute a user-specified job
file that contains, in Python, the user's script that defines both
facility geometry and gas-path details.

Usage::

    $ l_script.py -f <job>

The simulation control data is then organised via the classes:
GlobalData, GasSlug, Piston and Diaphragm.  These classes
provide places to store the configuration information and their
function/method names appear as commands in the user's
job description file.

When setting up a new simulation, first define the tube as a set
of (x,d) break-points and identify regions of head-loss and
regions where the wall-temperature varies from the nominal value.
Create the GasSlugs, Pistons, and Diaphragms that will make up the
gas path.  Note that places where two GasSlugs join will need a
GasInterface to be defined.  Once all of the components have been
created, assemble the gas path and then set any of the time-stepping
parameters for which you want values other than the default.

Here is an example script for the Sod shock-tube problem::

    # sod.py
    gdata.title = 'Sods ideal shock tube, 06-Jul-05'
    select_gas_model(model='ideal gas', species=['air'])
    
    Define the tube walls.
    add_break_point(0.0, 0.01)
    add_break_point(3.0, 0.01)
 
    # Create the gas-path.
    left_wall = VelocityEnd(x0=0.0, v=0.0)
    driver_gas = GasSlug(p=100.0e3, u=0.0, T=348.4, nn=100,
                         to_end_R=1, cluster_strength=1.1,
                         hcells=0)
    interface = GasInterface(x0=0.5)
    driven_gas = GasSlug(p=10.0e3, u=0.0, T=278.7, nn=100,
                         hcells=0)
    right_wall = VelocityEnd(x0=1.0, v=0.0)
    assemble_gas_path(left_wall, driver_gas, interface, driven_gas, right_wall)
    
    # Set some time-stepping parameters
    gdata.dt_init = 1.0e-7
    gdata.max_time = 0.6e-3
    gdata.max_step = 5000
    add_dt_plot(0.0, 10.0e-6, 5.0e-6)
    add_history_loc(0.7)

This script should define the gas path::

  .       |+----- driver-gas -----+|+----- driven-gas -----+|
  .       |                        |                        |
  .       |                        |                        |
  .   left-wall                interface               right-wall
    
and can be invoked with the command::

    $ l_script.py -f sod

Upon getting to the end of the user's script, this program should then
write a complete simulation parameter file (sod.Lp) in the INI format.  
Because this program (l_script) just gathers the data
in order to write the input parameter file, the old documentation for that
file is still (somewhat) relevant despite a few small name changes.

Note that Python is very picky about whitespace.  If you cut and paste the
example from above, make sure that the lines start in the first column and
that indentation is consistent with Python's syntax rules.

Globally-defined object
-----------------------

* gdata: Contains the GlobalData information describing the simulation.
  Note that there is one such variable set up by the main program and
  the user's script should directly set the attributes of this variable
  to adjust settings for the simulation.

.. Author: P.Jacobs

.. Versions: 
   June 2005
   24-Jul-2006 Ported to Rowan's new C++ chemistry.
   Mar-May 2010 More reworking for L1d3 and 
                the latest-greatest thermochemistry
   Sep-Oct-2012 Much cleaning up and Sphinx docs.
"""

# ----------------------------------------------------------------------
#
import sys
import os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from getopt import getopt
from gaspy import *

shortOptions = "hf:"
longOptions = ["help", "job="]

def printUsage():
    print ""
    print "Usage: scriptit.py" + \
          " [--help | -h]" + \
          " [--job=<jobName> | -f <jobName>]"
    print ""
    return

#----------------------------------------------------------------------
# This is where we store the core data for the simulation.

class GlobalData(object):
    """
    Contains the tube and simulation control parameters.

    The user's script should not create one of these
    but should specify the simulation parameters by
    altering the attributes of the global object "gdata"
    that already exists by the time the user's script executes.
    
    The following attributes are available:

    * case_id: (int) Specifies a special case that has custom C-code
      in the main simulation. See l1d.h for possible values.
      Use a value of 0 for a generic simulation; this is the usual case.
      This is not much used in the current code.

    * title: Short title string for embedding in the parameter and solution files.
    
    * gas_name: the name of the thermo-chemical model
      (which, in turn, determines the number of species).

    * reacting_flag: If set to 1, Rowan's finite-rate chemistry will
      be active.  (Default is 0)

    * dt_init: (float) The size of the time-step that will be used for the
      first few simulation steps.
      After a few steps, the cfl condition takes over the determination
      of a suitable time-step.

    * max_time: (float) The simulation will stop if it reaches this time.
      It is most usual to use this critereon to stop the simulation.

    * max_step: The simulation will be stopped if it reaches
      this number of steps.
      This is mostly used to catch the problem of the calculation taking
      a very long time (measured by one's patience), possibly because
      the time-step size has decreased to an extremely small value.

    * cfl: (float) Largest allowable CFL number.
      The time step is adjusted to ensure that this value is not exceeded
      in any particular cell.
      A typical value of 0.25 seems to work well for simulations with
      sudden events such as diaphragm bursting, while a value as high as
      0.5 should be considered only for well-behaved flows.

    * t_order: (int) 
      1=Euler time-stepping. This is generally cheap and nasty.
      2=predictor-corrector time-stepping, nominally second order.
      This is the default setting.
      It is, however, twice as CPU intensive as Euler time-stepping.

    * x_order: (int) 
      1=use cell averages without high-order reconstruction.
      Use this only if the second-order calculation is showing problems.
      2=use limited reconstruction (nominally second order).
      This is the default selection. 

    * dt_plot_list: (list of tuples) 
      Specifies the frequency of writing complete solutions
      (for later plotting, maybe) and also for the writing of data at
      history locations.
      It may be convenient to have different frequencies of writing such
      output at different stages of the simulation.
      For example, free-piston driven shock tunnels have a fairly long
      period during which the piston travels the length of the compression
      tube and then a relatively short period, following diaphragm rupture,
      when all the interesting things happen.
      It is good to have low-frequency output during most of the compression
      process and higher-frequency output starting just before diaphragm
      rupture.
      Arranging good values may require some trial and error.
      Add entries to this list via the add_dt_plot function.

    * xd_list: List of break-point tuples defining the tube wall.
      Add elements to the list via the function add_break_point.

    * n: (int) The number of small segments that will be used to describe
      the tube's area distribution internal to the simulation.
      To enable a fast lookup process for the area calculation,
      the area variation between equally-spaced x-positions is taken
      to be linear.
      The default value is 4000 and probably won't need to be changed
      except for geometries with rapidly changing cross-sections.

    * T_nominal: (float) The nominal wall temperature (in degrees K)
      in the absence of a patch of differing temperature.

    * T_patch_list: (list of tuples)
      Regions of the tube wall that have temperature different to the 
      nominal value can be specified via the function add_T_patch.

    * loss_region_list: (list of tuples)
      List of head-loss regions, usually associated
      with sudden changes in tube cross-section and diaphragm stations.
      Add regions via the function add_loss_region.

    * hloc_list: (list of floats)
      List of x-coordinates for the history locations.
      Add entries via the function add_history_loc.
    """
    count = 0

    # We want to prevent the user's script from introducing new attributes
    # via typographical errors.
    __slots__ = 'param_file', 'title', 'case_id', 'gas_model_file', 'gmodel', \
                'reaction_scheme_file', 'reacting_flag', \
                'dt_init', 'cfl', 'dt_plot_list', \
                'max_time', 'max_step', \
                'x_order', 't_order', 'thermal_damping', \
                'T_nominal', 'T_patch_list', 'loss_region_list', \
                'xd_list', 'n', 'hloc_list'
    
    def __init__(self):
        """Accepts user-specified data and sets defaults. Make one only."""
        if GlobalData.count >= 1:
            raise Exception, "Already have a GlobalData object defined."
        
        self.param_file = "simulation.Lp"
        self.title = "Another L1d3 Simulation."
        self.case_id = 0
        self.gas_model_file = "gas-model.lua"
        self.gmodel = None
        self.reaction_scheme_file = "None"
        self.reacting_flag = 0
        self.dt_init = 1.0e-6
        self.cfl = 0.5
        # If dt_plot_list is still an empty list when we write
        # the parameter file, just use the max_time value for
        # both dt_plot and dt_his.  Fill in later.
        self.dt_plot_list = []
        self.max_time = 1.0e-3
        self.max_step = 10
        self.x_order = 2
        self.t_order = 2
        self.thermal_damping = 0
        # Tube definition is a list of (x,diameter) tuples
        # defining the break-points of the tube wall.
        # The transition_flag with a value of 1
        # indicates linear transitions from break-point i to point i+1.
        # The alternative is a cubic transition (I think).
        self.xd_list = []
        self.n = 4000
        # Wall temperature is specified as a nominal value with
        # patches of other temperatures.
        self.T_nominal = 300.0  # Nominal tube-wall temperature in Kelvin
        self.T_patch_list = []
        # Head-losses are also spread over finite-length patches.
        self.loss_region_list = []
        # History locations are collected in a list.
        self.hloc_list = []
        #
        GlobalData.count += 1
        return

    def write_to_ini_file(self, fp, nslug, npiston, ndiaphragm):
        """
        Writes the configuration data to the specified file in .ini format.
        """
        fp.write("[global_data]\n")
        fp.write("    title = %s\n" % self.title)
        fp.write("    case_id = %d\n" % self.case_id)
        fp.write("    gas_model_file = %s\n" % self.gas_model_file)
        fp.write("    reaction_scheme_file = %s\n" % self.reaction_scheme_file)
        fp.write("    reacting_flag = %d\n" % self.reacting_flag)
        fp.write("    max_time = %e\n" % self.max_time)
        fp.write("    max_step = %d\n" % self.max_step)
        fp.write("    dt_init = %e\n" % self.dt_init)
        fp.write("    cfl = %e\n" % self.cfl)
        fp.write("    x_order = %d\n" % self.x_order)
        fp.write("    t_order = %d\n" % self.t_order)
        fp.write("    thermal_damping = %e\n" % self.thermal_damping)
        #
        if len(gdata.dt_plot_list) == 0:
            # Since the user did not specify any, default to the end.
            self.add_dt_plot(0.0, gdata.max_time, gdata.max_time)
        n_dt_plot = len(self.dt_plot_list)
        fp.write("    n_dt_plot = %d\n" % n_dt_plot)
        fp.write("    t_change =");
        for i in range(n_dt_plot):
            fp.write(" %e" % gdata.dt_plot_list[i][0])
        fp.write("\n")
        fp.write("    dt_plot =");
        for i in range(n_dt_plot):
            fp.write(" %e" % gdata.dt_plot_list[i][1])
        fp.write("\n")
        fp.write("    dt_his =");
        for i in range(n_dt_plot):
            fp.write(" %e" % gdata.dt_plot_list[i][2])
        fp.write("\n")
        #
        n_hloc = len(gdata.hloc_list)
        fp.write("    hloc_n = %d\n" % n_hloc)
        fp.write("    hloc_x =")
        for i in range(n_hloc):
            fp.write(" %e" % gdata.hloc_list[i])
        fp.write("\n")
        #
        fp.write("    tube_n = %d\n" % gdata.n)
        nseg = len(gdata.xd_list) - 1
        fp.write("    tube_nseg = %d\n" % nseg)
        fp.write("    tube_xb =")
        for i in range(nseg+1):
            fp.write(" %e" % gdata.xd_list[i][0])
        fp.write("\n")
        fp.write("    tube_d =")
        for i in range(nseg+1):
            fp.write(" %e" % gdata.xd_list[i][1])
        fp.write("\n")
        fp.write("    tube_linear =")
        for i in range(nseg+1):
            fp.write(" %d" % gdata.xd_list[i][2])
        fp.write("\n")
        #
        nKL = len(gdata.loss_region_list)
        fp.write("    KL_n = %d\n" % nKL)
        fp.write("    KL_xL =")
        for i in range(nKL):
            fp.write(" %e" % gdata.loss_region_list[i][0])
        fp.write("\n")
        fp.write("    KL_xR =")
        for i in range(nKL):
            fp.write(" %e" % gdata.loss_region_list[i][1])
        fp.write("\n")
        fp.write("    KL_K =")
        for i in range(nKL):
            fp.write(" %e" % gdata.loss_region_list[i][2])
        fp.write("\n")
        #
        fp.write("    T_nominal = %e\n" % gdata.T_nominal)
        nT = len(gdata.T_patch_list)
        fp.write("    Tpatch_n = %d\n" % nT)
        fp.write("    Tpatch_xL =")
        for i in range(nT):
            fp.write(" %e" % gdata.T_patch_list[i][0])
        fp.write("\n")
        fp.write("    Tpatch_xR =")
        for i in range(nT):
            fp.write(" %e" % gdata.T_patch_list[i][1])
        fp.write("\n")
        fp.write("    Tpatch_T =")
        for i in range(nT):
            fp.write(" %e" % gdata.T_patch_list[i][2])
        fp.write("\n")
        #
        fp.write("    nslug = %d\n" % nslug)
        fp.write("    npiston = %d\n" % npiston)
        fp.write("    ndiaphragm = %d\n" % ndiaphragm)
        return
    
# We will create just one GlobalData object that the user can alter.
gdata = GlobalData()

#----------------------------------------------------------------------
# copied the following couple of functions from e3prep.py 
# (and modified them slightly)

def select_gas_model(model=None, species=None, fname=None):
    """
    Selects a gas model for the simulation.

    :param model: (string) name of the gas model as shown in the list below.
    :param species: list of species names (strings).
    :param fname: (string) name of the gas-model file.

    :returns: a list of species names

    The new gas models are configured by stand-alone files.
    This function initialises a gas model for present use
    in the preparation program (ie. sets it in kernel code)
    and stores the name for later user at simulation time.

    If you already have a gas-model.lua file already set up,
    give its name as fname.
    
    If you have not already set up the gas-model.lua file,
    this function is provided as a simple but limited means to do so.

    Look-up-table (LUT) gas models and LUT_plus_composite cannot be
    created directly with this function.
    If you want a single LUT gas model, just set up the LUT table 
    file externally and supply the name of that file as fname.
    If you want a LUT-plus-composite gas model, set up the LUT table
    externally and then set up the rest of the composite gas model
    using create_gas_file() directly, then select the gas model
    by specifying the gas file name when calling this function.
    The create_gas_file() function has the capability of prepending
    the LUT gas species to the composite gas species list.
    """
    if fname is None:
        # Help the user to set up the gas-model file.
        fname = "gas-model.lua"
        if model == None:
            print "select_gas_model():"
            print "    A gas 'model' or 'fname' must be specified."
            print "    Bailing out!"
            sys.exit(1)
        if species == None:
            print "select_gas_model():"
            print "    When setting up a gas model, a list of species must be specified."
            print "    Bailing out!"
            sys.exit(1)
        create_gas_file(model, species, fname)
    # At this point, the gas model file exists as required.
    gdata.gmodel = create_gas_model(fname)
    gdata.gas_model_file = fname
    nsp = gdata.gmodel.get_number_of_species()
    return [ gdata.gmodel.species_name(isp) for isp in range(nsp) ]

def set_reaction_scheme(fname, reacting_flag=1):
    """
    Sets the reaction update model and specifies a reacting simulation.

    This function sets the name of the input file for the reaction update.
    It also sets the reacting flag.  Note that the user may later reset
    that flag, possibly to turn off reactions for a particular simulation.

    :param fname: (string) name of the file containing the reactions scheme
    :param reacting_flag: (int) =1 to activate the reaction scheme;
       =0 frozen flow.
    """
    gdata.reacting_flag = reacting_flag
    gdata.reaction_scheme_file = fname
    return

# We want to keep the old name, for a while.
set_reaction_update = set_reaction_scheme

# --------------------------------------------------------------------
# The following functions are to provide convenient ways of setting
# some of the GlobalData elements.

def add_dt_plot(t_change, dt_plot, dt_his):
    """
    Add a dt tuple to the dt_plot tuple list in GlobalData.

    :param t_change: (float) The time, in seconds, 
        at which this dt_plot and dt_his should take effect.
    :param dt_plot: (float) Time interval between writing whole solutions
        for later plotting.
    :param dt_his: (float) Time interval between writing data to history file.
    :returns: the new length of the list
    """
    global gdata
    if len(gdata.dt_plot_list) > 0:
        # Check that we are adding points monotonically in x.
        if t_change > gdata.dt_plot_list[-1][0]:
            gdata.dt_plot_list.append((t_change, dt_plot, dt_his))
        else:
            print "Warning: did not add dt_plot tuple (", \
                  t_change, dt_plot, dt_his, ")."
    else:
        gdata.dt_plot_list.append((t_change, dt_plot, dt_his))
    return len(gdata.dt_plot_list)


def add_break_point(x, d, transition_flag=0):
    """
    Add a break-point tuple to the tube-diameter description
    contained in GlobalData.

    The tube is described as a set of (x,d)-coordinate pairs that
    define break points in the profile of the tube wall.

    :param x: (float) x-coordinate, in metres, of the break point
    :param d: (float) diameter, in metres, of the tube wall at the break-point.
    :param transition_flag: (int) Indicates the variation in diameter between
       this break-point and the next. 1=linear, 0=Hermite-cubic.
    :returns: Number of break points defined so far.
    """
    global gdata
    if len(gdata.xd_list) > 0:
        # Check that we are adding points monotonically in x.
        if x > gdata.xd_list[-1][0]:
            gdata.xd_list.append((x, d, transition_flag))
        else:
            print "Warning: did not add new break-point (", x, d, ")."
    else:
        gdata.xd_list.append((x, d, transition_flag))
    return len(gdata.xd_list)


def add_loss_region(xL, xR, K):
    """
    Add a head-loss region to the tube description in L{GlobalData}.

    There is a momentum-sink term much like the so-called minor-loss terms
    in the fluid mechanics text books.
    The effect of the loss is spread over a finite region so that the cells
    are gradually affected as they pass through the region

    :param xL: (float) Left-end location, in metres, of the loss region.
    :param xR: (float) Right-end location, in metres, of the loss region.
    :param K: (float) Head-loss coefficient.  A value of 0.25 seems to be good for a
        reasonably smooth contraction such as the T4 main diaphragm station.
    :returns: Number of loss regions defined so far.
    """
    if xR < xL:
        # Keep x-values in increasing order
        xL, xR = xR, xL
    if abs(xR - xL) < 1.0e-3:
        print "Warning: loss region is very short: (", xL, xR, ")"
    gdata.loss_region_list.append((xL, xR, K))
    return len(gdata.loss_region_list)


def add_T_patch(xL, xR, T):
    """
    Add a temperature patch for a region where the wall temperature
    is different from the nominal value.

    :param xL: (float) Left-end location, in metres, of the loss region.
    :param xR: (float) Right-end location, in metres, of the loss region.
    :param T: (float) Wall temperature in degrees K.
    :returns: Number of temperature patches defined so far.
    """
    if xR < xL:
        # Keep x-values in increasing order
        xL, xR = xR, xL
    if abs(xR - xL) < 1.0e-3:
        print "Warning: temperature patch is very short: (", xL, xR, ")"
    gdata.T_patch_list.append((xL, xR, T))
    return len(gdata.T_patch_list)


def add_history_loc(x):
    """
    Add a location to the history-location list in L{GlobalData}.

    :param x: (float) x-coordinate, in metres, of the sample point.
    :returns: Number of sample points defined so far.
    """
    if isinstance(x, list):
        gdata.hloc_list.extend(x)
    else:
        gdata.hloc_list.append(float(x))
    return len(gdata.hloc_list)


#----------------------------------------------------------------------
# The following classes define the objects that will be assembled into
# the gas path.
# ---------------------------------------------------------------------

class GasSlug(object):
    """
    Contains the gas properties and discretisation for each gas slug.

    The user may create more than one gas slug to describe the initial
    gas properties throughout the facility.
    
    Note that a slug needs to have appropriate end-conditions.
    This is achieved by creating end-condition objects such as
    FreeEnd and VelocityEnd objects and then assembling 
    the gas-path via a call to the function assemble_gas_path.
    """
  
    # We will accumulate references to defined objects.
    slugList = []

    __slots__ = 'gas', 'u', 'indx', 'label', 'xL', 'xR', \
                'bcL', 'bcR', 'bcL_which_end', 'bcR_which_end', \
                'nn', 'to_end_L', 'to_end_R', 'cluster_strength', \
                'viscous_effects', 'adiabatic_flag', 'hcells', \
                'nnmax', 'adaptive', 'dxmin', 'dxmax'
        
    def __init__(self,
                 p = 100.0e3,
                 u = 0.0,
                 T = [300.0,],
                 massf = [1.0,],
                 label="",
                 nn = 10,
                 to_end_L=0,
                 to_end_R=0,
                 cluster_strength=0.0,
                 viscous_effects=0, # several options, see Lp-file documentation
                 adiabatic_flag=0,
                 hcells=[],
                 nnmax=None, # parameters for adaptive discretization
                 adaptive=0,
                 dxmin=0.01,
                 dxmax=0.05
                 ):
        """
        Creates a gas slug with user-specified properties.

        Most parameters have default properties so that only the user
        needs to override the ones that they wish to set differently.

        Note that the locations of the ends of the slug are communicated
        through end-condition objects that are attached during assembly
        of the gas path.
        
        :param p: (float) Pressure in Pa.
        :param u: (float) Velocity in m/s.
        :param T: (float or list of floats) Temperature in degrees K.
        :param massf: Mass fractions supplied as a list of floats 
            or a dictionary of species names and floats. 
            The number of mass fraction values should match the number 
            of species expected by the selected gas model.
        :param label: Optional (string) label for the gas slug.
        :param nn: (int) Number of cells within the gas slug.
        :param to_end_L: (int) Flag to indicate that cells should 
            be clustered to the left end.
        :param to_end_R: (int) Flag to indicate that cells should
            be clustered to the right end.
        :param cluster_strength: (float) As this value approaches 1.0 from above,
            the clustering gets stronger.
            A value of zero indicates no clustering.
        :param viscous_effects: (int) A nonzero value activates the viscous effects.
            0 = inviscid equations only;
            1 = include viscous source terms F_wall, loss, q,
            friction factor for pipe flow;
            2 = use Con Doolan's laminar mass-loss model if the mass within
            a cell is greater than MINIMUM_MASS as set in l1d.hh;
            3 = use Con Doolan's turbulent mass-loss model if the mass within
            a cell is greater than MINIMUM_MASS as set in l1d.hh;
            4 = include viscous source terms F_wall, loss, q,
            friction factor for flat plate;
            5 = use David Buttsworth's mass-loss model with
            pipe-flow friction factor;
            6 = use David Buttsworth's mass-loss model with
            flat-plate friction factor;
            7 = include viscous source terms F_wall, loss, q,
            friction factor for pipe flow; half heat flux.
        :param adiabatic_flag: (int) Flag to indicate that there should
            be no heat transfer at the tube wall.
        :param hcells: Either the index (int) of a single cell or 
            a list of indices of cells for which the data are 
            to be written every dt_his seconds, as set by add_dt_plot.
            Note that cells are indexed from 0 to nn-1.
        :param nnmax: (int) Maximum number of cells that can be carried in
            the adaptive discretization.
        :param adaptive: (int) Flag to indicate that adaptive discretization
            is to be used (=1).
        :param dxmin: (float) Minimum cell size, below which a cell will be
            fused into its neighbour.
        :param dxmax: (float) Maximum cell size, above which a cell will be
            split into 2.
        """
        self.indx = len(GasSlug.slugList) # next available index
        if gdata.gmodel == None:
            print "ERROR: The gas model has not yet been set."
            print "       A gas model should be set by calling"
            print "       set_type_of_gas(...) before declaring a GasSlug."
            print "       If unsure, read the documentation."
            sys.exit(-1)
        # Gas data related values
        self.gas = Gas_data(gdata.gmodel)
        nsp = gdata.gmodel.get_number_of_species()
        nmodes = gdata.gmodel.get_number_of_modes()
        self.gas.p = p
        if type(massf) is dict:
            set_massf(self.gas, gdata.gmodel, massf)
        else:
            # Assume that it is a full list.
            sum_mf = sum(massf)
            if abs(sum_mf - 1.0) > 1.0e-5:
                print "Warning: mass fractions sum to %g" % (sum_mf,)
            for isp in range(nsp): self.gas.massf[isp] = massf[isp]
        if type(T) is float:
            T = [T,] * nmodes
        if type(T) is int:
            T = [float(T),] * nmodes
        for imode in range(nmodes): self.gas.T[imode] = T[imode]
        gdata.gmodel.eval_thermo_state_pT(self.gas)
        gdata.gmodel.eval_transport_coefficients(self.gas, gdata.gmodel)
        self.u = u
        self.label = label
        #
        self.nn = nn
        self.to_end_L = to_end_L
        self.to_end_R = to_end_R
        self.cluster_strength = cluster_strength
        #
        # The adaptive functions have not been working well,
        # so we'll usually set some dummy values.
        if nnmax == None:
            self.nnmax = nn
        else:
            self.nnmax = nnmax
        self.adaptive = adaptive
        self.dxmin = dxmin
        self.dxmax = dxmax
        #
        self.viscous_effects = viscous_effects
        self.adiabatic_flag = adiabatic_flag
        if isinstance(hcells,int):
            self.hcells=[hcells,]
        elif isinstance(hcells,list):
            self.hcells=hcells
        else:
            print "Warning: hcells reset to empty list."
            hcells = []
        #
        # Boundary object at each end of the slug will be
        # attached later when the gas-path is assembled.
        self.bcL = None
        self.bcR = None
        # The spatial limits of the gas slug are determined from
        # the boundary-condition object.
        self.xL = None
        self.xR = None
        # We also want to know which end of the other object we
        # are attached to.
        self.bcL_which_end = 'R'
        self.bcR_which_end = 'L'
        #
        GasSlug.slugList.append(self)
        return

    def write_to_ini_file(self, fp):
        """
        Writes the flow state information to the specified file.
        """
        fp.write("[slug-%d]\n" % self.indx)
        fp.write("    label = %s\n" % self.label)
        fp.write("    nn = %d\n" % self.nn)
        fp.write("    cluster_to_end_L = %d\n" % self.to_end_L)
        fp.write("    cluster_to_end_R = %d\n" % self.to_end_R)
        fp.write("    cluster_strength = %e\n" % self.cluster_strength)
        fp.write("    nnmax = %d\n" % self.nnmax)
        fp.write("    adaptive = %d\n" % self.adaptive)
        fp.write("    dxmin = %e\n" % self.dxmin)
        fp.write("    dxmax = %e\n" % self.dxmax)
        fp.write("    viscous_effects = %d\n" % self.viscous_effects)
        fp.write("    adiabatic_flag = %d\n" % self.adiabatic_flag)
        #
        fp.write("    BC_L = %s\n" % boundary_control_string(self.bcL, self.bcL_which_end))
        fp.write("    BC_R = %s\n" % boundary_control_string(self.bcR, self.bcR_which_end))
        #
        hncell = len(self.hcells)
        fp.write("    hncell = %d\n" % hncell)
        fp.write("    hxcell =")
        for i in range(hncell):
            fp.write(" %d" % self.hcells[i])
        fp.write("\n")
        #
        fp.write("    initial_xL = %e\n" % self.xL)
        fp.write("    initial_xR = %e\n" % self.xR)
        fp.write("    initial_p = %e\n" % self.gas.p)
        fp.write("    initial_u = %e\n" % self.u)
        fp.write("    initial_T = %e\n" % self.gas.T[0])
        nsp = gdata.gmodel.get_number_of_species()
        fp.write("    massf =")
        for i in range(nsp):
            fp.write(" %e" % (self.gas.massf[i]))
        fp.write("\n")
        return

def boundary_control_string(other_object, other_object_which_end):
    """
    Assembles a boundary-condition control string for the supplied object.

    Helper function for the GasSlug class.
    """
    if isinstance(other_object, FreeEnd):
        bcs = "F  free-end"
    elif isinstance(other_object, VelocityEnd):
        bcs = "V %e  specified-velocity-end: velocity" % other_object.v
    elif isinstance(other_object, Piston):
        bcs = "P %d  piston: piston-id" % other_object.indx
    elif isinstance(other_object, GasSlug):
        bcs = "S %d %s  slug: slug-id, slug-end-id" % \
              (other_object.indx, other_object_which_end)
    elif isinstance(other_object, Diaphragm):
        # We need to get the details of the slug attached to
        # the other side of the diaphragm.
        if other_object_which_end == 'L':
            slug_id = other_object.slugR.indx
            slug_end_id = other_object.slugR_which_end
        elif other_object_which_end == 'R':
            slug_id = other_object.slugL.indx
            slug_end_id = other_object.slugL_which_end
        else:
            raise Exception, "boundary_control_string() is confused"
        bcs = "SD %d %s %d  diaphragm+slug: slug-id, slug-end-id, diaphragm-id" % \
              (slug_id, slug_end_id, other_object.indx)
    return bcs

#----------------------------------------------------------------------------

class Piston(object):
    """
    Contains the information for a piston.

    * The left- and right-end positions of the piston are
      also used to locate the ends of adjoining GasSlugs.
    * The basic piston model has inertia but no friction.
      To make accurate simulations of a particular facility,
      it is usually important to have some account of
      the friction caused by gas-seals and guide-rings that
      may be present on the piston.
    * The f_decay parameter is used to model secondary diaphragms 
      in expansion tubes as pistons which lose their mass over time.
    """

    __slots__ = 'indx', 'label', \
                'm', 'd', 'L', 'xL0', 'xR0', 'x0', 'v0', \
                'front_seal_f', 'front_seal_area', \
                'back_seal_f', 'back_seal_area', \
                'p_restrain', 'is_restrain', 'with_brakes', 'brakes_on', \
                'x_buffer', 'hit_buffer', \
                'slugL', 'slugL_which_end', \
                'slugR', 'slugR_which_end', \
                'f_decay', 'mass_limit'
    pistonList = []
    
    def __init__(self, m, d, xL0, xR0, v0,
                 front_seal_f=0.0, front_seal_area=0.0,
                 back_seal_f=0.0, back_seal_area=0.0,
                 p_restrain=0.0, is_restrain=0,    
                 with_brakes=0, brakes_on=0,      
                 x_buffer=10.e6, hit_buffer = 0, label="", f_decay=0.0,
                 mass_limit = 0.0
                 ):
        """
        Create a piston with specified properties.
            
        :param m: (float) Mass in kg.
        :param d: (float) Face diameter, metres.
        :param xL0: (float) Initial position of left-end, metres.
            The initial position of the piston centroid is set
            midway between xL0 and xR0 while piston length is the
            difference (xR0 - xL0).
        :param xR0: (float) Initial position of right-end, metres.
        :param v0: (float) Initial velocity (of the centroid), m/s.
        :param front_seal_f: (float) friction coefficient.
            Typical value might be 0.2.
        :param front_seal_area: (float) Seal area over which the front-side 
            pressure acts.
            This is the effective area over which the compressed gas pressed the 
            front-side seal against the tube wall.
            Friction force is this area multiplied by downstream-pressure by
            friction coefficient.
        :param back_seal_f: (float) friction coefficient. 
            A typical value might be 0.2.
        :param back_seal_area: (float) Seal area over which the back-side 
            pressure acts.
            Friction force is this area multiplied by downstream-pressure by
            friction coefficient.  This is for gun tunnel pistons that have
            flexible skirts that are pressed onto the tube wall by the pushing gas.
        :param p_restrain: (float) Pressure at which restraint will release.
            Some machines, such as two-stage light-gas guns, will
            hold the projectile in place with some form of mechanical
            restraint until the pressure behind the piston reaches
            a critical value.  The piston is then allowed to slide.
        :param is_restrain: (int) Status flag for restraint.
            0=free-to-move, 1=restrained, 2=predefined trajectory read from external file
        :param with_brakes: (int) Flag to indicate the presence of brakes.
            0=no-brakes, 1=piston-does-have-brakes.
            Such brakes, as on the T4 shock tunnel, allow forward
            motion of the piston but prevent backward motion by
            locking the piston against the tube wall.
        :param brakes_on: (int) Flag to indicate the state of the brakes.
            0=off, 1=on.
        :param x_buffer: (float) Position of the stopping buffer in metres.
            This is the location of the piston centroid at which the piston
            would strike the buffer (or brake, in HEG terminology).
            Note that it is different to the location of the front of
            the piston at strike.
        :param hit_buffer: (int) Flag to indicate state of buffer interaction.
            A value of 0 indicates that the piston has not (yet) hit the
            buffer.
            A value of 1 indicates that it has.
            Details of the time and velocity of the strike are recorded in
            the event file.
        :param label: (string) A bit of text for corresponding line in the Lp file.
        :param f_decay: (float) dm/dt = m * f_decay, thus a pseudo- time-constant
            for diaphragm mass decay. 
        :param mass_limit: (float) Mass limit for decaying diaphragm
        """
        self.indx = len(Piston.pistonList) # next available index
        if len(label) > 0:
            self.label = label
        else:
            # Construct a simple label.
            self.label = 'piston-' + str(self.indx)
        self.m = m
        self.d = d
        if xR0 < xL0:
            # We would like the x-values to be increasing to the right
            # but we really don't care if the piston length is zero.
            xL0, xR0 = xR0, xL0
        self.xL0 = xL0
        self.xR0 = xR0
        self.L = xR0 - xL0
        self.x0 = 0.5*(xL0 + xR0)
        self.v0 = v0
        self.front_seal_f = front_seal_f
        self.front_seal_area = front_seal_area
        self.back_seal_f = back_seal_f
        self.back_seal_area = back_seal_area
        self.p_restrain = p_restrain
        self.is_restrain = is_restrain
        self.with_brakes = with_brakes
        self.brakes_on = brakes_on
        self.x_buffer = x_buffer
        self.hit_buffer = hit_buffer
        self.f_decay = f_decay
        self.mass_limit = mass_limit
        #
        # The following will be assigned during assembly.
        self.slugL = None
        self.slugL_which_end = 'R'
        self.slugR = None
        self.slugR_which_end = 'L'
        #
        Piston.pistonList.append(self)
        return

    def write_to_ini_file(self, fp):
        """
        Writes the piston information to the specified file.
        """
        fp.write("[piston-%d]\n" % self.indx)
        fp.write("    label = %s\n" % self.label)
        fp.write("    front_seal_f = %e\n" % self.front_seal_f)
        fp.write("    front_seal_area = %e\n" % self.front_seal_area)
        fp.write("    back_seal_f = %e\n" % self.back_seal_f)
        fp.write("    back_seal_area = %e\n" % self.back_seal_area)
        fp.write("    mass = %e\n" % self.m)
        fp.write("    diameter = %e\n" % self.d)
        fp.write("    length = %e\n" % self.L)
        fp.write("    p_restrain = %e\n" % self.p_restrain)
        fp.write("    is_restrain = %d\n" % self.is_restrain)
        fp.write("    x_buffer = %e\n" % self.x_buffer)
        fp.write("    hit_buffer = %d\n" % self.hit_buffer)
        fp.write("    with_brakes = %d\n" % self.with_brakes)
        fp.write("    brakes_on = %d\n" % self.brakes_on)
        if self.slugL != None:
            indx = self.slugL.indx
        else:
            indx = -1
        fp.write("    left-slug-id = %d\n" % indx)
        fp.write("    left-slug-end-id = %s\n" % self.slugL_which_end)
        if self.slugR != None:
            indx = self.slugR.indx
        else:
            indx = -1
        fp.write("    right-slug-id = %d\n" % indx)
        fp.write("    right-slug-end-id = %s\n" % self.slugR_which_end)
        fp.write("    x0 = %e\n" % self.x0)
        fp.write("    v0 = %e\n" % self.v0)
        fp.write("    f_decay = %e\n" % self.f_decay)
        fp.write("    mass_limit = %e\n" % self.mass_limit)
        return
    
#----------------------------------------------------------------------------

class Diaphragm(object):
    """
    Contains the information for a diaphragm which controls the
    interaction of two GasSlugs.
    """

    __slots__ = 'indx', 'x0', 'p_burst', 'is_burst', \
                'slugL', 'slugR', \
                'slugL_which_end', 'slugR_which_end', \
                'dt_hold', 'dt_blend', 'dx_blend', \
                'RSP_dt', \
                'dxL', 'dxR', 'label'
    diaphragmList = []
    
    def __init__(self,
                 x0, p_burst, is_burst=0, dt_hold=0.0,
                 dt_blend=0.0, dx_blend=0.0, RSP_dt=0.0, 
                 dxL=0.0, dxR=0.0, label=""):
        """
        Creates a diaphragm with specified properties.

        The connections to GasSlugs are made later via the function
        assemble_gas_path.

        :param x0: (float) x-position in the tube, metres.
            This value is used to determine the end-points of the GasSlugs.
        :param p_burst: (float) Pressure, in Pa, at which rupture is triggered.
        :param is_burst: (int) Flag to indicate the state of diaphragm.
            A value of 0 indicates that the diaphragm is intact while
            a value of 1 indicates that the diaphragm is ruptured and the
            GasSlugs are interacting.
        :param dt_hold: (float) Time delay, in seconds, from rupture trigger
            to actual rupture.
        :param dt_blend: (float) Time delay, in seconds, from rupture to a
            blend event.
            This models the mixing of the two gas slugs some time after
            rupture of the diaphragm.
            Blending events are seldom used so this is usually set to 0.0.
        :param dx_blend: (float) Distance, in metres, over which blending occurs.
            Set to 0.0 to have no effective blending.
        :param dxL: (float) The distance over which p is averaged on left of
            the diaphragm.  The pressure difference between the left-
            and right-sided of the diaphragm is used to trigger rupture.
            The default value of 0.0 will cause the pressure in the
            gas cell immediately adjacent to the diaphragm to be used.
        :param dxR: (float) The distance, in metres, over which p is averaged
            on right-side of the diaphragm.
        :param label: A (string) label that will appear in the parameter file
            for this diaphragm.
        """
        self.indx = len(Diaphragm.diaphragmList) # next available index
        if len(label) > 0:
            self.label = label
        else:
            self.label = "diaphragm-" + str(self.indx)
        self.x0 = x0
        self.p_burst = p_burst
        self.is_burst = is_burst
        self.dt_hold = dt_hold + RSP_dt
        self.dt_blend = dt_blend
        self.dx_blend = dx_blend
        self.RSP_dt = RSP_dt
        self.dxL = dxL
        self.dxR = dxR
        #
        # The following will be assigned in assembly.
        self.slugL = None
        self.slugL_which_end = 'R'
        self.slugR = None
        self.slugR_which_end = 'L'
        #
        Diaphragm.diaphragmList.append(self)
        return

    def write_to_ini_file(self, fp):
        """
        Writes the diaphragm information to the specified file.
        """
        fp.write("[diaphragm-%d]\n" % self.indx)
        fp.write("    label = %s\n" % self.label)
        fp.write("    is_burst = %d\n" % self.is_burst)
        fp.write("    p_burst = %e\n" % self.p_burst)
        fp.write("    dt_hold = %e\n" % self.dt_hold)
        fp.write("    dt_blend = %e\n" % self.dt_blend)
        fp.write("    dx_blend = %e\n" % self.dx_blend)
        fp.write("    RSP_dt = %e\n" % self.RSP_dt)
        if self.slugL != None:
            indx = self.slugL.indx
        else:
            indx = -1
        fp.write("    left-slug-id = %d\n" % indx)
        fp.write("    left-slug-end-id = %s\n" % self.slugL_which_end)
        fp.write("    dxL = %e\n" % self.dxL)
        if self.slugR != None:
            indx = self.slugR.indx
        else:
            indx = -1
        fp.write("    right-slug-id = %d\n" % indx)
        fp.write("    right-slug-end-id = %s\n" % self.slugR_which_end)
        fp.write("    dxR = %e\n" % self.dxR)
        return
    
    
#----------------------------------------------------------------------------

class GasInterface(object):
    """
    Contains the information for an interface between two slugs.

    The primary use of this class is to locate the ends of
    the connected GasSlugs.
    Implicitly, the logical connections are also made via the
    function assemble_gas_path.
    """

    __slots__ = 'x0', 'slugL', 'slugL_which_end', \
                'slugR', 'slugR_which_end'
    interfaceList = []
    
    def __init__(self, x0):
        """
        Creates as interface between two L{GasSlug}s at specified location.

        :param x0: (float) Initial position, in metres.
        """
        self.x0 = x0
        self.slugL = None
        self.slugL_which_end = 'R'
        self.slugR = None
        self.slugR_which_end = 'L'
        GasInterface.interfaceList.append(self)
        return
    
#----------------------------------------------------------------------------

class FreeEnd(object):
    """
    Contains the information for a free-end condition.
    """

    __slots__ = 'x0'
    freeEndList = []
    
    def __init__(self, x0):
        """
        Creates a GasSlug end-condition with a specified location.

        :param x0: (float) Initial position, in metres.
        """
        self.x0 = x0
        FreeEnd.freeEndList.append(self)
        return

# --------------------------------------------------------------------

class VelocityEnd(object):
    """
    Contains the information for a fixed-velocity end condition
    for a GasSlug.
    """

    __slots__ = 'x0', 'v'
    velocityEndList = []
    
    def __init__(self, x0, v=0.0):
        """
        Creates a GasSlug end-condition with a specified location
        and velocity.

        :param x0: (float) Initial position, in metres.
        :param v: (float) Velocity, in m/s, of the end-point of the GasSlug.
        """
        self.x0 = x0
        self.v = v
        VelocityEnd.velocityEndList.append(self)
        return

# --------------------------------------------------------------------

def assemble_gas_path(*components):
    """
    Assembles a gas path by making the logical connections between
    adjacent components.

    The components are assembled left-to-right,
    as they are supplied to this function.

    :param components: An arbitrary number of arguments representing
        individual components or lists of components.
        Each component may be a GasSlug, Piston, or any
        other gas-path object, however, it doesn't always make sense
        to connect arbitrary components.
        For example, connecting a GasSlug to a Piston is reasonable
        but connecting a Piston to a Diaphragm without an intervening
        GasSlug does not make sense in the context of this simulation
        program.
    """
    print "Assemble gas path:"
    clist = []
    for c in components:
        if isinstance(c,tuple) or isinstance(c,list):
            clist.extend(c)
        else:
            clist.append(c)
    for i in range(len(clist)-1):
        connect_pair(clist[i], clist[i+1])
    #
    # We now need to go through the component list and,
    # for any GasInterface components, we need to connect
    # the slugs on either side.  Once this is done,
    # the GasInterface objects have done their job.
    for i in range(len(clist)):
        if isinstance(clist[i], GasInterface):
            connect_slugs(clist[i-1], clist[i+1])
    return


def connect_slugs(cL, cR):
    """
    Make the logical connection between a pair of gas slugs.
    
    :param cL: is left slug
    :param cR: is right slug

    Usually called by assemble_gas_path.
    """
    print "connect_slugs()"
    if isinstance(cL, GasSlug) and isinstance(cR, GasSlug):
        print "   Info: make slug <--> slug connection"
        cL.bcR = cR
        cL.bcR_which_end = 'L'
        cR.bcL = cL
        cR.bcL_which_end = 'R'
    else:
        raise Exception, "Error: Both objects must be GasSlugs."
    return


def connect_pair(cL, cR):
    """
    Make the logical connection between a pair of components.
    
    :param cL: is left object
    :param cR: is right object

    Usually called by assemble_gas_path.
    """
    # print "connect_pair()"
    # print "    left component", cL
    # print "    right component", cR

    if isinstance(cL,VelocityEnd) and isinstance(cR, GasSlug):
        cR.bcL = cL
        cR.xL = cL.x0
        print "    velocity-end <--> gas-slug is done"
    elif isinstance(cL,FreeEnd) and isinstance(cR, GasSlug):
        cR.bcL = cL
        cR.xL = cL.x0
        print "    free-end <--> gas-slug is done"
    elif isinstance(cL,GasInterface) and isinstance(cR, GasSlug):
        cR.bcL = cL
        cR.xL = cL.x0
        print "    gas-interface <--> gas-slug is done"
    elif isinstance(cL,Piston) and isinstance(cR, GasSlug):
        cL.slugR = cR
        cL.slugR_which_end = 'L'
        cR.bcL = cL
        cR.bcL_which_end = 'R'
        cR.xL = cL.xR0
        print "    piston <--> gas-slug is done"
    elif isinstance(cL,Diaphragm) and isinstance(cR, GasSlug):
        cL.slugR = cR
        cL.slugR_which_end = 'L'
        cR.bcL = cL
        cR.bcL_which_end = 'R'
        cR.xL = cL.x0
        print "    diaphragm <--> gas-slug is done"
    elif isinstance(cL,GasSlug) and isinstance(cR, VelocityEnd):
        cL.bcR = cR
        cL.xR = cR.x0
        print "    gas-slug <--> velocity-end is done"
    elif isinstance(cL,GasSlug) and isinstance(cR, FreeEnd):
        cL.bcR = cR
        cL.xR = cR.x0
        print "    gas-slug <--> free-end is done"
    elif isinstance(cL,GasSlug) and isinstance(cR, GasInterface):
        cL.bcR = cR
        cL.xR = cR.x0
        print "    gas-slug <--> gas-interface is done"
    elif isinstance(cL,GasSlug) and isinstance(cR, Piston):
        cL.bcR = cR
        cL.bcR_which_end = 'L'
        cL.xR = cR.xL0
        cR.slugL = cL
        cR.slugL_which_end = 'R'
        print "    gas-slug <--> piston is done"
    elif isinstance(cL,GasSlug) and isinstance(cR, Diaphragm):
        cL.bcR = cR
        cL.bcR_which_end = 'L'
        cL.xR = cR.x0
        cR.slugL = cL
        cR.slugL_which_end = 'R'
        print "    gas-slug <--> diaphragm is done"
    else:
        raise Exception, "    Invalid pair to connect."
    return


# --------------------------------------------------------------------

def write_parameter_file():
    """
    Writes the traditional (ugly) input-parameter file
    from the data that is stored in the GlobalData object and
    the other lists of objects that make up the gas-path.
    """
    global gdata
    print "Begin write parameter file."
    fp = open(gdata.param_file, "w")
    gdata.write_to_ini_file(fp, len(GasSlug.slugList),
                            len(Piston.pistonList),
                            len(Diaphragm.diaphragmList))
    for p in Piston.pistonList: p.write_to_ini_file(fp)
    for d in Diaphragm.diaphragmList: d.write_to_ini_file(fp)
    for s in GasSlug.slugList: s.write_to_ini_file(fp)
    print "End write parameter file."
    return


# --------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin l_script.py..."

    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or \
           uoDict.has_key("--help") or \
           uoDict.has_key("-h"):
        printUsage()
    else:
        if uoDict.has_key("--job"):
            jobName = uoDict.get("--job", "test")
        elif uoDict.has_key("-f"):
            jobName = uoDict.get("-f", "test")
        else:
            jobName = "test"
        rootName, ext = os.path.splitext(jobName)
        if os.path.exists(jobName):
            jobFileName = jobName
        else:
            jobFileName = rootName + ".py"
        print "Job file: %s" % jobFileName
        gdata.param_file = rootName + ".Lp"

        # The user-specified input comes in the form of Python code.
        # In a parallel calculation, all processes should see the same setup.
        execfile(jobFileName)
        print "Summary of components:"
        print "    gas slugs         :", len(GasSlug.slugList)
        print "    pistons           :", len(Piston.pistonList)
        print "    diaphragms        :", len(Diaphragm.diaphragmList)
        print "    free-ends         :", len(FreeEnd.freeEndList)
        print "    velocity-ends     :", len(VelocityEnd.velocityEndList)
        print "    gas-gas interfaces:", len(GasInterface.interfaceList)
        if len(GasSlug.slugList) < 1:
            print "Warning: no gas slugs defined; this is unusual."
        write_parameter_file()
    print "Done."
    sys.exit(0)



#! /usr/bin/env python
"""
e3prep.py: Python program to specify the flow problem.

It is intended for the user to define the flow simulation in terms of
the data objects defined in this program.
As part of its initialization, e3prep.py will execute a user-specified
job file that contains, in Python, the user's script that defines both
geometry and flow details.

The flow simulation definition is organised via the classes:
GlobalData2D, FlowCondition, HeatZone, SimplePiston, Block2D and Block3D.
These classes provide places to store the configuration information and
their function (or method) names appear as commands in the job
description file.
See the __init__() method for each class to determine what parameters
can be specified when creating an object of that class.

The user will define the particular geometry in terms
of the data objects defined in the geom and gpath modules.
This geometry definition is created in a bottom-up approach
by successively defining Nodes, simple path elements
(such as Line, Arc and Bezier elements) and, possibly,
compound path elements (such as Spline objects and Polyline objects).
Four paths are then used to define a patch over which
a block of finite-volume cells is defined.
The four bounding faces of each block also carry
boundary-condition information.

note: Physical quantities should be specified in MKS units.

gdata: Contains the GlobalData information describing the simulation.
    Note that there is one such variable set up by the main program and
    the user's script should directly set the attributes of this variable
    to adjust settings for the simulation.
sketch: A global variable holding the state of the SketchEnvironment.
    Note that there is one such variable set up by the main program and
    the user's script should directly set the attributes of this variable
    to adjust settings (such as scale and ranges of the axes)
    for the SVG output.

\author P.Jacobs
\version February 2005 as scriptit.py
\version 15-Jan-2006 renamed as mbcns_prep.py and modified to work with libgeom2.
\version 17-Jan-2006 reinstated Node rendering and eliminated option --do-grid-gen
                     because we will always want to generate grids (replacing mb_prep.c)
\version 02-Mar-2008 Elmer3 port started
\version 17-Mar-2008 Block3D merged.
\version 24-Feb-2009 changed to using folders for the grid and flow files.
"""

# ----------------------------------------------------------------------
#
import sys
import os
sys.path.append("/sw/lib/python2.3/site-packages/Numeric") # for Tim's MacOSX
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from getopt import getopt
from gzip import GzipFile
import copy
import math
import traceback

from libprep3 import *
from e3_defs import *
from e3_block import *
from e3_render import *
# from e3_gas import * # now in gaspy.i included in libprep3.i
# Make one instance to accumulate the settings for 2D SVG rendering.
sketch = SketchEnvironment()

shortOptions = ""
longOptions = ["help", "job=", "do-svg", "do-vrml", 
               "zip-files", "no-zip-files", "show-names"]

def printUsage():
    print ""
    print "Usage: e3prep.py [--help] [--job=<jobFileName>]"
    print "       [--do-svg] [--do-vrml] [--zip-files|--no-zip-files]"
    print "       [--show-names] [--split-input-file]"
    return

#----------------------------------------------------------------------

turbulence_model_list = ['k_omega', 'baldwin_lomax', 'spalart_allmaras']


def select_gas_model(model=None, species=None, fname=None):
    """
    Selects a gas model for the simulation.

    Input:
    model   : (string) name of the gas model as shown in the list below.
    species : list of species names (strings).
    fname   : (string) name of the gas-model file.

    The new gas models are configured by stand-alone files.
    This function initialises a gas model for present use
    in the preparation program (ie. sets it in kernel code)
    and stores the name for later user at simulation time.

    If you already have a gas-model.lua file already set up,
    give its name af fname.
    
    If you have not already set up the gas-model.lua file,
    this function is provided as a simple but limited means to do so.

    Look-up-table (LUT) gas models and LUT_plus_composite cannot be
    created directly with this function.
    If you want a single LUT gas model, just set up the LUT table 
    file externally and supply the name of that file as fname.
    If you want a LUT-plus-composite gas model, set up the LUT table
    externally and then set up the rest of the composite gas model
    using create_gas_file(), which has the capability of prepending 
    the LUT gas species.
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
    gmodel = set_gas_model_ptr(create_gas_model(fname))
    gdata.gas_model_file = fname
    nsp = gmodel.get_number_of_species()
    return [ gmodel.species_name(isp) for isp in range(nsp) ]

def set_reaction_scheme(fname, reacting_flag=1):
    """
    Sets the reaction update model and specifies a reacting simulation.

    This function sets the name of the input file for the reaction update.
    It also sets the reacting flag.  Note that the user may later reset
    that flag, possibly to turn off reactions for a particular simulation.
    """
    gdata.reacting_flag = reacting_flag
    gdata.reaction_update = fname
    return

# We want to keep the old name, for a while.
set_reaction_update = set_reaction_scheme

def set_energy_exchange_scheme(fname, energy_exchange_flag=1):
    """
    Sets the energy exchange update model and specifies a TNE simulation.

    This function sets the name of the input file for the energy exchange update.
    It also sets the energy exchange flag.  Note that the user may later reset
    that flag, possibly to turn off energy exchange for a particular simulation.
    """
    gdata.energy_exchange_flag = energy_exchange_flag
    gdata.energy_exchange_update = fname
    return

# We want to keep the old name, for a while.
set_energy_exchange_update = set_energy_exchange_scheme
    
def select_radiation_model( input_file="no_file", update_frequency=1 ):
    """
    Selects the radiation spectral and transport models 
    for the simulation.
    
    This function is provided to the user as a means
    to setup a radiation model.
    
    This function initialises a radiation model for present use
    in the preparation program (ie. sets it in kernel code)
    and stores the name for later user at simulation time.
    """
    
    if input_file!="no_file":
	# 1. store model and input file names for config file creation
	gdata.radiation_flag = 1
	gdata.radiation_input_file = input_file
	gdata.radiation_update_frequency = update_frequency
	# 2. set the models for immediate use in the simulation preparation
	set_radiation_transport_model( input_file )
	set_radiation_update_frequency( update_frequency ) 
    else:
	print "The field 'input_file' is required when selecting the"
	print "radiation model via the functionn'select_radiation_model'."
	print "Bailing out!"
	sys.exit(1)
	
    return

#----------------------------------------------------------------------

class GlobalData(object):
    """
    Organise, remember and write the global data.

    The user's script should not create one of these
    but should specify the simulation parameters by
    altering the attributes of the global object L{gdata}.

    title: (string) A piece of text that will be propagated through the
        solution files and subsequently generated plots.
    case_id: (int) An identifier for special cases in which pieces of
        specialised code have been embedded into the main simulation
        program.
        If you don't have such code to activate, the default value of 0 is fine.
    gas_model_file: (string) The input file required for configuring the gas model
    axisymmetric_flag: (0/1) A value of 0 sets two-dimensional, planar flow.
        A value of 1 sets axisymmetric flow with the x-axis being the axis
        of symmetry.
    viscous_flag: (0/1) Set to 1 to activate viscous transport terms.
        Set to 0 (the default) for inviscid flow simulation.
    viscous_delay: (float) Sometimes, the viscous terms make it difficult to
        start a calculation without encountering numerical instability.
        Set this parameter to the delay (in seconds) from simulation start
        to the time at which the viscous terms will be allowed to become
        active (if L{viscous_flag} was set to 1).
    viscous_factor_increment: (float) Increment for the viscous_factor
        after viscous_delay has passed.
    max_mu_t_factor: (float) Limiting factor for the turbulence viscosity
        relative to the laminar viscosity.
    transient_mu_t_factor: (float) A value of 1.0 indicates a fully-developed
        turbulent flow.  Set values less than 1.0 to model transient flows
        where the turbulence is not fully developed.  This may be good for
        Matt McGilvray's scramjet experiments in the expansion tube which have
        very short flow durations.
    diffusion_flag: (0/1) Set to 1 to compute multi-component diffusion of species
    diffusion_model: (string) set the type (by name) of multi-component diffusion model
    k_omega_flag: (0/1) Allow the k-omega turbulence model to be active.
    turbulence_prandtl_number: (float) default 0.89
    turbulence_schmidt_number: (float) default 0.75
    heat_time_start: (float) start time for heating zones to be adding energy
    heat_time_stop: (float) final time for heating zones to be adding heat
    heat_factor_increment: (float) the fraction of full heat load that will be
        added with each time step from t=heat_time_start
    reacting_flag: (0/1) A value of 1 will make Rowan Gollan's finite-rate
        chemistry active if the appropriate gas_name (e.g. 'perf_gas_mix')
        has been specified.
    reaction_time_start: (float) Time after which reactions are allowed to proceed.
    reaction_update: (string) A (file) name for the chemical kinetics scheme
    energy_exchange_flag: (0/1) A flag indicting finite-rate evolution of thermal state
    energy_exchange_update: (string) A (file) name for the thermal energy exchange scheme
    x_order: (int 1 or 2) Specifies the form of reconstruction from cell-average
        data to cell interface data.
        Select 1 for low-order (i.e. no) reconstruction.
        Select 2 for a higer-order (limited quadratic) reconstruction.
    t_order: (int 1 or 2) Specifies the form of time stepping scheme.
        Select 1 for Euler stepping.
        Select 2 for predictor-corrector stepping.
    stringent_cfl: (0/1) Set to 1 to get a very strict CFL check.
    flux_calc: (int or string) Specifies the form of flux calculation at cell interfaces.
        See module L{flux_dict.py} for options.
    interpolation_type: (string) Choose the set of thermo variables to use in interpolation.
        options: "rhoe", "rhop", "rhoT", "pT", default "rhoe"
    compression_tolerance: (float) relative velocity change for the shock detector.
        This value is expected to be a negative number (for compression)
        and not too large in magnitude. default -0.30
    shear_tolerance: (float) normalised tangential velocity change beyond which we
        we suppress the application of EFM in the ADAPTIVE flux calculator.
    apply_limiter_flag : (0/1) Set to 1 to have reconstruction limiter enabled (default)
        Set to 0 for no limiting
    dt: (float) Size of the initial time step.
        After a few steps, the solver will have enough information to select
        a suitable time step, based on the L{cfl} number.
    cfl: (float) The ratio of the actual time step to the allowed time
        step as determined by the flow condition and grid.
        Typically the default value of 0.5 is good but you may want
        to try smaller values if you are having the solution go unstable,
        especially for viscous simulations.
    print_count: (int) The number of time steps between printing a status message
        to the console.
    cfl_count: (int) The number of time steps between checking the CFL number.
    fixed_time_step: (boolean) Flag to indicate fixed time-stepping
    max_invalid_cells: (int) Number of cells that will be tolerated without too
        much complaint.
    dt_reduction_factor: (float) If a time step fails because of any bad cells,
        this is the factor to reduce the size of the next attempted time step.
    max_time: (float) The (simulation) time (in seconds) at which time stepping
        should be terminated.
    max_step: (int) Time stepping will be terminated if the simulation reached
        this number of steps.
    dt_plot: (float) Period between writing all of the flow field data to the
        solution file.
        Multiple instances can be written to the one file but be careful not to
        write too many and fill up your disk.
    dt_history: (float) Period (in seconds) between writing the data for the
        selected cells to the history file.
    """
    count = 0

    # We want to prevent the user's script from introducing new attributes
    # via typographical errors.
    __slots__ = 'dimensions', 'param_file', 'mpost_file', 'title', 'case_id', \
                'gas_model_file', 'reaction_update', \
                'axisymmetric_flag', \
                'radiation_input_file', 'radiation_update_frequency', 'radiation_flag', \
                'implicit_flag', 'control_count', \
                'radiation_update_frequency', \
                'viscous_flag', 'viscous_delay', 'viscous_factor_increment', \
                'max_mu_t_factor', 'transient_mu_t_factor', \
                'diffusion_flag', 'diffusion_model', \
                'turbulence_flag', 'turbulence_model', \
                'turbulence_prandtl_number', 'turbulence_schmidt_number', \
                'scalar_pdf_flag', 'reacting_flag', 'reaction_time_start', \
                'x_order', 'flux_calc', 'compression_tolerance', 'shear_tolerance', \
                't_order', 'stringent_cfl', \
                't0', 'dt', 'cfl', 'dt_chem', 'dt_therm', \
                'interpolation_type', 'sequence_blocks', \
                'print_count', 'cfl_count', 'max_invalid_cells', 'dt_reduction_factor', \
                'max_time', 'max_step', 'dt_plot', 'dt_history', \
                'displacement_thickness', 'time_average_flag', 'perturb_flag', \
                'perturb_frac', 'tav_0', 'tav_f', 'dt_av', \
                'fixed_time_step', 'apply_limiter_flag', \
                'energy_exchange_flag', 'energy_exchange_update', \
                'udf_file', 'udf_source_vector_flag', \
                'heat_time_start', 'heat_time_stop', 'heat_factor_increment'
    
    def __init__(self):
        """
        Accepts user-specified data and sets defaults. Make one only.
        """
        if GlobalData.count >= 1:
            raise Exception, "Already have a GlobalData object defined."
        self.dimensions = 2
        self.param_file = "temp.p"
        self.title = "Another Eilmer3 Simulation."
        self.case_id = 0
        self.udf_file = "" # default is no file
        self.udf_source_vector_flag = 0
        self.gas_model_file = "gas_model.lua"
        self.reaction_update = "dummy_scheme"
        self.reacting_flag = 0
        self.reaction_time_start = 0.0
        self.energy_exchange_flag = 0
        self.energy_exchange_update = "dummy_scheme"
        self.radiation_input_file = "no_file"
        self.radiation_update_frequency = 1
        self.radiation_flag = 0
        self.axisymmetric_flag = 0
        self.implicit_flag = 0
        self.radiation_update_frequency = 1
        self.viscous_flag = 0
        self.viscous_delay = 0.0
        self.viscous_factor_increment = 0.01
        self.max_mu_t_factor = 300.0
        self.transient_mu_t_factor = 1.0
        self.diffusion_flag = 0
        self.diffusion_model = "None"
        self.turbulence_flag = 0
        self.turbulence_model = turbulence_model_list[0]
        self.turbulence_prandtl_number = 0.89
        self.turbulence_schmidt_number = 0.75
        self.heat_time_start = 0.0
        self.heat_time_stop = 0.0 # nonzero indicates that we want heating some time
        self.heat_factor_increment = 0.01
        self.x_order = 2
        self.interpolation_type = "rhoe"
        self.apply_limiter_flag = 1
        self.flux_calc = ADAPTIVE
        self.compression_tolerance = -0.30
        self.shear_tolerance = 0.20
        self.t_order = 2
        self.cfl = 0.5
        self.stringent_cfl = 0
        self.fixed_time_step = False
        self.dt_reduction_factor = 0.2
        self.t0 = 0.0 
        # may be useful to change t0 if we are restarting from another job
        self.dt = 1.0e-6
        self.dt_chem = -1.0
        self.dt_therm = -1.0
        self.sequence_blocks = 0 
        # 0 = normal time iteration all blocks; 
        # 1 = iterate blocks one at a time (like space marching)
        self.print_count = 20
        self.control_count = 10
        self.cfl_count = 10
        self.max_invalid_cells = 10
        self.max_time = 1.0e-3
        self.max_step = 10
        self.dt_plot = 1.0e-3
        self.dt_history = 1.0e-3
        GlobalData.count += 1
        return

    def write_to_control_file(self, fp):
        """
        Writes the time-stepping control data to the specified file in .ini format.

        Since the main simulation loop checks this file every time step,
        it is possible to change these parameters as the simulation runs.
        One common use is to terminate the simulation cleanly by changing
        halt_now from 0 to 1.
        """
        fp.write("[control_data]\n")
        fp.write("x_order = %d\n" % self.x_order)
        fp.write("t_order = %d\n" % self.t_order)
        fp.write("dt = %e\n" % self.dt)
        fp.write("fixed_time_step = %s\n" % self.fixed_time_step)
        fp.write("dt_reduction_factor = %e\n" % self.dt_reduction_factor)
        fp.write("cfl = %e\n" % self.cfl)
        fp.write("stringent_cfl = %d\n" % self.stringent_cfl)
        fp.write("print_count = %d\n" % self.print_count)
        fp.write("cfl_count = %d\n" % self.cfl_count)
        fp.write("dt_plot = %e\n" % self.dt_plot)
        fp.write("dt_history = %e\n" % self.dt_history)
        fp.write("max_time = %e\n" % self.max_time)
        fp.write("max_step = %d\n" % self.max_step)
        fp.write("implicit_flag = %d\n" % self.implicit_flag)
        fp.write("radiation_update_frequency = %d\n" % self.radiation_update_frequency)
        fp.write("halt_now = 0\n"); # presumably, we want the simulation to proceed
        return

    def write_to_ini_file(self, fp):
        """
        Writes the configuration data to the specified file in .ini format.
        """
        fp.write("[global_data]\n")
        fp.write("title = %s\n" % self.title)
        fp.write("dimensions = %d\n" % self.dimensions)
        fp.write("case_id = %d\n" % self.case_id)
        fp.write("udf_file = %s\n" % self.udf_file)
        fp.write("udf_source_vector_flag = %d\n" % self.udf_source_vector_flag)
        fp.write("gas_model_file = %s\n" % self.gas_model_file)
        fp.write("reaction_update = %s\n" % self.reaction_update)
        fp.write("reacting_flag = %d\n" % self.reacting_flag)
        fp.write("reaction_time_start = %e\n"% self.reaction_time_start)
        fp.write("energy_exchange_flag = %d\n" % self.energy_exchange_flag)
        fp.write("energy_exchange_update = %s\n" % self.energy_exchange_update)
        fp.write("radiation_input_file = %s\n" % self.radiation_input_file)
        fp.write("radiation_update_frequency = %s\n" % self.radiation_update_frequency)
        fp.write("radiation_flag = %d\n" % self.radiation_flag)
        fp.write("viscous_flag = %d\n" % self.viscous_flag)
        fp.write("viscous_delay = %e\n"% self.viscous_delay)
        fp.write("viscous_factor_increment = %e\n"% self.viscous_factor_increment)
        fp.write("max_mu_t_factor = %e\n"% self.max_mu_t_factor)
        fp.write("transient_mu_t_factor = %e\n"% self.transient_mu_t_factor)
        fp.write("diffusion_flag = %d\n" % self.diffusion_flag)
        fp.write("diffusion_model = %s\n" % self.diffusion_model)
        fp.write("axisymmetric_flag = %d\n"% self.axisymmetric_flag)
        fp.write("turbulence_flag = %d\n" % self.turbulence_flag)
        fp.write("turbulence_model = %s\n" % self.turbulence_model.lower())
        fp.write("turbulence_prandtl_number = %g\n" % self.turbulence_prandtl_number)
        fp.write("turbulence_schmidt_number = %g\n" % self.turbulence_schmidt_number)
        fp.write("heat_time_start = %e\n"% self.heat_time_start)
        fp.write("heat_time_stop = %e\n"% self.heat_time_stop)
        fp.write("heat_factor_increment = %e\n"% self.heat_factor_increment)
        fp.write("flux_calc = %d\n" % self.flux_calc)
        fp.write("compression_tolerance = %e\n"% self.compression_tolerance)
        fp.write("shear_tolerance = %e\n"% self.shear_tolerance)
        fp.write("interpolation_type = %s\n" % self.interpolation_type)
        fp.write("apply_limiter_flag = %d\n" % self.apply_limiter_flag )
        fp.write("sequence_blocks = %d\n" % self.sequence_blocks)
        fp.write("max_invalid_cells = %d\n" % self.max_invalid_cells)
        fp.write("control_count = %d\n" % self.control_count)
        return

# We will create just one GlobalData object that the user can alter.
gdata = GlobalData()

#----------------------------------------------------------------------------

def lower_upper(p0, p1):
    """
    Returns the lower point and upper point to a volume, 2D-patch or line.

    Used as a helper function when setting a HeatZone, ReactionZone or
    TurbulenceZone.
    The user can specify any diagonally-opposite corners and the code will
    see the presumed lower and upper corners.
    """
    p_lower = Vector(min(p0.x,p1.x), min(p0.y,p1.y), min(p0.z,p1.z))
    p_upper = Vector(max(p0.x,p1.x), max(p0.y,p1.y), max(p0.z,p1.z))
    return p_lower, p_upper


class HeatZone(object):
    """
    Organise the description for heat-addition zones.

    Each zone is intended to model energy deposition into the flow domain
    maybe approximating a electric spark.
    """
    zoneList = []

    __slots__ = 'zoneId', 'qdot', 'point0', 'point1', 'label'

    def __init__(self, qdot, point0, point1, label=""):
        """
        Sets up a hexahedral heat-addition zone defined by its diagonal corner points.

        Note that it doesn't matter if the defined zone extends outside of a block
        because of the way the test is done within the source terms.

        qdot: (float) rate of heat addition in W/m**3
        point0: (Vector) lower corner of the zone
        point1: (Vector) upper corner of the zone
        
        """
        if not isinstance(qdot, float):
            raise TypeError, ("qdot should be a float but it is: %s" % type(qdot))
        if not isinstance(point0, Vector):
            raise TypeError, ("point0 should be a Vector but it is: %s" % type(point0))
        if not isinstance(point1, Vector):
            raise TypeError, ("point1 should be a Vector but it is: %s" % type(point1))
        self.qdot = qdot
        self.point0, self.point1 = lower_upper(Vector(point0), Vector(point1))
        self.zoneId = len(HeatZone.zoneList)    # next available index
        if len(label) == 0:
            label = "zone-" + str(self.zoneId)
        self.label = label
        HeatZone.zoneList.append(self)
        return

    def write_to_ini_file(self, fp):
        """
        Writes the HeatZone information to the specified file-object in .ini format.
        """
        fp.write("\n[heat_zone/%d]\n" % self.zoneId)
        fp.write("label = %s\n" % self.label)
        fp.write("qdot = %e\n" % self.qdot)
        fp.write("x0 = %e\n" % self.point0.x)
        fp.write("y0 = %e\n" % self.point0.y)
        fp.write("z0 = %e\n" % self.point0.z)
        fp.write("x1 = %e\n" % self.point1.x)
        fp.write("y1 = %e\n" % self.point1.y)
        fp.write("z1 = %e\n" % self.point1.z)
        return

#----------------------------------------------------------------------------

class ReactionZone(object):
    """
    Describe the zones in which finite-rate reactions are allowed.

    These zones will be used to set a flag contained within
    each finite-volume cell of the main simulation.
    """
    zoneList = []

    __slots__ = 'zoneId', 'qdot', 'point0', 'point1', 'label'

    def __init__(self, point0, point1, label=""):
        """
        Sets up a hexahedral heat-addition zone defined by its diagonal corner points.

        Note that it doesn't matter if the defined zone extends outside of a block.

        point0: (Vector) lower corner of the zone
        point1: (Vector) upper corner of the zone
        
        """
        if not isinstance(point0, Vector):
            raise TypeError, ("point0 should be a Vector but it is: %s" % type(point0))
        if not isinstance(point1, Vector):
            raise TypeError, ("point1 should be a Vector but it is: %s" % type(point1))
        self.point0, self.point1 = lower_upper(Vector(point0), Vector(point1))
        self.zoneId = len(ReactionZone.zoneList)    # next available index
        if len(label) == 0:
            label = "zone-" + str(self.zoneId)
        self.label = label
        ReactionZone.zoneList.append(self)
        return

    def write_to_ini_file(self, fp):
        """
        Writes the ReactionZone information to the specified file-object in .ini format.
        """
        fp.write("\n[reaction_zone/%d]\n" % self.zoneId)
        fp.write("label = %s\n" % self.label)
        fp.write("x0 = %e\n" % self.point0.x)
        fp.write("y0 = %e\n" % self.point0.y)
        fp.write("z0 = %e\n" % self.point0.z)
        fp.write("x1 = %e\n" % self.point1.x)
        fp.write("y1 = %e\n" % self.point1.y)
        fp.write("z1 = %e\n" % self.point1.z)
        return
        

#----------------------------------------------------------------------------

class TurbulenceZone(object):
    """
    Describe the zones in which the turbulence models are active.

    These zones will be used to set a flag contained within
    each finite-volume cell of the main simulation.
    """
    zoneList = []

    __slots__ = 'zoneId', 'qdot', 'point0', 'point1', 'label'

    def __init__(self, point0, point1, label=""):
        """
        Sets up a hexahedral turbulent-addition zone defined by its diagonal corner points.

        Note that it doesn't matter if the defined zone extends outside of a block.

        point0: (Vector) lower corner of the zone
        point1: (Vector) upper corner of the zone
        
        """
        if not isinstance(point0, Vector):
            raise TypeError, ("point0 should be a Vector but it is: %s" % type(point0))
        if not isinstance(point1, Vector):
            raise TypeError, ("point1 should be a Vector but it is: %s" % type(point1))
        self.point0, self.point1 = lower_upper(Vector(point0), Vector(point1))
        self.zoneId = len(TurbulenceZone.zoneList)    # next available index
        if len(label) == 0:
            label = "zone-" + str(self.zoneId)
        self.label = label
        TurbulenceZone.zoneList.append(self)
        return

    def write_to_ini_file(self, fp):
        """
        Writes the TurbulenceZone information to the specified file-object in .ini format.
        """
        fp.write("\n[turbulence_zone/%d]\n" % self.zoneId)
        fp.write("label = %s\n" % self.label)
        fp.write("x0 = %e\n" % self.point0.x)
        fp.write("y0 = %e\n" % self.point0.y)
        fp.write("z0 = %e\n" % self.point0.z)
        fp.write("x1 = %e\n" % self.point1.x)
        fp.write("y1 = %e\n" % self.point1.y)
        fp.write("z1 = %e\n" % self.point1.z)
        return
        
#----------------------------------------------------------------------------

class SimplePiston(object):
    """
    Contains the information for a simple piston model.

    The present model is restricted to a simple, frictionless piston that
    slides in a tube defined by the (parallel) NORTH and SOUTH boundaries
    of a set of blocks.
    We'll expand this description later to include friction, etc.

    x0 : x-coordinate position of the centre of mass
    v0 : x-coordinate velocity of the centre of mass
    D  : y-coordinate diameter, m
    L  : x-coordinate length, m
    m  : mass, kg
    """
    __slots__ = 'pistonId', 'label', 'd', 'L', 'm', 'x0', 'v0', 'f', 'const_v_flag', 'postv_v_flag'

    # Accumulate references to pistons
    pistonList = []

    def __init__(self, m, d, xL0, xR0, v0=0.0, f=0.0, label="", \
                 const_v_flag="False", postv_v_flag="False"):
        """
        Initialises a piston's geometric parameters and state.
        m : mass, kg
        d : effective diameter, m
        xL0 : x-position of left-face, m
        xR0 : x-position of right-face, m
        v0 : initial velocity, m/s
        f : friction coefficient
        label : (optional) string label
        const_v_flag: piston will not accelerate
        postv_v_flag: negative velocities forbidden
        """
        self.pistonId = len(SimplePiston.pistonList)    # next available index
        if len(label) == 0:
            label = "piston-" + str(self.pistonId)
        self.label = label
        self.d = d
        self.L = xR0 - xL0
        self.m = m
        self.x0 = 0.5 * (xL0 + xR0)
        self.v0 = v0
        self.f = f
        self.const_v_flag = const_v_flag
        self.postv_v_flag = postv_v_flag
        
        SimplePiston.pistonList.append(self)

        return

    def write_starting_solution(self, fp):
        """
        Writes the initial (t=0) piston state to an already open file.
        """
        print "Begin write initial piston state: label=", self.label
        fp.write("%d %20.12e %20.12e %20.12e\n" % (self.pistonId, 0.0, self.x0, self.v0) )
        return

    def write_to_ini_file(self, fp):
        """
        Writes the piston config information to the specified file-object in .ini format.
        """
        fp.write("\n[piston/%d]\n" % self.pistonId)
        fp.write("label = %s\n" % self.label)
        fp.write("D = %e\n" % self.d)
        fp.write("L = %e\n" % self.L)
        fp.write("m = %e\n" % self.m)
        fp.write("x0 = %e\n" % self.x0)
        fp.write("v0 = %e\n" % self.v0)
        fp.write("f = %e\n" % self.f)
        fp.write("const_v_flag = %s\n" % self.const_v_flag)
        fp.write("postv_v_flag = %s\n" % self.postv_v_flag)
        return

#----------------------------------------------------------------------------

class HistoryLocation(object):
    """
    Contains the Cartesian (x,y,z) location of a history object.
    """
    __slots__ = 'x', 'y', 'z', 'label', 'i', 'j', 'k', 'i_offset', 'j_offset', \
        'k_offset', 'blkId', 'cellId'

    # Accumulate references to history locations
    historyList = []

    def __init__(self, x, y, z=0.0, i_offset=0, j_offset=0, k_offset=0, label=""):
        """
        Initialises a history location.

        This location will later be tied to a cell in the grid.
        """
        if not isinstance(x, float):
            raise TypeError, ("x should be a float but it is: %s" % type(x))
        if not isinstance(y, float):
            raise TypeError, ("y should be a float but it is: %s" % type(y))
        if not isinstance(z, float):
            raise TypeError, ("z should be a float but it is: %s" % type(z))
        self.x = x
        self.y = y
        self.z = z
        self.label = label
        self.i = 0
        self.j = 0
        self.k = 0
        self.i_offset = i_offset
        self.j_offset = j_offset
        self.k_offset = k_offset
        self.blkId = 0
        self.cellId = 0
        HistoryLocation.historyList.append(self)
        return

# --------------------------------------------------------------------

def locate_history_cells():
    """
    Given the Cartesian coordinates, locate appropriate cells.
    """
    global gdata
    #
    for h in HistoryLocation.historyList:
        x_target = h.x
        y_target = h.y
        z_target = h.z
        #
        # Create an initial guess
        x, y, z, vol = Block.blockList[0].cell_centre_location(0,0,0,gdata)
        min_dist = math.sqrt( (x-x_target)**2 + (y-y_target)**2 + (z-z_target)**2 )
        best_block = 0; best_i = 0; best_j = 0; best_k = 0
        #
        for b in Block.blockList:
            # print "BLOCK-%d" % b.blkId
            if gdata.dimensions == 3:
                for k in range(b.grid.nk-1):
                    for j in range(b.grid.nj-1):
                        for i in range(b.grid.ni-1):
                            x,y,z,vol = b.cell_centre_location(i,j,k,gdata)
                            dist = math.sqrt( (x-x_target)**2 + (y-y_target)**2 + (z-z_target)**2 )
                            if dist < min_dist:
                                min_dist = dist
                                best_block = b.blkId
                                best_i = i; best_j = j; best_k = k
            else:
                k = 0
                for j in range(b.grid.nj-1):
                    for i in range(b.grid.ni-1):
                        x,y,z,vol = b.cell_centre_location(i,j,k,gdata)
                        dist = math.sqrt( (x-x_target)**2 + (y-y_target)**2 )
                        if dist < min_dist:
                            min_dist = dist
                            best_block = b.blkId
                            best_i = i; best_j = j
        #
        b = Block.blockList[best_block]
        print "For history location: ", x_target, y_target, z_target
        print "    Closest grid cell is at: block= ", best_block, \
              "i=", best_i, "j=", best_j, "k=", best_k
        best_i += h.i_offset
        if best_i < 0: best_i = 0
        if best_i > b.nni-1: best_i = b.nni-1
        best_j += h.j_offset
        if best_j < 0: best_j = 0
        if best_j > b.nnj-1: best_j = b.nnj-1
        best_k += h.k_offset
        if best_k < 0: best_k = 0
        if best_k > b.nnk-1: best_k = b.nnk-1
        print "    After offsets: i=", best_i, "j=", best_j, "k=", best_k
        Block.blockList[best_block].hcell_list.append( (best_i, best_j, best_k) )
        print "    For block", best_block ,"this becomes history cell index=", \
              len(Block.blockList[best_block].hcell_list) - 1
        h.i = best_i; h.j = best_j; h.k = best_k;
        h.blkId = best_block
        h.cellId = len(Block.blockList[best_block].hcell_list) - 1
    #
    # FIX-ME update the following table for 3D also.
    if len(HistoryLocation.historyList) > 0:
        hfp = open("history_cells.list", 'w')
        hfp.write("List of history cells associated with %s \n" % gdata.title )
        hfp.write("\n")
        hfp.write("========================================================")
        hfp.write("========================================================\n")
        hfp.write("|   history cell   |    x    |    y    |    z    |  block  ")
        hfp.write("|   i  |   j  |   k  |  cell id  |  cells in block  |\n")
        hfp.write("========================================================")
        hfp.write("========================================================\n")
    for h in HistoryLocation.historyList:
        hfp.write("|   %-11s    | %6.4g  | %6.4g  | %6.4g  |  %3d    " % 
                  (h.label, h.x, h.y, h.z, h.blkId))
        hfp.write("| %3d  | %3d  | %3d  |   %3d     |       %3d        |\n" %
                  (h.i, h.j, h.k, h.cellId, len(Block2D.blockList[h.blkId].hcell_list)))
        hfp.write("--------------------------------------------------------")
        hfp.write("--------------------------------------------------------\n")
    #
    if len(HistoryLocation.historyList) > 0:
        hfp.close()
    return

#----------------------------------------------------------------------------

def write_times_file(rootName):
    print "Begin write zero-entry into times file."
    fp = open(rootName+".times", "w")
    fp.write("# tindx sim_time dt_global\n")
    fp.write("%04d %e %e\n" % (0, gdata.t0, gdata.dt))
    fp.close()
    return


def write_parameter_file(rootName):
    print "Begin write configuration file (INI format)."
    fp = open(rootName+".config", "w")
    gdata.write_to_ini_file(fp)
    fp.write("npiston = %d\n" %len(SimplePiston.pistonList) )
    fp.write("nflow = %d\n" % len(FlowCondition.flowList) )
    fp.write("nheatzone = %d\n" % len(HeatZone.zoneList))
    fp.write("nreactionzone = %d\n" % len(ReactionZone.zoneList))
    fp.write("nturbulencezone = %d\n" % len(TurbulenceZone.zoneList))
    fp.write("nblock = %d\n" % len(Block.blockList) )
    for piston in SimplePiston.pistonList:
        piston.write_to_ini_file(fp)
    for flow in FlowCondition.flowList:
        flow.write_to_ini_file(fp)
    for zone in HeatZone.zoneList:
        zone.write_to_ini_file(fp)
    for zone in ReactionZone.zoneList:
        zone.write_to_ini_file(fp)
    for zone in TurbulenceZone.zoneList:
        zone.write_to_ini_file(fp)
    for block in Block.blockList:
        block.write_to_ini_file(fp, gdata.dimensions)
    fp.write("\n# end file\n")
    fp.close()
    print "End write config file."
    return


def write_control_file(rootName):
    print "Begin write control file (INI format)."
    fp = open(rootName+".control", "w")
    gdata.write_to_control_file(fp)
    fp.close()
    return


def write_grid_files(rootName, blockList, zipFiles=0):
    print "Begin write grid file(s)."
    # Grid already created in main loop
    # Write one file per block.
    gridPath = os.path.join("grid", "t0000")
    if not os.access(gridPath, os.F_OK):
        os.makedirs(gridPath)
    for b in blockList:
        fileName = rootName+(".grid.b%04d.t0000" % b.blkId)
        fileName = os.path.join(gridPath, fileName)
        if zipFiles:
            fp = GzipFile(fileName+".gz", "wb")
        else:
            fp = open(fileName, "w")
        b.grid.write(fp)
        fp.close()
    print "End write grid file(s)."
    return

def split_input_file(blockList):
    print "Begin split input file(s)."
    # Write one input file per block.
    ncell = 0
    start_indx = 1
    end_indx = 1
    for b in blockList:
        if gdata.dimensions == 3:
            faceList = faceList3D
        else:
            faceList = faceList2D
        for iface in faceList:
            split_file_flag = 0
            bc = b.bc_list[iface]
            if bc.type_of_BC == 20 or bc.type_of_BC == 18 and os.path.exists(bc.filename):
                split_file_flag = 1
            if split_file_flag == 1:
                if iface == 0 or iface == 2:
                    ncell += b.nni
                    end_indx += b.nni
                    block_ncell = b.nni
                else:
                    ncell += b.nnj
                    end_indx += b.nnj
                    block_ncell = b.nnj
                input_file = bc.filename
                bc.filename = input_file+(".blk%04d" % b.blkId)
                inputData = open(input_file, 'r').read().split('\n')
                outputData = inputData[start_indx:end_indx]
                fp = open(bc.filename, 'w')
                fp.write("%d\n" % block_ncell)
                fp.write('\n'.join(outputData))
                fp.close()
                start_indx = end_indx
                ncell_for_profile = int(open(input_file, 'r').readline())
                if ncell_for_profile != ncell and iface == faceList[-1]:
                    print "split_input_file():"
                    print "    Inconsistent numbers of cells: ncell = %d" % ncell
                    print "    ncell_for_profile = %d" % ncell_for_profile
                    print "    Bailing out!"
                    sys.exit(1)

    print "End split input file(s)."
    return

def write_starting_solution_files(rootName, blockList, pistonList, zipFiles=0):
    print "Begin write starting solution file(s)."
    # Write one file per block.
    flowPath = os.path.join("flow", "t0000")
    if not os.access(flowPath, os.F_OK):
        os.makedirs(flowPath)
    for b in blockList:
        fileName = rootName+(".flow.b%04d.t0000" % b.blkId)
        fileName = os.path.join(flowPath, fileName)
        if zipFiles:
            fp = GzipFile(fileName+".gz", "wb")
        else:
            fp = open(fileName, "w")
        b.write_starting_solution(fp, gdata)
        fp.close()
    if len(pistonList) > 0:
        fp = open(rootName + ".piston.t0000", "w")
        for p in pistonList:
            p.write_starting_solution(fp)
        fp.close()
    print "End write starting solution file(s)."
    return

# --------------------------------------------------------------------
RenderList = []

def add_item_to_render_list(item):
    """
    Add the given item to the RenderList
    """
    RenderList.append(item)
    return

def render_to_vrml(rootName):
    """
    Render the items referenced in RenderList to the wrl file.
 
    Since the 3D rendering is a bit more complicated, we will start simple
    and provide a way to let the user specify what will be rendered to VRML.
    Items may be surfaces or volumes (rendered as 6 surfaces).

    Maybe we'll work out a way to add annotations and
    boundary-condition information at a later date.
    """
    fileName = rootName + ".wrl"
    print "Render to VRML in file", fileName
    fp = open(fileName, "w")
    fp.write("#VRML V2.0 utf8\n")
    for item in RenderList:
        fp.write(item.vrml_str() + "\n")
    fp.close()
    return

# --------------------------------------------------------------------

def main(uoDict):
    """
    Top-level function for the e3prep application.

    It may be handy to be able to embed most of the functions in 
    this file into a custom preprocessing script.
    """
    jobName = uoDict.get("--job", "test")
    rootName, ext = os.path.splitext(jobName)
    sketch.root_file_name = rootName
    if os.path.exists(jobName):
        jobFileName = jobName
    else:
        jobFileName = rootName + ".py"
    print "Job file: %s" % jobFileName
    zipFiles = 1  # Default: use zip file format for grid and flow data files.
    if uoDict.has_key("--zip-files"): zipFiles = 1
    if uoDict.has_key("--no-zip-files"): zipFiles = 0
    #
    # The user-specified input comes in the form of Python code.
    # In a parallel calculation, all processes should see the same setup.
    execfile(jobFileName, globals())
    #
    if uoDict.has_key("--split-input-file"):
        # Split input file
        split_input_file(Block.blockList)
        sys.exit(0)
    #
    locate_history_cells()
    write_times_file(rootName)
    write_parameter_file(rootName)
    write_control_file(rootName)
    #
    if len(FlowCondition.flowList) < 1:
        print "Warning: no flow conditions defined."
    if len(Block.blockList) < 1:
        print "Warning: no blocks defined."
    #
    bll = open("block_labels.list", "w")
    bll.write("# indx label\n");
    for b in Block.blockList:
        bll.write("%d %s\n" % (b.blkId,b.label) )
    bll.close()
    #
    write_grid_files(rootName, Block.blockList, zipFiles)
    write_starting_solution_files(rootName, Block.blockList, SimplePiston.pistonList, zipFiles)
    if gdata.viscous_flag and (not os.access("heat", os.F_OK)):
        os.makedirs("heat")
    #
    if uoDict.has_key("--do-svg") and gdata.dimensions == 2:
        sketch.write_svg_file(gdata, FlowCondition.flowList, Block.blockList, faceList2D)
    #
    if uoDict.has_key("--do-vrml") and gdata.dimensions == 3:
        render_to_vrml(rootName)
    return

def pretty_print_names(nameList):
    """
    Display the list of names, a few per line.
    """
    nameList.sort()
    print "There are", len(nameList), "names defined."
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for letter in letters:
        localList = []
        for name in nameList:
            if name[0].upper() == letter:
                localList.append(name)
        if len(localList) > 0:
            print "Names beginning with %s or %s (%d):" % (letter, letter.lower(), len(localList))
            count = 0
            print "    ",
            for name in localList:
                if count == len(localList)-1:
                    print name
                elif (count+1)%3 == 0:
                    print name
                    print "    ",
                else:
                    print name,
                count += 1
    return

# --------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin e3prep.py..."
    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
        sys.exit(1) # abnormal exit
    elif uoDict.has_key("--show-names"):
        pretty_print_names(dir())
        sys.exit(0) # normal exit, presuming this is the desired response
    else:
        try:
            main(uoDict)
        except:
            print "This run of e3prep.py has gone bad."
            traceback.print_exc(file=sys.stdout)
            sys.exit(2) # abnormal exit; we tried and failed
    print "Done."
    sys.exit(0) # Finally, a normal exit.

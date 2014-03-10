## \file bc_defs.py
## \ingroup eilmer3
##
## \author P.Jacobs
## \version 31-Jan-2005 extracted from e_model_spec.py
## \version 2013 clean up object definitions.  Introduce symbols for the BCs.
##

import copy

"""
Historically, within the C/C++ simulation code, 
the boundary conditions were identified by integer constants
for which macro names were defined.
These definitions are now gone and symbols are used instead:

* ADJACENT: This boundary joins that of another block.
    Normally, this boundary condition would be set implicitly
    when making block connections.
* COMMON: Synonym for ADJACENT.
* SUP_IN: Fully-prescribed inflow (e.g. supersonic inflow).
* SUP_OUT: Synonym for EXTRAPOLATE_OUT.
* SHOCK_FITTING_IN: Moving boundary with fully-prescribed inflow flux.
* EXTRAPOLATE_OUT: Extrapolate all flow properties from
   just inside the boundary into the ghost-cells outside
   the boundary.  This works fine for a strong supersonic outflow.
* SLIP: Synonym for SLIP_WALL
* SLIP_WALL: A solid but inviscid wall.
    Effectively, this boundary condition copies and reflects the
    properties just inside the boundary into the ghost cells.
* ADIABATIC: A solid, no-slip wall without heat transfer.
    (i.e. the near-wall temperature is reflected in the ghost cells)
* FIXED_T: A solid, no-slip wall with a user-specified
    temperature.
* SLIDING_T: A solid, no-slip wall with a sliding user-specified
    temperature that linearly varies with time.
* SUBSONIC_IN: An inflow boundary for which the total pressure
    and temperature have been specified and the pressure from
    just inside the boundary is copied into the ghost cells
    along with a velocity that is consistent with the stagnation conditions.
    If you specify a non-zero mass-flux per unit area, the pressure
    will be adjusted to achieve this mass-flux.
* SUBSONIC_OUT: An outflow boundary which will try to prevent
    wave reflection at the boundary in the presence of subsonic flow.
    (Doesn't work so well at present.)
* SUB_OUT: Synonym for SUBSONIC_OUT
* TRANSIENT_UNI: An transient inflow boundary which has
    a uniform flow condition applied across the full boundary.
* TRANSIENT_PROF
* STATIC_PROF: A steady inflow boundary with a variable set of
    flow conditions across the boundary.
* FIXED_P_OUT: Something like EXTRAPOLATE_OUT but with the
    pressure set to some user-specified value.
    It is probably best to set this pressure at the same value as
    the initial fill pressure so that this boundary condition will
    be passive until a wave arrives at the boundary.
* RRM: Andrew Denman's recycled and renormalised boundary
    condition.
* SEB: Surface energy balance.
* USER_DEFINED
* ADJACENT_PLUS_UDF
* ABLATING
* SLIDING_T
* FSTC
* SHOCK_FITTING_IN
* USER_DEFINED_MASS_FLUX
* CONJUGATE_HT
* MASS_FLUX_OUT
* MAPPED_CELL
"""
#
# The following symbol definitions are for use in the user's
# input script if they wish to refer to BC types.
#
ADJACENT = object()
SUP_IN = object()
SUP_OUT = object()
EXTRAPOLATE_OUT = object()
SLIP_WALL = object()
ADIABATIC = object()
FIXED_T = object()
SUBSONIC_IN = object()
SUBSONIC_OUT = object()
TRANSIENT_UNI = object()
TRANSIENT_PROF = object()
STATIC_PROF = object()
FIXED_P_OUT = object()
RRM = object()
TRANSIENT_T_WALL = object()
SEB = object()
USER_DEFINED = object()
ADJACENT_PLUS_UDF = object()
ABLATING = object()
SLIDING_T = object()
FSTC = object()
SHOCK_FITTING_IN = object()
USER_DEFINED_MASS_FLUX = object()
CONJUGATE_HT = object()
MOVING_WALL = object()
MASS_FLUX_OUT = object()
MAPPED_CELL = object()

#
# When we ust the set_BC method for a Block object, we will want to look up
# the correct boundary condition object by name.
# The integer values in the following dictionary are a reminder of the old
# macro definitions in the C code.  They are retained here for reference.
# Newer boundary conditions will just have a name, in upper case.
#
bcSymbolFromName = {
     0: ADJACENT, "0": ADJACENT, "ADJACENT": ADJACENT, "COMMON": ADJACENT,
     1: SUP_IN, "1": SUP_IN, "SUP_IN": SUP_IN,
     2: EXTRAPOLATE_OUT, "2":  EXTRAPOLATE_OUT, "SUP_OUT": EXTRAPOLATE_OUT,
         "EXTRAPOLATE_OUT": EXTRAPOLATE_OUT,
     3: SLIP_WALL, "3": SLIP_WALL, "SLIP": SLIP_WALL, "SLIP_WALL": SLIP_WALL,
     4: ADIABATIC, "4": ADIABATIC, "ADIABATIC": ADIABATIC,
     5: FIXED_T, "5": FIXED_T, "FIXED_T": FIXED_T,
     6: SUBSONIC_IN, "6": SUBSONIC_IN, "SUBSONIC_IN": SUBSONIC_IN,
     7: SUBSONIC_OUT, "7": SUBSONIC_OUT, "SUBSONIC_OUT": SUBSONIC_OUT,
         "SUB_OUT": SUBSONIC_OUT,
     8: TRANSIENT_UNI, "8": TRANSIENT_UNI, "TRANSIENT_UNI": TRANSIENT_UNI,
     9: TRANSIENT_PROF, "9":  TRANSIENT_PROF,  "TRANSIENT_PROF": TRANSIENT_PROF,
    10: STATIC_PROF, "10": STATIC_PROF,  "STATIC_PROF": STATIC_PROF,
    11: FIXED_P_OUT, "11": FIXED_P_OUT, "FIXED_P_OUT": FIXED_P_OUT, 
    12: RRM, "12": RRM, "RRM": RRM,    
    13: TRANSIENT_T_WALL, "13" : TRANSIENT_T_WALL, "TRANSIENT_T_WALL": TRANSIENT_T_WALL,
    15: SEB, "15" : SEB, "SEB" : SEB, "SURFACE_ENERGY_BALANCE" : SEB,
    16: USER_DEFINED, "16": USER_DEFINED, "USER_DEFINED": USER_DEFINED,
    17: ADJACENT_PLUS_UDF, "17": ADJACENT_PLUS_UDF, "ADJACENT_PLUS_UDF": ADJACENT_PLUS_UDF,
    18: ABLATING, "18" : ABLATING, "ABLATING": ABLATING,
    19: SLIDING_T, "19" : SLIDING_T, "SLIDING_T": SLIDING_T,
    20: FSTC, "20" : FSTC, "FSTC": FSTC,
    21: SHOCK_FITTING_IN, "21" : SHOCK_FITTING_IN, "SHOCK_FITTING_IN": SHOCK_FITTING_IN,
    "USER_DEFINED_MASS_FLUX": USER_DEFINED_MASS_FLUX,
    "CONJUGATE_HT": CONJUGATE_HT,
    "MOVING_WALL": MOVING_WALL,
    "MASS_FLUX_OUT": MASS_FLUX_OUT,
    "MAPPED_CELL": MAPPED_CELL
}
bcName = {
    ADJACENT: "ADJACENT",
    SUP_IN: "SUP_IN",
    EXTRAPOLATE_OUT: "EXTRAPOLATE_OUT",
    SLIP_WALL: "SLIP_WALL",
    ADIABATIC: "ADIABATIC",
    FIXED_T: "FIXED_T",
    SUBSONIC_IN: "SUBSONIC_IN",
    SUBSONIC_OUT: "SUBSONIC_OUT",
    TRANSIENT_UNI: "TRANSIENT_UNI",
    TRANSIENT_PROF: "TRANSIENT_PROF",
    STATIC_PROF: "STATIC_PROF",
    FIXED_P_OUT: "FIXED_P_OUT",
    RRM: "RRM",
    TRANSIENT_T_WALL: "TRANSIENT_T_WALL",
    SEB: "SEB",
    USER_DEFINED: "USER_DEFINED",
    ADJACENT_PLUS_UDF: "ADJACENT_PLUS_UDF",
    ABLATING: "ABLATING",
    SLIDING_T: "SLIDING_T",
    FSTC: "FSTC",
    SHOCK_FITTING_IN: "SHOCK_FITTING_IN",
    USER_DEFINED_MASS_FLUX: "USER_DEFINED_MASS_FLUX",
    CONJUGATE_HT: "CONJUGATE_HT",
    MOVING_WALL: "MOVING_WALL",
    MASS_FLUX_OUT: "MASS_FLUX_OUT",
    MAPPED_CELL: "MAPPED_CELL"
    }

class BoundaryCondition(object):
    """
    Base class for boundary condition specifications.
    """
    __slots__ = 'type_of_BC', 'Twall', 'Pout', 'inflow_condition', \
                'x_order', 'sponge_flag', 'other_block', 'other_face', 'orientation', \
                'filename', 'n_profile', 'is_wall', 'sets_conv_flux', 'sets_visc_flux', \
                'assume_ideal', 'mdot', 'Twall_i', 'Twall_f', 't_i', 't_f', 'emissivity', \
                'r_omega', 'centre', 'v_trans', 'Twall_flag', \
                'reorient_vector_quantities', 'Rmatrix', \
                'mass_flux', 'p_init', 'relax_factor', \
                'direction_type', 'direction_vector', 'direction_alpha', 'direction_beta', \
                'ghost_cell_trans_fn', 'label'
    def __init__(self,
                 type_of_BC=SLIP_WALL,
                 Twall=300.0,
                 Pout=100.0e3,
                 inflow_condition=None,
                 x_order=0,
                 sponge_flag=0,
                 other_block=-1,
                 other_face=-1,
                 orientation=0,
                 filename="",
                 n_profile=1,
                 is_wall=0,
                 sets_conv_flux=0,
                 sets_visc_flux=0,
                 assume_ideal=0,
                 mdot = [],
                 Twall_i=300.0,
                 Twall_f=300.0,
                 t_i=0.0,
                 t_f=0.0,
                 emissivity=1.0,
                 r_omega=None,
                 centre=None,
                 v_trans=None,
                 Twall_flag=False,
                 reorient_vector_quantities=False,
                 Rmatrix=[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                 mass_flux=0.0,
                 p_init=100.0e3,
                 relax_factor=0.05,
                 direction_type="normal",
                 direction_vector=[1.0, 0.0, 0.0],
                 direction_alpha=0.0,
                 direction_beta=0.0,
                 ghost_cell_trans_fn=lambda x, y, z: (x, y, z),
                 label=""):
        """
        Construct a generic boundary condition object.

        This constructor is usually invoked by one of the more specific
        boundary-condition constructors. 
        It is a catch-all for the bits of data that might be required by
        any particular boundary condition.

        :param type_of_BC: specifies the boundary condition (symbol value)
        :param Twall: fixed wall temperature (in degrees K) that will be used if
            the boundary conditions needs such a value.
        :param Pout: fixed outside pressure (in Pascals) that will be used if
            the boundary conditions needs such a value.
        :param inflow_condition: the flow condition that will be applied if the
            specified boundary condition needs it.
        :param x_order: Extrapolation order of the boundary conduition.
            0=just copy the nearest cell data into both ghost cells 
            (zero-order extrapolation).
            1=linear extrapolation of the interior data into the ghost cells. 
        :param sponge_flag: A value of 1 will activate Andrew Denman's damping
            terms near the boundary.
        :param other_block: index to an adjacent block, if any. 
            A value of -1 will indicate that there is no adjacent block.
        :param other_face: index of the adjacent face of the other block, if any.
        :param orientation: for 3D connections the other block face can have one of
            4 rotational orientations.
        :param filename: Name of the UDF source file (in Lua)
            or of the profile data, if relevant.
        :param n_profile: Number of profiles to be found in the input data file
            for the StaticProfileBC.
        :param is_wall: Flag to indicate that various parts of the simulation code 
            should treat this boundary as a solid wall.
        :param sets_conv_flux: For this boundary, the fluxes are computed directly.
            Typically, this relates to using a user-supplied Lus script which
            provides a convective_flux() function. This pretty much ignores the ghost-cell
            data, however, it does not relieve the user of supplying a suitable 
            function for setting that data.
        :param sets_visc_flux: As for sets_conv_flux except that this relates to
            setting the viscous component of flux due to the effect of the boundary.
        :param assume_ideal:
        :param mdot: species ablation rate (list)
        :param emissivity: surface radiative emissivity (between 0 and 1)
        :param Twall_i: initial temperature for sliding temperature BC
        :param Twall_f: final temperature for sliding temperature BC
        :param t_i: initial time for sliding temperature BC
        :param t_f: final time for sliding temperature BC
        :param r_omega: angular velocity for Jason Qin's moving-wall boundary
        :param centre: a point on the axis of rotation for the moving-wall boundary
        :param v_tran: a translational velocity to superimpost on the moving-wall boundary 
        :param Twall_flag: a boolean parameter to select fixed Twall for moving-wall boundary
        :param reorient_vector_quantities: for exchange of vector quantities between adjacent boundaries
        :param Rmatrix: the 9 elements of the rotation matrix
        :param mass_flux: mass flux per unit area (in kg/s/m**2) across the block boundary
        :param p_init: initial pressure (in Pa) at the mass-flux-out boundary
        :param relax_factor: relaxation factor for adjustment of the actual pressure applied
            to the ghost-cells of the mass-flux-out boundary or the subsonic-in boundary
        :param direction_type: "normal" (default) is to have the inflow velocity
            locally-normal to the subsonic-in boundary.
            "uniform" has the inflow velocity aligned with direction_vector
            "radial" radial-inflow (turbine) through a cylindrical surface
            with flow angles direction_alpha and direction_beta.
            "axial" axial-flow (turbine) through a circular surface 
            with flow angles direction_alpha and direction_beta.
        :param direction_vector: List of x,y,z-components specifying the direction that the
            inflow velocity follows if direction_type=="uniform"
        :param direction_alpha: flow angle splitting the in-plane velocity into radial
            and tangential components for direction_type=="radial" or "axial" (radians)
        :param direction_beta: flow angle determining axial-velocity component
            for direction_type=="radial" or "axial" (radians)
        :param ghost_cell_trans_fn: User-supplied transform function mapping from 
            ghost-cell position to source-cell position.
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        self.type_of_BC = type_of_BC
        if Twall is None:
            self.Twall = 300.0
        else:
            self.Twall = Twall
        self.Pout = Pout
        self.inflow_condition = inflow_condition
        self.x_order = x_order
        self.sponge_flag = sponge_flag
        self.other_block = other_block
        self.other_face = other_face
        self.orientation = orientation
        self.filename = filename
        self.n_profile = n_profile
        self.is_wall = is_wall
        self.sets_conv_flux = sets_conv_flux
        self.sets_visc_flux = sets_visc_flux
        self.assume_ideal = assume_ideal
        self.mdot = mdot
        self.Twall_i = Twall_i
        self.Twall_f = Twall_f
        self.t_i = t_i
        self.t_f = t_f
        self.emissivity = emissivity
        if r_omega is None:
            self.r_omega = [0.0, 0.0, 0.0]
        else:
            self.r_omega = [r_omega[0], r_omega[1], r_omega[2]]
        if centre is None:
            self.centre = [0.0, 0.0, 0.0]
        else:
            self.centre = [centre[0], centre[1], centre[2]]
        if v_trans is None:
            self.v_trans = [0.0, 0.0, 0.0]
        else:
            self.v_trans = [v_trans[0], v_trans[1], v_trans[2]]
        self.reorient_vector_quantities = reorient_vector_quantities
        self.Twall_flag = Twall_flag
        assert (type(Rmatrix) is list) and (len(Rmatrix) == 9)
        self.Rmatrix = Rmatrix
        self.mass_flux = mass_flux
        self.p_init = p_init
        self.relax_factor = relax_factor
        self.direction_type = direction_type
        assert (type(direction_vector) is list) and (len(direction_vector) == 3)
        self.direction_vector = [direction_vector[0], direction_vector[1], direction_vector[2]]
        self.direction_alpha = direction_alpha
        self.direction_beta = direction_beta
        self.ghost_cell_trans_fn = ghost_cell_trans_fn
        self.label = label
            
        return
    def __str__(self):
        str_rep = "BoundaryCondition("
        str_rep += "type_of_BC=%s" % bcName[self.type_of_BC]
        str_rep += ", Twall=%g" % self.Twall
        str_rep += ", Pout=%g" % self.Pout
        str_rep += ", inflow_condition=%d" % self.inflow_condition
        str_rep += ", x_order=%d" % self.x_order
        str_rep += ", sponge_flag=%d" % self.sponge_flag
        str_rep += ", other_block=%d" % self.other_block
        str_rep += ", other_face=%d" % self.other_face
        str_rep += ", orientation=%d" % self.orientation
        str_rep += ", filename=\"%s\"" % self.filename
        str_rep += ", n_profile=%d" % self.n_profile
        str_rep += ", is_wall=%d" % self.is_wall
        str_rep += ", sets_conv_flux=%d" % self.sets_conv_flux
        str_rep += ", sets_visc_flux=%d" % self.sets_visc_flux
        str_rep += ", assume_ideal=%d" % self.assume_ideal
        str_rep += ", mdot=["
        for mdi in mdot: str_rep += "%g," % mdi
        str_rep += "]"
        str_rep += ", emissivity=%g" % self.emissivity
        str_rep += ", r_omega=[%g, %g, %g]" % (self.r_omega[0], self.r_omega[1], self.r_omega[2])
        str_rep += ", centre=[%g, %g, %g]" % (self.centre[0], self.centre[1], self.centre[2])
        str_rep += ", v_trans=[%g, %g, %g]" % (self.v_trans[0], self.v_trans[1], self.v_trans[2])
        str_rep += ", reorient_vector_quantities=%d" % self.reorient_vector_quantities
        str_rep += ", Twall_flag=%s" % self.Twall_flag
        str_rep += ", Rmatrix=["
        for elem in Rmatrix: str_rep += "%g, " % elem
        str_rep += "]"
        str_rep += ", mass_flux=%g" % self.mass_flux
        str_rep += ", p_init=%g" % self.p_init
        str_rep += ", relax_factor=%g" % self.relax_factor
        str_rep += ", direction_type=\"%s\"" % self.direction_type
        str_rep += ", direction_vector=[%g, %g, %g]" % \
            (self.direction_vector[0], self.direction_vector[1], self.direction_vector[2])
        str_rep += ", direction_alpha=%g" % self.direction_alpha
        str_rep += ", direction_beta=%g" % self.direction_beta
        str_rep += ", ghost_cell_trans_fn=%s" % self.ghost_cell_trans_fn
        str_rep += ", label=\"%s\")" % self.label
        return str_rep
    def __copy__(self):
        return BoundaryCondition(type_of_BC=self.type_of_BC,
                                 Twall=self.Twall,
                                 Pout=self.Pout,
                                 inflow_condition=self.inflow_condition,
                                 x_order=self.x_order,
                                 sponge_flag=self.sponge_flag,
                                 other_block=self.other_block,
                                 other_face=self.other_face,
                                 orientation=self.orientation,
                                 filename=self.filename,
                                 n_profile=self.n_profile,
                                 is_wall=self.is_wall,
                                 sets_conv_flux=self.sets_conv_flux,
                                 sets_visc_flux=self.sets_visc_flux,
                                 assume_ideal=self.assume_ideal,
                                 mdot=copy.copy(self.mdot),
                                 Twall_i=self.Twall_i,
                                 Twall_f=self.Twall_f,
                                 t_i=self.t_i,
                                 t_f=self.t_f,
                                 emissivity=self.emissivity,
                                 r_omega=copy.copy(self.r_omega),
                                 centre=copy.copy(self.centre),
                                 v_trans=copy.copy(self.v_trans),
                                 reorient_vector_quantities=self.reorient_vector_quantities,
                                 Twall_flag=self.Twall_flag,
                                 Rmatrix=copy.copy(self.Rmatrix),
                                 mass_flux=self.mass_flux,
                                 p_init=self.p_init,
                                 relax_factor=self.relax_factor,
                                 direction_type=self.direction_type,
                                 direction_vector=copy.copy(self.direction_vector),
                                 direction_alpha=self.direction_alpha,
                                 direction_beta=self.direction_beta,
                                 ghost_cell_trans_fn=self.ghost_cell_trans_fn,
                                 label=self.label)
    
class AdjacentBC(BoundaryCondition):
    """
    This boundary joins (i.e. is adjacent to) a boundary of another block.
    """
    def __init__(self, other_block=-1, other_face=-1, orientation=0, 
                 reorient_vector_quantities=False,
                 Rmatrix=[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], label=""):
        """
        Join the boundary face to a boundary-face of another block.

        This condition is usually not set manually but is set as part of the
        connect_blocks() function.
        There should be a corresponding AdjacentBC on the other block and
        the connect_blocks() function ensures this.

        :param other_block: index to an adjacent block, if any. 
            A value of -1 will indicate that there is no adjacent block.
        :param other_face: index of the adjacent face of the other block, if any.
        :param orientation: for 3D connections the other block face can have one of
            4 rotational orientations.
        :param reorient_vector_quantities: for exchange of vector quantities between adjacent boundaries
        :param Rmatrix: the 9 elements of the rotation matrix
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        if reorient_vector_quantities:
            assert (type(Rmatrix) is list) and (len(Rmatrix) == 9)
        BoundaryCondition.__init__(self, type_of_BC=ADJACENT, other_block=other_block,
                                   other_face=other_face, orientation=orientation,
                                   reorient_vector_quantities=reorient_vector_quantities,
                                   Rmatrix=Rmatrix,
                                   label=label)
        return
    def __str__(self):
        return ("AdjacentBC(other_block=%d, other_face=%d, orientation=%d, " +
                "reorient_vector_quantities=%d, Rmatrix=%s, " +
                "label=\"%s\")" %
                (self.other_block, self.other_face, self.orientation, 
                 self.reorient_vector_quantities, self.Rmatrix,
                 self.label))
    def __copy__(self):
        return AdjacentBC(other_block=self.other_block,
                          other_face=self.other_face,
                          orientation=self.orientation,
                          reorient_vector_quantities=self.reorient_vector_quantities,
                          Rmatrix=self.Rmatrix,
                          label=self.label)
    
class SupInBC(BoundaryCondition):
    """
    Apply a (presumably) supersonic inflow condition to the boundary.

    The inflow_condition data is copied into the ghost cells.
    """
    def __init__(self, inflow_condition, label=""):
        """
        Construct a supersonic-inflow boundary condition.

        :param inflow_condition: A reference to a previously-constructed
            FlowCondition object.
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=SUP_IN,
                                   inflow_condition=inflow_condition,
                                   label=label)
        return
    def __str__(self):
        return "SupInBC(inflow_condition=%s, label=\"%s\")" % \
            (self.inflow_condition, self.label)
    def __copy__(self):
        return SupInBC(self.inflow_condition, self.label)

class ExtrapolateOutBC(BoundaryCondition):
    """
    Fill the ghost cells with data from just inside the boundary.

    This boundary condition will work best if the flow is supersonic,
    directed out of the flow domain.
    """
    def __init__(self, x_order=0, sponge_flag=0, label=""):
        """
        Construct an outflow BC that extrapolates the interior flow data.

        :param x_order: Extrapolation order of the boundary conduition.
            0=just copy the nearest cell data into both ghost cells 
            (zero-order extrapolation).
            1=linear extrapolation of the interior data into the ghost cells. 
        :param sponge_flag: A value of 1 will activate Andrew Denman's damping
            terms near the boundary.
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=EXTRAPOLATE_OUT,
                                   x_order=x_order,
                                   sponge_flag=sponge_flag,
                                   label=label)
        return
    def __str__(self):
        return "ExtrapolateOutBC(x_order=%d, sponge_flag=%d, label=\"%s\")" % \
            (self.x_order, self.sponge_flag, self.label)
    def __copy__(self):
        return ExtrapolateOutBC(x_order=self.x_order, sponge_flag=self.sponge_flag,
                                label=self.label)

class ShockFittingInBC(BoundaryCondition):
    """
    Apply a (presumably) supersonic inflow condition to the boundary.
    
    The inflow_condition data is copied into the ghost cells, and unlike a 
    SupInBC, the flux from the ghost cells to the domain is specified.
    """
    def __init__(self, inflow_condition, label=""):
        """
        Construct a shock-inflow boundary condition.

        :param inflow_condition: A reference to a previously-constructed
            FlowCondition object.
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=SHOCK_FITTING_IN,
                                   inflow_condition=inflow_condition,
                                   label=label)
        return
    def __str__(self):
        return "ShockFittingInBC(inflow_condition=%s, label=\"%s\")" % \
            (self.inflow_condition, self.label)
    def __copy__(self):
        return ShockFittingInBC(self.inflow_condition, self.label)
        
class SlipWallBC(BoundaryCondition):
    """
    An inviscid-flow solid-boundary.

    Effectively, this boundary condition copies and reflects the
    properties just inside the boundary into the ghost cells.
    This is the default boundary condition applied to a block face
    if you don't otherwise specify a boundary condition.
    """
    def __init__(self, emissivity=1.0, label=""):
        """
        Construct a no-friction, solid-wall boundary.

        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=SLIP_WALL, is_wall=1,
                                   emissivity=emissivity, label=label)
        return
    def __str__(self):
        return "SlipWallBC(label=\"%s\")" % self.label
    def __copy__(self):
        return SlipWallBC(emissivity=self.emissivity,label=self.label)
    
class AdiabaticBC(BoundaryCondition):
    """
    A solid, no-slip wall without heat transfer.
    
    The near-wall temperature is reflected in the ghost cells.
    This BC is only effective if viscous effects are active
    else it acts as another solid (slip) wall.
    """
    def __init__(self, label=""):
        """
        Construct a no-slip, no-heat-transfer, solid-wall boundary.

        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=ADIABATIC, is_wall=1, label=label)
        return
    def __str__(self):
        return "AdiabaticBC(label=\"%s\")" % self.label
    def __copy__(self):
        return AdiabaticBC(label=self.label)

class FixedTBC(BoundaryCondition):
    """
    A solid boundary with no-slip and a user specified temperature.

    Like the AdiabaticBC, this is completey effective only when viscous
    effects are active.  Else, it is just like another solid (slip) wall.
    """
    def __init__(self, Twall, emissivity=1.0, label=""):
        """
        Construct a no-slip, fixed-temperature, solid-wall boundary.

        :param Twall: fixed wall temperature (in degrees K) 
        :param emissivity: surface radiative emissivity (between 0 and 1) 
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=FIXED_T, Twall=Twall, 
                                   is_wall=1, emissivity=emissivity, label=label)
        return
    def __str__(self):
        return "FixedTBC(Twall=%g, label=\"%s\")" % (self.Twall, self.label)
    def __copy__(self):
        return FixedTBC(Twall=self.Twall, emissivity=self.emissivity, label=self.label)

class fstcBC(BoundaryCondition):
    """
    A solid boundary with no-slip and a user-specified temperature.

    Like the AdiabaticBC, this is completey effective only when viscous
    effects are active.  Else, it is just like another solid (slip) wall.
    """
    def __init__(self, filename="fstc_temperatures.dat", label=""):
        """
        Construct a no-slip, fixed-temperature, solid-wall boundary.

        :param filename: containing the specified temperatures
            for each cell along the boundary. 
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=FSTC, filename=filename, is_wall=1, label=label)
        return
    def __str__(self):
        return "fstcBC(filename=\"%s\", label=\"%s\")" % (self.filename, self.label)
    def __copy__(self):
        return fstcBC(filename=self.filename, label=self.label)

class SlidingTBC(BoundaryCondition):
    """
    A solid boundary with no-slip and a user specified sliding temperature range.

    Like the AdiabaticBC, this is completey effective only when viscous
    effects are active.  Else, it is just like another solid (slip) wall.
    """
    def __init__(self, Twall_i, Twall_f, t_i, t_f, label=""):
        """
        Construct a no-slip, sliding-temperature, solid-wall boundary.

        :param Twall_i:
        :param Twall_f:
        :param t_i:
        :param t_f:
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=SLIDING_T, 
            Twall_i=Twall_i, Twall_f=Twall_f, 
            t_i=t_i, t_f=t_f, is_wall=1, label=label)
        return
    def __str__(self):
        return "SlidingTBC(%g, %g, %g, %g, label=\"%s\")" % \
            (self.Twall_i, self.Twall_f, self.t_i, self.t_f, self.label)
    def __copy__(self):
        return SlidingTBC(Twall_i=self.Twall_i, Twall_f=self.Twall_f, 
            t_i=self.t_i, t_f=self.t_f, label=self.label)

class SubsonicInBC(BoundaryCondition):
    """
    Apply a possibly subsonic inflow condition to the boundary.

    assume_ideal==0: use generalized stepping, down from stagnation to get conditions.
    """
    def __init__(self, inflow_condition, mass_flux=0.0, relax_factor=0.05, 
                 direction_type="normal", direction_vector=[1.0,0.0,0.0],
                 direction_alpha=0.0, direction_beta=0.0,
                 assume_ideal=False, label=""):
        """
        Construct a subsonic-inflow boundary.

        :param inflow_condition: Refers to a FlowCondition object that represents
            the total conditions for flow at the boundary.
        :param mass_flux: required inflow mass-flux per unit area (in kg/s/m**2)
            Set to 0.0 (default) if you don't wish to specify a value.
        :param relax_factor: under-relaxation is advised.
        :param direction_type: "normal" (default) is to have the inflow velocity
            locally-normal to the boundary.
            "uniform" has the inflow velocity aligned with direction_vector
            "radial" radial-inflow (turbine) through a cylindrical surface
            with flow angles direction_alpha and direction_beta.
            "axial" axial-flow (turbine) through a circular surface 
            with flow angles direction_alpha and direction_beta.
        :param direction_vector: List of x,y,z-components specifying the direction that the
            inflow velocity follows if direction_type=="uniform"
        :param direction_alpha: flow angle splitting the in-plane velocity into radial
            and tangential components for direction_type=="radial" or "axial" (radians)
        :param direction_beta: flow angle determining axial-velocity component
            for direction_type=="radial" or "axial" (radians)
        :param assume_ideal: (not working) A value of True allows the code to use 
            ideal gas relations to get an estimate of the ghost-cell conditions.
            A value of False, causes the code to compute the ghost-cell flow conditions
            from the total-conditions by making a number of finite steps through
            the isentropic expansion while allowing a general equation of state.
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.

        The flow is assumed to enter the domain in a direction the is locally-normal
        to the boundary.
        """
        BoundaryCondition.__init__(self, type_of_BC=SUBSONIC_IN,
            inflow_condition=inflow_condition, mass_flux=mass_flux, relax_factor=relax_factor,
            direction_type=direction_type, direction_vector=direction_vector,
            direction_alpha=direction_alpha, direction_beta=direction_beta,
            assume_ideal=assume_ideal, label=label)
        return
    def __str__(self):
        return "SubsonicInBC(inflow_condition=%s, mass_flux=%g, relax_factor=%g, " \
            "direction_type=\"%s\", direction_vector=[%g,%g,%g], " \
            "direction_alpha=%g, direction_beta=%g, " \
            "assume_ideal=%s, label=\"%s\")" % \
            (self.inflow_condition, self.mass_flux, self.relax_factor,
             self.direction_type, self.direction_vector[0], self.direction_vector[1],
             self.direction_vector[2], self.direction_alpha, self.direction_beta,
             self.assume_ideal, self.label)
    def __copy__(self):
        return SubsonicInBC(inflow_condition=self.inflow_condition,
                            mass_flux=self.mass_flux, relax_factor=self.relax_factor,
                            direction_type=self.direction_type,
                            direction_vector=self.direction_vector,
                            direction_alpha=self.direction_alpha,
                            direction_beta=self.direction_beta,
                            assume_ideal=self.assume_ideal, label=self.label)

class SubsonicOutBC(BoundaryCondition):
    """
    An outflow boundary which will try to prevent wave reflection.
    
    (Doesn't work so well at present.)
    """
    def __init__(self, sponge_flag=0, label=""):
        BoundaryCondition.__init__(self, type_of_BC=SUBSONIC_OUT,
            sponge_flag=sponge_flag, label=label)
        return
    def __str__(self):
        return "SubsonicOutBC(sponge_flag=%d, label=\"%s\")" % (self.sponge_flag, self.label)
    def __copy__(self):
        return SubsonicOutBC(sponge_flag=self.sponge_flag, label=self.label)

class TransientUniBC(BoundaryCondition):
    """
    Transient but uniform inflow is applied to the ghost cells at the boundary.
    
    The actual flow data is read (at run time) from the file transient_uni.dat.
    Mostly, this boundary condition is used to get L1d3 simulation data to drive
    an axisymmetric simulation here.
    """
    def __init__(self, filename="transient_uniform.dat", label=""):
        """
        Construct a uniform-inflow boundary with user-specified, transient properties.

        :param filename: containing the specified flow conditions that will be
            applied uniformly along the boundary. 
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=TRANSIENT_UNI, filename=filename, label=label)
        return
    def __str__(self):
        return "TransientUniBC(filename=\"%s\", label=\"%s\")" % (self.filename, self.label)
    def __copy__(self):
        return TransientUniBC(filename=self.filename, label=self.label)

class TransientProfBC(BoundaryCondition):
    """
    Transient non-uniform inflow is applied to the ghost cells at the boundary.
    
    The actual flow data is read (at run time) from the specified file.

    Actually, this BC is vapourware.
    """
    def __init__(self, filename="transient_profile.dat", label=""):
        BoundaryCondition.__init__(self, type_of_BC=TRANSIENT_PROF, filename=filename, label=label)
        return
    def __str__(self):
        return "TransientProfBC(filename=\"%s\", label=\"%s\")" % (self.filename, self.label)
    def __copy__(self):
        return TransientProfBC(filename=self.filename, label=self.label)
    
class StaticProfBC(BoundaryCondition):
    """
    Static, non-uniform inflow is applied to the ghost cells at the boundary.
    
    The actual flow data is read (at run time) from the specified file.
    """
    def __init__(self, filename="profile.dat", n_profile=1, label=""):
        """
        Construct an inflow boundary with user-specified properties.

        :param filename: containing the specified flow conditions that may
            vary along the boundary. 
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=STATIC_PROF, filename=filename,
                                   n_profile=n_profile, label=label)
        return
    def __str__(self):
        return "StaticProfBC(filename=\"%s\", n_profile=%d, label=\"%s\")" % \
            (self.filename, self.n_profile, self.label)
    def __copy__(self):
        return StaticProfBC(filename=self.filename, n_profile=self.n_profile, label=self.label)
 
class FixedPOutBC(BoundaryCondition):
    """
    Something like ExtrapolateOutBC but with the pressure set
    to some user-specified value.
    
    It is probably best to set this pressure at the same value as
    the initial fill pressure so that this boundary condition will
    be passive until a wave arrives at the boundary.
    """
    def __init__(self, Pout, x_order=0, label=""):
        """
        Construct an outflow BC that extrapolates the interior flow data but
        specifies pressure directly.

        :param Pout: fixed outside pressure (in Pascals)
        :param x_order: Extrapolation order of the boundary condition.
            0=just copy the nearest cell data into both ghost cells 
            (zero-order extrapolation).
            1=linear extrapolation of the interior data into the ghost cells. 
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=FIXED_P_OUT, Pout=Pout, 
                                   x_order=x_order, label=label)
        return
    def __str__(self):
        return "FixedPOutBC(Pout=%g, x_order=%d, label=\"%s\")" % \
            (self.Pout, self.x_order, self.label)
    def __copy__(self):
        return FixedPOutBC(Pout=self.Pout, x_order=self.x_order, label=self.label)

class RRMBC(BoundaryCondition):
    """
    Andrew Denman's recycled and renormalised boundary condition.
    """
    def __init__(self, sponge_flag=0, label=""):
        BoundaryCondition.__init__(self, type_of_BC=RRM, sponge_flag=sponge_flag, label=label)
        return
    def __str__(self):
        return "RRMBC(sponge_flag=%d, label=\"%s\")" % (self.sponge_flag, self.label)
    def __copy__(self):
        return RRMBC(sponge_flag=self.sponge_flag, label=self.label)

class UserDefinedBC(BoundaryCondition):
    """
    The user defines the flow properties to use the ghost cells via a Lua function.
    
    The actual flow data is computed (at run time) from the functions defined in that file.
    For details, please see the Appendix in the User Guide and Example Book.
    """
    def __init__(self, filename="udf.lua", is_wall=0,
                 sets_conv_flux=0, sets_visc_flux=0, label=""):
        """
        Construct a user-defined boundary condition.

        :param filename: Name of the file containing the Lua functions.
        :param is_wall: Flag to indicate that various parts of the simulation code 
            should treat this boundary as a solid wall.
        :param sets_conv_flux: For this boundary, the fluxes are computed directly.
            Typically, this relates to using a user-supplied Lus script which
            provides a convective_flux() function. This pretty much ignores the ghost-cell
            data, however, it does not relieve the user of supplying a suitable 
            function for setting that data.
        :param sets_visc_flux: As for sets_conv_flux except that this relates to
            setting the viscous component of flux due to the effect of the boundary.
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=USER_DEFINED, filename=filename,
                                   is_wall=is_wall,
                                   sets_conv_flux=sets_conv_flux, sets_visc_flux=sets_visc_flux,
                                   label=label)
        return
    def __str__(self):
        return ("UserDefinedBC(filename=\"%s\", is_wall=%d, sets_conv_flux=%d,  " +
                "sets_visc_flux=%d, label=\"%s\")") % \
            (self.filename, self.iswall, self.sets_conv_flux,
             self.sets_visc_flux, self.label)
    def __copy__(self):
        return UserDefinedBC(filename=self.filename, is_wall=self.is_wall, 
                             sets_conv_flux=self.sets_conv_flux, sets_visc_flux=self.sets_visc_flux,
                             label=self.label)
    
class AdjacentPlusUDFBC(BoundaryCondition):
    """
    This boundary joins (i.e. is adjacent to) a boundary of another block 
    and a user-defined (Lua) function is used.

    Usually, this condition is not set manually but is set as part of the
    connect_blocks() function.
    """
    def __init__(self, other_block=-1, other_face=-1, orientation=0,
                 filename="udf.lua", is_wall=0, sets_conv_flux=0, sets_visc_flux=0, 
                 reorient_vector_quantities=False, Rmatrix=None, 
                 label=""):
        """
        Construct a connecting boundary condition that also has some user-defined behaviour.

        :param other_block: index to an adjacent block, if any. 
            A value of -1 will indicate that there is no adjacent block.
        :param other_face: index of the adjacent face of the other block, if any.
        :param orientation: for 3D connections the other block face can have one of
            4 rotational orientations.
        :param filename: Name of the file containing the Lua functions.
        :param is_wall: Flag to indicate that various parts of the simulation code 
            should treat this boundary as a solid wall.
        :param sets_conv_flux: For this boundary, the fluxes are computed directly.
            Typically, this relates to using a user-supplied Lus script which
            provides a convective_flux() function. This pretty much ignores the ghost-cell
            data, however, it does not relieve the user of supplying a suitable 
            function for setting that data.
        :param sets_visc_flux: As for sets_conv_flux except that this relates to
            setting the viscous component of flux due to the effect of the boundary.
        :param reorient_vector_quantities: for exchange of vector quantities between adjacent boundaries
        :param Rmatrix: the 9 elements of the rotation matrix
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=ADJACENT_PLUS_UDF, other_block=other_block,
                                   other_face=other_face, orientation=orientation,
                                   filename=filename, is_wall=is_wall,
                                   sets_conv_flux=sets_conv_flux, sets_visc_flux=sets_visc_flux, 
                                   reorient_vector_quantities=reorient_vector_quantities,
                                   Rmatrix=Rmatrix,
                                   label=label)
        return
    def __str__(self):
        return ("AdjacentPlusUDFBC(other_block=%d, other_face=%d, orientation=%d, " + 
                "filename=\"%s\", is_wall=%d, sets_conv_flux=%d, sets_visc_flux=%d, " +
                "reorient_vector_quantities=%d, Rmatrix=%s, " +
                "label=\"%s\")" %
                (self.other_block, self.other_face, self.orientation, 
                 self.filename, self.is_wall, self.sets_conv_flux, self.sets_visc_flux,
                 self.reorient_vector_quantities, self.Rmatrix,
                 self.label))
    def __copy__(self):
        return AdjacentPlusUDFBC(other_block=self.other_block, other_face=self.other_face,
                                 orientation=self.orientation, filename=self.filename,
                                 is_wall=self.is_wall,
                                 sets_conv_flux=self.sets_conv_flux, sets_visc_flux=self.sets_visc_flux,
                                 reorient_vector_quantities=self.reorient_vector_quantities,
                                 Rmatrix=self.Rmatrix,
                                 label=self.label)
     
class SurfaceEnergyBalanceBC(BoundaryCondition):
    """
    Apply the surface energy balance boundary condition, for radiating flows.
    """
    def __init__(self, emissivity, label=""):
        BoundaryCondition.__init__(self, type_of_BC=SEB, is_wall=1, 
                                   emissivity=emissivity, label=label)
        return
    def __str__(self):
        return "SurfaceEnergyBalanceBC(emissivity=%g, label=\"%s\")" % (self.emissivity, self.label)
    def __copy__(self):
        return SurfaceEnergyBalanceBC(emissivity=self.emissivity, label=self.label)

class AblatingBC(BoundaryCondition):
    """
    A reacting wall boundary where the flowfield interacts with the surface
    char layer (carbon), and injects mass into the flowfield.
    
    At present, the temperature is set by the user (c.f. fixedTBC) and a
    reaction scheme for the surface reactions ONLY is an input.

    Like the AdiabaticBC, this is completey effective only when viscous
    effects are active.  Else, it is just like another solid (slip) wall.
    """
    def __init__(self, Twall, filename="reaction-scheme.lua", emissivity=1.0, label=""):
        BoundaryCondition.__init__(self, type_of_BC=ABLATING, Twall=Twall, 
                                   filename=filename, is_wall=1, 
                                   emissivity=emissivity, label=label)
        return
    def __str__(self):
        return "AblatingBC(Twall=%g, label=\"%s\")" % (self.Twall, self.label)
    def __copy__(self):
        return AblatingBC(Twall=self.Twall, filename=self.filename,
                          emissivity=self.emissivity, label=self.label)

class UserDefinedMassFluxBC(BoundaryCondition):
    """
    The user defines the diffusive flux of species mass.
    This is intended to model the effect of species production/destruction
    at a surface with finite-rate chemisitry.
    """
    def __init__(self, filename="udf.lua", label=""):
        """
        Construct a user-defined boundary condition.

        :param filename: Name of the file containing the Lua functions.
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=USER_DEFINED_MASS_FLUX,
                                   filename=filename,
                                   label=label)
        return
    def __str__(self):
        return ("UserDefinedMassFluxBC(filename=\"%s\", label=\"%s\")") % \
            (self.filename, self.label)
    def __copy__(self):
        return UserDefinedMassFluxBC(filename=self.filename, label=self.label)

class ConjugateHTBC(BoundaryCondition):
    """
    A wall with conjugate heat transfer coupling in wall and flow domain.

    This boundary condition is only effective when a wall model for the
    conjugate heat transfer at the boundary between fluid and solid
    is active. See gdata.conjugate_ht_active and gdata.conjugate_ht_file.
    Aside from the heat transfer, this wall acts as a no-slip wall.
    """
    def __init__(self, emissivity=1.0, label=""):
        """
        Construct a conjugate heat transfer, solid-wall boundary.

        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=CONJUGATE_HT, is_wall=1,
                                   emissivity=emissivity, label=label)
        return
    def __str__(self):
        return "ConjugateHTBC(label=\"%s\")" % self.label
    def __copy__(self):
        return ConjugateHTBC(emissivity=self.emissivity,label=self.label)

class MovingWallBC(BoundaryCondition):
    """
    A solid boundary with no-slip but non-zero surface velocity.

    Like the AdiabaticBC, this is completey effective only when viscous
    effects are active.  Else, it is just like another solid (slip) wall.
    """
    def __init__(self, r_omega=None, centre=None, v_trans=None, Twall_flag=False, Twall=None, label=""):
        """
        Construct a no-slip, solid-wall boundary that has a non-zero surface velocity.

        :param r_omega: angular velocity vector (list or Vector) in rad/s
        :param centre: point on axis of rotation (list or Vector)
        :param v_trans: translational velocity to superimpose (list or Vector3)
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        import numpy
        if r_omega is None:
            my_r_omega = [0.0, 0.0, 0.0]
        elif type(r_omega) is list and len(r_omega) >= 3:
            my_r_omega = [r_omega[0], r_omega[1], r_omega[2]]
        elif type(r_omega) is numpy.ndarray and r_omega.shape == (3,):
            my_r_omega = [r_omega[0], r_omega[1], r_omega[2]]
        elif type(r_omega) is Vector:
            my_r_omega = [r_omega.x, r_omega.y, r_omega.z]
        else:
            raise RuntimeError("Invalid input for r_omega: " + str(r_omega))
        if centre is None:
            my_centre = [0.0, 0.0, 0.0]
        elif type(centre) is list and len(centre) >= 3:
            my_centre = [centre[0], centre[1], centre[2]]
        elif type(centre) is numpy.ndarray and centre.shape == (3,):
            my_centre = [centre[0], centre[1], centre[2]]
        elif type(centre) is Vector:
            my_centre = [centre.x, centre.y, centre.z]
        else:
            raise RuntimeError("Invalid input for centre: " + str(centre))
        if v_trans is None:
            my_v_trans = [0.0, 0.0, 0.0]
        elif type(v_trans) is list and len(v_trans) >= 3:
            my_v_trans = [v_trans[0], v_trans[1], v_trans[2]]
        elif type(v_trans) is numpy.ndarray and v_trans.shape == (3,):
            my_v_trans = [v_trans[0], v_trans[1], v_trans[2]]
        elif type(v_trans) is Vector:
            my_v_trans = [v_trans.x, v_trans.y, v_trans.z]
        else:
            raise RuntimeError("Invalid input for v_trans: " + str(v_trans))
        BoundaryCondition.__init__(self, type_of_BC=MOVING_WALL, 
                                   r_omega=my_r_omega, centre=my_centre, v_trans=my_v_trans,
                                   Twall_flag=Twall_flag, Twall=Twall, label=label)
        return
    def __str__(self):
        return "MovingWallBC(r_omega=[%g,%g,%g], centre=[%g,%g,%g], " \
            "v_trans=[%g,%g,%g], Twall_flag=%s, Twall=%g, label=\"%s\")" % \
            (self.r_omega, self.centre, self.Twall_flag, self.Twall, self.label)
    def __copy__(self):
        return MovingWallBC(r_omega=self.r_omega, centre=self.centre, v_trans=self.v_trans, 
                            Twall_flag=self.Twall_flag, Twall=self.Twall, label=self.label)
 
class MassFluxOutBC(BoundaryCondition):
    """
    Something like FixedPOutBC but with the pressure computed to
    achieve some user-specified mass-flux outflow for the boundary.
    
    It is probably best to set the initial pressure at the same value
    as the initial fill pressure so that this boundary condition will
    start changing things gradually.
    """
    def __init__(self, mass_flux, p_init, relax_factor=0.05, label=""):
        """
        Construct an outflow BC that extrapolates the interior flow data but
        applies a computed pressure to achieve a required mass flux.

        :param mass_flux: required outflow mass-flux per unit area (in kg/s/m**2)
        :param p_init: initial outside pressure (in Pascals)
        :param relax_factor: under-relaxation is advised. 
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=MASS_FLUX_OUT,
                                   mass_flux=mass_flux, p_init=p_init, relax_factor=relax_factor, 
                                   label=label)
        return
    def __str__(self):
        return "MassFluxOutBC(mass_flux=%g, p_init=%g, relax_factor=%g, label=\"%s\")" % \
            (self.mass_flux, self.p_init, self.relax_factor, self.label)
    def __copy__(self):
        return MassFluxOutBC(mass_flux=self.mass_flux, p_init=self.p_init,
                             relax_factor=self.relax_factor, label=self.label)
 
class MappedCellBC(BoundaryCondition):
    """
    Ghost-cell flow data is obtained from source-cells.

    The source-cell for each ghost-cell is identified as being the closest one
    to the ghost-cell's mapped position.
    """
    def __init__(self, ghost_cell_trans_fn=lambda x, y, z: (x, y, z),
                 reorient_vector_quantities=False, 
                 Rmatrix=[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                 mapped_cell_list=[],
                 label=""):
        """
        Construct an outflow BC that extrapolates the interior flow data but
        applies a computed pressure to achieve a required mass flux.

        :param ghost_cell_trans_fn: User-supplied transform function mapping from 
            ghost-cell position to source-cell position.
        :param reorient_vector_quantities: for exchange of vector quantities between adjacent boundaries
        :param Rmatrix: the 9 elements of the rotation matrix
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.

        Note that, because the MappedCellBCs store boundary-specific information,
        we must assign a unique MappedCellBC object on each block boundary.
        We cannot reuse the one object for many boundaries without great confusion.
        """
        BoundaryCondition.__init__(self, type_of_BC=MAPPED_CELL,
                                   ghost_cell_trans_fn=ghost_cell_trans_fn,
                                   reorient_vector_quantities=reorient_vector_quantities,
                                   Rmatrix=Rmatrix,
                                   label=label)
        assert callable(ghost_cell_trans_fn)
        self.mapped_cell_list = copy.copy(mapped_cell_list)
        return
    def __str__(self):
        return "MappedCellBC(ghost_cell_trans_fn=%s, reorient_vector_quantities=%s, " \
            "Rmatrix=%s, label=\"%s\")" % \
            (self.ghost_cell_trans_fn, self.reorient_vector_quantities, self.Rmatrix, self.label)
    def __copy__(self):
        return MappedCellBC(ghost_cell_trans_fn=self.ghost_cell_trans_fn, 
                            reorient_vector_quantities=self.reorient_vector_quantities,
                            Rmatrix=copy.copy(self.Rmatrix),
                            mapped_cell_list=copy.copy(self.mapped_cell_list),
                            label=self.label)
    def write_mapped_cells_to_file(self, fileName):
        """
        Write the ghost-cell to source-cell mapping
        in a format suitable for the main simulation program.
        """
        fp = open(fileName, "w")
        fp.write("# blkId   i   j   k\n")
        for c in self.mapped_cell_list:
            fp.write("%d %d %d %d\n" % c)
        fp.close()
        return


#####################################################################################
# FIX-ME -- should we merge the catalycity bcs with the main boundary-condition list?
#####################################################################################

"""
Dictionary to look up wall catalycity boundary-condition index from name or number.

Boundary conditions are implemented within the simulation by setting
flow data in ghost cells to suitable values.
This is done once per time step, before evaluating the fluxes.:

* NON_CATALYTIC: The wall has no chemical influence on the
    neighbouring gas particles.
* PARTIALLY_CATALYTIC: The wall induces some finite recombination
    on the neighbouring gas particles.
* CATALYTIC: The wall induces complete chemical equilibrium on the
    neighbouring gas particles a the local temperature (T_wall) and pressure.
* FULLY_CATALYTIC: A synonym for CATALYTIC
* SUPER_CATALYTIC: The wall causes the gas to rcombine to freestream values.
"""

wc_bcNames = ['non-catalytic', 'equil-catalytic', 'super-catalytic', 'partially-catalytic' ]

NON_CATALYTIC        = 22
EQUIL_CATALYTIC      = 23
SUPER_CATALYTIC      = 24
PARTIALLY_CATALYTIC  = 25

wc_bcIndexFromName = {
    22: NON_CATALYTIC, "22": NON_CATALYTIC, "NON_CATALYTIC": NON_CATALYTIC,
    23: EQUIL_CATALYTIC, "23": EQUIL_CATALYTIC, "EQUIL_CATALYTIC": EQUIL_CATALYTIC,
    24: SUPER_CATALYTIC, "24": SUPER_CATALYTIC, "SUPER_CATALYTIC": SUPER_CATALYTIC,
    25: PARTIALLY_CATALYTIC, "25": PARTIALLY_CATALYTIC, "PARTIALLY_CATALYTIC": PARTIALLY_CATALYTIC
    }


class WallCatalycityBoundaryCondition(object):
    """
    Base class for wall catalycity boundary condition specifications.

    type_of_WCBC: specifies the boundary condition
    Twall: fixed wall temperature (in degrees K) that will be used if
        the boundary conditions needs such a value.
    f_wall: fixed mass fractions at the wall that will be used if
        the boundary conditions needs such a value.
    input_file: name of input file if required for implementation of
        of the specified boundray condition
    """
    __slots__ = 'type_of_WCBC', 'f_wall', 'input_file'
    def __init__(self,
                 type_of_WCBC=NON_CATALYTIC,
                 f_wall=[1.0, ],
                 input_file=None,
                 label=""):
        self.type_of_WCBC = type_of_WCBC
        self.f_wall = copy.copy(f_wall)
        self.input_file = input_file
        self.label = label
        return
    def __str__(self):
        str = "WallCatalyticBoundaryCondition(%d, \n" % (self.type_of_WCBC)
        str += "  [ "
        for val in self.f_wall:
            str += "%g, " % val
        str += "]\n"
        str += ", input_file=\"%s\", label=\"%s\" )\n" % (self.input_file, self.label)
        return str
    def __copy__(self):
        return WallCatalycityBoundaryCondition(type_of_WCBC=self.type_of_WCBC,
                                               f_wall=self.f_wall,
                                               input_file=self.input_file,
                                               label=self.label)

class NonCatalyticWBC(WallCatalycityBoundaryCondition):
    """
    A non-catalytic wall boundary condition which has no effect on the
    chemical composition near the wall.
    """
    def __init__(self, label=""):
        WallCatalycityBoundaryCondition.__init__(self, type_of_WCBC=NON_CATALYTIC, label=label)
        return
    def __str__(self):
        return "NonCatalyticWBC(label=\"%s\")" % (self.label)
    def __copy__(self):
        return NonCatalyticWBC(label=self.label)


class PartiallyCatalyticWBC(WallCatalycityBoundaryCondition):
    """
    A partially catalytic wall boundary condition sets the chemical composition
    at the wall corresponding to a finite rate recombination coefficient, and a
    set of finite-rate reactions. The local temperature is fixed at the wall
    temperature and this should be used to evaluate the rate constant for the
    specified set of chemical reactions at the wall. These reactions can be set
    in a lua file, similar to setting the reactions for the overall flowfield.
    """
    def __init__(self, input_file, label=""):
        WallCatalycityBoundaryCondition.__init__(self, type_of_WCBC=PARTIALLY_CATALYTIC,
                                                 input_file=input_file, label=label)
        return
    def __str__(self):
        return "PartiallyCatalyticWBC(input_file=\"%s\", label=\"%s\")" % \
            (self.input_file, self.label)
    def __copy__(self):
        return PartiallyCatalyticWBC(self.input_file,self.label)


class EquilCatalyticWBC(WallCatalycityBoundaryCondition):
    """
    An equilibrium catalytic wall boundary condition sets the chemical composition
    at the wall in chemical equilibrium at the local temperature
    and pressure.  The local temperature is fixed at the wall temperature
    and this value should be used to create the LUT of chemical compositions.
    The chemical composition at the wall changes with pressure.  This composition
    is found from a pre-computed look-up table at various pressures.
    """
    def __init__(self, input_file, label=""):
        WallCatalycityBoundaryCondition.__init__(self, type_of_WCBC=EQUIL_CATALYTIC,
                                                 input_file=input_file, label=label)
        return
    def __str__(self):
        return "EquilCatalyticWBC(input_file=\"%s\", label=\"%s\")" % (self.input_file, self.label)
    def __copy__(self):
        return EquilCatalyticWBC()

class SuperCatalyticWBC(WallCatalycityBoundaryCondition):
    """
    A super-catalytic wall boundary condition sets the chemical composition
    at the wall to the freestream chemical composition.
    """
    def __init__(self, f_wall, label=""):
        WallCatalycityBoundaryCondition.__init__(self, type_of_WCBC=SUPER_CATALYTIC,
                                                 f_wall=f_wall, label=label)
        return
    def __str__(self):
        str = "SuperCatalyticWBC(%d\n" % (self.type_of_WCBC)
        str += " fwall=[ "
        for val in self.f_wall:
            str += "%g, " % val
        str += "], label=\"%s\" )\n" % self.label
        return str
    def __copy__(self):
        return SuperCatalyticWBC(f_wall=self.f_wall, label=self.label)
        




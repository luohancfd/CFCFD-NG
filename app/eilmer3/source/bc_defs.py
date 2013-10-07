## \file cns_bc_defs.py
## \ingroup mb_cns
##
## \author P.Jacobs
## \version 31-Jan-2005 extracted from e_model_spec.py
##
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
    and temperature have been specified and the velocity from
    just inside the boundary is copied into the ghost cells.
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
#
# The integer values in the following dictionary are a reminder of the old
# macro definitions in the C code.  They are retained here for fererence.
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
    SHOCK_FITTING_IN: "SHOCK_FITTING_IN"
    }

class BoundaryCondition(object):
    """
    Base class for boundary condition specifications.
    """
    __slots__ = 'type_of_BC', 'Twall', 'Pout', 'inflow_condition', \
                'x_order', 'sponge_flag', 'other_block', 'other_face', 'orientation', \
                'filename', 'n_profile', 'is_wall', 'sets_conv_flux', 'sets_visc_flux', 'assume_ideal', \
                'mdot', 'Twall_i', 'Twall_f', 't_i', 't_f', 'emissivity', 'label'
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
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        self.type_of_BC = type_of_BC
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
                                 label=self.label)
    
class AdjacentBC(BoundaryCondition):
    """
    This boundary joins (i.e. is adjacent to) a boundary of another block.
    """
    def __init__(self, other_block=-1, other_face=-1, orientation=0, label=""):
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
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=ADJACENT, other_block=other_block,
                                   other_face=other_face, orientation=orientation,
                                   label=label)
        return
    def __str__(self):
        return "AdjacentBC(other_block=%d, other_face=%d, orientation=%d, label=\"%s\")" % \
            (self.other_block, self.other_face, self.orientation, self.label)
    def __copy__(self):
        return AdjacentBC(other_block=self.other_block,
                          other_face=self.other_face,
                          orientation=self.orientation,
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
    def __init__(self, inflow_condition, assume_ideal=0, label=""):
        """
        Construct a subsonic-inflow boundary.

        :param inflow_condition: Refers to a FlowCondition object that represents
            the total conditions for flow at the boundary.
        :param assume_ideal: A value of 1 allows the code to use ideal gas relations
            to get an estimate of the ghost-cell conditions.
            A value of 0, causes the code to compute the ghost-cell flow conditions
            from the total-conditions by making a number of finite steps through
            the isentropic expansion while allowing a general equation of state.
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=SUBSONIC_IN,
            inflow_condition=inflow_condition, assume_ideal=assume_ideal, label=label)
        return
    def __str__(self):
        return "SubsonicInBC(inflow_condition=%s, assume_ideal=%d, label=\"%s\")" % \
            (self.inflow_condition, self.assume_ideal, self.label)
    def __copy__(self):
        return SubsonicInBC(inflow_condition=self.inflow_condition, 
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
        :param x_order: Extrapolation order of the boundary conduition.
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
                 filename="udf.lua", is_wall=0, sets_conv_flux=0,
                 sets_visc_flux=0, label=""):
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
        :param label: A string that may be used to assist in identifying the boundary
            in the post-processing phase of a simulation.
        """
        BoundaryCondition.__init__(self, type_of_BC=ADJACENT_PLUS_UDF, other_block=other_block,
                                   other_face=other_face, orientation=orientation,
                                   filename=filename, is_wall=is_wall,
                                   sets_conv_flux=sets_conv_flux, sets_visc_flux=sets_visc_flux, 
                                   label=label)
        return
    def __str__(self):
        return ("AdjacentPlusUDFBC(other_block=%d, other_face=%d, orientation=%d, " + 
                "filename=\"%s\", is_wall=%d, sets_conv_flux=%d, sets_visc_flux=%d, label=\"%s\")") % \
               (self.other_block, self.other_face, self.orientation, 
                self.filename, self.is_wall, self.sets_conv_flux, self.sets_visc_flux,
                self.label)
    def __copy__(self):
        return AdjacentPlusUDFBC(other_block=self.other_block, other_face=self.other_face,
                                 orientation=self.orientation, filename=self.filename,
                                 is_wall=self.is_wall,
                                 sets_conv_flux=self.sets_conv_flux, sets_visc_flux=self.sets_visc_flux,
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

import copy

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
        




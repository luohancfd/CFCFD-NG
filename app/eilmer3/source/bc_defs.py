## \file cns_bc_defs.py
## \ingroup mb_cns
##
## \author P.Jacobs
## \version 31-Jan-2005 extracted from e_model_spec.py
##
"""
Dictionary to look up boundary-condition index from name or number.

Boundary conditions are implemented within the simulation by setting
flow data in ghost cells to suitable values.
This is done once per time step, before evaluating the fluxes.

@var ADJACENT: This boundary joins that of another block.
    Normally, this boundary condition would be set implicitly
    when making block connections.
@type ADJACENT: int
@var COMMON: Synonym for L{ADJACENT}.
@type COMMON: int
@var SUP_IN: Fully-prescribed inflow (e.g. supersonic inflow).
@type SUP_IN: int
@var SUP_OUT: Synonym for L{EXTRAPOLATE_OUT}.
@type SUP_OUT: int
@var EXTRAPOLATE_OUT: Extrapolate all flow properties from
   just inside the boundary into the ghost-cells outside
   the boundary.  This works fine for a strong supersonic outflow.
@type EXTRAPOLATE_OUT: int
@var SLIP: Synonym for L{SLIP_WALL}
@type SLIP: int
@var SLIP_WALL: A solid but inviscid wall.
    Effectively, this boundary condition copies and reflects the
    properties just inside the boundary into the ghost cells.
@type SLIP_WALL: int
@var ADIABATIC: A solid, no-slip wall without heat transfer.
    (i.e. the near-wall temperature is reflected in the ghost cells)
@type ADIABATIC: int
@var FIXED_T: A solid, no-slip wall with a user-specified
    temperature.
@type FIXED_T: int
@var SLIDING_T: A solid, no-slip wall with a sliding user-specified
    temperature that linearly varies with time.
@type SLIDING_T: int
@var SUBSONIC_IN: An inflow boundary for which the total pressure
    and temperature have been specified and the velocity from
    just inside the boundary is copied into the ghost cells.
@type SUBSONIC_IN: int
@var SUBSONIC_OUT: An outflow boundary which will try to prevent
    wave reflection at the boundary in the presence of subsonic flow.
    (Doesn't work so well at present.)
@type SUBSONIC_OUT: int
@var SUB_OUT: Synonym for L{SUBSONIC_OUT}
@type SUB_OUT: int
@var TRANSIENT_UNI: An transient inflow boundary which has
    a uniform flow condition applied across the full boundary.
@type TRANSIENT_UNI: int
@undocumented: TRANSIENT_PROF
@var STATIC_PROF: A steady inflow boundary with a variable set of
    flow conditions across the boundary.
@type STATIC_PROF: int
@var FIXED_P_OUT: Something like L{EXTRAPOLATE_OUT} but with the
    pressure set to some user-specified value.
    It is probably best to set this pressure at the same value as
    the initial fill pressure so that this boundary condition will
    be passive until a wave arrives at the boundary.
@type FIXED_P_OUT: int
@var RRM: Andrew Denman's recycled and renormalised boundary
    condition.
@type RRM: int
@var SEB: Surface energy balance.
@type SEB: int
@undocumented: SPECIAL

@undocumented: bcIndexFromName
"""
ADJACENT        = 0
COMMON          = 0
SUP_IN          = 1
SUP_OUT         = 2
EXTRAPOLATE_OUT = 2
SLIP            = 3
SLIP_WALL       = 3
ADIABATIC       = 4
FIXED_T         = 5
SUBSONIC_IN     = 6
SUBSONIC_OUT    = 7
SUB_OUT         = 7
TRANSIENT_UNI   = 8
TRANSIENT_PROF  = 9
STATIC_PROF     = 10
FIXED_P_OUT     = 11
RRM             = 12
TRANSIENT_T_WALL = 13
# EULER_MANUFACTURED was 14
# now deprecated: user-defined B.C
# is used instead.
SEB             = 15
USER_DEFINED    = 16
ADJACENT_PLUS_UDF = 17
ABLATING        = 18
SLIDING_T       = 19
FSTC            = 20
SPECIAL         = -1
bcIndexFromName = {
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
    12: RRM, "12": RRM,  "RRM": RRM,
    13: TRANSIENT_T_WALL, "13" : TRANSIENT_T_WALL, "TRANSIENT_T_WALL": TRANSIENT_T_WALL,
    15: SEB, "15" : SEB, "SEB" : SEB,
    16: USER_DEFINED, "16": USER_DEFINED, "USER_DEFINED": USER_DEFINED,
    17: ADJACENT_PLUS_UDF, "17": ADJACENT_PLUS_UDF, "ADJACENT_PLUS_UDF": ADJACENT_PLUS_UDF,
    18: ABLATING, "18" : ABLATING, "ABLATING": ABLATING,
    19: SLIDING_T, "19" : SLIDING_T, "SLIDING_T": SLIDING_T,
    20: FSTC, "20" : FSTC, "FSTC": FSTC,
    -1: SPECIAL, "-1": SPECIAL,  "SPECIAL": SPECIAL,
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
    SEB: "SEB",
    USER_DEFINED: "USER_DEFINED",
    ADJACENT_PLUS_UDF: "ADJACENT_PLUS_UDF",
    ABLATING: "ABLATING",
    SLIDING_T: "SLIDING_T",
    SPECIAL: "SPECIAL",
    FSTC: "FSTC"
    }

class BoundaryCondition(object):
    """
    Base class for boundary condition specifications.

    @ivar type_of_BC: specifies the boundary condition
    @ivar Twall: fixed wall temperature (in degrees K) that will be used if
        the boundary conditions needs such a value.
    @ivar Pout: fixed outside pressure (in Pascals) that will be used if
        the boundary conditions needs such a value.
    @ivar inflow_condition: the flow condition that will be applied if the
        specified boundary condition needs it
    @ivar sponge_flag: A value of 1 will activate Andrew Denman's damping
        terms near the boundary.
    @type sponge_flag: int
    """
    __slots__ = 'type_of_BC', 'Twall', 'Pout', 'inflow_condition', \
                'x_order', 'sponge_flag', 'other_block', 'other_face', 'orientation', \
                'filename', 'n_profile', 'is_wall', 'use_udf_flux', 'assume_ideal', \
                'mdot', 'epsilon', 'Twall_i', 'Twall_f', 't_i', 't_f', 'label'
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
                 use_udf_flux=0,
                 assume_ideal=0,
                 mdot = [],
                 epsilon=0.9,
                 Twall_i=300.0,
                 Twall_f=300.0,
                 t_i=0.0,
                 t_f=0.0,
                 label=""):
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
        self.use_udf_flux = use_udf_flux
        self.assume_ideal = assume_ideal
        self.mdot = mdot
        self.epsilon = epsilon
        self.Twall_i = Twall_i
        self.Twall_f = Twall_f
        self.t_i = t_i
        self.t_f = t_f
        self.label = label
            
        return
    def __str__(self):
        str_rep = "BoundaryCondition("
        str_rep += "type_of_BC=%d" % self.type_of_BC
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
        str_rep += ", use_udf_flux=%d" % self.use_udf_flux
        str_rep += ", assume_ideal=%d" % self.assume_ideal
        str_rep += ", mdot=["
        for mdi in mdot: str_rep += "%g," % mdi
        str_rep += "]"
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
                                 use_udf_flux=self.use_udf_flux,
                                 assume_ideal=self.assume_ideal,
                                 mdot=copy.copy(self.mdot),
                                 Twall_i=self.Twall_i,
                                 Twall_f=self.Twall_f,
                                 t_i=self.t_i,
                                 t_f=self.t_f,
                                 label=self.label)
    
class AdjacentBC(BoundaryCondition):
    """
    This boundary joins (i.e. is adjacent to) a boundary of another block.

    This condition is usually not set manually but is set as part of the
    connect_blocks() function.
    """
    def __init__(self, other_block=-1, other_face=-1, orientation=0, label=""):
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

class SlipWallBC(BoundaryCondition):
    """
    An inviscid-flow solid-boundary.

    Effectively, this boundary condition copies and reflects the
    properties just inside the boundary into the ghost cells.
    """
    def __init__(self, label=""):
        BoundaryCondition.__init__(self, type_of_BC=SLIP_WALL, is_wall=1, label=label)
        return
    def __str__(self):
        return "SlipWallBC(label=\"%s\")" % self.label
    def __copy__(self):
        return SlipWallBC(label=self.label)
    
class AdiabaticBC(BoundaryCondition):
    """
    A solid, no-slip wall without heat transfer.
    
    The near-wall temperature is reflected in the ghost cells.
    This BC is really only effective if viscous effects are active
    else it acts as another solid (slip) wall.
    """
    def __init__(self, label=""):
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
    def __init__(self, Twall, label=""):
        BoundaryCondition.__init__(self, type_of_BC=FIXED_T, Twall=Twall, is_wall=1, label=label)
        return
    def __str__(self):
        return "FixedTBC(Twall=%g, label=\"%s\")" % (self.Twall, self.label)
    def __copy__(self):
        return FixedTBC(Twall=self.Twall, label=self.label)

class fstcBC(BoundaryCondition):
    """
    A solid boundary with no-slip and a user specified temperature.

    Like the AdiabaticBC, this is completey effective only when viscous
    effects are active.  Else, it is just like another solid (slip) wall.
    """
    def __init__(self, filename="fstc_temperatures.dat", label=""):
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
    
    The actual flow data is read (at run time) from the file transient_uni.dat
    """
    def __init__(self, filename="transient_uniform.dat", label=""):
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
    Something like L{ExtrapolateOutBC} but with the pressure set
    to some user-specified value.
    
    It is probably best to set this pressure at the same value as
    the initial fill pressure so that this boundary condition will
    be passive until a wave arrives at the boundary.
    """
    def __init__(self, Pout, label=""):
        BoundaryCondition.__init__(self, type_of_BC=FIXED_P_OUT, Pout=Pout, x_order=0, label=label)
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
    
    The actual flow data is computed (at run time) from the specified file.
    """
    def __init__(self, filename="udf.lua", is_wall=0, use_udf_flux=0, label=""):
        BoundaryCondition.__init__(self, type_of_BC=USER_DEFINED, filename=filename,
                                   is_wall=is_wall, use_udf_flux=use_udf_flux, label=label)
        return
    def __str__(self):
        return "UserDefinedBC(filename=\"%s\", is_wall=%d, use_udf_flux=%d, label=\"%s\")" % \
            (self.filename, self.iswall, self.use_udf_flux, self.label)
    def __copy__(self):
        return UserDefinedBC(filename=self.filename, is_wall=self.is_wall, 
                             use_udf_flux=self.use_udf_flux, label=self.label)
    
class AdjacentPlusUDFBC(BoundaryCondition):
    """
    This boundary joins (i.e. is adjacent to) a boundary of another block 
    and a user-defined (Lua) function is used.

    This condition is usually not set manually but is set as part of the
    connect_blocks() function.
    """
    def __init__(self, other_block=-1, other_face=-1, orientation=0,
                 filename="udf.lua", is_wall=0, use_udf_flux=0, label=""):
        BoundaryCondition.__init__(self, type_of_BC=ADJACENT_PLUS_UDF, other_block=other_block,
                                   other_face=other_face, orientation=orientation,
                                   filename=filename, is_wall=is_wall, use_udf_flux=use_udf_flux,
                                   label=label)
        return
    def __str__(self):
        return ("AdjacentPlusUDFBC(other_block=%d, other_face=%d, orientation=%d, " + 
                "filename=\"%s\", is_wall=%d, use_udf_flux=%d, label=\"%s\")") % \
               (self.other_block, self.other_face, self.orientation, 
                self.filename, self.is_wall, self.use_udf_flux, self.label)
    def __copy__(self):
        return AdjacentPlusUDFBC(other_block=self.other_block, other_face=self.other_face,
                                 orientation=self.orientation, filename=self.filename,
                                 is_wall=self.is_wall, use_udf_flux=self.use_udf_flux,
                                 label=self.label)
     
class SurfaceEnergyBalanceBC(BoundaryCondition):
    """
    Apply the surface energy balance boundary condition, for radiating flows.
    """
    def __init__(self, epsilon, label=""):
        BoundaryCondition.__init__(self, type_of_BC=SEB, is_wall=1, label=label)
        self.epsilon = epsilon
        return
    def __str__(self):
        return "SurfaceEnergyBalanceBC(epsilon=%g, label=\"%s\")" % (self.epsilon, self.label)
    def __copy__(self):
        return SurfaceEnergyBalanceBC(epsilon=self.epsilon, label=self.label)

class AblatingBC(BoundaryCondition):
    """
    A solid boundary with no-slip and user specified wall mass-flux and temperature.

    Like the AdiabaticBC, this is completey effective only when viscous
    effects are active.  Else, it is just like another solid (slip) wall.
    """
    def __init__(self, Twall, mdot, filename="T_profile.dat", label=""):
        BoundaryCondition.__init__(self, type_of_BC=ABLATING, Twall=Twall, 
                                   filename=filename, is_wall=1, label=label)
        self.mdot = copy.copy(mdot)
        return
    def __str__(self):
    	mdot_str = "  [ "
        for val in self.mdot:
            mdot_str += "%g, " % val
        mdot_str += "]"
        return "AblatingBC(Twall=%g, mdot=%s, label=\"%s\")" % (self.Twall, mdot_str, self.label)
    def __copy__(self):
        return AblatingBC(Twall=self.Twall, mdot=self.mdot, filename=self.filename, label=self.label)

"""
Dictionary to look up wall catalycity boundary-condition index from name or number.

Boundary conditions are implemented within the simulation by setting
flow data in ghost cells to suitable values.
This is done once per time step, before evaluating the fluxes.

@var NON_CATALYTIC: The wall has no chemical influence on the
    neighbouring gas particles.
@type NON_CATALYTIC: int
@var PARTIALLY_CATALYTIC: The wall induces some finite recombination
    on the neighbouring gas particles.
@type PARTIALLY_CATALYTIC: int
@var CATALYTIC: The wall induces complete chemical equilibrium on the
    neighbouring gas particles a the local temperature (T_wall) and pressure.
@type CATALYTIC: int
@var FULLY_CATALYTIC: A synonym for CATALYTIC
@type FULLY_CATALYTIC: int
@var SUPER_CATALYTIC: The wall causes the gas to rcombine to freestream
    values.
@type SUPER_CATALYTIC: int

@undocumented: wc_bcIndexFromName
"""

wc_bcNames = ['non-catalytic', 'equil-catalytic', 'super-catalytic' ]

NON_CATALYTIC        = 21
EQUIL_CATALYTIC      = 22
SUPER_CATALYTIC      = 23

wc_bcIndexFromName = {
    21: NON_CATALYTIC, "21": NON_CATALYTIC, "NON_CATALYTIC": NON_CATALYTIC,
    22: EQUIL_CATALYTIC, "22": EQUIL_CATALYTIC, "EQUIL_CATALYTIC": EQUIL_CATALYTIC,
    23: SUPER_CATALYTIC, "23": SUPER_CATALYTIC, "SUPER_CATALYTIC": SUPER_CATALYTIC
    }

import copy

class WallCatalycityBoundaryCondition(object):
    """
    Base class for wall catalycity boundary condition specifications.

    @ivar type_of_WCBC: specifies the boundary condition
    @ivar Twall: fixed wall temperature (in degrees K) that will be used if
        the boundary conditions needs such a value.
    @ivar f_wall: fixed mass fractions at the wall that will be used if
        the boundary conditions needs such a value.
    @ivar input_file: name of input file if required for implementation of
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

#
# PartiallyCatalyticWBC() --- NOT IMPLEMENTED
#

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



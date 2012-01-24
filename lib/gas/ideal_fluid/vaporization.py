# vaporization.py
#
# Helper functions for defining the fluid's vapour dome.
# PJ, 29-Nov-2009
#
# For an explanation of what's going on, see
# W.C. Reynolds
# Thermodynamic properties in SI
# graphs, tables and computational equations for 40 substances
# Section 2C Saturation Data.

from libfluid import State


def vaporization_jump(fluid, T):
    """
    Returns p_sat and s_sat for the saturated vapour line, plus
    enthalpy-of-vaporization and entropy-of-vaporization.
    """
    st = State()
    p_sat = fluid.p_sat(T)
    dpsatdT = fluid.dpsatdT(T)
    st.p = p_sat
    st.T = T
    fluid.eval_state(st, "pt")
    s_sat = st.s
    v_g = 1.0/st.rho # saturated gas density
    v_f = 1.0/fluid.sat_liquid_density(T)
    # enthalpy from Clapeyron equation (eq 11 in Reynolds' text)
    h_fg = T * (v_g - v_f) * dpsatdT
    s_fg = h_fg / T
    return p_sat, s_sat, h_fg, s_fg


def vaporization_dome_points(fluid, T_low, T_steps=10):
    """
    Computes pairs of T,s points defining the fluid vapour dome.
    """
    # Start at the critical point...
    st = State()
    st.p = fluid.Pc
    st.T = fluid.Tc
    fluid.eval_state(st, "pt")
    Ts_points = [(st.s,st.T),]
    # ...and work down in temperature.
    dT = (fluid.Tc - T_low) / T_steps
    for i in range(T_steps):
        st.T = fluid.Tc - (i+1)*dT
        st.p = fluid.p_sat(st.T)
        fluid.eval_state(st, "pt")
        # points on the saturated vapour line go on the end
        Ts_points.append((st.s,st.T))
        # points on the saturated liquid line on the front
        p_sat, s_sat, h_fg, s_fg = vaporization_jump(fluid, st.T)
        Ts_points.insert(0, (st.s-s_fg,st.T))
    return Ts_points

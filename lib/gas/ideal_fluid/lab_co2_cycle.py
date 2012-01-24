# lab_co2_cycle.py
# Supercritical CO2 Brayton cycle for QEGC lab.
#
# Modelled on Hal's MATLAB code (co2_co2.m)
# PJ Apr,May 2008, Aug 2009

from libfluid import CarbonDioxide, State

Treservoir = 273.0 + 235.0 # hot rock temperature
Preservoir = 25.01e6       # pressure at the bottom of the bore
                           # as computed by Hal's matlab code
Tcooler = 320.0            # ambient air temperature
Pcooler = 8.0e6            # pressure in heat exchanger on surface

fluid = CarbonDioxide()
st = [State() for i in range(16)] # set of CO2 states
eta_c = 0.70; # compressor efficiency
eta_t = 0.70; # turbine efficiency

print "# Supercritical CO2 Brayton Cycle."
print "# (using Hal's pressures and temperatures to anchor cycle)"
print "# Pressure in MPa, T in degrees K, enthalpy in kJ/kg, entropy in kJ/kg.K"
print ""

print "# [0] Inlet to compressor as used by Hal."
st[0].p = Pcooler; st[0].T = Tcooler
fluid.eval_state(st[0], "pt")
print "# [0]", st[0].str()

print "# [1] Outlet from compressor."
st[1].p = 13.39e6; st[1].s = st[0].s # specified
st[1].rho = st[0].rho * 1.5; st[1].T = st[0].T * 1.5 # guess
fluid.eval_state(st[1], "ps" ,1)
print "# [1]", st[1].str()

print "# [11] imperfect compressor: efficiency=", eta_c
st[11].h = 1.0/eta_c * (st[1].h - st[0].h) + st[0].h
st[11].p = st[1].p
st[11].rho = st[1].rho; st[11].T = st[1].T # guess
fluid.eval_state(st[11], "ph" ,1)
print "# [11]", st[11].str()

print "# [2] Bottom of bore (no change in elevation in lab)."
st[2].p = st[1].p; st[2].s = st[0].s # specified
st[2].rho = st[1].rho; st[2].T = st[1].T # must be...
fluid.eval_state(st[2], "ps" ,1)
print "# [2]", st[2].str()

print "# [3] End of hot reservoir, entry to upward bore."
st[3].p = 13.39e6; st[3].T = 444.0 # specified from thermosiphon cycle
st[3].rho = st[2].rho*st[2].T/st[3].T # guess
fluid.eval_state(st[3], "pt" ,1)
print "# [3]", st[3].str()

print "# [4] Exit of bore at surface, inlet to turbine."
st[4].p = 13.39e6; st[4].s = st[3].s # specified, isentropic
st[4].rho = st[3].rho; st[4].T = st[3].T # must be...
fluid.eval_state(st[4], "ps" ,1)
print "# [4]", st[4].str()

print "# [5] Exit of turbine."
st[5].p = Pcooler; st[5].s = st[4].s # specified, isentropic
st[5].rho = st[4].rho/2.0; st[5].T = st[4].T/1.5 # guess
fluid.eval_state(st[5], "ps" ,1)
print "# [5]", st[5].str()

print "# [15] imperfect turbine: efficiency=", eta_t
st[15].h = st[4].h - eta_t * (st[4].h - st[5].h)
st[15].p = Pcooler
st[15].rho = st[5].rho; st[15].T = st[5].T # guess
fluid.eval_state(st[15], "ph" ,1)
print "# [15]", st[15].str()

print "------------------------------"
print "With ideal turbomachinery:"
work_from_turbine = st[4].h - st[5].h;
heat_added = st[3].h - st[2].h;
# heat_added = (st[3].s - st[2].s) * 0.5 * (st[3].T + st[2].T); # Hal's approach.
thermal_eff = work_from_turbine / heat_added;
heat_rejected = st[5].h - st[0].h;
work_by_pump = st[1].h - st[0].h;
cycle_eff = (work_from_turbine - work_by_pump) / heat_added;
print "# work-from-turbine=", work_from_turbine/1000.0, " kJ/kg"
print "# heat-added=", heat_added/1000.0, " kJ/kg"
print "# geothermal-efficiency=", thermal_eff
print "# cycle-thermal-efficiency=", cycle_eff
print "# work-done-by-pump=", work_by_pump/1000.0, " kJ/kg"
print "# work-done-down-well=", (st[2].h-st[1].h)/1000.0, " kJ/kg"
print "# work-done-up-well=", (st[4].h-st[3].h)/1000.0, " kJ/kg"
print "# heat-rejected=", heat_rejected/1000.0, " kJ/kg"

print "------------------------------"
print "With imperfect turbomachinery:"
work_from_turbine = st[4].h - st[15].h;
heat_added = st[3].h - st[11].h;
thermal_eff = work_from_turbine / heat_added;
heat_rejected = st[15].h - st[0].h;
work_by_pump = st[11].h - st[0].h;
cycle_eff = (work_from_turbine - work_by_pump) / heat_added;
print "# work-from-turbine=", work_from_turbine/1000.0, " kJ/kg"
print "# heat-added=", heat_added/1000.0, " kJ/kg"
print "# geothermal-efficiency=", thermal_eff
print "# cycle-thermal-efficiency=", cycle_eff
print "# work-done-by-pump=", work_by_pump/1000.0, " kJ/kg"
print "# work-done-down-well=", (st[2].h-st[1].h)/1000.0, " kJ/kg"
print "# work-done-up-well=", (st[4].h-st[3].h)/1000.0, " kJ/kg"
print "# heat-rejected=", heat_rejected/1000.0, " kJ/kg"

print "# Done."

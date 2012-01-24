# single_stage_cycle.py
# Supercritical CO2 Brayton cycle for QEGC lab.
#
# Modelled on Hal's MATLAB code (co2_co2.m)
# PJ Apr,May 2008, Aug 2009, Nov 2009

from libfluid import CarbonDioxide, State
from vaporization import vaporization_dome_points
from pylab import plot, show, xlabel, ylabel, axis

T_brine = 225.0 + 273.15              # degrees K
T_ambient = 30.0 + 273.15             # degrees K
T_heatex_approach = 8.0               # degrees C or K
Thot  = T_brine - T_heatex_approach   # degree K
Tcold = T_ambient + T_heatex_approach # degree K
Plow  = 8.0e6                         # Pa
if 0:
    print "Hal's original high pressure value."
    Phigh = 22.0e6                    # Pa
    dh_regen = 9.0e3                  # guess for heat regenerated, J/kg
    eta_c = 0.90 # compressor efficiency
    eta_t = 0.90 # turbine efficiency
if 0:
    print "Hal's original high pressure value with eta=0.85."
    Phigh = 22.0e6                    # Pa
    dh_regen = 10.5e3                 # guess for heat regenerated, J/kg
    eta_c = 0.85 # compressor efficiency
    eta_t = 0.85 # turbine efficiency
if 0:
    print "Carlos and Emilie's limit."
    Phigh = Plow * 2.5                # Pa
    dh_regen = 22.4e3                 # guess for heat regenerated, J/kg
    eta_c = 0.85 # compressor efficiency
    eta_t = 0.85 # turbine efficiency
if 1:
    print "Carlos and Emilie's preferred values."
    Phigh = Plow * 1.9                # Pa
    dh_regen = 57.3e3                 # guess for heat regenerated, J/kg
    eta_c = 0.85 # compressor efficiency
    eta_t = 0.85 # turbine efficiency

fluid = CarbonDioxide()
st = [State() for i in range(16)] # set of CO2 states

print "# Supercritical CO2 Brayton Cycle."
print "# (using Hal's pressures and temperatures to anchor cycle)"
print "# Pressure in MPa, T in degrees K, enthalpy in kJ/kg, entropy in kJ/kg.K"
print ""

print "# T_brine=", T_brine
print "# T_ambient=", T_ambient
print "# T_heatex_approach=", T_heatex_approach
print "# Thot=", Thot
print "# Tcold=", Tcold
print "# Phigh=", Phigh
print "# Plow=", Plow
print "# dh_regen=", dh_regen
print "# eta_c=", eta_c
print "# eta_t=", eta_t
print ""

print "# [1] Inlet to compressor."
st[1].p = Plow; st[1].T = Tcold
fluid.eval_state(st[1], "pt")
print "# [1]", st[1].str()

# Accumulate points in a pair of lists so that we may plot them
# on a Ts diagram.
s_cycle = [st[1].s/1000,]
T_cycle = [st[1].T,]
s_key_points = [st[1].s/1000,]
T_key_points = [st[1].T,]

print "# [2] Outlet from compressor."
st[2].p = Phigh; st[2].s = st[1].s # specified
st[2].rho = st[1].rho * 1.5; st[2].T = st[1].T * 1.5 # guess
fluid.eval_state(st[2], "ps", 1)
print "# [2]", st[2].str()

print "# [12] imperfect compressor: efficiency=", eta_c
st[12].h = 1.0/eta_c * (st[2].h - st[1].h) + st[1].h
st[12].p = st[2].p
st[12].rho = st[2].rho; st[12].T = st[2].T # guess
fluid.eval_state(st[12], "ph", 1)
print "# [12]", st[12].str()
s_cycle.append(st[12].s/1000)
T_cycle.append(st[12].T)
s_key_points.append(st[12].s/1000)
T_key_points.append(st[12].T)

print "# [3] Exit of regenerator."
st[3].p = Phigh; st[3].h = st[12].h + dh_regen
st[3].rho = st[12].rho; st[3].T = st[12].T # guess
fluid.eval_state(st[3], "ph", 1)
print "# [3]", st[3].str()
s_cycle.append(st[3].s/1000)
T_cycle.append(st[3].T)
s_key_points.append(st[3].s/1000)
T_key_points.append(st[3].T)

print "# [4] Exit of brine heatex, inlet to turbine."
st[4].p = st[3].p; st[4].T = Thot # specified
st[4].rho = st[3].rho*st[3].T/st[4].T # guess
fluid.eval_state(st[4], "pt", 1)
print "# [4]", st[4].str()
s_key_points.append(st[4].s/1000)
T_key_points.append(st[4].T)

# Step up through heat exchanger to get a better plotted curve
nstep = 20
dT = (st[4].T - st[3].T)/nstep
for i in range(nstep):
    st[0].p = st[3].p; st[0].T = st[3].T + (i+1)*dT
    st[0].rho = st[3].rho * st[3].T/st[0].T # guess
    fluid.eval_state(st[0], "pt", 1)
    s_cycle.append(st[0].s/1000)
    T_cycle.append(st[0].T)

print "# [5] Exit of turbine."
st[5].p = Plow; st[5].s = st[4].s # specified, isentropic
st[5].rho = st[4].rho/2.0; st[5].T = st[4].T/1.5 # guess
fluid.eval_state(st[5], "ps", 1)
print "# [5]", st[5].str()

print "# [15] imperfect turbine: efficiency=", eta_t
st[15].h = st[4].h - eta_t * (st[4].h - st[5].h)
st[15].p = Plow
st[15].rho = st[5].rho; st[15].T = st[5].T # guess
fluid.eval_state(st[15], "ph", 1)
print "# [15]", st[15].str()
s_cycle.append(st[15].s/1000)
T_cycle.append(st[15].T)
s_key_points.append(st[15].s/1000)
T_key_points.append(st[15].T)

print "# [6] After regenerator."
st[6].p = Plow; st[6].h = st[15].h - dh_regen
st[6].rho = st[15].rho; st[6].T = st[15].T # guess
fluid.eval_state(st[6], "ph", 1)
print "# [6]", st[6].str()
s_cycle.append(st[6].s/1000)
T_cycle.append(st[6].T)
s_key_points.append(st[6].s/1000)
T_key_points.append(st[6].T)

print "# temperature difference across regenerator:", \
    st[6].T - st[3].T, "degrees K"

print "# [7] After cooling tower."
st[7].p = Plow; st[7].T = Tcold # suck the rest of the heat out
fluid.eval_state(st[7], "pt" ,1)
print "# [7]", st[7].str()

# Step down through heat exchanger to get a better plotted curve
nstep = 20
dT = (st[7].T - st[6].T)/nstep
for i in range(nstep):
    st[0].p = st[6].p; st[0].T = st[6].T + (i+1)*dT
    st[0].rho = st[6].rho * st[6].T/st[0].T # guess
    fluid.eval_state(st[0], "pt", 1)
    s_cycle.append(st[0].s/1000)
    T_cycle.append(st[0].T)


print "------------------------------"
print "With ideal turbomachinery:"
work_from_turbine = st[4].h - st[5].h
heat_added = st[4].h - st[3].h
heat_rejected = st[6].h - st[7].h
work_by_pump = st[2].h - st[1].h
cycle_eff = (work_from_turbine - work_by_pump) / heat_added
print "# work-from-turbine=", work_from_turbine/1000.0, " kJ/kg"
print "# heat-added=", heat_added/1000.0, " kJ/kg"
print "# cycle-thermal-efficiency=", cycle_eff * 100.0, "%"
print "# work-done-by-pump=", work_by_pump/1000.0, " kJ/kg"
print "# heat-rejected=", heat_rejected/1000.0, " kJ/kg"
print "------------------------------"
print "With non-ideal turbomachinery:"
work_from_turbine = st[4].h - st[15].h
heat_added = st[4].h - st[3].h
heat_rejected = st[6].h - st[7].h
work_by_pump = st[12].h - st[1].h
cycle_eff = (work_from_turbine - work_by_pump) / heat_added
print "# work-from-turbine=", work_from_turbine/1000.0, " kJ/kg"
print "# heat-added=", heat_added/1000.0, " kJ/kg"
print "# cycle-thermal-efficiency=", cycle_eff * 100.0, "%"
print "# work-done-by-pump=", work_by_pump/1000.0, " kJ/kg"
print "# heat-rejected=", heat_rejected/1000.0, " kJ/kg"

Ts_points = vaporization_dome_points(fluid, 220.0, 150)
# print "Ts_points=", Ts_points
s_kJK = [sT[0]/1000.0 for sT in Ts_points]
T_K = [sT[1] for sT in Ts_points]
plot(s_kJK, T_K, '-', s_cycle, T_cycle, '-', s_key_points, T_key_points, 'o')
xlabel('s, kJ/kg.K')
ylabel('T, degree K')
axis([0.8,2.0,250.0,500.0])
show()

print "# Done."

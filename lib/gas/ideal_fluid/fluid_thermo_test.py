#! /usr/bin/env python
# fluid_thermo_test.py
from libfluid import CarbonDioxide, State

print "fluid_thermo_test:" 

co2 = CarbonDioxide()
print co2.name 

print "Pressure from density and temperature:" 
rho = 1.0/0.00779
T = 300.0
dT = 0.01
print "rho=", rho, " T=", T, " p=", co2.p_rhoT(rho,T)
print " expected p=5.0e6" 
print "    dpdT=", co2.dpdT_rhoT(rho,T) \
    , "finite-diff=", 0.5*(co2.p_rhoT(rho,T+dT) - co2.p_rhoT(rho,T-dT))/dT \
    , "rhoR=", rho*co2.R 

rho = 1.0/0.01824
T = 500.0
print "rho=", rho, " T=", T, " p=", co2.p_rhoT(rho,T)
print " expected p=5.0e6" 
print "    dpdT=", co2.dpdT_rhoT(rho,T) \
    , "finite-diff=", 0.5*(co2.p_rhoT(rho,T+dT) - co2.p_rhoT(rho,T-dT))/dT \
    , "rhoR=", rho*co2.R 

rho = 1.0/0.03039
T = 800.0
print "rho=", rho, " T=", T, " p=", co2.p_rhoT(rho,T)
print " expected p=5.0e6" 
print "    dpdT=", co2.dpdT_rhoT(rho,T) \
    , "finite-diff=", 0.5*(co2.p_rhoT(rho,T+dT) - co2.p_rhoT(rho,T-dT))/dT \
    , "rhoR=", rho*co2.R 

rho = 1.0/0.1893
T = 1000.0
print "rho=", rho, " T=", T, " p=", co2.p_rhoT(rho,T)
print " expected p=1.0e6" 
print "    dpdT=", co2.dpdT_rhoT(rho,T) \
    , "finite-diff=", 0.5*(co2.p_rhoT(rho,T+dT) - co2.p_rhoT(rho,T-dT))/dT \
    , "rhoR=", rho*co2.R 

print "Try utility function..." 
st = State()
st.p = 4.9993e6
st.T = 800.0
co2.eval_state(st, "pt")
print "PT: rho=", st.rho, " T=", st.T, " st.p=", st.p 
st.s = 2472.44
st.T = 800.0
co2.eval_state(st, "st")
print "ST: rho=", st.rho, " T=", st.T, " s=", st.s 
st.rho = 32.9056
st.u = 789170.0
co2.eval_state(st, "du")
print "DU: rho=", st.rho, " T=", st.T, " u=", st.u 
st.p = 4.9993e6
st.h = 941098.0
co2.eval_state(st, "ph")
print "PH: rho=", st.rho, " T=", st.T, " p=", st.p, " h=", st.h 
st.p = 4.9993e6
st.s = 2472.44
co2.eval_state(st, "ps")
print "PS: rho=", st.rho, " T=", st.T, " p=", st.p, " s=", st.s 

print "Saturated liquid density:" 
T = 280.0
print "T=", T, " rho_f=", co2.sat_liquid_density(T)
print " expected rho_f=895.48" 

print "Saturated pressure:" 
T = 287.5
rho = 1.0/0.00641
print "T=", T, "p_sat=", co2.p_sat(T)
print " expected p_sat=5.0e6" 

from vaporization import vaporization_jump, vaporization_dome_points
p_sat, s_sat, h_fg, s_fg = vaporization_jump(co2, T)
print "At T=", T, "p_sat=", p_sat, "s_sat=", s_sat, "h_fg=", h_fg, "s_fg=", s_fg

Ts_points = vaporization_dome_points(co2, 220.0)
print "Ts_points=", Ts_points

print "Specific heat, constant volume, ideal gas:" 
T = 287.5
print "T=", T, "cv0=", co2.cv0(T) 

print "Done." 

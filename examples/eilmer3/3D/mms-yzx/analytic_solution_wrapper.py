# analytic_solution_wrapper.py
# 07-Nov-2012: yzx-version
from analytic_solution import AnalyticSolution
import sys

R_air = 287.0

fp = open('case.txt', 'r'); case_str = fp.readline().strip(); fp.close()
case = int(case_str)

if case == 1:
    ev = AnalyticSolution(L=1.0,
         rho0=1.0, rhoy=0.15, rhoz=-0.1, rhoyz=0.0, arhoy=1.0, arhoz=0.5, arhoyz=0.0,
         u0=800.0, uy=50.0, uz=-30.0, uyz=0.0, auy=1.5, auz=0.6, auyz=0.0,
         v0=800.0, vy=-75.0, vz=40.0, vyz=0.0, avy=0.5, avz=2.0/3, avyz=0.0,
         p0=1.0e5, py=0.2e5, pz=0.5e5, pyz=0.0, apy=2.0, apz=1.0, apyz=0.0, case=1)
elif case == 2:
    ev = AnalyticSolution(L=1.0,
         rho0=1.0, rhoy=0.1, rhoz=0.15, rhoyz=0.08, arhoy=0.75, arhoz=1.0, arhoyz=1.25,
         u0=70.0, uy=4.0, uz=-12.0, uyz=7.0, auy=5.0/3, auz=1.5, auyz=0.6,
         v0=90.0, vy=-20.0, vz=4.0, vyz=-11.0, avy=1.5, avz=1.0, avyz=0.9,
         p0=1.0e5, py=-0.3e5, pz=0.2e5, pyz=-0.25e5, apy=1.0, apz=1.25, apyz=0.75, case=2)
elif case == 3:
    ev = AnalyticSolution(L=1.0,
         rho0=1.0, rhoy=0.15, rhoz=-0.1, rhoyz=0.0, arhoy=1.0, arhoz=0.5, arhoyz=0.0,
         u0=800.0, uy=50.0, uz=-30.0, uyz=0.0, auy=1.5, auz=0.6, auyz=0.0,
         v0=800.0, vy=-75.0, vz=40.0, vyz=0.0, avy=0.5, avz=2.0/3, avyz=0.0,
         p0=1.0e5, py=0.2e5, pz=0.5e5, pyz=0.0, apy=2.0, apz=1.0, apyz=0.0, case=3)
elif case == 4:
    ev = AnalyticSolution(L=1.0,
         rho0=1.0, rhoy=0.1, rhoz=0.15, rhoyz=0.08, arhoy=0.75, arhoz=1.0, arhoyz=1.25,
         u0=70.0, uy=4.0, uz=-12.0, uyz=7.0, auy=5.0/3, auz=1.5, auyz=0.6,
         v0=90.0, vy=-20.0, vz=4.0, vyz=-11.0, avy=1.5, avz=1.0, avyz=0.9,
         p0=1.0e5, py=-0.3e5, pz=0.2e5, pyz=-0.25e5, apy=1.0, apz=1.25, apyz=0.75, case=4)
else:
    print "UNKNOWN CASE"
    sys.exit()

def ref_function(x, y, z, t):
    rho = ev.rho(y, z)
    p = ev.p(y, z)
    T = p / (rho*R_air)
    u = ev.u(y, z)
    v = ev.v(y, z)
    return {"rho":rho, "p":p, "T":T, "vel.y":u, "vel.z":v}



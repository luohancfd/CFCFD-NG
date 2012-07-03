# analytic_solution_wrapper.py
from analytic_solution import AnalyticSolution
import sys

R_air = 287.0

fp = open('case.txt', 'r'); case_str = fp.readline().strip(); fp.close()
case = int(case_str)

if case == 1:
    ev = AnalyticSolution(L=1.0,
         rho0=1.0, rhox=0.15, rhoy=-0.1, rhoxy=0.0, arhox=1.0, arhoy=0.5, arhoxy=0.0,
         u0=800.0, ux=50.0, uy=-30.0, uxy=0.0, aux=1.5, auy=0.6, auxy=0.0,
         v0=800.0, vx=-75.0, vy=40.0, vxy=0.0, avx=0.5, avy=2.0/3, avxy=0.0,
         p0=1.0e5, px=0.2e5, py=0.5e5, pxy=0.0, apx=2.0, apy=1.0, apxy=0.0, case=1)
elif case == 2:
    ev = AnalyticSolution(L=1.0,
         rho0=1.0, rhox=0.1, rhoy=0.15, rhoxy=0.08, arhox=0.75, arhoy=1.0, arhoxy=1.25,
         u0=70.0, ux=4.0, uy=-12.0, uxy=7.0, aux=5.0/3, auy=1.5, auxy=0.6,
         v0=90.0, vx=-20.0, vy=4.0, vxy=-11.0, avx=1.5, avy=1.0, avxy=0.9,
         p0=1.0e5, px=-0.3e5, py=0.2e5, pxy=-0.25e5, apx=1.0, apy=1.25, apxy=0.75, case=2)
elif case == 3:
    ev = AnalyticSolution(L=1.0,
         rho0=1.0, rhox=0.15, rhoy=-0.1, rhoxy=0.0, arhox=1.0, arhoy=0.5, arhoxy=0.0,
         u0=800.0, ux=50.0, uy=-30.0, uxy=0.0, aux=1.5, auy=0.6, auxy=0.0,
         v0=800.0, vx=-75.0, vy=40.0, vxy=0.0, avx=0.5, avy=2.0/3, avxy=0.0,
         p0=1.0e5, px=0.2e5, py=0.5e5, pxy=0.0, apx=2.0, apy=1.0, apxy=0.0, case=3)
elif case == 4:
    ev = AnalyticSolution(L=1.0,
         rho0=1.0, rhox=0.1, rhoy=0.15, rhoxy=0.08, arhox=0.75, arhoy=1.0, arhoxy=1.25,
         u0=70.0, ux=4.0, uy=-12.0, uxy=7.0, aux=5.0/3, auy=1.5, auxy=0.6,
         v0=90.0, vx=-20.0, vy=4.0, vxy=-11.0, avx=1.5, avy=1.0, avxy=0.9,
         p0=1.0e5, px=-0.3e5, py=0.2e5, pxy=-0.25e5, apx=1.0, apy=1.25, apxy=0.75, case=4)
else:
    print "UNKNOWN CASE"
    sys.exit()

def ref_function(x, y, z, t):
    rho = ev.rho(x, y)
    p = ev.p(x, y)
    T = p / (rho*R_air)
    u = ev.u(x, y)
    v = ev.v(x, y)
    return {"rho":rho, "p":p, "T":T, "vel.x":u, "vel.y":v}



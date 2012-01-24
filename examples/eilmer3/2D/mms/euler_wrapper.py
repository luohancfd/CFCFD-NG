# euler_wrapper.py
from analytic_solution import AnalyticSolution

R_air = 287.1

ev = AnalyticSolution(L=1.0,
        rho0=1.0, rhox=0.15, rhoy=-0.1, rhoxy=0.0, arhox=1.0, arhoy=0.5, arhoxy=0.0,
        u0=800.0, ux=50.0, uy=-30.0, uxy=0.0, aux=1.5, auy=0.6, auxy=0.0,
        v0=800.0, vx=-75.0, vy=40.0, vxy=0.0, avx=0.5, avy=2.0/3, avxy=0.0,
        p0=1.0e5, px=0.2e5, py=0.5e5, pxy=0.0, apx=2.0, apy=1.0, apxy=0.0)

def ref_function(x, y, z, t):
    rho = ev.rho(x, y)
    p = ev.p(x, y)
    T = p / (rho*R_air)
    u = ev.u(x, y)
    v = ev.v(x, y)
    return {"rho":rho, "p":p, "T":T, "vel.x":u, "vel.y":v}



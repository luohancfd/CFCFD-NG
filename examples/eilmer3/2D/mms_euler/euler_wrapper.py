# euler_wrapper.py
from euler_verify import *

ev = EulerManufacturedSolution( 1.0,
                               1.0, 0.15, -0.1, 1.0, 0.5,
                               1.0e5, 0.2e5, 0.5e5, 2.0, 1.0,
                               800.0, 50.0, -30.0, 1.5, 0.6,
                               800.0, -75.0, 40.0, 0.5, 2.0/3.0)

def ref_function(x, y, z, t):
    rho = ev.calculate_rho(x, y)
    p = ev.calculate_p(x, y)
    T = p / (rho*R_air)
    u = ev.calculate_u(x, y)
    v = ev.calculate_v(x, y)
    return {"rho":rho, "p":p, "T":T, "vel.x":u, "vel.y":v}



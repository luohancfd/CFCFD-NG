# odw_analytical.py
#
# Small script to help the mbcns_verify.py find
# the correct solution function.

from oblique_detonation import *
from math import pi

od = ObliqueDetonation( pi/4.0, 300.0, 3.0, 1.0)

def ref_function(x, y, z, t):
    x1, y1, rho, p, T, f, u, v, X, Y = od.solution(x, y)
    return {"rho":rho, "T[0]":T,
            "vel.x":u, "vel.y":v,
            "massf[0]":f[0], "massf[1]":f[1]}



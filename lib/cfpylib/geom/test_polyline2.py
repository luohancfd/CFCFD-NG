# test_polyline2.py

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from libprep3 import *
a = Vector(0.0,0.0)
b = Vector(0.0,1.0)
c = Vector(1.0,1.0)
d = Vector(1.0,0.0)
ab = Line(a,b, "ab")
bc = Line(b,c, "bc")
cd = Line(c,d, "cd")

e = Polyline2([ab,bc],cd,a)
N = 8
for i in range(N+1):
    t = i * 1.0/N
    print t, e.eval(t)

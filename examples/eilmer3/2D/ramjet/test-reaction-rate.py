# test-reaction-rate.py
from pylab import *

print "Begin..."
def k(T, A=1000000.0, n=0.0, Ta=2500.0):
    return A * T**n * exp(-Ta/T)

ts = arange(100, 1000, 10)
plot(ts, k(ts))
show()
print "Done."

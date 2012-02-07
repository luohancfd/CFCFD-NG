"""
stats.py: Simple statistics for arrays of values.

To replace those in scipy, just in case scipy is not installed.

PJ, 27-Mar-2007
"""

import math
try:
    from Numeric import array
except:
    try:
        from numpy import array
    except:
        print "Could not import Numeric nor numpy."

def mean(a):
    """
    Return the mean value of an array.
    """
    return sum(a)/len(a)

def std(a):
    """
    Return the standard deviation of an array of values.
    """
    if len(a) <= 1: return 0.0
    m = mean(a)
    dev = a - m
    return math.sqrt(sum(dev*dev)/(len(a)-1))

if __name__ == '__main__':
    values = array([1,2.0,3.0])
    print "values=", values
    print "mean=", mean(values)
    print "std-dev=", std(values)
    print "Done."



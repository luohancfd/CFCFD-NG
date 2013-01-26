#!/usr/bin/env python
"""
collect_run_stats.py: Run an executable a number of times and report.

This is a small application program rather than a library function.
Usage::

    $ ${E3BIN}/cfpylib/nm/collect_run_stats.py executable no_times
"""

import sys, os

from time import time
from numpy import zeros

try:
    from scipy.stats import mean, std
except:
    from stats import mean, std
    
def printUsage():
    print "collect_run_stats.py"
    print "Usage:"
    print "collect_run_stats.py executable no_times"
    sys.exit(1)

def main():
    if( len(sys.argv) != 3 ):
        printUsage()

    executable = sys.argv[1]
    no_times = int(sys.argv[2])

    cmd = "time ./%s" % executable

    times = zeros(no_times, 'f');

    for i in range(no_times):
        t = time(); os.system(cmd); t_store = time() - t;
        print "t_store= ", t_store
        times[i] = t_store;

    print "Mean time for operation = ", mean(times)
    print "Standard deviation for operation = ", std(times)


if __name__ == '__main__': main()
        

    

    






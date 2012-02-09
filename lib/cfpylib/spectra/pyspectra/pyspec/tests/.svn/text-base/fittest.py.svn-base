from pylab import *
from pyspec import spec, fit, fitfuncs
import time

def fittest():
    data = loadtxt('testdata.dat')

    times = []
    for optimizer in ['leastsq', 'levmar', 'mpfit']:
        f = fit.fit(x = data[:,0], y = data[:,1], 
                    funcs = [fitfuncs.linear, fitfuncs.gauss],
                    optimizer = optimizer)
        t1 = time.time()
        f.go()
        t2 = time.time()
        times.append('Optimizer %s took %0.3f s' % (optimizer, (t2-t1)))
        
    for t in times:
        print t
        
if __name__ == '__main__':
    fittest()

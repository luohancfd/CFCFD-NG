from pylab import *
from pyspec import spec, fit, fitfuncs
from pyspec.ccd.PrincetonSPE import PrincetonSPEFile
import time

def twodgauss(x, p, mode='eval'):
    if mode == 'eval':
        out = p[0] * exp( -( (x[0] - p[1])**2 / (2 * p[2])**2 + (x[1] - p[3])**2 / (2 * p[4])**2) )
    elif mode == 'params':
        out = ['A', 'CentX', 'SigmaX', 'CentY', 'SigmaY']
    elif mode == 'name':
        out = "2D Gauss"
    return out

def twodlin(x, p, mode='eval'):
    if mode == 'eval':
        out = p[0] + p[1] * x[0] + p[2] * x[1]
    elif mode == 'params':
        out = ['a00', 'a10', 'a01']
    elif mode == 'name':
        out = "2D Lin"
    return out

def test2Dfit(fitType = None, verbose = False):
    testDebug = 0
    testOptimizer = 'levmar'
    roi = [1, 325, 1, 335]
    conIm = PrincetonSPEFile('testimage.spe')[0].astype(float32)
    conZ = array(meshgrid(arange(roi[0], roi[0] + roi[1]), arange(roi[2], roi[2] + roi[3])))
    if fitType == 'twodgauss':    
        guess = array([ conIm.max(), 162.5, 5, 167.5, 6])
        func = [twodgauss]
    elif fitType == 'twodgausslin':
        guess = array([ 0.01*conIm.max(), 1e-8, 1e-8, conIm.max(), 162.5, 5, 167.5, 6])
        func = [twodlin, twodgauss]

    times = []
    for optimizer in ['leastsq', 'levmar', 'mpfit']:
        f = fit.fit(x = conZ, y = conIm, funcs = func , guess = guess, 
                    debug = testDebug, optimizer = optimizer)
        t1 = time.time()
        f.go()
        t2 = time.time()

        times.append('Optimizer %s took %0.3f s' % (optimizer, (t2-t1)))
   
    for t in times:
        print t

if __name__ == "__main__":
    test2Dfit(fitType = 'twodgauss', verbose = False)
    print 'ready'

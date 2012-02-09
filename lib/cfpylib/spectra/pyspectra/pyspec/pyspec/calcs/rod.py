from scipy import *

import scipy.fftpack as fft
from scipy.interpolate import splrep, splev
import numpy.ma as masked
import pylab, numpy

try:
    import psyco
    psyco.profile()
except ImportError:
    pass

class BraggRod(object):
    def __init__(self, hkl, structure, *args, **kwargs):
        """Initialize structure

        Possibe kwargs are:

        'pdepth' : Penetration depth of x-rays 
        'alpha'  : Alpha incidence angle

        """
        self.hkl = hkl
        self.structure = structure
        self.params = kwargs

        # Start by calculating some useful constants

        self.q = self.structure[0].getRLattice() * self.hkl
        
        if not self.params.has_key('pdepth'):
            self.pdepth = self.structure[0].calcPenetrationDepth(kwargs['alpha'])
    
    def go(self):
        """Start Calculation"""
        print "---- Calculating Fctr"
        self._calcFctr()
        print "---- Calculating F"
        self._calcF()
        print "---- Calculating Intensity"
        self._calcIntensity()
        print "---- Correcting Intensity"
        self._calcCorrections()

    def _calcFctr(self):
        """Calculate CTR intensity for given momentum transfer
        
        Calculate CTR Intensity.
        Uses:
        F_CTR = 1 / (1 - exp(-2 * pi * i * l) exp(- c / mu))
        
        Where:
        c    = "c" lattice parameter
        l    = "l" l in reciprocal lattice units
        mu   = penetration depth in angstroms

        """
        c = self.structure[0].getLattice()[2]
        self.Fctr = exp(2 * pi * complex(0, -1) * self.hkl[:,2])
        self.Fctr = self.Fctr * exp(-1.0 * c / self.pdepth)
        self.Fctr = 1. / (1 - self.Fctr)
        
    def _calcF(self, *args, **kwargs):
        """Calculate structure factors for all crystal objects

        This routine calls <crystal>.calcF(q, *args, **kwargs) for 
        all crystal objects to calculate the total amplitude. 

        A = F_0 * F_CTR + \sum_n [ F_n * exp(i*dot(h, q))]
        
        where h is the position from the origin of the new cell 

        Phase: 
        The phase of A is such that after the calculation, the origin of
        direct space is at the top of the crystal stack and the -ve direction
        moves into the crystal stack, i.e. +ve is in the vacuum. 

        """
        # Initialize array
        self.F = zeros((size(self.structure), self.q.shape[0]), numpy.complex)

        height = array([0, 0, self.structure[0].getLattice()[2]]) # c lattice parameter
        self.F[0] = self.structure[0].calcF(self.q, *args, **kwargs)
        for s in range(1, len(self.structure)):
            self.F[s] = self.structure[s].calcF(self.q, *args, **kwargs)
            self.F[s] = self.F[s] * exp(complex(0, 1) * inner(self.q, height))
            height = height + array([0, 0, self.structure[s].getLattice()[2]])  
            
        phasechange = exp(complex(0, -1) * inner(self.q, height))
        self.F = self.F * phasechange

    def _calcIntensity(self):
        self.A = self.F[0] * self.Fctr
        for i in range(1, self.F.shape[0]):
            self.A = self.A + self.F[i]

        # Remove points from forbidden peaks
        if False:
            self.I = (self.A * self.A.conj()).real
            points = where(abs(self.A) < 1e-5)[0]
            for p in points:
                sr = splrep(concatenate((self.q[p-5:p-1][:,2], self.q[p+1:p+5][:,2])),
                            concatenate((self.A[p-5:p-1].real, self.A[p+1:p+5].real)))
                si = splrep(concatenate((self.q[p-5:p-1][:,2], self.q[p+1:p+5][:,2])),
                            concatenate((self.A[p-5:p-1].imag, self.A[p+1:p+5].imag)))
                self.A[p] = complex(splev(self.q[p][2], sr), splev(self.q[p][2], si))
                
        self.I = (self.A * self.A.conj()).real

    def _calcCorrections(self):
        self.footprintCorrection = ones(self.q.shape[0])
        if self.params.has_key('footprint'):
            p = self.params['footprint']
            if p:
                l = self.structure[0].getLambda() * p[0] / (4 * pi * p[1]) 
                for i in range(self.q.shape[0]):
                    qs = sqrt(vdot(self.q[i], self.q[i]))
                    self.footprintCorrection[i] = qs * l
                self.footprintCorrection[self.footprintCorrection > 1] = 1.
                self.I = self.I * self.footprintCorrection
                print "Corrected for footprint"

        self.lorentzCorrection = ones(self.q.shape[0])
        if self.params.has_key('lorentz'):
            if self.params['lorentz']:
                l = self.structure[0].getLambda()
                for i in range(self.q.shape[0]):
                    qs = sqrt(vdot(self.q[i], self.q[i]))
                    theta = arcsin(l * qs / (4 * pi))
                    self.lorentzCorrection[i] = 1 / sin(2 * theta)
                self.I = self.I * self.lorentzCorrection
                print "Corrected for Lorentz factor."

def IBraggRoughness(l, theta):
    """Calculate Bragg component of CTR"""
    I = 1 - (2 * theta * (1 - theta) * (1 - cos(2. * pi * l)))
    return I

def IDiffRoughness(l, theta, L, a1, a2, deltah):
    """Calculate Diffuse component of CTR"""
    I = 2 * theta * (1 - theta) * (1 - cos(2 * pi * l))
    I = I * a2 * L * sqrt(pi) * exp(-1. * pow(pi * L * deltah / a1, 2))
    #I = I * pi * pow(L, 2) * exp(-1. * pow(pi * L ,2)   
    return I 

def gauss(x, c, sigma):
    """Return the gaussian distribution

    y = [ 1 / (sigma * sqrt(2 * pi) ] exp(-(x - c)^2 / 2 * sigma^2)

    """
    a = 1. / (sigma * sqrt(2. * pi))
    y = a * exp(-1. * pow(x - c, 2) / (2. * pow(sigma, 2)))

    return y

def responceGaussFunction(len, sigma):
    """Returns the responce function for integration into an FFT"""

    x1 = arange(0, (len / 2) + 1)
    x2 = arange(-1, -1 * (len / 2), -1)

    x = hstack((x1, x2[ :: -1]))
    print x

    a = exp(-1. * pow(x * sigma, 2))

    return a

def convolve(y, sigma):
    ly = concatenate((y, zeros(len(y) * 20)))

    f = fft.fft(ly)
    r = responceGaussFunction(len(ly), sigma)
    pylab.figure()
    pylab.plot(f*r)
    cy = fft.ifft(f * r)
    pylab.figure()
    pylab.semilogy(abs(cy[:2000]))
    return cy

## Some helper routines

def hklmesh(h = 0., k =0., lmin = 0., lmax = 2, lstep = 0.01):
    """Create hkl mesh for passing to rod routines"""
    c = arange(lmin, lmax, lstep)
    a = ones(len(c)) * h
    b = ones(len(c)) * k
    return vstack((a, b, c)).T



    

# analytic_solution.py
# Python version of the analytic solution described in Appendix A of
# C.J. Roy, C.C. Nelson, T.M. Smith and C.C. Ober
# Verification of Euler/Navier-Stokes codes using the method
# of manufactured solutions.
# Int J for Numerical Methods in Fluids 2004; 44:599-620
#
# PJ, 28-May-2011
# It essentially Rowan's code with more and renamed variables
# to bring it closer to the original paper.
# PJ, 30-June-2012
# Scale the disturbance to reduce its magnitude away from the centre.

from math import sin, cos, pi, exp

class AnalyticSolution:
    def __init__(self, L,
                 rho0, rhox, rhoy, rhoxy, arhox, arhoy, arhoxy,
                 u0, ux, uy, uxy, aux, auy, auxy,
                 v0, vx, vy, vxy, avx, avy, avxy,
                 p0, px, py, pxy, apx, apy, apxy,
                 case):
        self.L = L
        self.rho0=rho0; self.rhox=rhox; self.rhoy=rhoy; self.rhoxy=rhoxy;
        self.arhox=arhox; self.arhoy=arhoy; self.arhoxy=arhoxy
        self.u0=u0; self.ux=ux; self.uy=uy; self.uxy=uxy;
        self.aux=aux; self.auy=auy; self.auxy=auxy
        self.v0=v0; self.vx=vx; self.vy=vy; self.vxy=vxy;
        self.avx=avx; self.avy=avy; self.avxy=avxy
        self.p0=p0; self.px=px; self.py=py; self.pxy=pxy;
        self.apx=apx; self.apy=apy; self.apxy=apxy; self.case=case
        return

    def S(self, x, y):
        if self.case == 1 or self.case == 2:
            return 1.0
        else:
            # Reduce the disturbance away from the centre of the domain.
            rsq = (x - self.L/2)**2 + (y - self.L/2)**2
            return exp(-16.0*rsq/(self.L**2))

    def rho(self, x, y):
        return self.rho0 + self.S(x,y)*self.rhox*sin(self.arhox*pi*x/self.L) \
            + self.S(x,y)*self.rhoy*cos(self.arhoy*pi*y/self.L) \
            + self.S(x,y)*self.rhoxy*cos(self.arhoxy*pi*x*y/(self.L*self.L))

    def u(self, x, y):
        return self.u0 + self.S(x,y)*self.ux*sin(self.aux*pi*x/self.L) \
            + self.S(x,y)*self.uy*cos(self.auy*pi*y/self.L) \
            + self.S(x,y)*self.uxy*cos(self.auxy*pi*x*y/(self.L*self.L))

    def v(self, x, y):
        return self.v0 + self.S(x,y)*self.vx*cos(self.avx*pi*x/self.L) \
            + self.S(x,y)*self.vy*sin((self.avy*pi*y)/self.L) \
            + self.S(x,y)*self.vxy*cos(self.avxy*pi*x*y/(self.L*self.L))
    
    def p(self, x, y):
        return self.p0 + self.S(x,y)*self.px*cos((self.apx*pi*x)/self.L) \
        + self.S(x,y)*self.py*sin(self.apy*pi*y/self.L) \
        + self.S(x,y)*self.pxy*sin(self.apxy*pi*x*y/(self.L*self.L))
    
if __name__ == "__main__":
    print "Display the analytic solution graphically."
    from numpy import array, linspace
    from pylab import contour, show
    xs = linspace(0.0, 1.0, 0.1)
    ys = linspace(0.0, 1.0, 0.1)
    nx = len(xs)
    ny = len(ys)
    supersonic = AnalyticSolution(L=1.0,
        rho0=1.0, rhox=0.15, rhoy=-0.1, rhoxy=0.0, arhox=1.0, arhoy=0.5, arhoxy=0.0,
        u0=800.0, ux=50.0, uy=-30.0, uxy=0.0, aux=1.5, auy=0.6, auxy=0.0,
        v0=800.0, vx=-75.0, vy=40.0, vxy=0.0, avx=0.5, avy=2.0/3, avxy=0.0,
        p0=1.0e5, px=0.2e5, py=0.5e5, pxy=0.0, apx=2.0, apy=1.0, apxy=0.0, case=1)
    x = 0.5
    y = 0.5
    print 'rho=', supersonic.rho(x,y),\
        'u=', supersonic.u(x,y), \
        'v=', supersonic.v(x,y), \
        'p=', supersonic.p(x,y)
        

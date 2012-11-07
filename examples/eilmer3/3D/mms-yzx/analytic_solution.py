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
# PJ, 07-Nov-2012
# yzx version

from math import sin, cos, pi, exp

class AnalyticSolution:
    def __init__(self, L,
                 rho0, rhoy, rhoz, rhoyz, arhoy, arhoz, arhoyz,
                 u0, uy, uz, uyz, auy, auz, auyz,
                 v0, vy, vz, vyz, avy, avz, avyz,
                 p0, py, pz, pyz, apy, apz, apyz,
                 case):
        self.L = L
        self.rho0=rho0; self.rhoy=rhoy; self.rhoz=rhoz; self.rhoyz=rhoyz;
        self.arhoy=arhoy; self.arhoz=arhoz; self.arhoyz=arhoyz
        self.u0=u0; self.uy=uy; self.uz=uz; self.uyz=uyz;
        self.auy=auy; self.auz=auz; self.auyz=auyz
        self.v0=v0; self.vy=vy; self.vz=vz; self.vyz=vyz;
        self.avy=avy; self.avz=avz; self.avyz=avyz
        self.p0=p0; self.py=py; self.pz=pz; self.pyz=pyz;
        self.apy=apy; self.apz=apz; self.apyz=apyz; self.case=case
        return

    def S(self, y, z):
        if self.case == 1 or self.case == 2:
            return 1.0
        else:
            # Reduce the disturbance away from the centre of the domain.
            rsq = (y - self.L/2)**2 + (z - self.L/2)**2
            return exp(-16.0*rsq/(self.L**2))

    def rho(self, y, z):
        return self.rho0 + self.S(y,z)*self.rhoy*sin(self.arhoy*pi*y/self.L) \
            + self.S(y,z)*self.rhoz*cos(self.arhoz*pi*z/self.L) \
            + self.S(y,z)*self.rhoyz*cos(self.arhoyz*pi*y*z/(self.L*self.L))

    def u(self, y, z):
        return self.u0 + self.S(y,z)*self.uy*sin(self.auy*pi*y/self.L) \
            + self.S(y,z)*self.uz*cos(self.auz*pi*z/self.L) \
            + self.S(y,z)*self.uyz*cos(self.auyz*pi*y*z/(self.L*self.L))

    def v(self, y, z):
        return self.v0 + self.S(y,z)*self.vy*cos(self.avy*pi*y/self.L) \
            + self.S(y,z)*self.vz*sin((self.avz*pi*z)/self.L) \
            + self.S(y,z)*self.vyz*cos(self.avyz*pi*y*z/(self.L*self.L))
    
    def p(self, y, z):
        return self.p0 + self.S(y,z)*self.py*cos((self.apy*pi*y)/self.L) \
        + self.S(y,z)*self.pz*sin(self.apz*pi*z/self.L) \
        + self.S(y,z)*self.pyz*sin(self.apyz*pi*y*z/(self.L*self.L))
    
if __name__ == "__main__":
    print "Display the analytic solution graphically."
    from numpy import array, linspace
    from pylab import contour, show
    ys = linspace(0.0, 1.0, 0.1)
    zs = linspace(0.0, 1.0, 0.1)
    ny = len(ys)
    nz = len(zs)
    supersonic = AnalyticSolution(L=1.0,
        rho0=1.0, rhoy=0.15, rhoz=-0.1, rhoyz=0.0, arhoy=1.0, arhoz=0.5, arhoyz=0.0,
        u0=800.0, uy=50.0, uz=-30.0, uyz=0.0, auy=1.5, auz=0.6, auyz=0.0,
        v0=800.0, vy=-75.0, vz=40.0, vyz=0.0, avy=0.5, avz=2.0/3, avyz=0.0,
        p0=1.0e5, py=0.2e5, pz=0.5e5, pyz=0.0, apy=2.0, apz=1.0, apyz=0.0, case=1)
    y = 0.5
    z = 0.5
    print 'rho=', supersonic.rho(y,z),\
        'u=', supersonic.u(y,z), \
        'v=', supersonic.v(y,z), \
        'p=', supersonic.p(y,z)
        

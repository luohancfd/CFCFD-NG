#
# diffractometer.py (c) Stuart B. Wilkins 2010 and (c) Sven Partzsch 2010
#
# $Id: diffractometer.py 176 2011-02-13 19:08:16Z stuwilkins $
# $HeadURL: https://pyspec.svn.sourceforge.net/svnroot/pyspec/trunk/pyspec/diffractometer.py $
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Part of the "pyspec" package
#
import numpy as np
import exceptions


__version__   = "$Revision: 176 $"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>" +\
                "Sven Partzsch <SvenPartzsch@gmx.de>"
__date__      = "$LastChangedDate: 2011-02-14 05:08:16 +1000 (Mon, 14 Feb 2011) $"
__id__        = "$Id: diffractometer.py 176 2011-02-13 19:08:16Z stuwilkins $"

class Diffractometer():
    """Diffractometer class

    This class provides various functions to perform calculations to
    and from the sample and instrument frames

    used lab frame:
    Z up, Y along X-ray beam, X = Y x Z;
    sample rotations:
    mu    : along +Z     -> S'    (mu-frame),
    theta : along +X'    -> S''   (theta-frame),
    chi   : along +Y''   -> S'''  (chi-frame),
    phi   : along +X'''  -> S'''' (phi-frame);
    detector rotations:
    mu    : along +Z     -> S'    (mu-frame),
    delta : along +X'    -> S*    (delta-frame),
    gamma : along +Z*    -> S**   (gamma-frame)"""

    # HC_OVER_E to convert from E to Lambda
    hc_over_e = 12398.4
    # transformation from sixc to fourc (tardis)
    sixcToFourc = np.array([ [ 0, 0, -1],
                             [ 0, 1,  0],
                             [ 1, 0,  0] ])

    def __init__(self, mode = 'sixc'):
        """Initialize the class
        'mode' defines the type of diffractometer"""
        self.mode = mode
        self.Qtheta = None

    #
    # set part
    #

    def setAllAngles(self, angles, mode = 'deg'):
        """Sets angles for calculation.
        Angles are expected in spec 'sixc' order:
        Delta     Theta       Chi       Phi        Mu     Gamma"""
        self._settingAngles = angles
        if mode == 'deg':
            self._settingAngles = self._settingAngles / 180.0 * np.pi

    def setAngles(self, delta = None, theta = None, chi = None,
                  phi = None, mu = None, gamma = None, mode = 'deg'):
        """Set the angles for calculation"""
        
        # Internally angles are stored in sixc (spec) order
        # Delta     Theta       Chi       Phi        Mu     Gamma

        maxlen = 1
        #for aa in zip([delta, theta, chi, phi, mu, gamma], range(6)):
        for aa in [delta, theta, chi, phi, mu, gamma]:
            if type(aa) == np.ndarray:
                if maxlen == 1:
                    maxlen = aa.size
                elif maxlen != aa.size:
                    raise exceptions.ValueError("Values must be numpy array of same size or scalar (float).")
                
        #maxlen = np.array([theta.size, delta.size, gamma.size, mu.size, phi.size, chi.size]).max()
        self._settingAngles = np.zeros((maxlen, 6))
        for aa, i in zip([delta, theta, chi, phi, mu, gamma], range(6)):
            if type(aa) == float:
                self._settingAngles[:,i] = np.ones(maxlen) * aa
            else:
                self._settingAngles[:,i] = aa
               
            """
            if aa.size == maxlen:
                self._settingAngles[:,i] = aa
            elif aa.size == 1:
                self._settingAngles[:,i] = ones(maxlen) * aa
            elif aa is not None:
                raise exceptions.ValueError("Values must be numpy array or scalar (float).")
            """
        if mode == 'deg':
            self._settingAngles = self._settingAngles / 180.0 * np.pi

    def setEnergy(self, energy):
        """Set the energy (in eV) for calculations"""
        self._waveLen = self.hc_over_e / energy

    def setLambda(self, waveLen):
        """Set the wavelength (in Angstroms) for calculations"""
        self._waveLen = waveLen

    def setUbMatrix(self, UBmat):
        """Sets the UB matrix"""
        self.UBmat = np.matrix(UBmat)
        self.UBinv = self.UBmat.I

    #
    # calc part
    #
    
    def _calc_QTheta(self):
        """Calculate (Qx, Qy, Qz) set in theta-frame from angles

        k0      = (0, kl, 0), kl: wave vector length,
        ki''    = rotX(-theta)*rotZ(-mu)*k0,
        kf''    = rotX(-theta+delta)*rotZ(gamma)*k0,
        QTheta = Q'' = kf'' - ki'' """

        # wave vector length in 1/A
        kl = 2 * np.pi / self._waveLen
        # wave vector for all diffractometer angles zero
        # k0 = [0, kl, 0]

        # alias for used angles
        mu    = self._settingAngles[:,4]
        theta = self._settingAngles[:,1]
        delta = self._settingAngles[:,0]
        gamma = self._settingAngles[:,5]

        ki = np.zeros((self._settingAngles.shape[0], 3))
        kf = np.zeros(ki.shape)
        
        # initial wave vector in theta-frame
        # ki'' = rotX(-theta)*rotZ(-mu)*k0
        ki[:,0] =  np.sin(mu)              *kl
        ki[:,1] =  np.cos(theta)*np.cos(mu)*kl
        ki[:,2] = -np.sin(theta)*np.cos(mu)*kl
        
        # final   wave vector in theta-frame
        # kf'' = rotX(-theta+delta)*rotZ(gamma)*k0
        kf[:,0] = -np.sin(gamma)                    *kl
        kf[:,1] =  np.cos(delta-theta)*np.cos(gamma)*kl
        kf[:,2] =  np.sin(delta-theta)*np.cos(gamma)*kl
        
        #   scattering vector in theta-frame
        q  = kf - ki

        self.QTheta = q

    def calc(self):
        self._calc_QTheta()
        
    #
    # get part
    #
    
    def getQTheta(self):
        """Return transformed coordinates"""
        return self.QTheta
    
    def getQPhi(self):
        """Calculate (Qx, Qy, Qz) set in phi-frame from (Qx, Qy, Qz) set in theta-frame

        QPhi = rotZ(-phi) rotY(-chi) QTheta """

        # alias for used angles
        chi = self._settingAngles[:,2]
        phi = self._settingAngles[:,3]
        # alias for q-vector in theta-frame
        QTh = self.QTheta
        
        # QPhi = rotZ(-phi) rotY(-chi) QTheta
        # matrix coefficients are calculated by hand for convenience
        r11 =               np.cos(chi)
        r12 =  0.0
        r13 =              -np.sin(chi)
        r21 =  np.sin(phi)* np.sin(chi)
        r22 =  np.cos(phi)
        r23 =  np.sin(phi)* np.cos(chi)
        r31 =  np.cos(phi)* np.sin(chi)
        r32 = -np.sin(phi)
        r33 =  np.cos(phi)* np.cos(chi)

        QPhi = np.zeros(QTh.shape)
        QPhi[:,0] = r11*QTh[:,0] + r12*QTh[:,1] + r13*QTh[:,2]
        QPhi[:,1] = r21*QTh[:,0] + r22*QTh[:,1] + r23*QTh[:,2]
        QPhi[:,2] = r31*QTh[:,0] + r32*QTh[:,1] + r33*QTh[:,2]
        
        return QPhi

    def getQCart(self):
        """Calculate (Qx, Qy, Qz) set in cartesian reciprocal space from (Qx, Qy, Qz) set in theta-frame

        still under construction """
        
        # still under construction

        #return (np.dot( self.sixcToFourc.T, self.getQPhi().T ) ).T

    def getQHKL(self):
        """Calc HKL values from (Qx, Qy, Qz) set in theta-frame with UB-matrix

        QHKL = UB^-1 sixcToFourc^T QPhi"""
 
        return ( self.UBinv * np.dot( self.sixcToFourc.T, self.getQPhi().T )  ).T



####################################
#
# main program, for test
#
####################################


if __name__ == "__main__":
                
    testDiff = Diffractometer()
    testDiff.setAngles(delta=40, theta=15, chi = 30, phi = 25, mu = 10.0, gamma=5.0)
    #testDiff.setLambda(2*np.pi)
    testDiff.setEnergy(640)
    testDiff.calc()
    testDiff.setUbMatrix([[-0.01231028454, 0.7405370482 , 0.06323870032], 
                          [ 0.4450897473 , 0.04166852402,-0.9509449389 ],
                          [-0.7449130975 , 0.01265920962,-0.5692399963 ]])
    #print 'QTheta = \n%s' % (testDiff.getQTheta())
    #print 'QPhiSixC = \n%s' % (testDiff.getQCart()  ) # under construction !!! #
    print 'HKL = \n%s' % (testDiff.getQHKL())
    print 'Ready'

    ###########
    # comparison with sixc simulation mode, checked all six angles
    # first with old spangleq.py and again for diffractometer class
    ###########

    """
    4.SIXC> setlat 6.28319 6.28319 6.28319 90 90 90
    (UB recalculated from orientation reflections and lattice.)
    """
    #14.SIXC> p UB
    #UB["0"] = 0.999999253114992
    #UB["1"] = -1.53360282755493e-16
    #UB["2"] = -1.60814167124132e-16
    #UB["3"] = -7.4538843686392e-18
    #UB["4"] = 0.999999253114992
    #UB["5"] = -6.12302719591125e-17
    #UB["6"] = 0
    #UB["7"] = 0
    #UB["8"] = 0.999999253114992
    """
    24.SIXC> LAMBDA = hc_over_e / 640

    34.SIXC> pa

    Six-Circle Geometry, Omega fixed (four circle, Mu = Gamma = 0) (mode 0)
    Sector 0

    Primary Reflection (at lambda 1.54):
     del th chi phi mu gam = 60 30 0 0 0 0 
                     H K L = 1 0 0

    Secondary Reflection (at lambda 1.54):
     del th chi phi mu gam = 60 30 0 90 0 0 
                     H K L = 0 1 0

    Lattice Constants (lengths / angles):
              real space = 6.283 6.283 6.283 / 90 90 90
        reciprocal space = 1 1 1 / 90 90 90

    Azimuthal Reference:
                   H K L = 0 0 1
               sigma tau = 0 0

    Monochromator:
                  Lambda = 19.3725 

    Cut Points:
      del   th  chi  phi   mu  gam
     -180 -180 -180 -180 -180 -180

    """

    #Qsam = angle2Qsam(mu, theta, chi, phi, delta, gamma, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    """
    Qsam = angle2Qsam(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    python: In [81]: Qsam
            Out[81]: array([ 0.,  0.,  0.])
    38.SIXC> wh

    H K L =  0  0  0
    Alpha = 0  Beta = 0  Azimuth = 180
    Two Theta = 0  Omega = 0  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000

    
    Qsam = angle2Qsam(0.0, 15.0, 0.0, 0.0, 30.0, 0.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    python: In [83]: Qsam
            Out[83]: array([ 0.16788546,  0.        ,  0.        ])
    40.SIXC> wh

    H K L =  0.16789  3.9558e-09  0
    Alpha = 0  Beta = 0  Azimuth = -90
    Two Theta = 30  Omega = 1.35e-06  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    30.0000   15.0000    0.0000    0.0000    0.0000    0.0000         


    Qsam = angle2Qsam(0.0, 15.0, 0.0, 0.0, 40.0, 0.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [85]: Qsam
    Out[85]: array([ 0.22101043, -0.01933591,  0.        ])
    42.SIXC> wh

    H K L =  0.22101  -0.019336  0
    Alpha = 0  Beta = 0  Azimuth = -90
    Two Theta = 40  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000    0.0000    0.0000    0.0000    0.0000

    
    Qsam = angle2Qsam(0.0, 15.0, 0.0, 0.0, 40.0, 5.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [87]: Qsam
    Out[87]: array([ 0.22048884, -0.02045445,  0.0282672 ])
    44.SIXC> wh

    H K L =  0.22049  -0.020455  0.028268
    Alpha = 0  Beta = 5  Azimuth = -87.318
    Two Theta = 40.259  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000    0.0000    0.0000    0.0000    5.0000

    
    Qsam = angle2Qsam(10.0, 15.0, 0.0, 0.0, 40.0, 5.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [89]: Qsam
    Out[89]: array([ 0.21921356, -0.01569504,  0.08458648])
    46.SIXC> wh

    H K L =  0.21922  -0.015695  0.084588
    Alpha = 10  Beta = 5  Azimuth = -92.851
    Two Theta = 42.574  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000    0.0000    0.0000   10.0000    5.0000

    
    Qsam = angle2Qsam(10.0, 15.0, 30.0, 0.0, 40.0, 5.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [91]: Qsam
    Out[91]: array([ 0.14755127, -0.01569504,  0.18286083])
    48.SIXC> wh

    H K L =  0.14755  -0.015695  0.18286
    Alpha = 16.131  Beta = 16.618  Azimuth = -89.602
    Two Theta = 42.574  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000   30.0000    0.0000   10.0000    5.0000

    
    Qsam = angle2Qsam(10.0, 15.0, 30.0, 25.0, 40.0, 5.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [93]: Qsam
    Out[93]: array([ 0.14035988,  0.04813332,  0.18286083])
    50.SIXC> wh

    H K L =  0.14036  0.048134  0.18286
    Alpha = 16.131  Beta = 16.618  Azimuth = -89.602
    Two Theta = 42.574  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000   30.0000   25.0000   10.0000    5.0000
    """


"""




Notes : We define the refractive index the following way

n = 1 - \delta + i\Beta

"""

import scipy
import sfact
import rod
import spacegroup
import structures
import copy as _copy
import numpy as np
import numpy.ma as masked

try:
    import psyco
except ImportError:
    pass
else:
    psyco.profile()

import pylab
from pylab import *

__verbose__ = True

class Atom(object):
    """Class defining a single scatterer in a crystal (atom)

    This class reprosents a single scatterer. It will have the atomic
    species, phase factor and scattering length.

    """
    def __init__(self, 
                 name = '', Z = '', charge = '',
                 x = None, y = None, z = None,  
                 occ = None, U = None, sg = False, r = None):
        """Initialize the Atom by defining parameters

        Parameters
        ----------
        name      : [string] Name to define the atom 
        Z         : [string] Atomic species e.g. 'Fe'
        charge    : [string] Atomic charge e.g. '2+'
        x,y,z     : [float]  Position in unit cell (0..1)
        occ       : [float]  Occupancy (0..1)
        U         : [float]  Debye-Waller factor
        sg        : [string] Spacegroup e.g. 'I4/mmm'
        r         : [array]  Postion of atom in unit cell

        Returns
        -------

        Atom object

        """

        self.name = name
        self.Z = Z
        self.charge = charge

        self.f = 0.0
        self.p = 0.0

        self.sgGenerated = sg

        if r is not None:
            self.r = r
        else:
            if x is not None:
                self.r = scipy.array([x, y, z])
            else:
                self.r = scipy.array([0, 0, 0])

        if U is not None:
            self.U = U
        else:
            self.U = array([0.0])

        if occ is not None:
            self.occ = occ
        else:
            self.occ = 1.0
        
    def setName(self, name = ""):
        """Sets the name of the atom"""
        self.name = name

    def setZ(self, Z):
        """Sets the atomic species"""
        self.Z = Z

    def getZ(self):
        """Returns the atomic species""" 
        return self.Z

    def getU(self):
        """Returns the Debye-Waller factor""" 
        return self.U

    def getName(self):
        """Returns the name of the atom"""
        return self.name

    def getCharge(self):
        """Returns the charge of the atom"""
        return self.charge

    def getSgGenerated(self):
        """Returns the name of the space group"""
        return self.sgGenerated

    def setSgGenerated(self, sg):
        """Sets the name of the spacegroup"""
        self.sgGenerated = sg

    def getPositionVector(self):
        """Set the position vector of the atom"""
        return self.r

    def setPositionVector(self, r):
        """Returns the position vector of the atom"""
        self.r = r

    def setPosition(self, 
                    x = 0.0, y = 0.0, z = 0.0, 
                    o = 1.0):
        """Set the position and occupancy of the scatterer"""
        self.r = scipy.array([x, y, z])
        self.o = o

    def getOccupancy(self):
        """Returns the occupancy of the atom"""
        return self.occ

    def setOccupancy(self, o):
        """Sets the occupancy of the atom"""
        self.occ = o

    def calc(self, energy = None):
        """Force calculation of atom parameters

        Parameters
        ----------
        
        energy : [float] Incident photon energy in keV

        This function will read both interpolate objects for
        f1 and f2 and also the parameters for f0. If energy is
        None type then no calulation of f1,f2 will be performed

        """
        self.f0p = sfact.getF0Params(self.Z + self.charge)
        self.f1, self.f2 = sfact.getF1F2(self.Z)

        if energy is not None:
            self.f1f2 = complex(self.f1(energy), self.f2(energy))
        else:
            self.f1f2 = None
    
    def getScatLen(self, q = 0, energy = None, anom = True):
        """Caculates and returns the scattering length of the atom

        Parameters
        ----------

        q        : [float] Magnitude of the scattering vector.
        energy   : [float] Energy in keV
        anom     : [bool]  If false do not include f1,f2

        """
        
        self.f = sfact.calcF0(self.f0p, q)

        if self.U is not None:
            if self.U.shape == (1,):
                if q.ndim == 1:
                    self.dw = exp(-1.0 * self.U * pow(q, 2))
                else:
                    self.dw = exp(-1.0 * self.U * pow(q, 2).sum(1))
            elif self.U.shape == (3,):
                self.dw = exp(-1.0 * inner(self.U, pow(q, 2)))
            else:
                raise Exception('Invalid shape of Debye-Waller Factor')
            self.f = self.dw * self.f

        if anom == True:
            if energy is not None:
                f = complex(self.f1(energy), self.f2(energy))
                self.f = self.f + f
            elif self.f1f2 is not None:
                self.f = self.f + self.f1f2
            
        return self.f

    def getElectronDensity(self, x = 0.):
        """Calculates and returns the electron density dist. for the atom"""
        return sfact.getRealF0(self.Z + self.charge, x)

class Crystal(object):
    """Class defining an x-ray diffraction crystal

    This class reprosents a crystal. It is made up of a collection of
    atoms that are defined by the Atom class. It holds information such
    as the energy (wavelength) of the incident x-rays and the current 
    wavevector transfer. It can calculate structure factors for x-ray 
    diffraction experiments.

    This module is linked to the sginfo c-code to do spacegroup calcs
    and will generate all atoms from a given space group.

    """
    N_A = 6.02214179e-23 
    R_0 = 2.8179402894e-5 # Angstroms 
    def __init__(self, 
                 a = None, b = None, c = None,
                 alpha = 90.0, beta = 90.0, gamma = 90.0,
                 lattice = None):
        """Initialize the Crystal class

        Parameters
        ----------
        a,b,c              : [float] Lattice parameters in Angstroms
        alpha, beta, gamma : [float] Cell angles in degrees
        lattice            : [array] array of the above in thr format
                             array([a, b, c, alpha, beta, gamma])
        
        This initialization will calculate the reciprocal lattice
        a further call to calc() should be made to calculate crystal
        parameters.

        """

        if lattice is None:
            if a is None:
                self.lattice = scipy.array([5., 5., 5., 90., 90., 90.])
            else:
                self.lattice = scipy.array([a, b, c, alpha, beta, gamma])
        else:
            self.lattice = lattice

        self.calcRLattice()

        self.n = 1
        self.alpha_c = 0.
        self.n = 0.
        self.delta = 0.
        self.beta = 0.
        self.energy = None
        self.k = None
        self.atoms = []
        
        self.sg = None

    def calcRLattice(self):
        """Calculate reciprocal lattice parameters"""
        self.rlattice = 2 * scipy.pi / self.lattice[0:3]

    def getLattice(self):
        """Returns array of direct-space lattice""" 
        return self.lattice

    def getRLattice(self):
        """Returns array of reciprocal lattice"""
        return self.rlattice

    def setLattice(self, lattice):
        """Sets lattice parameters
        
        Parameters
        ----------
        lattice     : [array] array of the lattice parameters in the format
                              array([a, b, c, alpha, beta, gamma])

        """
        self.lattice = lattice

    def setEnergy(self, e):
        """Set the energy of incident x-rays

        Parameters
        ----------
        e   : [float] energy in keV

        Calcuates |k| along with energy.
        Note: The magnitude of k is defined as 
        |k| = 2 * pi * e / 12.39
        """
        self.k = 2.0 * scipy.pi * e / 12.390
        self.energy = e * 1000
        
    def setLambda(self, l):
        """Set the wavelength of incident x-rays

        Parameters
        ----------
        l   : [float] wavelength in Angstroms

        """
        self.k = 2 * scipy.pi / l
        self.energy = 12390 / l
    
    def getLambda(self):
        """Get the wavelength of the incident x-rays"""
        return 12390./self.energy

    def setAtoms(self, a):
        """ Set the atom names and positions as a tuple.

        Atoms should be a list of Atom class objects.

        """
        self.atoms = a

    def addAtom(self, *args, **kwargs):
        """Add an atom to the crystal"""
        newa = Atom(*args, **kwargs)
        self.atoms.append(newa)

    def truncate(self, value):
        """Truncate the unit cell

        Here we truncate the unit cell by a vector value.
        The atoms which are truncated have their occupancy set
        to zero. For example:

        """
        
        for atom in self.atoms:
            r = atom.getPositionVector()
            if (r >= value).sum():
                print "Removing atom", atom.getName(), "at", r
                atom.setOccupancy(0.)
           

        #self.lattice[:3] = value * self.lattice[:3]
        #self.calcRLattice()

    def relax(self, value):
        """Relax the unit cell

        Relaxes the unit cell by the value in the array passed.
        For example: if value = [1, 1, 1.1] there will be a relaxation
        of 10% along the c-direction

        """

        l = self.lattice[:3] * value
        self.lattice[:3] = l
        self.calcRLattice()

    def getCellPhase(self, q):
        """Returns the phase change over the unit cell"""
        return exp(complex(0, inner(self.lattice[:3], q)))

    def setOrigin(self, origin, wrap = True):
        """Move the origin of the unit cell

        Here the origin is displaced by "origin" along the [001] axis.
        All atoms which have negative positions along the [001] direction
        are moved by one unit along the [001] axis.

        If wrap is True and if the resulting position is -ve then 1 is added
        to the position vector

        """

        o = array([0, 0, origin])

        for a in self.atoms:
            r = a.getPositionVector()
            r = r - o
            if wrap:
                if r[2] < 0:
                    r[2] = r[2] + 1
            
            a.setPositionVector(r)

    def calcF(self, q, anom = True):
        # Calclate Structure Factor

        amplitude = complex(0., 0.)
        lattice = self.lattice[:3]
        for a in self.atoms:
            r = lattice * a.r
            phase = exp(complex(0, 1.) * inner(r, q))
            aa = (a.getOccupancy() * a.getScatLen(q, anom = anom) * phase)
            amplitude += aa
        
        self.xrayAmplitude = amplitude
        return self.xrayAmplitude

    def calc(self, energy = None):
        """Calculate all parameters for crystal independant of lambda"""

        # Unit cell volume

        self.cellVolume = (self.lattice[0] * self.lattice[1] * self.lattice[2])

        # Generate spacegroup atoms

        self.generateAtoms()

        if energy is not None:
            self.calcE(energy)
            return 

        if self.energy is not None:
            self.calcE()
            return

    def calcE(self, energy = None):
        """Calculate parameters dependant on lambda"""
        
        if energy is not None:
            self.setEnergy(energy)
        
        # Calculate Refractive Index and mu
        
        refindex = complex(0,0)
        for a in self.atoms:
            a.calc(energy = self.energy)
            
            r = a.getScatLen(q = scipy.array([0]))
            refindex = refindex + complex(r.real, -1.0 * r.imag)
        refindex = (refindex * self.R_0 * 2 * pi / (pow(self.k,2) * self.cellVolume))

        self.delta = refindex.real
        self.beta = refindex.imag
        self.n = 1 - refindex

        self.alpha_c = pow(self.delta * 2, 0.5)
        self.mu = 2.0 * self.beta * self.k

    def getTransmittivity(self, alpha):
        """Returns the ampltiude transmittivity"""
        alpha_prime = sqrt(pow(alpha, 2) - pow(self.alpha_c, 2) + complex(0, 2* self.beta))
        t = 2 * alpha / (alpha + alpha_prime)

        return t

    def getReflectivity(self, alpha):
        """Returns the amplitude reflectivity"""
        alpha_prime = sqrt(pow(alpha, 2) - pow(self.alpha_c, 2) + complex(0, 2* self.beta))
        r = (alpha - alpha_prime) / (alpha + alpha_prime)
        return r

    def getMu(self):
        return self.mu

    def setHKL(self, hkl):
        self.hkl = hkl

    def calcPenetrationDepth(self, alpha):
        q = alpha / self.alpha_c
        b_mu = pow(2 * self.k / (2 * self.k * self.alpha_c), 2) * self.beta 
        q_prime = pow(pow(q,2) - 1 + complex(0,2 * b_mu), 0.5)
        self.pdepth = 1 / (2 * self.k * self.alpha_c * q_prime.imag)
        self.pdepth = self.pdepth / (2 * self.k * self.alpha_c)
        return self.pdepth

    def setSpaceGroup(self, name):
        """Sets the spacegroup of the crystal"""
        self.sg = spacegroup.Spacegroup(name)

    def generateAtoms(self):

        # Check if we have a spacegroup
        if self.sg is None:
            return 0
        
        generators = self.sg.getGenerators()
        trvectors = self.sg.getTranslationVectors()

        while True:
            atomAdded = False
            for atom in self.atoms:
                for generator in generators:
                    if self.applyGenerator(atom, generator):
                        attomAdded = True

                    if self.sg.isCentric():
                        centric_generator = generator * -1.
                        centric_generator[3, 3] = 1
                        if self.applyGenerator(atom, centric_generator):
                            attomAdded = True

            if not atomAdded:
                break

    def applyGenerator(self, atom, generator):
        # Now apply the symmetry equivalent position
        r = atom.getPositionVector()
        r = hstack((r, array([1.])))
        r_eq = dot(generator, r.T)

        # Normalize this vector from 0 to 1
        while ((r_eq[:3] < 0.).sum()) or ((r_eq[:3] >= 1.).sum()):
            r_eq = r_eq + hstack(((r_eq[:3] < 0.) * 1., array([1])))
            r_eq = r_eq - hstack(((r_eq[:3] >= 1.) * 1., array([1])))
                
        # Now check if we have an atom at that location
        if not self.isAtomAt(r_eq[0:3], atom.getOccupancy()):
            occ = atom.getOccupancy()
            Z = atom.getZ()
            name = atom.getName()
            charge = atom.getCharge()
            U = atom.getU()
            self.addAtom(r = r_eq[0:3], 
                         occ = occ, name = name, 
                         Z=Z, U=U, charge = charge, sg = True)
            return True
        
        return False
            
    def isAtomAt(self, r, occ = 1., tol = 5e-5):
        """Looks to see if there is an atom at a given location

        Returns a list of the atoms which are at the location specified

        """
        atoms_at = []
        for a in self.atoms:
            e = abs(a.getPositionVector() - r).sum()
            if (e < tol) and (occ == a.getOccupancy()):
                atoms_at.append(a)

        return atoms_at

    def __str__(self):
        out = self.showCrystal()
        out = out + '\n'
        out = out + self.showAtoms()
        out = out + '\n'
        out = out + self.showOccupancy()
        return out

    def showOccupancy(self):
        out =  "--------------------------------------------------------------------------------\n"
        occsum = 0.0
        
        for a in self.atoms:
            occsum += a.getOccupancy()

        out += "No of atoms in crystal : %5d\n" % len(self.atoms)
        out += "Sum of occupancy : %5.3f\n" % occsum
        return out

    def showAtoms(self):
        out = "%-10s %-10s %-6s    %-6s    %-6s    %-10s    %-5s    %s\n" % (
            "Name", "Z", "rx", "ry", "rz", "U", "Occ.", "SGG")
        out = out + "--------------------------------------------------------------------------------\n"
        out = out + '\n'

        for a in self.atoms:
            out = out + "%-10s %-10s" % (a.getName(), a.getZ() + a.getCharge())
            r = a.getPositionVector()
            out = out + " %5.4f    %5.4f    %5.4f" % (r[0], r[1], r[2])
            out = out + "    %5.4e" % a.getU()
            out = out + "    %4.3f" % a.getOccupancy()
            out = out + "    %s" % a.getSgGenerated()
            out = out + "\n"
        
        return out

    def showCrystal(self):
        fmte = "%-20s = %5.4e %s\n"
        fmtf = "%-20s = %5.4f %s\n"
        sep = "--------------------------------------------------------------------------------\n"
        
        out = ""

        out = out + "Lattice Parameters:\n"
        out = out + sep
        out = out + fmtf % ("a", self.lattice[0], "Angstroms")
        out = out + fmtf % ("b", self.lattice[1], "Angstroms")
        out = out + fmtf % ("c", self.lattice[2], "Angstroms")
        out = out + fmtf % ("alpha", self.lattice[3], "Angstroms")
        out = out + fmtf % ("beta", self.lattice[4], "Angstroms")
        out = out + fmtf % ("gamma", self.lattice[5], "Angstroms")
        out = out + "\n"

        out = out + "Spacegroup:\n"
        out = out + sep
        out = out + str(self.sg)
        out = out + '\n'

        out = out + "Calculated Parameters:\n"
        out = out + sep
        out = out + fmte % ("Critical Angle", 180 * self.alpha_c / pi, "Deg")
        out = out + fmte % ("Delta", self.delta, "")
        out = out + fmte % ("Beta", self.beta, "")
        #out = out + fmte % ("Ref. Index", self.n, "")
        out = out + fmte % ("mu", self.mu, "Angstroms^-1")
        
        return out

    def showDensity(self, x, d = array([0, 0, 1])):
        """Return electron density along a given direction"""
        
        abc = self.getLattice()[:3]
        rho = zeros(len(x))
        for a in self.atoms:
            r = dot(d, a.getPositionVector() * abc)
            rho += (a.getElectronDensity(x - r) * a.getOccupancy())
            
        return rho


def _test():
    pass

if __name__ == "__main__":
    _test()

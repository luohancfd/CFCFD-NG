"""

    pyxray (c) Stuart B. Wilkins 2008
    spacegroup module

    
    
"""
#
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
#  $Id: graph.py 88 2009-01-07 02:20:19Z swilkins $
#  $HeadURL: svn+ssh://solids/home/swilkins/.svnrepos/pyVacuum/pyVacuum/graph.py $
#
#
from numpy import identity, concatenate, array
import sginfo as sginfo

class Spacegroup():
    """Spacegroup object"""
    def __init__(self, name = None):
        self.si = None
        if name is not None:
            self.setSpacegroup(name)

    def __str__(self):
        fmts = "%-20s = %s\n"

        if self.si is not None:
            out = "Spacegroup %s (No. %d)\n" % (self.sg_name_info['Label'].replace('_',' '),
                                                self.sg_name_info['Number'])
            out = out + fmts % ("Hall name", self.sg_name_info['Hall'])
            out = out + fmts % ("System", self.sg_stats['System'])
            if self.sg_name_info['Setting'] is not '':
                out = out + fmts % ("Setting", self.sg_name_info['Setting'])
            out = out + fmts % ("Point Group", self.sg_stats['PointGroup'])
            out = out + fmts % ("Laue Group", self.sg_stats['LaueGroup'])

            out = out + "\n"
            out = out + "Generators:\n\n"
            for i in range(self.sg_generators.shape[0]):
                out = out + self._generator_to_str(self.sg_generators[i])

            out = out + "\n"
            out = out + "Lattice Translations:\n\n"
            for t in self.sg_translations:
                out = out + "%3.2f x + %3.2f y + %3.2f z\n" % (t[0], t[1], t[2])
            return out
        else:
            return "None"

    def _generator_to_str(self, generator):
        xyz = 'xyz'
        out = ""
        for y in range(3):
            for x in range(0,3):
               if generator[y, x] > 0:
                   if out is not "":
                       out = out + "  +  " + xyz[x]
                   else:
                       out = out + "     " + xyz[x]

               elif generator[y, x] < 0:
                   out = out + "  -  " + xyz[x]
        
        if generator[y, 3]:
            out = out + "%s\n" % generator[y,3] 
        return out + "\n"

    def setSpacegroup(self, name):
        self.si = sginfo.getsg(name)
        self.sg_name_info = self.si[0]
        self.sg_stats = self.si[1]
        
        # Make translation vectors into seitz 
        # matricies, then they become generators

        self.sg_translations = self.si[2]
        self.sg_generators = self.si[3]
        self.generators = self.sg_generators

        for t in self.sg_translations:
            if sum(t): 
                seitz = identity(4)
                seitz[0,3] = t[0]
                seitz[1,3] = t[1]
                seitz[2,3] = t[2]
                self.generators = concatenate((self.generators, array([seitz])))

    def getHallName(self):
        return self.sg_name_info['Hall']

    def stats(self):
        out = ''
        return out

    def getGenerators(self):
        return self.generators
    
    def getTranslationVectors(self):
        return self.sg_translations

    def isCentric(self):
        return self.sg_stats['Centric'] == -1

    def isInversionOffOrigin(self):
        return self.sg_stats['InversionOffOrigin']

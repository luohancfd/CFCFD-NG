#
# utils.py (c) Stuart B. Wilkins 2011
#
# $Id: transformations.py 209 2011-05-03 19:32:58Z stuwilkins $
# $HeadURL: https://pyspec.svn.sourceforge.net/svnroot/pyspec/trunk/pyspec/ccd/transformations.py $
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

def binArray(a, bins):
    """Bin 2D array

    Bin 
    """
    if (a.size[0] % bins[0]):
        raise Exception("To apply bin data x dimension should be divisible by no of bins")
    if (a.size[1] % bins[1]):
        raise Exception("To apply bin data y dimension should be divisible by no of bins")
    a = a.resize((img.size[0],-1,self._binOnRead[0])).sum(2)
    a = a.resize((-1, self._binOnRead[1], img.size[1])).sum(1)
    return a

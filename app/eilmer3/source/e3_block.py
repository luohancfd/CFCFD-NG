"""
e3_block.py -- Classes for defining 2D and 3D blocks.

It is expected that this module be imported by the application program e3prep.py.
The classes and functions will then be available for use in the user's input scripts.

.. Author: PJ
.. Versions:
   16-Mar-2008: extracted from e3prep.py (formerly mbcns_prep.py)
   29-Nov-2010: moved face and vertex definitions to separate file.
"""

import sys
from libprep3 import *
# Dictionaries to look up face index values from name or number.
from e3_defs import *
from e3_grid import *
from e3_flow import *
# Dictionaries to look up boundary-condition,
from bc_defs import *
from flux_dict import *
import string

# The folowing dictionary provides a map from inter-block connections
# that are defined by (A,B) vertex pairs to inter-block connections that
# are defined by (A-face, B-face, orientation, axis-map) tuples.
# The orientation is used in the C-code to decide where to transfer data
# at the block boundaries and the axis-map is used to interpret GridPro
# connectivity data.  The i,j,k axes of *this* block are aligned with the
# specified axes of the *other* block. 
connectionDict3D = {}
# North-North
vpairs = [(3,2),(7,6),(6,7),(2,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, NORTH, 0, '-i-j+k')
vpairs = [(3,3),(7,2),(6,6),(2,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, NORTH, 1, '+k-j+i')
vpairs = [(3,7),(7,3),(6,2),(2,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, NORTH, 2, '+i-j-k')
vpairs = [(3,6),(7,7),(6,3),(2,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, NORTH, 3, '-k-j-i')
# North-South
vpairs = [(3,0),(7,4),(6,5),(2,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, SOUTH, 0, '+i+j+k')
vpairs = [(3,1),(7,0),(6,4),(2,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, SOUTH, 1, '+k+j-i')
vpairs = [(3,5),(7,1),(6,0),(2,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, SOUTH, 2, '-i+j-k')
vpairs = [(3,4),(7,5),(6,1),(2,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, SOUTH, 3, '-k+j+i')
# North-East
vpairs = [(3,1),(7,5),(6,6),(2,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, EAST, 0, '+j-i+k')
vpairs = [(3,2),(7,1),(6,5),(2,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, EAST, 1, '+k-i-j')
vpairs = [(3,6),(7,2),(6,1),(2,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, EAST, 2, '-j-i-k')
vpairs = [(3,5),(7,6),(6,2),(2,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, EAST, 3, '-k-i+j')
# North-West
vpairs = [(3,3),(7,7),(6,4),(2,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, WEST, 0, '-j+i+k')
vpairs = [(3,0),(7,3),(6,7),(2,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, WEST, 1, '+k+i+j')
vpairs = [(3,4),(7,0),(6,3),(2,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, WEST, 2, '+j+i-k')
vpairs = [(3,7),(7,4),(6,0),(2,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, WEST, 3, '-k+i-j')
# North-Top
vpairs = [(3,4),(7,7),(6,6),(2,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, TOP, 0, '+i-k+j')
vpairs = [(3,5),(7,4),(6,7),(2,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, TOP, 1, '+j-k-i')
vpairs = [(3,6),(7,5),(6,4),(2,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, TOP, 2, '-i-k-j')
vpairs = [(3,7),(7,6),(6,5),(2,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, TOP, 3, '-j-k+i')
# North-Bottom
vpairs = [(3,1),(7,2),(6,3),(2,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, BOTTOM, 0, '-i+k+j')
vpairs = [(3,0),(7,1),(6,2),(2,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, BOTTOM, 1, '+j+k+i')
vpairs = [(3,3),(7,0),(6,1),(2,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, BOTTOM, 2, '+i+k-j')
vpairs = [(3,2),(7,3),(6,0),(2,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (NORTH, BOTTOM, 3, '-j+k-i')

# South-North
vpairs = [(1,2),(5,6),(4,7),(0,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, NORTH, 0, '+i+j+k')
vpairs = [(1,3),(5,2),(4,6),(0,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, NORTH, 1, '-k+j+i')
vpairs = [(1,7),(5,3),(4,2),(0,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, NORTH, 2, '-i+j-k')
vpairs = [(1,6),(5,7),(4,3),(0,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, NORTH, 3, '+k+j-i')
# South-South
vpairs = [(1,0),(5,4),(4,5),(0,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, SOUTH, 0, '-i-j+k')
vpairs = [(1,1),(5,0),(4,4),(0,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, SOUTH, 1, '-k-j-i')
vpairs = [(1,5),(5,1),(4,0),(0,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, SOUTH, 2, '+i-j-k')
vpairs = [(1,4),(5,5),(4,1),(0,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, SOUTH, 3, '+k-j+i')
# South-East
vpairs = [(1,1),(5,5),(4,6),(0,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, EAST, 0, '-j+i+k')
vpairs = [(1,2),(5,1),(4,5),(0,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, EAST, 1, '-k+i-j')
vpairs = [(1,6),(5,2),(4,1),(0,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, EAST, 2, '+j+i-k')
vpairs = [(1,5),(5,6),(4,2),(0,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, EAST, 3, '+k+i+j')
# South-West
vpairs = [(1,3),(5,7),(4,4),(0,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, WEST, 0, '+j-i+k')
vpairs = [(1,0),(5,3),(4,7),(0,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, WEST, 1, '-k-i+j')
vpairs = [(1,4),(5,0),(4,3),(0,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, WEST, 2, '-j-i-k')
vpairs = [(1,7),(5,4),(4,0),(0,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, WEST, 3, '+k-i-j')
# South-Top
vpairs = [(1,4),(5,7),(4,6),(0,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, TOP, 0, '-i+k+j')
vpairs = [(1,5),(5,4),(4,7),(0,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, TOP, 1, '-j+k-i')
vpairs = [(1,6),(5,5),(4,4),(0,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, TOP, 2, '+i+k-j')
vpairs = [(1,7),(5,6),(4,5),(0,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, TOP, 3, '+j+k+i')
# South-Bottom
vpairs = [(1,1),(5,2),(4,3),(0,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, BOTTOM, 0, '+i-k+j')
vpairs = [(1,0),(5,1),(4,2),(0,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, BOTTOM, 1, '-j-k+i')
vpairs = [(1,3),(5,0),(4,1),(0,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, BOTTOM, 2, '-i-k-j')
vpairs = [(1,2),(5,3),(4,0),(0,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (SOUTH, BOTTOM , 3, '+j-k-i')

# East-North
vpairs = [(2,2),(6,6),(5,7),(1,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, NORTH, 0, '-j+i+k')
vpairs = [(2,3),(6,2),(5,6),(1,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, NORTH, 1, '-j-k+i')
vpairs = [(2,7),(6,3),(5,2),(1,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, NORTH, 2, '-j-i-k')
vpairs = [(2,6),(6,7),(5,3),(1,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, NORTH, 3, '-j+k-i')
# East-South
vpairs = [(2,0),(6,4),(5,5),(1,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, SOUTH, 0, '+j-i+k')
vpairs = [(2,1),(6,0),(5,4),(1,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, SOUTH, 1, '+j-k-i')
vpairs = [(2,5),(6,1),(5,0),(1,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, SOUTH, 2, '+j+i-k')
vpairs = [(2,4),(6,5),(5,1),(1,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, SOUTH, 3, '+j+k+i')
# East-East
vpairs = [(2,1),(6,5),(5,6),(1,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, EAST, 0, '-i-j+k')
vpairs = [(2,2),(6,1),(5,5),(1,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, EAST, 1, '-i-k-j')
vpairs = [(2,6),(6,2),(5,1),(1,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, EAST, 2, '-i+j-k')
vpairs = [(2,5),(6,6),(5,2),(1,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, EAST, 3, '-i+k+j')
# East-West
vpairs = [(2,3),(6,7),(5,4),(1,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, WEST, 0, '+i+j+k')
vpairs = [(2,0),(6,3),(5,7),(1,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, WEST, 1, '+i-k+j')
vpairs = [(2,4),(6,0),(5,3),(1,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, WEST, 2, '+i-j-k')
vpairs = [(2,7),(6,4),(5,0),(1,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, WEST, 3, '+i+k-j')
# East-Top
vpairs = [(2,4),(6,7),(5,6),(1,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, TOP, 0, '-k-i+j')
vpairs = [(2,5),(6,4),(5,7),(1,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, TOP, 1, '-k-j-i')
vpairs = [(2,6),(6,5),(5,4),(1,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, TOP, 2, '-k+i-j')
vpairs = [(2,7),(6,6),(5,5),(1,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, TOP, 3, '-k+j+i')
# East-Bottom
vpairs = [(2,1),(6,2),(5,3),(1,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, BOTTOM, 0, '+k+i+j')
vpairs = [(2,0),(6,1),(5,2),(1,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, BOTTOM, 1, '+k-j+i')
vpairs = [(2,3),(6,0),(5,1),(1,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, BOTTOM, 2, '+k-i-j')
vpairs = [(2,2),(6,3),(5,0),(1,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (EAST, BOTTOM, 3, '+k+j-i')

# West-North
vpairs = [(0,2),(4,6),(7,7),(3,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, NORTH, 0, '+j-i+k')
vpairs = [(0,3),(4,2),(7,6),(3,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, NORTH, 1, '+j+k+i')
vpairs = [(0,7),(4,3),(7,2),(3,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, NORTH, 2, '+j+i-k')
vpairs = [(0,6),(4,7),(7,3),(3,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, NORTH, 3, '+j-k-i')
# West-South
vpairs = [(0,0),(4,4),(7,5),(3,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, SOUTH, 0, '-j+i+k')
vpairs = [(0,1),(4,0),(7,4),(3,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, SOUTH, 1, '-j+k-i')
vpairs = [(0,5),(4,1),(7,0),(3,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, SOUTH, 2, '-j-i-k')
vpairs = [(0,4),(4,5),(7,1),(3,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, SOUTH, 3, '-j-k+i')
# West-East
vpairs = [(0,1),(4,5),(7,6),(3,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, EAST, 0, '+i+j+k')
vpairs = [(0,2),(4,1),(7,5),(3,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, EAST, 1, '+i+k-j')
vpairs = [(0,6),(4,2),(7,1),(3,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, EAST, 2, '+i-j-k')
vpairs = [(0,5),(4,6),(7,2),(3,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, EAST, 3, '+i-k+j')
# West-West
vpairs = [(0,3),(4,7),(7,4),(3,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, WEST, 0, '-i-j+k')
vpairs = [(0,0),(4,3),(7,7),(3,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, WEST, 1, '-i+k+j')
vpairs = [(0,4),(4,0),(7,3),(3,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, WEST, 2, '-i+j-k')
vpairs = [(0,7),(4,4),(7,0),(3,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, WEST, 3, '-i-k-j')
# West-Top
vpairs = [(0,4),(4,7),(7,6),(3,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, TOP, 0, '+k+i+j')
vpairs = [(0,5),(4,4),(7,7),(3,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, TOP, 1, '+k+j-i')
vpairs = [(0,6),(4,5),(7,4),(3,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, TOP, 2, '+k-i-j')
vpairs = [(0,7),(4,6),(7,5),(3,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, TOP, 3, '+k-j+i')
# West-Bottom
vpairs = [(0,1),(4,2),(7,3),(3,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, BOTTOM, 0, '-k-i+j')
vpairs = [(0,0),(4,1),(7,2),(3,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, BOTTOM, 1, '-k+j+i')
vpairs = [(0,3),(4,0),(7,1),(3,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, BOTTOM, 2, '-k+i-j')
vpairs = [(0,2),(4,3),(7,0),(3,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (WEST, BOTTOM, 3, '-k-j-i')

# Top-North
vpairs = [(5,2),(6,6),(7,7),(4,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, NORTH, 0, '+i+k-j')
vpairs = [(5,3),(6,2),(7,6),(4,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, NORTH, 1, '-k+i-j')
vpairs = [(5,7),(6,3),(7,2),(4,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, NORTH, 2, '-i-k-j')
vpairs = [(5,6),(6,7),(7,3),(4,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, NORTH, 3, '+k-i-j')
# Top-South
vpairs = [(5,0),(6,4),(7,5),(4,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, SOUTH, 0, '-i+k+j')
vpairs = [(5,1),(6,0),(7,4),(4,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, SOUTH, 1, '-k-i+j')
vpairs = [(5,5),(6,1),(7,0),(4,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, SOUTH, 2, '+i-k+j')
vpairs = [(5,4),(6,5),(7,1),(4,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, SOUTH, 3, '+k+i+j')
# Top-East
vpairs = [(5,1),(6,5),(7,6),(4,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, EAST, 0, '-j+k-i')
vpairs = [(5,2),(6,1),(7,5),(4,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, EAST, 1, '-k-j-i')
vpairs = [(5,6),(6,2),(7,1),(4,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, EAST, 2, '+j-k-i')
vpairs = [(5,5),(6,6),(7,2),(4,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, EAST, 3, '+k+j-i')
# Top-West
vpairs = [(5,3),(6,7),(7,4),(4,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, WEST, 0, '+j+k+i')
vpairs = [(5,0),(6,3),(7,7),(4,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, WEST, 1, '-k+j+i')
vpairs = [(5,4),(6,0),(7,3),(4,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, WEST, 2, '-j-k+i')
vpairs = [(5,7),(6,4),(7,0),(4,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, WEST, 3, '+k-j+i')
# Top-Top
vpairs = [(5,4),(6,7),(7,6),(4,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, TOP, 0, '-i+j-k')
vpairs = [(5,5),(6,4),(7,7),(4,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, TOP, 1, '-j-i-k')
vpairs = [(5,6),(6,5),(7,4),(4,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, TOP, 2, '+i-j-k')
vpairs = [(5,7),(6,6),(7,5),(4,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, TOP, 3, '+j+i-k')
# Top-Bottom
vpairs = [(5,1),(6,2),(7,3),(4,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, BOTTOM, 0, '+i+j+k')
vpairs = [(5,0),(6,1),(7,2),(4,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, BOTTOM, 1, '-j+i+k')
vpairs = [(5,3),(6,0),(7,1),(4,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, BOTTOM, 2, '-i-j+k')
vpairs = [(5,2),(6,3),(7,0),(4,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (TOP, BOTTOM, 3, '+j-i+k')

# Bottom-North
vpairs = [(0,2),(3,6),(2,7),(1,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, NORTH, 0, '-i+k+j')
vpairs = [(0,3),(3,2),(2,6),(1,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, NORTH, 1, '+k+i+j')
vpairs = [(0,7),(3,3),(2,2),(1,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, NORTH, 2, '+i-k+j')
vpairs = [(0,6),(3,7),(2,3),(1,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, NORTH, 3, '-k-i+j')
# Bottom-South
vpairs = [(0,0),(3,4),(2,5),(1,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, SOUTH, 0, '+i+k-j')
vpairs = [(0,1),(3,0),(2,4),(1,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, SOUTH, 1, '+k-i-j')
vpairs = [(0,5),(3,1),(2,0),(1,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, SOUTH, 2, '-i-k-j')
vpairs = [(0,4),(3,5),(2,1),(1,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, SOUTH, 3, '-k+i-j')
# Bottom-East
vpairs = [(0,1),(3,5),(2,6),(1,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, EAST, 0, '+j+k+i')
vpairs = [(0,2),(3,1),(2,5),(1,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, EAST, 1, '+k-j+i')
vpairs = [(0,6),(3,2),(2,1),(1,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, EAST, 2, '-j-k+i')
vpairs = [(0,5),(3,6),(2,2),(1,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, EAST, 3, '-k+j+i')
# Bottom-West
vpairs = [(0,3),(3,7),(2,4),(1,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, WEST, 0, '-j+k-i')
vpairs = [(0,0),(3,3),(2,7),(1,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, WEST, 1, '+k+j-i')
vpairs = [(0,4),(3,0),(2,3),(1,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, WEST, 2, '+j-k-i')
vpairs = [(0,7),(3,4),(2,0),(1,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, WEST, 3, '-k-j-i')
# Bottom-Top
vpairs = [(0,4),(3,7),(2,6),(1,5)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, TOP, 0, '+i+j+k')
vpairs = [(0,5),(3,4),(2,7),(1,6)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, TOP, 1, '+j-i+k')
vpairs = [(0,6),(3,5),(2,4),(1,7)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, TOP, 2, '-i-j+k')
vpairs = [(0,7),(3,6),(2,5),(1,4)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, TOP, 3, '-j+i+k')
# Bottom-Bottom
vpairs = [(0,1),(3,2),(2,3),(1,0)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, BOTTOM, 0, '-i+j-k')
vpairs = [(0,0),(3,1),(2,2),(1,3)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, BOTTOM, 1, '+j+i-k')
vpairs = [(0,3),(3,0),(2,1),(1,2)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, BOTTOM, 2, '+i-j-k')
vpairs = [(0,2),(3,3),(2,0),(1,1)]; vpairs.sort()
connectionDict3D[tuple(vpairs)] = (BOTTOM, BOTTOM, 3, '-j-i-k')

# When reading GridPro block connectivity file, 
# we need to look up Eilmer notation for connection orientations.
eilmer_orientation = {}
for vpairs in connectionDict3D.keys():
    this_face, other_face, orientation, axis_map = connectionDict3D[vpairs]
    eilmer_orientation[this_face,other_face,axis_map] = orientation 

# print "eilmer_orientation=", eilmer_orientation

vpairsDict3D = {}
for vpairs in connectionDict3D.keys():
    this_face, other_face, orientation, axis_map = connectionDict3D[vpairs]
    vpairsDict3D[this_face, other_face, orientation] = vpairs

def to_eilmer_axis_map(this_ijk):
    """
    Convert from GridPro axis_map string to Eilmer3 axis_map string.

    From GridPro manual, Section 7.3.2 Connectivity Information.
    Example, 123 --> '+i+j+k'
    """
    gridpro_other_axis = {0:'xx', 1:'+i', 2:'+j', 3:'+k', 4:'-i', 5:'-j', 6:'-k'}
    if type(this_ijk) == type('a'):
        pass
    elif type(this_ijk) == type(0):
        this_ijk = "%03d" % this_ijk
    else:
        print "Expected a string or integer of three digits but got this_ijk=", this_ijk
        sys.exit(-1)
    other_ijk = string.join([gridpro_other_axis[int(d)] for d in this_ijk], '')
    return other_ijk

def get_eilmer_orientation(this_face, other_face, gridpro_axis_map):
    """
    Returns an orientation integer, 0..3, for the orientation of the block connection.

    Note that, when accessing the dictionary of orientations, 
    we convert to the form of keys actually present in the dictionary.
    """
    return eilmer_orientation[faceDict[this_face], faceDict[other_face], 
                              to_eilmer_axis_map(gridpro_axis_map)]

#----------------------------------------------------------------------

def close_enough(vA, vB, tolerance=1.0e-4):
    """
    Decide if two Vector quantities are close enough to being equal.

    This will be used to test that the block corners coincide.
    """
    return (vabs(vA - vB)/(vabs(vA + vB)+1.0)) <= tolerance

#----------------------------------------------------------------------------

def make_patch(north, east, south, west, grid_type="TFI", tol=1.0e-6):
    """
    Defines a 2D patch (or region) by its boundary paths (in order NESW).

    :param north: bounding path on the NORTH side
    :param east: bounding path on the EAST side
    :param south: bounding path on the SOUTH side
    :param west: bounding path on the WEST side
    :param grid_type: indicates the type of interpolation within the patch:

        * TFI, COONS: transfinite interpolation or Coons patch
        * AO: interpolation via an area-orthogonality grid, used as background

    :param tol: relative tolerance for testing coincidence of the corner points

    A patch is defined by its four bounding Path objects
    with assumed positive directions as shown::

                       NORTH
         1       +------->-------+
         |       |               |
         s  WEST |               | EAST
         |       |               |
         0       +------->-------+
                       SOUTH
        
                 0-------r------>1

    NORTH and SOUTH boundaries progress WEST to EAST while
    EAST and WEST boundaries progress SOUTH to NORTH.

    To reuse a Path object when building multiple blocks,
    you will need to pay attention to the orientation of the blocks
    and the defined positive direction of the Path object.
    """
    if not isinstance(north, Path):
        raise TypeError, ("north should be a Path but it is: %s" % type(north))
    if not isinstance(east, Path):
        raise TypeError, ("east should be a Path but it is: %s" % type(east))
    if not isinstance(south, Path):
        raise TypeError, ("south should be a Path but it is: %s" % type(south))
    if not isinstance(west, Path):
        raise TypeError, ("west should be a Path but it is: %s" % type(west))
    # Check that the corner points match.
    p00 = south.eval(0.0); p10 = south.eval(1.0)
    p01 = north.eval(0.0); p11 = north.eval(1.0)
    p00_alt = west.eval(0.0); p01_alt = west.eval(1.0)
    p10_alt = east.eval(0.0); p11_alt = east.eval(1.0)
    average_length = 0.25*(vabs(p10-p00)+vabs(p01-p00)+vabs(p11-p10)+vabs(p11-p01))
    tolerance = tol * (average_length + 1.0)
    corners_OK = True
    if not close_enough(p00, p00_alt, tolerance):
        print "Error: south and west boundaries do not coincide at corner p00."
        print "   ", str(p00), str(p00_alt)
        corners_OK = False
    if not close_enough(p10, p10_alt, tolerance):
        print "Error: south and east boundaries do not coincide at corner p10."
        print "   ", str(p10), str(p10_alt)
        corners_OK = False
    if not close_enough(p01, p01_alt, tolerance):
        print "Error: north and west boundaries do not coincide at corner p01."
        print "   ", str(p01), str(p01_alt)
        corners_OK = False
    if not close_enough(p11, p11_alt, tolerance):
        print "Error: north and east boundaries do not coincide at corner p11."
        print "   ", str(p11), str(p11_alt)
        corners_OK = False
    if not corners_OK:
        print "----------------------------------------------------------------"
        print "Preparation is stopping at this point because there is a problem"
        print "with the geomtric data presented to mesh_patch."
        print "The following stack trace should tell you where this occured"
        print "in your input script. The preceeding data should give you a hint"
        print "about where the problem points are in modelling space."
        print "----------------------------------------------------------------"
        raise RuntimeError("Corner points do not coincide")
    # Finally, proceed with mesh generation.
    grid_type.upper()
    if grid_type == "AO":
        return AOPatch(south, north, west, east)
    else:
        return CoonsPatch(south, north, west, east)

#----------------------------------------------------------------------------

class Block(object):
    """
    Python base class to organise the setting of block parameters.

    We will create instances of its subclasses: Block2D or Block3D.
    """
    # We will accumulate references to defined objects.
    blockList = []
    # Minimum number of cells in any direction
    # The boundary conditions affect two cells into the mesh.
    nmin = 2

    def set_BC(self, face_name, type_of_BC,
               inflow_condition=None, x_order=0, sponge_flag=None,
               Twall=None, Pout=None, filename=None, n_profile=1,
               is_wall=0, sets_conv_flux=0, sets_visc_flux=0,
               assume_ideal=0, mdot=None, epsilon=None,
               Twall_i=None, Twall_f=None, t_i=None, t_f=None,
               label=''):
        """
        Sets a boundary condition on a particular face of the block.

        :param face_name: int identifier to select the appropriate boundary
             within the block.
        :param type_of_BC: Name or index value of the requested boundary
            condition.  See module bc_defs.py for the available options.
        :param inflow_condition: If the type of boundary requires the user to
            specify the inflow FlowCondition object, this is the parameter to do so.
        :param sponge_flag: Set to 1 to activate Andrew Denman's damping layer
            near the boundary.
        :param Twall: If appropriate, specify the boundary-wall temperature 
            in degrees Kelvin.
        :param Pout: If appropriate, specify the value of static pressure
            (in Pascals) just outside the boundary.

        Sometimes it is good to be able to adjust properties after
        block creation; this function provides that capability.
        """
        iface = faceDict[face_name]
        if type(type_of_BC) == type(object()) and type_of_BC in bcName.keys():
            pass # We already have a valid symbol.
        elif type(type_of_BC) == type(1):
            type_of_BC = bcSymbolFromName[type_of_BC]
        elif type(type_of_BC) == type('string'):
            type_of_BC = bcSymbolFromName[str(type_of_BC).upper()]
        else:
            print "Error in setting BC:"
            print "    We've got a type_of_BC value that we don't know."
            print "    type()=", type(type_of_BC)
            print "    value=", str(type_of_BC)
            return
        print "Set block:", self.blkId, "face:", faceName[iface], \
              "BC:", bcName[type_of_BC]
        if inflow_condition != None and not isinstance(inflow_condition, FlowCondition):
            print "Error in setting BC:"
            print "    The inflow_condition argument is not a FlowCondition object."
            return
        if sponge_flag != None:
            sponge_flag = int(sponge_flag)
        else:
            sponge_flag = 0
        if Twall != None:
            Twall = float(Twall)
        if Pout != None:
            Pout = float(Pout)
        if epsilon != None:
            epsilon = float(epsilon)
        if Twall_i != None:
            Twall_i = float(Twall_i)
        if Twall_f != None:
            Twall_f = float(Twall_f)
        if t_i != None:
            t_i = float(t_i)
        if t_f != None:
            t_f = float(t_f)
        # Now, create a new boundary condition object.
        if type_of_BC == ADJACENT:
            print "Error in setting BC:"
            print "    Should not be setting an ADJACENT BC on a single block."
            return
        if type_of_BC == SUP_IN:
            newbc = SupInBC(inflow_condition, label=label)
        if type_of_BC == EXTRAPOLATE_OUT:
            newbc = ExtrapolateOutBC(x_order, sponge_flag, label=label)
        if type_of_BC == SHOCK_FITTING_IN:
            newbc = ShockFittingInBC(inflow_condition, label=label)
        if type_of_BC == SLIP_WALL:
            newbc = SlipWallBC(label=label)
        if type_of_BC == ADIABATIC:
            newbc = AdiabaticBC(label=label)
        if type_of_BC == FIXED_T:
            newbc = FixedTBC(Twall, label=label)
        if type_of_BC == SUBSONIC_IN:
            newbc = SubsonicInBC(inflow_condition, assume_ideal=assume_ideal, label=label)
        if type_of_BC == SUBSONIC_OUT:
            newbc = SubsonicOutBC(sponge_flag, label=label)
        if type_of_BC == TRANSIENT_UNI:
            if not filename: filename = "transient_uniform.dat"
            newbc = TransientUniBC(filename, label=label)
        if type_of_BC == TRANSIENT_PROF:
            newbc = TransientProfBC(label=label)
        if type_of_BC == STATIC_PROF:
            if not filename: filename = "profile.dat"
            newbc = StaticProfBC(filename, n_profile, label=label)
        if type_of_BC == FIXED_P_OUT:
            newbc = FixedPOutBC(Pout, x_order, label=label)
        if type_of_BC == RRM:
            newbc = RRMBC(sponge_flag, label=label)
        if type_of_BC == USER_DEFINED:
            if not filename: filename = "udf.lua"
            newbc = UserDefinedBC(filename, is_wall, sets_conv_flux, sets_visc_flux, label=label)
        if type_of_BC == SEB:
            newbc = SurfaceEnergyBalanceBC(epsilon, label=label)
        if type_of_BC == ABLATING:
            if not filename: filename = ""
            if mdot != None:
            	print "Error in setting BC:"
            	print "    A mass-flux list (mdot) is required to create an AblatingBC"
            	return
            newbc = AblatingBC(Twall, mdot, filename, label=label)
        if type_of_BC == FSTC:
            if not filename: filename = "fstc_temperatures.dat"
            newbc = fstcBC(filename, label=label)
        if type_of_BC == SLIDING_T:
            newbc = SlidingTBC(Twall_i, Twall_f, t_i, t_f, label=label)
        #
        try:
            self.bc_list[iface] = newbc
        except:
            print "Boundary condition not set correctly, type_of_BC=", bcName[type_of_BC]
            sys.exit()
        return

    def set_WBC(self, face_name, type_of_WBC, f_wall=[1.0,], input_file=None, label=''):
        """
        Sets a wall catalytic boundary condition on a particular face of the block.

        :param face_name: String or int identifier to select the appropriate Face2D 
            within the block.
        :param type_of_WBC: Name or index value of the requested wall catalycity 
            boundary condition.  See module cns_bc_defs for the available options.
        :param f_wall: If the user is required to set the chemical composition
            at the wall, then a list of mass fractions (floats) should be 
            supplied as this parameter.
        :param input_file: If the boundary condition requires an input file,
            then its name is supplied as this parameter.

        Sometimes it is good to be able to adjust properties after
        block creation; this function provides that capability.
        """
        iface = faceDict[face_name]
        type_of_WBC = str(type_of_WBC).upper()
        fbc = self.wc_bc_list[iface]
        fbc.type_of_WCBC = wc_bcIndexFromName[type_of_WCBC]
        fbc.label = label
        print "Set block:", self.blkId, "face:", faceName[iface], "WCBC:", type_of_WBC
        fbc.f_wall = copy.copy(f_wall)
        if input_file != None:
            fbc.input_file = input_file
        return

    def write_to_ini_file(self, fp, dimensions):
        """
        Writes the information to the specified file-object in .ini format.
        """
        fp.write("\n[block/%d]\n" % self.blkId)
        fp.write("label = %s\n" % self.label)
        fp.write("active = %d\n" % self.active)
        fp.write("nni = %d\n" % self.nni)
        fp.write("nnj = %d\n" % self.nnj)
        fp.write("nnk = %d\n" % self.nnk)
        fp.write("omegaz = %f\n" % self.omegaz)
        if self.hcell_list != None:
            fp.write("nhcell = %d\n" % len(self.hcell_list))
            count = 0
            for hcell in self.hcell_list:
                if dimensions == 3:
                    fp.write("history-cell-%d = %d %d %d\n" % 
                             (count,hcell[0],hcell[1],hcell[2]))
                else:
                    fp.write("history-cell-%d = %d %d\n" % 
                             (count,hcell[0],hcell[1]))
                count += 1
        if dimensions == 3:
            faceList = faceList3D
        else:
            faceList = faceList2D
        for iface in faceList:
            fp.write("[block/%d/face/%s]\n" % (self.blkId, faceName[iface]))
            bc = self.bc_list[iface]
            fp.write("label = %s\n" % bc.label)
            fp.write("bc = %s\n" % bcName[bc.type_of_BC])
            if bc.inflow_condition == None:
                inflow_indx = 0
            else:
                inflow_indx = bc.inflow_condition.indx
            fp.write("inflow_condition = %d\n" % inflow_indx)
            fp.write("filename = %s\n" % bc.filename)
            fp.write("n_profile = %d\n" % bc.n_profile)
            fp.write("x_order = %d\n" % bc.x_order)
            fp.write("sponge_flag = %d\n" % bc.sponge_flag)
            fp.write("xforce_flag = %d\n" % self.xforce_list[iface])
            fp.write("Twall = %e\n" % bc.Twall)
            fp.write("Pout = %e\n" % bc.Pout)
            fp.write("is_wall = %d\n" % bc.is_wall)
            fp.write("sets_conv_flux = %d\n" % bc.sets_conv_flux)
            fp.write("sets_visc_flux = %d\n" % bc.sets_visc_flux)
            fp.write("assume_ideal = %d\n" % bc.assume_ideal)
            fp.write("mdot = ")
            for val in bc.mdot:
                fp.write("%.6e " % val )
            fp.write("\n")
            fp.write("epsilon = %e\n" % bc.epsilon)
            #
            fp.write("other_block = %d\n" % bc.other_block)
            if bc.other_face >= 0:
                fp.write("other_face = %s\n" % faceName[bc.other_face])
            else:
                fp.write("other_face = none\n")
            fp.write("neighbour_orientation = %d\n" % bc.orientation)
            if dimensions == 3:
                # I think that this neighbour-vtx list is redundant, PJ.
                fp.write("neighbour-vtx = ")
                for ivtx in range(8):
                    fp.write(" %d" % self.neighbour_vertex_list[iface][ivtx])
                fp.write("\n");
            #
            wbc = self.wc_bc_list[iface]
            fp.write("wcbc = %d\n" % wbc.type_of_WCBC )
            if wbc.input_file != None:
                fp.write("wcbc_input_file = %s \n" % wbc.input_file )
            fp.write("f_wall = ")
            for val in wbc.f_wall:
                fp.write("%.6e " % val )
            fp.write("\n")
            fp.write("Twall_i = %e\n" % bc.Twall_i)
            fp.write("Twall_f = %e\n" % bc.Twall_f)
            fp.write("t_i = %e\n" % bc.t_i)
            fp.write("t_f = %e\n" % bc.t_f)
        return

    def write_starting_solution(self, fp, gdata, fb=None):
        """
        Writes the flow-solution in Eilmer3-native format to a file.

        This is almost-Tecplot format.
        """
        print "Write initial flow solution: label=", self.label
        if isinstance(self.fill_condition, FlowCondition):
            label = self.fill_condition.flow.label
            if len(label) == 0: label = "<no-label>"
            print "   using FlowCondition:", label
        elif isinstance(self.fill_condition, ExistingSolution):
            print "   using ExistingSolution:", self.fill_condition.rootName
        elif callable(self.fill_condition):
            print "   using user-defined Python function:", self.fill_condition.__name__
        fp.write("%20.12e\n" % gdata.t0)
        fp.write("%s\n" % quoted_string(variable_list_for_cell(gdata)))
        fp.write("%d %d %d\n" % (self.nni, self.nnj, self.nnk)) # number of cells in each dir
        
        if fb:
            fb.write("%20.12e\n" % gdata.t0)
            fb.write("%s\n" % quoted_string(bgk_list_for_cell(gdata)))
            fb.write("%d %d %d\n" % (self.nni, self.nnj, self.nnk)) # number of cells in each dir
            
        # 
        # Initialise for instances when ExistingSolution() is used, so that 
        # we can start the search from the last located block and cell.
        i_found = 0; j_found = 0; k_found = 0; jb_found = 0
        for k in range(self.nnk):
            for j in range(self.nnj):
                for i in range(self.nni):
                    x, y, z, vol = self.cell_centre_location(i, j, k, gdata)
                    if isinstance(self.fill_condition, FlowCondition):
                        # A simple FlowCondition object knows how to write its
                        # properties to a dictionary.
                        cell_properties = self.fill_condition.to_dict()
                        cell_properties['pos.x'] = x
                        cell_properties['pos.y'] = y
                        cell_properties['pos.z'] = z
                        cell_properties['volume'] = vol
                    elif isinstance(self.fill_condition, ExistingSolution):
                        # An ExistingSolution object also has similar capabilities.
                        if self.fill_condition.assume_same_grid:
                            interp_values = self.fill_condition.interpolate_flow_condition(x, y, z, vol, i, j, k, self.blkId )
                        else:
                            interp_values = self.fill_condition.interpolate_flow_condition(x, y, z, vol, i_found, j_found, k_found, jb_found ) 
                        
                        if fb:
                            cell_properties, bgk_properties, i_found, j_found, k_found, jb_found = interp_values
                        else:
                            cell_properties, i_found, j_found, k_found, jb_found = interp_values                 

                    elif callable(self.fill_condition):
                        # We expect that the user-supplied function will accept
                        # x,y,z coordinate values and return a dictionary of
                        # flow-condition values.
                        # See FlowCondition.to_dict() method for the expected elements
                        # in this dictionary.
                        cell_properties = self.fill_condition(x, y, z)
                        cell_properties['pos.x'] = x
                        cell_properties['pos.y'] = y
                        cell_properties['pos.z'] = z
                        cell_properties['volume'] = vol
                    else:
                        print "We have a problem in Block.write_starting_solution()."
                        print "    fill_condition is of type", type(self.fill_condition)
                        print "    and we don't know what to do."
                        raise RuntimeError, "Problem with fill_condition."
                    write_cell_data(fp, cell_properties, gdata)
                    if fb:
                        write_bgk_data(fb, bgk_properties, gdata)
                print ".",
                sys.stdout.flush()
        print
        # print "End write_initial_solution for block."
        return

#----------------------------------------------------------------------------

class Block2D(Block):
    """
    Python class to organise the setting of block parameters for 2D flow.
    """

    __slots__ = 'blkId', 'label', 'psurf', 'nni', 'nnj', 'nnk', 'grid', \
                'hcell_list', 'fill_condition', 'active', \
                'cf_list', 'bc_list', 'vtx', \
                'xforce_list', 'wc_bc_list', 'omegaz'
    
    def __init__(self,
                 psurf=None,
                 grid=None,
                 import_grid_file_name=None,
                 nni=2,
                 nnj=2,
                 cf_list=[None,]*4,
                 bc_list=[SlipWallBC(),]*4,
                 wc_bc_list=[NonCatalyticWBC(),]*4,
                 fill_condition=None,
                 hcell_list=[],
                 xforce_list = [0,]*4,
                 label="",
                 active=1
                 ):
        """
        Create a block from a parametric-surface object.

        You should specify one of the following three:

        :param psurf: The ParametricSurface object defining the block in space.
            Typically, this will be a CoonsPatch or an AOPatch.
        :param grid: A StructuredGrid object may be supplied directly.
        :param import_grid_file_name: name of a VTK file containing the grid

        :param nni: number of cells in the i-direction (west to east)
        :param nnj: number of cells in the j-direction (south to north)
        :param cf_list: List of the cluster_function objects, one for each boundary.
            The order within the list is NORTH, EAST, SOUTH and WEST.
        :param bc_list: List of BoundaryCondition objects, one for each face.
        :param fill_condition: Either a single FlowCondition or user-defined function.
        :param hcell_list: List of (i,j) tuples specifying the cells
            (for this block)
            whose flow data is to be recorded in the history file.
            For an MPI simulation, there is one history file for each
            block but, for a shared-memory simulation, the history cells
            for all blocks are written to a single history file.
        :param xforce_list: list of int flags to indicate that we want 
            boundary forces calculated
        :param label: Optional string label that will appear 
            in the generated parameter file.
        :param active: flag that indicates if the block is active:

            * =1 (default) the time integration operates for this block
            * =0 time integration for this block is suppressed

        Note: 
            The blocks are given their identity (counting from zero)
            according to the order in which they are created by the user's script.
            This identity is stored as blkId and is used internally
            in the preprocessing, simulation and postprocessing stages.
            For the MPI simulations, it also the same as the rank of the process.
        """
        if not isinstance(nni, int):
            raise TypeError, ("nni should be an int but it is: %s" % type(nni))
        if not isinstance(nnj, int):
            raise TypeError, ("nnj should be an int but it is: %s" % type(nnj))
        #
        self.blkId = len(Block.blockList)    # next available block index
        if len(label) == 0:
            label = "blk-" + str(self.blkId)
        self.label = label
        self.fill_condition_function = None
        self.active = active
        #
        # The grid may come from one of several sources, in order of priority:
        # 1. discretization of a parametric surface
        # 2. supplied directly as a StructuredGrid object
        # 3. read from a VTK file.
        print "Block2D.__init__: blkId=%d, label=%s" % (self.blkId, self.label)
        if psurf != None:
            print "Generate a grid from a user-specified parametric surface."
            self.psurf = psurf.clone() # needed for later rendering
            self.nni = nni
            self.nnj = nnj
            assert len(cf_list) == 4
            self.grid = StructuredGrid((self.nni+1, self.nnj+1))
            self.grid.make_grid_from_surface(self.psurf, cf_list)
        elif grid != None:
            print "Accept grid as given."
            # Should assert grid is StructuredGrid or e3_grid.StructuredGrid
            self.psurf = None # will be a problem for later rendering
            self.grid = grid
            self.nni = self.grid.ni - 1
            self.nnj = self.grid.nj - 1
        elif import_grid_file_name != None:
            print "Import a grid from a VTK data file:", import_grid_file_name
            self.grid = StructuredGrid()
            fin = open(import_grid_file_name, "r")
            self.grid.read_block_in_VTK_format(fin)
            fin.close()
            self.nni = self.grid.ni - 1
            self.nnj = self.grid.nj - 1
        else:
            raise ValueError("Block2D constructor was expecting one of psurf or grid.")
        self.vtx = [self.grid.get_vertex_coords(ivtx) for ivtx in range(4)]
        #
        assert self.nni >= Block.nmin
        assert self.nnj >= Block.nmin
        self.nnk = 1
        assert len(bc_list) == 4
        #
        # In 2D simulations, the rotational frame is not allowed.
        self.omegaz = 0.0
        #
        # Make copies of supplied lists in case we are given the same
        # (default) empty list for each block
        self.bc_list = copy.copy(bc_list)
        self.wc_bc_list = copy.copy(wc_bc_list)
        self.cf_list = copy.copy(cf_list)
        self.hcell_list = copy.copy(hcell_list)
        self.xforce_list = copy.copy(xforce_list)
        self.fill_condition = fill_condition
        #
        Block.blockList.append(self)
        return

    def cell_centre_location(self, i, j, k, gdata):
        """
        Return the cell geometry.

        Geometry of cell::

           ^ j
           |
           |
           NW-----NE
           |       |
           |       |
           SW-----SE  --> i
        """
        # Should match the code in C++ function calc_volumes_2D() in block.cxx.
        k = 0
        xSE = self.grid.x[i+1,j,k];   ySE = self.grid.y[i+1,j,k]
        xNE = self.grid.x[i+1,j+1,k]; yNE = self.grid.y[i+1,j+1,k]
        xNW = self.grid.x[i,j+1,k];   yNW = self.grid.y[i,j+1,k]
        xSW = self.grid.x[i,j,k];     ySW = self.grid.y[i,j,k]
        # Cell area in the (x,y)-plane.
        xyarea = 0.5 * ((xNE + xSE) * (yNE - ySE) + (xNW + xNE) * (yNW - yNE) +
                        (xSW + xNW) * (ySW - yNW) + (xSE + xSW) * (ySE - ySW))
        # Cell Centroid.
        centre_x = 1.0 / (xyarea * 6.0) * \
            ((yNE - ySE) * (xSE * xSE + xSE * xNE + xNE * xNE) + 
             (yNW - yNE) * (xNE * xNE + xNE * xNW + xNW * xNW) +
             (ySW - yNW) * (xNW * xNW + xNW * xSW + xSW * xSW) + 
             (ySE - ySW) * (xSW * xSW + xSW * xSE + xSE * xSE))
        centre_y = -1.0 / (xyarea * 6.0) * \
            ((xNE - xSE) * (ySE * ySE + ySE * yNE + yNE * yNE) + 
             (xNW - xNE) * (yNE * yNE + yNE * yNW + yNW * yNW) +
             (xSW - xNW) * (yNW * yNW + yNW * ySW + ySW * ySW) + 
             (xSE - xSW) * (ySW * ySW + ySW * ySE + ySE * ySE))
        #
        if gdata.axisymmetric_flag == 1:
            # volume per radian
            vol = xyarea * centre_y
        else:
            # volume per unit depth in z
            vol = xyarea
        #
        return (centre_x, centre_y, 0.0, vol)
    
def connect_blocks_2D(A, faceA, B, faceB, with_udf=0, 
                      filename=None, is_wall=0,
                      sets_conv_flux=0, sets_visc_flux=0):
    """
    Make the face-to-face connection between neighbouring blocks.

    :param A: first Block2D object
    :param faceA: indicates which face of block A is to be connected.
        The constants NORTH, EAST, SOUTH, and WEST may be convenient to use.
    :param B: second Block2D object
    :param faceB: indicates which face of block B is to be connected.
        The constants NORTH, EAST, SOUTH, and WEST may be convenient to use.
    """
    assert isinstance(A, Block2D)
    assert isinstance(B, Block2D)
    assert faceA in faceList2D
    assert faceB in faceList2D
    print "connect block", A.blkId, "face", faceName[faceA], \
          "to block", B.blkId, "face", faceName[faceB]
    if with_udf:
        # Exchange connection with user-defined function.
        A.bc_list[faceA] = AdjacentPlusUDFBC(B.blkId, faceB, filename=filename, 
                                             is_wall=is_wall,
                                             sets_conv_flux=sets_conv_flux, sets_visc_flux=sets_visc_flux)
        B.bc_list[faceB] = AdjacentPlusUDFBC(A.blkId, faceA, filename=filename, 
                                             is_wall=is_wall,
                                             sets_conv_flux=sets_conv_flux, sets_visc_flux=sets_visc_flux)
    else:
        # Classic exchange connection.
        A.bc_list[faceA] = AdjacentBC(B.blkId, faceB)
        B.bc_list[faceB] = AdjacentBC(A.blkId, faceA)
    return


def identify_block_connections_2D(block_list=None, exclude_list=[], tolerance=1.0e-6):
    """
    Identifies and makes block connections based on block-vertex locations.

    :param block_list: list of Block2D objects that are to be included in the search.
        If none is supplied, the whole collection is searched.
        This allows one to specify a limited selection of blocks
        to be connected.
    :param exclude_list: list of pairs of Block2D objects that should not be connected
    :param tolerance: spatial tolerance for colocation of vertices
    """
    if block_list == None:
        # Default to searching all defined blocks.
        block_list = Block.blockList
    print "Begin searching for block connections..."
    for thisBlock in block_list:
        for otherBlock in block_list:
            if thisBlock == otherBlock: continue
            inExcludeList = (exclude_list.count((thisBlock,otherBlock)) + \
                             exclude_list.count((otherBlock,thisBlock))) > 0
            if not inExcludeList:
                connections = 0
                #
                if (vabs(thisBlock.vtx[NE] - otherBlock.vtx[SW]) < tolerance) and \
                   (vabs(thisBlock.vtx[NW] - otherBlock.vtx[NW]) < tolerance) :
                    connect_blocks_2D(thisBlock, NORTH, otherBlock, WEST)
                    connections += 1
                if (vabs(thisBlock.vtx[NE] - otherBlock.vtx[NW]) < tolerance) and \
                   (vabs(thisBlock.vtx[NW] - otherBlock.vtx[NE]) < tolerance) :
                    connect_blocks_2D(thisBlock, NORTH, otherBlock, NORTH)
                    connections += 1
                if (vabs(thisBlock.vtx[NE] - otherBlock.vtx[NE]) < tolerance) and \
                   (vabs(thisBlock.vtx[NW] - otherBlock.vtx[SE]) < tolerance) :
                    connect_blocks_2D(thisBlock, NORTH, otherBlock, EAST)
                    connections += 1
                if (vabs(thisBlock.vtx[NE] - otherBlock.vtx[SE]) < tolerance) and \
                   (vabs(thisBlock.vtx[NW] - otherBlock.vtx[SW]) < tolerance) :
                    connect_blocks_2D(thisBlock, NORTH, otherBlock, SOUTH)
                    connections += 1
                #
                if (vabs(thisBlock.vtx[SE] - otherBlock.vtx[SW]) < tolerance) and \
                   (vabs(thisBlock.vtx[NE] - otherBlock.vtx[NW]) < tolerance) :
                    connect_blocks_2D(thisBlock, EAST, otherBlock, WEST)
                    connections += 1
                if (vabs(thisBlock.vtx[SE] - otherBlock.vtx[NW]) < tolerance) and \
                   (vabs(thisBlock.vtx[NE] - otherBlock.vtx[NE]) < tolerance) :
                    connect_blocks_2D(thisBlock, EAST, otherBlock, NORTH)
                    connections += 1
                if (vabs(thisBlock.vtx[SE] - otherBlock.vtx[NE]) < tolerance) and \
                   (vabs(thisBlock.vtx[NE] - otherBlock.vtx[SE]) < tolerance) :
                    connect_blocks_2D(thisBlock, EAST, otherBlock, EAST)
                    connections += 1
                if (vabs(thisBlock.vtx[SE] - otherBlock.vtx[SE]) < tolerance) and \
                   (vabs(thisBlock.vtx[NE] - otherBlock.vtx[SW]) < tolerance) :
                    connect_blocks_2D(thisBlock, EAST, otherBlock, SOUTH)
                    connections += 1
                #
                if (vabs(thisBlock.vtx[SW] - otherBlock.vtx[SW]) < tolerance) and \
                   (vabs(thisBlock.vtx[SE] - otherBlock.vtx[NW]) < tolerance) :
                    connect_blocks_2D(thisBlock, SOUTH, otherBlock, WEST)
                    connections += 1
                if (vabs(thisBlock.vtx[SW] - otherBlock.vtx[NW]) < tolerance) and \
                   (vabs(thisBlock.vtx[SE] - otherBlock.vtx[NE]) < tolerance) :
                    connect_blocks_2D(thisBlock, SOUTH, otherBlock, NORTH)
                    connections += 1
                if (vabs(thisBlock.vtx[SW] - otherBlock.vtx[NE]) < tolerance) and \
                   (vabs(thisBlock.vtx[SE] - otherBlock.vtx[SE]) < tolerance) :
                    connect_blocks_2D(thisBlock, SOUTH, otherBlock, EAST)
                    connections += 1
                if (vabs(thisBlock.vtx[SW] - otherBlock.vtx[SE]) < tolerance) and \
                   (vabs(thisBlock.vtx[SE] - otherBlock.vtx[SW]) < tolerance) :
                    connect_blocks_2D(thisBlock, SOUTH, otherBlock, SOUTH)
                    connections += 1
                #
                if (vabs(thisBlock.vtx[NW] - otherBlock.vtx[SW]) < tolerance) and \
                   (vabs(thisBlock.vtx[SW] - otherBlock.vtx[NW]) < tolerance) :
                    connect_blocks_2D(thisBlock, WEST, otherBlock, WEST)
                    connections += 1
                if (vabs(thisBlock.vtx[NW] - otherBlock.vtx[NW]) < tolerance) and \
                   (vabs(thisBlock.vtx[SW] - otherBlock.vtx[NE]) < tolerance) :
                    connect_blocks_2D(thisBlock, WEST, otherBlock, NORTH)
                    connections += 1
                if (vabs(thisBlock.vtx[NW] - otherBlock.vtx[NE]) < tolerance) and \
                   (vabs(thisBlock.vtx[SW] - otherBlock.vtx[SE]) < tolerance) :
                    connect_blocks_2D(thisBlock, WEST, otherBlock, EAST)
                    connections += 1
                if (vabs(thisBlock.vtx[NW] - otherBlock.vtx[SE]) < tolerance) and \
                   (vabs(thisBlock.vtx[SW] - otherBlock.vtx[SW]) < tolerance) :
                    connect_blocks_2D(thisBlock, WEST, otherBlock, SOUTH)
                    connections += 1
                #
                if connections > 0:
                    # Avoid doubling-up with the reverse connections.
                    exclude_list.append((thisBlock,otherBlock))
    print "Finished searching for block connections."
    return

# --------------------------------------------------------------------

class MultiBlock2D(object):
    """
    Allows us to specify a block of sub-blocks.

    A number of internally-connected Block2D objects will be created when
    one MultiBlock2D object is created.
    Individual blocks occupy subsections of the original parametric surface.
   
    Note that the collection of Block2D objects will be stored in
    a list of lists with each inner-list storing a j-column of blocks::
    
        .                       North
        .   1       +-------------+-------------+
        .   |       |   [0][1]    |   [1][1]    |
        .   s  West +-------------+-------------+ East
        .   |       |   [0][0]    |   [1][0]    |
        .   0       +-------------+-------------+
        .                       South
        .           0           --r-->          1

    The user script may access an individual block within the MultiBlock2D
    object as object.blks[i][j].
    This will be useful for connecting blocks within the MultiBlock cluster
    to other blocks as defined in the user's script.

    Some properties, such as fill_conditions and grid_type, will be propagated
    to all sub-blocks.  Individual sub-blocks can be later customised.
    """

    __slots__ = 'psurf', 'bc_list', 'nb_w2e', 'nb_s2n', 'nn_w2e',\
                'nn_s2n', 'cluster_w2e', 'cluster_s2n', 'fill_condition', \
                'grid_type', 'split_single_grid', 'label', 'blks', \
                'wc_bc_list'
    
    def __init__(self,
                 psurf=None,
                 nni=None,
                 nnj=None,
                 bc_list=[SlipWallBC(), SlipWallBC(), SlipWallBC(), SlipWallBC()],
                 wc_bc_list=[NonCatalyticWBC(), NonCatalyticWBC(),
                             NonCatalyticWBC(), NonCatalyticWBC()],
                 nb_w2e=1,
                 nb_s2n=1,
                 nn_w2e=None,
                 nn_s2n=None,
                 cluster_w2e=None,
                 cluster_s2n=None,
                 fill_condition=None,
                 label="blk",
                 active=1):
        """
        Create a cluster of blocks within an original parametric surface.

        :param psurf: ParametricSurface which defines the block.
        :param bc_list: List of boundary condition objects
            The order within the list is NORTH, EAST, SOUTH and WEST.
        :param nb_w2e: Number of sub-blocks from west to east.
        :param nb_s2n: Number of sub-blocks from south to north.
        :param nn_w2e: List of discretisation values for north and south
            boundaries of the sub-blocks.
            If a list is not supplied, the original number of cells for the
            outer boundary is divided over the individual sub-block boundaries.
        :param nn_s2n: List of discretisation values for west and east
            boundaries of the sub-blocks.
            If a list is not supplied, the original number of cells for the
            outer boundary is divided over the individual sub-block boundaries.
        :param cluster_w2e: If a list of cluster functions is supplied,
            individual clustering will be applied to the corresponding
            south and north boundaries of each sub-block.
            If not supplied, a default of no clustering will be applied.
        :param cluster_s2n: If a list of cluster functions is supplied,
            individual clustering will be applied to the corresponding
            west and east boundaries of each sub-block.
            If not supplied, a default of no clustering will be applied.
        :param fill_condition: A single FlowCondition object that is to be
            used for all sub-blocks
        :param grid_type: Select the type of grid generator from TFI or AO.
        :param split_single_grid: If this boolean flag is true, a single grid is
            generated which is then subdivided into the required blocks.
        :param label: A label that will be augmented with the sub-block index
            and then used to label the individual Block2D objects.
        :param active: (int) flag indicating if the block is active:

            *  =1 (default) the time integration operates for this block
            *  =0 time integration for this block is suppressed
        """
        default_cluster_function = LinearFunction()
        self.blks = []
        dt_s2n = 1.0 / nb_s2n
        dt_w2e = 1.0 / nb_w2e
        for i in range(nb_w2e):
            self.blks.append([])  # for a new j-column of blocks
            for j in range(nb_s2n):
                sublabel = label + "-" + str(i) + "-" + str(j)
                new_psurf = psurf.clone()
                #
                new_psurf.s0 = j * dt_s2n
                new_psurf.s1 = (j+1) * dt_s2n
                new_psurf.r0 = i * dt_w2e
                new_psurf.r1 = (i+1) * dt_w2e
                #
                # Transfer the boundary conditions.
                bc_newlist=[SlipWallBC(), SlipWallBC(), SlipWallBC(), SlipWallBC()]
                wc_bc_newlist = [NonCatalyticWBC(), NonCatalyticWBC(),
                                 NonCatalyticWBC(), NonCatalyticWBC()]
                if i == 0:
                    bc_newlist[WEST] = copy.copy(bc_list[WEST])
                    wc_bc_newlist[WEST] = copy.copy(wc_bc_list[WEST])
                if i == nb_w2e - 1:
                    bc_newlist[EAST] = copy.copy(bc_list[EAST])
                    wc_bc_newlist[EAST] = copy.copy(wc_bc_list[EAST])
                if j == 0:
                    bc_newlist[SOUTH] = copy.copy(bc_list[SOUTH])
                    wc_bc_newlist[SOUTH] = copy.copy(wc_bc_list[SOUTH])
                if j == nb_s2n - 1:
                    bc_newlist[NORTH] = copy.copy(bc_list[NORTH])
                    wc_bc_newlist[NORTH] = copy.copy(wc_bc_list[NORTH])
                #
                # For discretization, take the list of cell numbers
                # if it available, else divide overall number of nodes.
                if nn_s2n:
                    new_nnj = nn_s2n[j]
                else:
                    new_nnj = nnj / nb_s2n
                if nn_w2e:
                    new_nni = nn_w2e[i]
                else:
                    new_nni = nni / nb_w2e
                #
                cf_newlist = []
                try:
                    cf_newlist.append(cluster_w2e[i])
                except:
                    cf_newlist.append(default_cluster_function)
                try:
                    cf_newlist.append(cluster_s2n[j])
                except:
                    cf_newlist.append(default_cluster_function)
                try:
                    cf_newlist.append(cluster_w2e[i])
                except:
                    cf_newlist.append(default_cluster_function)
                try:
                    cf_newlist.append(cluster_s2n[j])
                except:
                    cf_newlist.append(default_cluster_function)
                #
                new_blk = Block2D(psurf=new_psurf, nni=new_nni, nnj=new_nnj,
                                  cf_list=cf_newlist, bc_list=bc_newlist,
                                  wc_bc_list=wc_bc_newlist,
                                  fill_condition=fill_condition,
                                  label=sublabel, active=active)
                self.blks[i].append(new_blk)
        #
        # print "blocks:", self.blks
        # Make the internal, sub-block to sub-block connections.
        if nb_w2e > 1:
            for i in range(1,nb_w2e):
                for j in range(nb_s2n):
                    connect_blocks_2D(self.blks[i-1][j], EAST, self.blks[i][j], WEST)
        if nb_s2n > 1:
            for i in range(nb_w2e):
                for j in range(1,nb_s2n):
                    connect_blocks_2D(self.blks[i][j-1], NORTH, self.blks[i][j], SOUTH)
        #
        return
    
# --------------------------------------------------------------------

def subdivide_vertex_range(ncells, nblocks):
    """
    Subdivide one index-direction into a number of sub-blocks within a SuperBlock.

    :param ncells  : number of cells in one index-direction
    :param nblocks : number of sub-blocks in the same direction
    
    :returns: a list of tuples specifying the subranges of the vertices 
        that form the sub-blocks.
    """
    vtx_subranges = []
    vtx_min = 0
    for iblk in range(nblocks):
        n_subblock = int(ncells / (nblocks - iblk))  # try to divide remaining cells uniformly
        vtx_max = vtx_min + n_subblock
        vtx_subranges.append((vtx_min, vtx_max))
        vtx_min = vtx_max
        ncells -= n_subblock # remaining cells to be sub-divided
    return vtx_subranges


class SuperBlock2D(object):
    """
    Creates a single grid over the region and then subdivides that grid.

    .. Original implementation by Rowan; refactored (lightly) by PJ Nov-2010.
    """

    __slots__ = 'psurf', 'grid', 'bc_list', 'nni', 'nnj', 'nbi', 'nbj', 'cf_list',\
                'fill_condition', 'label', 'blks', 'wc_bc_list'
    
    def __init__(self,
                 psurf=None,
                 grid=None,
                 nni=2,
                 nnj=2,
                 nbi=1,
                 nbj=1,
                 cf_list=[None, None, None, None],
                 bc_list=[SlipWallBC(), SlipWallBC(), SlipWallBC(), SlipWallBC()],
                 wc_bc_list=[NonCatalyticWBC(), NonCatalyticWBC(),
                             NonCatalyticWBC(), NonCatalyticWBC()],
                 fill_condition=None,
                 hcell_list=[],
                 label="sblk",
                 active=1
                 ):
        """
        Creates a single grid over the region and then subdivides that grid.

        On return, self.blks is a (nested) list-of-lists of subblock references.
        """
        self.blks = []
        # 1. Create the large grid for the super-block
        if psurf != None:
            print "Generate a grid from a user-specified parametric surface."
            self.grid = StructuredGrid((nni+1, nnj+1))
            self.grid.make_grid_from_surface(psurf, cf_list)
        elif grid != None:
            print "Accept grid as given, but subdivide."
            self.grid = grid
            nni = self.grid.ni-1
            nnj = self.grid.nj-1
        # 2. Create lists of the subgrid indices
        si_list = subdivide_vertex_range(nni, nbi)
        sj_list = subdivide_vertex_range(nnj, nbj)
        # 3. Do the actual subdivision.
        for i in range(nbi):
            self.blks.append([])  # for a new j-column of blocks
            for j in range(nbj):
                sublabel = label + "-" + str(i) + "-" + str(j)
                imin = si_list[i][0]; imax = si_list[i][1]
                jmin = sj_list[j][0]; jmax = sj_list[j][1]
                subgrid = self.grid.create_subgrid(imin, imax, jmin, jmax)
                #
                # Transfer the boundary conditions.
                bc_newlist=[SlipWallBC(), SlipWallBC(), SlipWallBC(), SlipWallBC()]
                wc_bc_newlist = [NonCatalyticWBC(), NonCatalyticWBC(),
                                 NonCatalyticWBC(), NonCatalyticWBC()]
                if i == 0:
                    bc_newlist[WEST] = copy.copy(bc_list[WEST])
                    wc_bc_newlist[WEST] = copy.copy(wc_bc_list[WEST])
                if i == nbi - 1:
                    bc_newlist[EAST] = copy.copy(bc_list[EAST])
                    wc_bc_newlist[EAST] = copy.copy(wc_bc_list[EAST])
                if j == 0:
                    bc_newlist[SOUTH] = copy.copy(bc_list[SOUTH])
                    wc_bc_newlist[SOUTH] = copy.copy(wc_bc_list[SOUTH])
                if j == nbj - 1:
                    bc_newlist[NORTH] = copy.copy(bc_list[NORTH])
                    wc_bc_newlist[NORTH] = copy.copy(wc_bc_list[NORTH])
                #
                new_blk = Block2D(grid=subgrid, bc_list=bc_newlist,
                                  nni=imax-imin, nnj=jmax-jmin, 
                                  wc_bc_list=wc_bc_newlist,
                                  fill_condition=fill_condition,
                                  label=sublabel, active=active)
                self.blks[i].append(new_blk)
        # print "blocks:", self.blks
        #
        # 4. Make the internal, sub-block to sub-block connections.
        if nbi > 1:
            for i in range(1,nbi):
                for j in range(nbj):
                    connect_blocks_2D(self.blks[i-1][j], EAST, self.blks[i][j], WEST)
        if nbj > 1:
            for i in range(nbi):
                for j in range(1,nbj):
                    connect_blocks_2D(self.blks[i][j-1], NORTH, self.blks[i][j], SOUTH)
        #
        return
    
# --------------------------------------------------------------------

class Block3D(Block):
    """
    Organises the setting of block parameters for 3D geometry.
    """
    __slots__ = 'blkId', 'nni', 'nnj', 'nnk', 'label', \
                'parametric_volume', 'grid', 'cf_list', 'hcell_list', \
                'fill_condition', 'omegaz', 'active', \
                'bc_list', 'Twall_list', 'Pout_list', \
                'xforce_list', 'wc_bc_list', 'sponge_flag_list', \
                'vertex_location_list', 'neighbour_vertex_list'
              
    def __init__(self,
                 parametric_volume=None,
                 grid=None,
                 import_grid_file_name=None,
                 nni=None,
                 nnj=None,
                 nnk=None,
                 cf_list=[None,]*12,
                 bc_list=[SlipWallBC(),]*6,
                 wc_bc_list=[NonCatalyticWBC(),]*6,
                 fill_condition=None,
                 hcell_list=None,
                 xforce_list=[0,]*6,
                 label="",
                 active=1,
                 omegaz=0.0
                 ):
        """
        Basic initialisation for a block.

        | The order of the cluster_function_list elements are
        | edge 0 is from p0 -> p1  (i-direction, bottom surface)
        | edge 1         p1 -> p2  (j-direction, bottom surface)
        |      2         p3 -> p2  (i-direction, bottom surface)
        |      3         p0 -> p3  (j-direction, bottom surface)
        |      4         p4 -> p5  (i-direction, top surface)
        |      5         p5 -> p6  (j-direction, top surface)
        |      6         p7 -> p6  (i-direction, top surface)
        |      7         p4 -> p7  (j-direction, top surface)
        |      8         p0 -> p4  (k-direction)
        |      9         p1 -> p5  (k-direction)
        |      10        p2 -> p6  (k-direction)
        |      11        p3 -> p7  (k-direction)

        The definitive list is in the function make_TFI_grid_from_volume() in
        the Python module e3_grid.py.
        
        See below for other methods to provide details of the boundary conditions
        and the inter-block connections.
        """
        #
        # The blocks are given their identity according to the
        # order in which they are created by the user's script.
        self.blkId = len(Block3D.blockList)    # next available block index
        if len(label) == 0:
            label = "blk-" + str(self.blkId)
        self.label = label
        print "Defining blkId=%d, label=%s" % (self.blkId, self.label)
        #
        self.active = active
        self.wc_bc_list = copy.copy(wc_bc_list)
        self.xforce_list = copy.copy(xforce_list)
        self.fill_condition=fill_condition
        #
        # Rotating frame of reference.
        self.omegaz = omegaz
        #
        # The grid may come from one of several sources, in order of priority:
        # 1. discretization of a parametric volume
        # 2. supplied directly as a StructuredGrid object
        # 3. read from a VTK file.
        if parametric_volume != None:
            print "Generate a grid from a user-specified parametric volume."
            print "    parametric subrange: r0=", parametric_volume.r0, \
                  "r1=", parametric_volume.r1,
            print "s0=", parametric_volume.s0, "s1=", parametric_volume.s1,
            print "t0=", parametric_volume.t0, "t1=", parametric_volume.t1
            self.parametric_volume = parametric_volume.copy()
            self.cf_list = copy.copy(cf_list)
            if nni != None: 
                self.nni = nni
            else:
                self.nni = Block.nmin
            if nnj != None: 
                self.nnj = nnj
            else:
                self.nnj = Block.nmin
            if nnk != None: 
                self.nnk = nnk
            else:
                self.nnk = Block.nmin
            self.grid = StructuredGrid((self.nni+1, self.nnj+1, self.nnk+1))
            print "Generate TFI grid for block."
            self.grid.make_TFI_grid_from_volume(parametric_volume, cf_list)
        elif grid != None:
            print "Use the grid as given."
            # Should assert grid is StructuredGrid or e3_grid.StructuredGrid
            self.grid = grid
            self.nni = self.grid.ni-1
            self.nnj = self.grid.nj-1
            self.nnk = self.grid.nk-1
            self.parametric_volume = None
        elif import_grid_file_name != None:
            print "Import a grid from a VTK data file:", import_grid_file_name
            self.grid = StructuredGrid()
            fin = open(import_grid_file_name, "r")
            self.grid.read_block_in_VTK_format(fin)
            fin.close()
            self.nni = self.grid.ni-1
            self.nnj = self.grid.nj-1
            self.nnk = self.grid.nk-1
            self.parametric_volume = None
        else:
            raise ValueError("Block3D did not receive a suitable data to get/create a grid.")
        # Because of the way reflection boundary conditions are implemented,
        # we need a minimum number of cells across each dimension.
        assert self.nni >= Block.nmin
        assert self.nnj >= Block.nmin
        assert self.nnk >= Block.nmin
        # Gather the vertex coordinates for later use in connecting blocks.
        self.vtx = [self.grid.get_vertex_coords(ivtx) for ivtx in range(8)]
        #
        if hcell_list != None:
            for hcell in hcell_list:
                assert isinstance(hcell, tuple) and len(hcell) == 3,\
                    "History cells should be specified as a list of tuples (i,j,k)"
        self.hcell_list = copy.copy(hcell_list)
        #
        # Dummy values for the boundary conditions.
        # Other values may be selected via the set_BC() method.
        self.bc_list = copy.copy(bc_list)
        self.Twall_list = [300.0,] * 6
        self.Pout_list = [100.0e3,] * 6
        self.sponge_flag_list = [0,] * 6
        #
        # Connections to other blocks are initialized to None
        # self.neighbour_block_list = [None,] * 6
        # self.neighbour_faceId_list = [-1,] * 6
        # self.neighbour_orientation_list = [0,] * 6
        # Each vertex is shared by three faces of the block and
        # each of these faces may be connected to another block.
        # For the moment, it is easier to over-do the storage
        # and only index into part of it.  That way, a vertex always
        # has the same index, no matter which face we are considering.
        self.neighbour_vertex_list = [[-1,] * 8,
                                      [-1,] * 8,
                                      [-1,] * 8,
                                      [-1,] * 8,
                                      [-1,] * 8,
                                      [-1,] * 8] # ensure separate elements
        #
        Block.blockList.append(self)
        return

    def cell_centre_location(self, i, j, k, gdata):
        """
        Cell geometry in 3D.
        """
        p0 = Vector(self.grid.x[i,j,k],       self.grid.y[i,j,k],       self.grid.z[i,j,k])
        p1 = Vector(self.grid.x[i+1,j,k],     self.grid.y[i+1,j,k],     self.grid.z[i+1,j,k])
        p2 = Vector(self.grid.x[i+1,j+1,k],   self.grid.y[i+1,j+1,k],   self.grid.z[i+1,j+1,k])
        p3 = Vector(self.grid.x[i,j+1,k],     self.grid.y[i,j+1,k],     self.grid.z[i,j+1,k])
        p4 = Vector(self.grid.x[i,j,k+1],     self.grid.y[i,j,k+1],     self.grid.z[i,j,k+1])
        p5 = Vector(self.grid.x[i+1,j,k+1],   self.grid.y[i+1,j,k+1],   self.grid.z[i+1,j,k+1])
        p6 = Vector(self.grid.x[i+1,j+1,k+1], self.grid.y[i+1,j+1,k+1], self.grid.z[i+1,j+1,k+1])
        p7 = Vector(self.grid.x[i,j+1,k+1],   self.grid.y[i,j+1,k+1],   self.grid.z[i,j+1,k+1])
        centre = hexahedron_centroid(p0, p1, p2, p3, p4, p5, p6, p7)
        vol = hexahedron_volume(p0, p1, p2, p3, p4, p5, p6, p7)
        return (centre.x, centre.y, centre.z, vol)

# --------------------------------------------------------------------

class MultiBlock3D(object):
    """
    Allows us to specify a block of sub-blocks.

    A number of internally-connected Block3D objects will be created when
    one MultiBlock2D object is created.
    Individual blocks occupy subsections of the original parametric volume.
    
    Note that the collection of Block3D objects will be stored in
    a list of lists of lists with each inner-most list storing
    a k-column of blocks.  For example, the k=constant layer will have
    a plan view as shown here::
    
        .                       North
        .   1       +-------------+-------------+
        .   |       |  [0][1][k]  |  [1][1][k]  |
        .   s  West +-------------+-------------+ East
        .   |       |  [0][0][k]  |  [1][0][k]  |
        .   0       +-------------+-------------+
        .                       South
        .           0           --r-->          1

    The user script may access an individual block within the MultiBlock3D
    object as object.blks[i][j][k].
    This will be useful for connecting blocks within the MultiBlock cluster
    to other blocks as defined in the user's script.

    Some properties, such as fill_condition, will be propagated to all sub-blocks.
    Individual sub-blocks can be later customised.
    """

    __slots__ = 'pvolume', 'nbi', 'nbj', 'nbk', 'nni', 'nnj', 'nnk', \
                'clusteri', 'clusterj', 'clusterk', 'fill_condition', \
                'label', 'blks'
    
    def __init__(self,
                 parametric_volume=None,
                 nbi=1, nbj=1, nbk=1,
                 nni=None, nnj=None, nnk=None,
                 clusteri=None, clusterj=None, clusterk=None,
                 fill_condition=None,
                 label="blk",
                 omegaz=0.0,
                 active=1
                 ):
        """
        Create a cluster of blocks within a given parametric volume.

        :param parametric_volume: The total volume that will be subdivided.
        :param nbi: integer number of sub-blocks in the i-direction (i.e. from west to east).
        :param nbj: integer number of sub-blocks in the j-direction (i.e. from south to north).
        :param nbk: integer number of sub-blocks in the k-direction (i.e. from south to north).
        :param nni: List of integer discretisation values along the i-direction edges
            of the sub-blocks.
        :param nnj: List of integer discretisation values for j-direction edges.
        :param nnk: List of integer discretisation values for k-direction edges.
        :param clusteri: If a list of cluster function objects is supplied,
            individual clustering will be applied to the corresponding
            i-direction edges of each sub-block.
            If not supplied, a default of no clustering will be applied.
        :param clusterj: If a list of cluster function objects is supplied,
            individual clustering will be applied to the corresponding
            j-direction edges of each sub-block.
            If not supplied, a default of no clustering will be applied.
        :param clusterk: If a list of cluster function objects is supplied,
            individual clustering will be applied to the corresponding
            k-direction edges of each sub-block.
            If not supplied, a default of no clustering will be applied.
        :param fill_condition: A single FlowCondition|ExistingSolution|function 
            that is to be used for all sub-blocks
        :param label: A string label that will be augmented with the sub-block index
            and then used to label the individual Block3D objects.
        :param active: flag to indicate if block is active:

            * =1 (default) the time integration operates for this block
            * =0 time integration for this block is suppressed
        """
        #
        # Number of blocks along each parametric direction.
        assert isinstance(nbi, int) and (nbi >= 1), \
            "MultiBlock3D.__init__: Expected nbi as int >= 1."
        assert isinstance(nbj, int) and (nbj >= 1), \
            "MultiBlock3D.__init__: Expected nbj as int >= 1."
        assert isinstance(nbk, int) and (nbk >= 1), \
            "MultiBlock3D.__init__: Expected nbk as int >= 1."
        #
        dr = 1.0 / nbi
        ds = 1.0 / nbj
        dt = 1.0 / nbk
        #
        # Discretization along each parametric direction.
        assert (isinstance(nni, list) and len(nni) == nbi) \
               or isinstance(nni, int), \
            "MultiBlock3D.__init__: Expected a list of int or a single int for nni."
        if isinstance(nni, int):
            # Split the number of cells equally.
            nni_per_block = int(math.ceil(float(nni) / nbi))
            nni = [nni_per_block,] * nbi
        #
        assert (isinstance(nnj, list) and len(nnj) == nbj) \
               or isinstance(nnj, int), \
            "MultiBlock3D.__init__: Expected a list of int or a single int for nnj."
        if isinstance(nnj, int):
            # Split the number of cells equally.
            nnj_per_block = int(math.ceil(float(nnj) / nbj))
            nnj = [nnj_per_block,] * nbj
        #
        assert (isinstance(nnk, list) and len(nnk) == nbk) \
               or isinstance(nnk, int), \
            "MultiBlock3D.__init__: Expected a list of int or a single int for nnk."
        if isinstance(nnk, int):
            # Split the number of cells equally.
            nnk_per_block = int(math.ceil(float(nnk) / nbk))
            nnk = [nnk_per_block,] * nbk
        #
        # Assemble the collection of blocks.
        self.blks = []
        flat_list = []
        for i in range(nbi):
            r_min = 0.0 + i * dr
            r_max = r_min + dr
            slab_of_blocks = []
            for j in range(nbj):
                s_min = 0.0 + j * ds
                s_max = s_min + ds
                column_of_blocks = []
                for k in range(nbk):
                    sublabel = label + "-" + str(i) + "-" + str(j) + "-" + str(k)
                    print "MultiBlock3D.__init__: assembling block:", sublabel
                    t_min = 0.0 + k * dt
                    t_max = t_min + dt
                    pv = parametric_volume.copy()
                    pv.r0 = r_min; pv.r1 = r_max
                    pv.s0 = s_min; pv.s1 = s_max
                    pv.t0 = t_min; pv.t1 = t_max
                    if isinstance(clusteri, list):
                        cfi = clusteri[i]
                    else:
                        cfi = clusteri
                    if isinstance(clusterj, list):
                        cfj = clusterj[j]
                    else:
                        cfj = clusterj
                    if isinstance(clusterk, list):
                        cfk = clusterk[k]
                    else:
                        cfk = clusterk
                    cf_list = [cfi, cfj, cfi, cfj, cfi, cfj, cfi, cfj, cfk, cfk, cfk, cfk]
                    new_block = Block3D(pv,
                                        nni=nni[i], nnj=nnj[j], nnk=nnk[k],
                                        cf_list=cf_list,
                                        fill_condition=fill_condition,
                                        label=sublabel, omegaz=omegaz, active=active)
                    column_of_blocks.append(new_block)
                    flat_list.append(new_block)
                slab_of_blocks.append(column_of_blocks)
            self.blks.append(slab_of_blocks)
        #
        ## print "blocks:", self.blks
        # Make the internal, sub-block to sub-block connections.
        identify_block_connections_3D(flat_list)
        return

# --------------------------------------------------------------------

class SuperBlock3D(object):
    """
    Creates a single grid over the region and then subdivides that grid.

    Implemented as an extension of SuperBlock2D to the third dimension.
    
    .. Wilson Chan 08-Dec-2008; Refactored (lightly) by PJ Nov-2010.
    """

    __slots__ = 'parametric_volume', 'bc_list', 'nni', 'nnj', 'nnk',\
                'nbi', 'nbj', 'nbk', 'cf_list', 'fill_condition',\
                'label', 'blks', 'wc_bc_list'

    def __init__(self,
                 parametric_volume=None,
                 nni=2,
                 nnj=2,
                 nnk=2,
                 nbi=1,
                 nbj=1,
                 nbk=1,
                 cf_list=[None,]*12,
                 bc_list=[SlipWallBC(), SlipWallBC(), SlipWallBC(),
                          SlipWallBC(), SlipWallBC(), SlipWallBC()],
                 wc_bc_list=[NonCatalyticWBC(), NonCatalyticWBC(),
                             NonCatalyticWBC(), NonCatalyticWBC(),
                             NonCatalyticWBC(), NonCatalyticWBC()],
                 fill_condition=None,
                 hcell_list=[],
                 label="sblk",
                 omegaz = 0.0,
                 active=1
                 ):
        """
        Creates a single grid over the region and then subdivides that grid.
        
        On return self.blks holds a list-of-lists collection of sub-blocks
        """
        # 1. Create the large grid for the super-block
        grid = StructuredGrid((nni+1, nnj+1, nnk+1))
        grid.make_TFI_grid_from_volume(parametric_volume, cf_list)
        # 2. Create lists of the subgrid indices
        si_list = subdivide_vertex_range(nni, nbi)
        sj_list = subdivide_vertex_range(nnj, nbj)
        sk_list = subdivide_vertex_range(nnk, nbk)
        # 3. Do the actual subdivision.
        self.blks = []
        flat_list = []
        for i in range(nbi):
            slab_of_blocks = []
            for j in range(nbj):
                column_of_blocks = []
                for k in range(nbk):
                    sublabel = label + "-" + str(i) + "-" + str(j) + "-" + str(k)
                    imin = si_list[i][0]; imax = si_list[i][1]
                    jmin = sj_list[j][0]; jmax = sj_list[j][1]
                    kmin = sk_list[k][0]; kmax = sk_list[k][1]
                    subgrid = grid.create_subgrid(imin, imax, jmin, jmax, kmin, kmax)
                    # Transfer the boundary conditions.
                    bc_newlist=[SlipWallBC(), SlipWallBC(), SlipWallBC(),
                                SlipWallBC(), SlipWallBC(), SlipWallBC()]
                    wc_bc_newlist = [NonCatalyticWBC(), NonCatalyticWBC(),
                                     NonCatalyticWBC(), NonCatalyticWBC(),
                                     NonCatalyticWBC(), NonCatalyticWBC()]
                    if i == 0:
                        bc_newlist[WEST] = copy.copy(bc_list[WEST])
                        wc_bc_newlist[WEST] = copy.copy(wc_bc_list[WEST])
                    if i == nbi - 1:
                        bc_newlist[EAST] = copy.copy(bc_list[EAST])
                        wc_bc_newlist[EAST] = copy.copy(wc_bc_list[EAST])
                    if j == 0:
                        bc_newlist[SOUTH] = copy.copy(bc_list[SOUTH])
                        wc_bc_newlist[SOUTH] = copy.copy(wc_bc_list[SOUTH])
                    if j == nbj - 1:
                        bc_newlist[NORTH] = copy.copy(bc_list[NORTH])
                        wc_bc_newlist[NORTH] = copy.copy(wc_bc_list[NORTH])
                    if k == 0:
                        bc_newlist[BOTTOM] = copy.copy(bc_list[BOTTOM])
                        wc_bc_newlist[BOTTOM] = copy.copy(wc_bc_list[BOTTOM])
                    if k == nbk - 1:
                        bc_newlist[TOP] = copy.copy(bc_list[TOP])
                        wc_bc_newlist[TOP] = copy.copy(wc_bc_list[TOP])
                    # Create the actual subblock.
                    new_blk = Block3D(grid=subgrid, bc_list=bc_newlist,
                                      nni=imax-imin, nnj=jmax-jmin, nnk=kmax-kmin,
                                      wc_bc_list=wc_bc_newlist,
                                      fill_condition=fill_condition,
                                      label=sublabel, omegaz=omegaz, active=active)
                    column_of_blocks.append(new_blk)
                    flat_list.append(new_blk)
                slab_of_blocks.append(column_of_blocks)
            self.blks.append(slab_of_blocks)
        # Exclicitly connect blocks together within this SuperBlock because the
        # automatic search breaks when surfaces are allowed to collapse to lines.
        # identify_block_connections(flat_list)
        if nbi > 1:
            # Make connections of faces EAST to WEST
            for i in range(nbi-1):
                for j in range(nbj):
                    for k in range(nbk):
                        connect_blocks_3D(self.blks[i][j][k], self.blks[i+1][j][k],
                                          ((2,3),(6,7),(5,4),(1,0)))
        if nbj > 1:
            # Make connections of faces NORTH to SOUTH
            for i in range(nbi):
                for j in range(nbj-1):
                    for k in range(nbk):
                        connect_blocks_3D(self.blks[i][j][k], self.blks[i][j+1][k],
                                          ((3,0),(7,4),(6,5),(2,1)))
        if nbk > 1:
            # Make connections of faces TOP to BOTTOM
            for i in range(nbi):
                for j in range(nbj):
                    for k in range(nbk-1):
                        connect_blocks_3D(self.blks[i][j][k], self.blks[i][j][k+1],
                                          ((5,1),(6,2),(7,3),(4,0)))
        return

#----------------------------------------------------------------------

def identify_block_connections_3D(block_list=None, exclude_list=[], tolerance=1.0e-6):
    """
    Identifies and makes block connections based on vertex locations.

    :param block_list: list of Block3D objects that are to be included in the search.
       If none is supplied, the whole collection is searched.
       This allows one to specify a limited selection of blocks to be connected.
    :param exclude_list: list of pairs of Block3D objects that should not be connected
    :param tolerance: spatial tolerance for colocation of vertices
    """
    if block_list == None:
        # Default to searching all defined blocks.
        block_list = Block.blockList
    print "Begin searching for block connections..."
    for thisBlock in block_list:
        for otherBlock in block_list:
            if thisBlock == otherBlock: continue
            inExcludeList = (exclude_list.count((thisBlock,otherBlock)) + \
                             exclude_list.count((otherBlock,thisBlock))) > 0
            if not inExcludeList:
                vtxPairList = identify_colocated_vertices(thisBlock, otherBlock, tolerance)
                if len(vtxPairList) == 4:
                    # We only make simple single-face-to-single-face connections.
                    connect_blocks_3D(thisBlock, otherBlock, vtxPairList)
                    exclude_list.append((thisBlock,otherBlock))
    print "Finished searching for block connections."
    return

def identify_colocated_vertices(A, B, tolerance):
    """
    Identify colocated vertices by looking at their position is 3D space.

    :param A: Block3D object
    :param B: Block3D object
    :param tolerance: Vertices are considered to be colocated if their Euclidian distance
        is less than tolerance.
    """
    from math import sqrt
    vtxPairList = []
    for iA in range(8):
        for iB in range(8):
            if vtxPairList.count((iA,iB)) > 0: continue
            if vabs(A.vtx[iA] - B.vtx[iB]) <= tolerance: vtxPairList.append((iA,iB))
    return vtxPairList

def connect_blocks_3D(A, B, vtx_pairs, with_udf=0, 
                      filename=None, is_wall=0,
                      sets_conv_flux=0, sets_visc_flux=0):
    """
    Make the specified vertex-to-vertex connection between neighbouring blocks.

    :param A: Block3D object
    :param B: Block3D object
    :param vtxPairs: list of 4 pairs of vertex indices specifying the corresponding corners. 
    """
    print "connect_blocks(): begin..."
    if not isinstance(vtx_pairs, list):
        vtx_pairs = list(vtx_pairs)
    if len(vtx_pairs) != 4:
        raise Exception, "connect_blocks(): Incorrect pairs: %s" % str(vtx_pairs)
    vtx_pairs.sort()  # want pairs in a standard order for dictionary lookup
    try:
        faceA, faceB, orientation, axis_map = connectionDict3D[tuple(vtx_pairs)]
        print "connectionDict3D lookup: faceA=", faceName[faceA], \
              "faceB=", faceName[faceB], "orientation=", orientation, "axis_map=", axis_map
    except KeyError:
        print "connect_blocks_3D(): error"
        print "   The following set of vertex pairs is not in"
        print "   the dictionary of known interblock connections."
        print "   vtx_pairs=", vtx_pairs
        faceA = -1; faceB = -1; orientation = 0
        sys.exit(-1)
    print "connect block", A.blkId, "face", faceName[faceA], \
          "to block", B.blkId, "face", faceName[faceB]
    if with_udf:
        # Exchange connection with user-defined function
        A.bc_list[faceA] = AdjacentPlusUDFBC(B.blkId, faceB, orientation, filename=filename, 
                                             is_wall=is_wall,
                                             sets_conv_flux=sets_conv_flux, sets_visc_flux=sets_visc_flux)
        B.bc_list[faceB] = AdjacentPlusUDFBC(A.blkId, faceA, orientation, filename=filename, 
                                             is_wall=is_wall,
                                             sets_conv_flux=sets_conv_flux, sets_visc_flux=sets_visc_flux)
    else:
        # Classic exchange connection.
        A.bc_list[faceA] = AdjacentBC(B.blkId, faceB, orientation)
        B.bc_list[faceB] = AdjacentBC(A.blkId, faceA, orientation)
    #
    if not check_block_connection_3D(A, faceA, B, faceB, orientation):
        print "connect_blocks(): Block vertex locations not consistent for this connection."
        print "   vtx_pairs=", vtx_pairs
    if not cell_count_consistent_3D(A, faceA, B, faceB, orientation):
        print "connect_blocks(): Cell counts not consistent for this connection."
        print "   vtx_pairs=", vtx_pairs
    #
    for vtxA, vtxB in vtx_pairs:
        A.neighbour_vertex_list[faceA][vtxA] = vtxB
        B.neighbour_vertex_list[faceB][vtxB] = vtxA
    print "connect_blocks(): done."
    return

def check_block_connection_3D(blkA, faceA, blkB, faceB, orientation, tolerance=1.0e-6):
    """
    Returns True if blocks connect at specified faces.
    """
    vtx_pairs = vpairsDict3D[faceA, faceB, orientation]
    for iA, iB in vtx_pairs:
        if vabs(blkA.vtx[iA] - blkB.vtx[iB]) <= tolerance:
            continue
        else:
            print "Corners of grid do not match."
            print "blkA= ", blkA.blkId, " vtxA= ", iA, " pos= ", blkA.vtxA[iA]
            print "blkB= ", blkB.blkId, " vtxB= ", iB, " pos= ", blkB.vtxB[iB]
            return False
    #
    # If we get this far, then all is ok.
    return True

def cell_count_consistent_3D(blkA, faceA, blkB, faceB, orientation):
    """
    Returns logical 1 if the cell-discretization matches for the faces.
    """
    result_flag = 0
    #
    if (faceA == NORTH or faceA == SOUTH) and (faceB == NORTH or faceB == SOUTH):
        if orientation == 0 or orientation == 2:
            result_flag = (blkA.nni == blkB.nni) and (blkA.nnk == blkB.nnk)
        elif orientation == 1 or orientation == 3:
            result_flag = (blkA.nnk == blkB.nni) and (blkA.nni == blkB.nnk)
    if (faceA == NORTH or faceA == SOUTH) and (faceB == EAST or faceB == WEST):
        if orientation == 0 or orientation == 2:
            result_flag = (blkA.nni == blkB.nnj) and (blkA.nnk == blkB.nnk)
        elif orientation == 1 or orientation == 3:
            result_flag = (blkA.nnk == blkB.nnj) and (blkA.nni == blkB.nnk)
    if (faceA == NORTH or faceA == SOUTH) and (faceB == TOP or faceB == BOTTOM):
        if orientation == 0 or orientation == 2:
            result_flag = (blkA.nni == blkB.nni) and (blkA.nnk == blkB.nnj)
        elif orientation == 1 or orientation == 3:
            result_flag = (blkA.nnk == blkB.nni) and (blkA.nni == blkB.nnj)
    #
    if (faceA == EAST or faceA == WEST) and (faceB == NORTH or faceB == SOUTH):
        if orientation == 0 or orientation == 2:
            result_flag = (blkA.nnj == blkB.nni) and (blkA.nnk == blkB.nnk)
        elif orientation == 1 or orientation == 3:
            result_flag = (blkA.nnk == blkB.nni) and (blkA.nnj == blkB.nnk)
    if (faceA == EAST or faceA == WEST) and (faceB == WEST or faceB == EAST):
        if orientation == 0 or orientation == 2:
            result_flag = (blkA.nnj == blkB.nnj) and (blkA.nnk == blkB.nnk)
        elif orientation == 1 or orientation == 3:
            result_flag = (blkA.nnk == blkB.nnj) and (blkA.nnj == blkB.nnk)
    if (faceA == EAST or faceA == WEST) and (faceB == TOP or faceB == BOTTOM):
        if orientation == 0 or orientation == 2:
            result_flag = (blkA.nnj == blkB.nni) and (blkA.nnk == blkB.nnj)
        elif orientation == 1 or orientation == 3:
            result_flag = (blkA.nnk == blkB.nni) and (blkA.nnj == blkB.nnj)
    #
    if (faceA == TOP or faceA == BOTTOM) and (faceB == NORTH or faceB == SOUTH):
        if orientation == 0 or orientation == 2:
            result_flag = (blkA.nni == blkB.nni) and (blkA.nnj == blkB.nnk)
        elif orientation == 1 or orientation == 3:
            result_flag = (blkA.nni == blkB.nnk) and (blkA.nnj == blkB.nni)
    if (faceA == TOP or faceA == BOTTOM) and (faceB == EAST or faceB == WEST):
        if orientation == 0 or orientation == 2:
            result_flag = (blkA.nni == blkB.nnj) and (blkA.nnj == blkB.nnk)
        elif orientation == 1 or orientation == 3:
            result_flag = (blkA.nni == blkB.nnk) and (blkA.nnj == blkB.nnj)
    if (faceA == TOP or faceA == BOTTOM) and (faceB == BOTTOM or faceB == TOP):
        if orientation == 0 or orientation == 2:
            result_flag = (blkA.nni == blkB.nni) and (blkA.nnj == blkB.nnj)
        elif orientation == 1 or orientation == 3:
            result_flag = (blkA.nni == blkB.nnj) and (blkA.nnj == blkB.nni)
    #
    return result_flag

#-------------------------------------------------------------------------------------

def identify_block_connections(block_list=None, exclude_list=[], tolerance=1.0e-6):
    """
    Identifies and makes block connections based on vertex locations.

    :param block_list: list of Block3D or Block2D objects that are to be included in the search.
       If none is supplied, the whole collection is searched.
       This allows one to specify a limited selection of blocks to be connected.
    :param exclude_list: list of pairs of Block3D objects that should not be connected
    :param tolerance: spatial tolerance for colocation of vertices

    Note that this function is just a proxy for the specialized 2D and 3D functions.
    """
    if isinstance(Block.blockList[0], Block3D):
        identify_block_connections_3D(block_list, exclude_list, tolerance)
    else:
        identify_block_connections_2D(block_list, exclude_list, tolerance)
    return

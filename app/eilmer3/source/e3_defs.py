"""
e3_defs.py: Definitions of the faces and vertex names for 2D and 3D blocks and cells.

This set of definitions is used by various other e3_xxxx.py modules
and may also be used by the user's input script.

.. Author: PJ
.. Versions:
   29-Nov-2010: extracted from e3_block.py
"""

# Face names and Ids
NORTH  = 0
EAST   = 1
SOUTH  = 2
WEST   = 3
TOP    = 4
BOTTOM = 5

# We sometimes want to look up the face identity based on the vertex ids.
northVtx  = [2, 3, 7, 6]; northVtx.sort()
eastVtx   = [1, 2, 6, 5]; eastVtx.sort()
southVtx  = [1, 0, 4, 5]; southVtx.sort()
westVtx   = [0, 3, 7, 4]; westVtx.sort()
topVtx    = [4, 5, 6, 7]; topVtx.sort()
bottomVtx = [0, 1, 2, 3]; bottomVtx.sort()

faceDict = { "NORTH": NORTH, "north": NORTH, NORTH: NORTH, tuple(northVtx): NORTH,
             "EAST": EAST, "east": EAST, EAST: EAST, tuple(eastVtx): EAST,
             "SOUTH": SOUTH, "south": SOUTH, SOUTH: SOUTH, tuple(southVtx): SOUTH,
             "WEST": WEST, "west": WEST, WEST: WEST, tuple(westVtx): WEST,
             "TOP": TOP, "top": TOP, TOP: TOP, tuple(topVtx): TOP,
             "BOTTOM": BOTTOM, "bottom": BOTTOM, BOTTOM: BOTTOM, tuple(bottomVtx): BOTTOM }

# The following list specifies the order for lists of information
# about the faces (and should remain consistent with mbcns2).
faceList2D = [NORTH, EAST, SOUTH, WEST]
faceList3D = [NORTH, EAST, SOUTH, WEST, TOP, BOTTOM]
faceName = ["north", "east", "south", "west", "top", "bottom"]

# For 2D geometry, it has been convenient to use compass points to label vertics.
SW = 0
SE = 1
NE = 2
NW = 3

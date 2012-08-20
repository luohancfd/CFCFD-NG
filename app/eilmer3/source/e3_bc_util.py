"""
e3_bc_util.py : Python utility functions to help with setting boundary conditions.

This module is a catch-all place for functions that support converting
boundary condition information from some external source.

.. Author: Rowan J. Gollan

.. Versions:
   16-Aug-2012 initial coding
   19-Aug-2012 Add support for parsing Gridpro .conn file.
"""
import sys

from e3_defs import *
from e3_flow import FlowCondition
from bc_defs import AdjacentBC
from e3_block import eilmer_orientation

#-----------------------------------------------------------------------

def apply_gridpro_connectivity(fname, blks):
    """
    Apply block connections from Gridpro connectivity file.
    
    :param fname: File name of Gridprop connectivity file.
    :param blks: a list of Block object(s)
    """
    fname = open(fname, 'r')
    while True:
        line = f.readline()
        if line.startswith("#"): continue
        else: break
    nb = int(line.split()[0])
    conns = []
    for ib in range(nb):
        cons.append({})
        while True:
            line = f.readline()
            if line.startswith("#"): continue
            else: break
        tks = line.split()
        # Work on faces in order.
        # Gridpro imin ==> Eilmer WEST face
        other_blk = int(tks[4])
        if other_blk > 0: # there is a connection
            conns[ib]["WEST"] = (other_blk-1, tks[5])
        # Gridpro imax ==> Eilmer EAST face
        other_blk = int(tks[8])
        if other_blk > 0:
            conns[ib]["EAST"] = (other_blk-1, tks[9])
        # Gridpro jmin ==> Eilmer SOUTH face
        other_blk = int(tks[12])
        if other_blk > 0:
            conns[ib]["SOUTH"] = (other_blk-1, tks[13])
        # Gridpro jmax ==> Eilmer NORTH face
        other_blk = int(tks[16])
        if other_blk > 0:
            conns[ib]["NORTH"] = (other_blk-1, tks[17])
        # Gridpro kmin ==> Eilmer BOTTOM face
        other_blk = int(tks[20])
        if other_blk > 0:
            conns[ib]["BOTTOM"] = (other_blk-1, tks[21])
        # Gridpro kmax ==> Eilmer TOP face
        other_blk = int(tks[24])
        if other_blk > 0:
            conns[ib]["TOP"] = (other_blk-1, tks[25])
    #
    for ib in range(nb):
        for faceA, (oblk, axis_map) in conns[ib].items():
            A = blks[ib]
            B = blks[oblk]
            # Retrieve connection info about other block, then
            # delete other info (as its no longer needed)
            for face, (oblk2, axis_map2) in conns[oblk].items():
                if oblk2 == ib:
                    faceB = face
                    break
            # delete entry for that face from other block
            conns[oblk].pop(faceB)
            orientation = eilmer_orientation(faceA, faceB, axis_map)
            # Connect blocks
            A.bc_list[faceA] = AdjacentBC(B.blkId, faceB, orientation)
            B.bc_list[faceB] = AdjacentBC(A.blkId, faceA, orientation)
    #
    return


def apply_gridpro_bcs(fname, blks, bc_map):
    """
    Apply Gridpro boundary conditions from file to group of blocks.
    
    :param fname: File name of Gridpro property file
    :param blks: A list of Block object(s)
    :param bc_map: A dict containing the extra information required
       to set certain boundary conditions. Keys are strings which
       correspond to the Eilmer names for boundary conditions, and
       the values are the extra parameter information, which is 
       dependent on boundary condition type. For example, to supply
       information about a supersonic inflow and a fixed temperature
       wall, the dict would look like:
           {'SUP_IN':inflow, 'FIXED_T':450.0} 
       Note that the keys must conform to the Eilmer names (in bc_defs.py)
       and the values vary depending on the boundary condition.
    """
    f = open(fname, 'r')
    nb = int(f.readline().split()[0])
    if nb != len(blks):
        print "Error in apply_gridpro_bcs(): mismatch in numbers of blocks."
        print "The number of blocks described in the Gridpro property file is", nb
        print "But the number of blocks supplied in the block list is", len(blks)
        print "Bailing out!"
        sys.exit(1)
    #
    bcs = []
    for ib in range(nb):
        while True:
            line = f.readline()
            if line.startswith("#"):
                continue
            else:
                break
        tks = line.split()
        bcs.append({'WEST' : int(tks[4]), 'EAST' : int(tks[6]),
                    'SOUTH' : int(tks[8]), 'NORTH' : int(tks[10]),
                    'BOTTOM' : int(tks[12]), 'TOP' : int(tks[14])})
    #
    nlabels = int(f.readline().split()[0])
    for il in range(nlabels):
        f.readline()
    #
    nbc_types = int(f.readline().split()[0])
    bc_type_map = {}
    for ibc in range(nbc_types):
        line = f.readline()
        if line.startswith("#"):
            continue
        tks = line.split()
        bc_type_map[int(tks[0])] = tks[1]
    #
    f.close()
    # All the information is gathered.
    # So now loop over the blocks and apply bcs as appropriate.
    for ib, blk in enumerate(blks):
        for face, bc_id in bcs[ib].items():
            # Now for some ugly switching based on boundary
            # condition type. There is probably some metaprogramming
            # way to do this more nicely, but I don't presently
            # know what that is in Python.
            bc_type = bc_type_map[bc_id]
            if bc_type == 'SUP_IN':
                # Need to check that a flow condition was supplied
                assert(isinstance(bc_map['SUP_IN'], FlowCondition))
                blk.set_BC(face, 'SUP_IN', inflow_condition=bc_map['SUP_IN'])
            elif bc_type == 'FIXED_T':
                assert('FIXED_T' in bc_map)
                blk.set_BC(face, 'FIXED_T', Twall=bc_map['FIXED_T'])
            elif bc_type == 'FIXED_P_OUT':
                assert('FIXED_P_OUT' in bc_map)
                blk.set_BC(face, 'FIXED_P_OUT', Pout=bc_map['FIXED_P_OUT'])
            elif bc_type == 'USER_DEFINED':
                assert('USER_DEFINED' in bc_map)
                blk.set_BC(face, 'USER_DEFINED',
                           filename=bc_map['USER_DEFINED']['filename'],
                           is_wall=bc_map['USER_DEFINED'].get('is_wall', 0),
                           use_udf_flux=bc_map['USER_DEFINED'].get('use_udf_flux', 0))
            elif bc_type == 'INTERBLK':
                # This should already be set when block connections were
                # identified.
                pass
            else:
                # For all other cases, no special information is required.
                # FIX ME: Need to work out what to do for Adjacent+UDF
                blk.set_BC(face, bc_type)
    return
            
            
    

    

    

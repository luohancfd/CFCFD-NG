/** \file exch3d.hh
 * \ingroup eilmer3
 * \brief Prototypes for functions to copy boundary data from one 3D block to another.
 *
 * \author PJ
 * \version August 2004.
 * \version July 2008 Elmer3 port.
 */

#ifndef EXCH3D_HH
#define EXCH3D_HH

#include "block.hh"

/* Direct-copy exchange functions */

int copy_boundary_data_3D( size_t jb, int type_of_copy );

int copy_into_east_boundary_3D(Block *bp, Block *bp_src, int type_of_copy);
int copy_into_west_boundary_3D(Block *bp, Block *bp_src, int type_of_copy);
int copy_into_north_boundary_3D(Block *bp, Block *bp_src, int type_of_copy);
int copy_into_south_boundary_3D(Block *bp, Block *bp_src, int type_of_copy);
int copy_into_top_boundary_3D(Block *bp, Block *bp_src, int type_of_copy);
int copy_into_bottom_boundary_3D(Block *bp, Block *bp_src, int type_of_copy);

/* Copy-via-buffer exchange functions */

int copy_into_send_buffer_3D(Block *bp, int bndry, int type_of_copy, double *send_buffer);
int copy_from_receive_buffer_3D(Block *bp, int bndry, int type_of_copy, double *receive_buffer);
int copy_from_receive_buffer_to_north(Block *bp, int type_of_copy, double *receive_buffer);
int copy_from_receive_buffer_to_south(Block *bp, int type_of_copy, double *receive_buffer);
int copy_from_receive_buffer_to_east(Block *bp, int type_of_copy, double *receive_buffer);
int copy_from_receive_buffer_to_west(Block *bp, int type_of_copy, double *receive_buffer);
int copy_from_receive_buffer_to_top(Block *bp, int type_of_copy, double *receive_buffer);
int copy_from_receive_buffer_to_bottom(Block *bp, int type_of_copy, double *receive_buffer);

#endif

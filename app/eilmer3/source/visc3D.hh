/** \file visc3D.hh
 * \ingroup eilmer3
 * \brief Prototypes for functions to compute viscous fluxes in 3D.
 *
 * \author Andrew Denman
 * \version August 2004 bring code over from mb_cns.
 * \version November 2008 port to Elmer3, PJ
 */

#ifndef VISC3D_HH

#include "cell.hh"
#include "block.hh"

int viscous_flux_3D(Block *A);
int viscous_derivatives_3D(Block *A, size_t time_level);
int viscous_derivatives_edge_3D(Block *A, size_t time_level);
int viscous_derivatives_corners_3D(Block *A, size_t time_level);
int copy_derivatives( struct cell_vertex *a, struct cell_vertex *b );
int copy_2_derivatives( struct cell_vertex *a, struct cell_vertex *b, 
			struct cell_vertex *c );
int copy_3_derivatives( struct cell_vertex *a, struct cell_vertex *b, 
			struct cell_vertex *c, struct cell_vertex *d );

#define VISC3D_HH
#endif

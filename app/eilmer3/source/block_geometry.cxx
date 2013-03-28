/// \file block_geometry.cxx
/// \ingroup eilmer3
/// \brief Geometric functions that need to operate across the block.
///
/// \version 23-Mar-2013 extracted from block.cxx
///

#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <unistd.h>
#include "cell.hh"
#include "kernel.hh"
#include "block.hh"
#include "bc.hh"

//-----------------------------------------------------------------------------

int Block::compute_primary_cell_geometric_data(size_t dimensions, size_t time_level)
// Compute cell and interface geometric properties.
{
    size_t i, j, k;
    FV_Cell *cell, *cell_1, *cell_2, *ghost_cell;
    Vector3 dummy;
    FV_Interface *iface;
    Vector3 *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7;
    Vector3 ds;

    if ( dimensions == 2 ) {
	calc_volumes_2D(time_level);
	calc_faces_2D(time_level);
	calc_ghost_cell_geom_2D(time_level);
	return SUCCESS;
    }

    // Cell properties of volume and position.
    // Estimates of cross-cell distances for use in high-order reconstruction.
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		cell = get_cell(i,j,k);
		p0 = &(get_vtx(i,j,k)->pos[time_level]);
		p1 = &(get_vtx(i+1,j,k)->pos[time_level]);
		p2 = &(get_vtx(i+1,j+1,k)->pos[time_level]);
		p3 = &(get_vtx(i,j+1,k)->pos[time_level]);
		p4 = &(get_vtx(i,j,k+1)->pos[time_level]);
		p5 = &(get_vtx(i+1,j,k+1)->pos[time_level]);
		p6 = &(get_vtx(i+1,j+1,k+1)->pos[time_level]);
		p7 = &(get_vtx(i,j+1,k+1)->pos[time_level]);
		hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				     cell->pos[time_level], cell->volume[time_level], cell->iLength,
				     cell->jLength, cell->kLength);
		cell->L_min = cell->iLength;
		if ( cell->jLength < cell->L_min ) cell->L_min = cell->jLength;
		if ( cell->kLength < cell->L_min ) cell->L_min = cell->kLength;
	    }
	}
    }

    // work on ifi face as a WEST face
    // t1 in the j-ordinate direction
    // t2 in the k-ordinate direction
    for ( i = imin; i <= imax + 1; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		iface = get_ifi(i,j,k);
		p0 = &(get_vtx(i,j,k)->pos[time_level]);
		p3 = &(get_vtx(i,j+1,k)->pos[time_level]);
		p7 = &(get_vtx(i,j+1,k+1)->pos[time_level]);
		p4 = &(get_vtx(i,j,k+1)->pos[time_level]);
		quad_properties( *p0, *p3, *p7, *p4,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area[time_level] );
	    }
	}
    }

    // work on ifj face as a SOUTH face
    // t1 in the k-ordinate direction
    // t2 in the i-ordinate direction
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax + 1; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		iface = get_ifj(i,j,k);
		p0 = &(get_vtx(i,j,k)->pos[time_level]);
		p4 = &(get_vtx(i,j,k+1)->pos[time_level]);
		p5 = &(get_vtx(i+1,j,k+1)->pos[time_level]);
		p1 = &(get_vtx(i+1,j,k)->pos[time_level]);
		quad_properties( *p0, *p4, *p5, *p1,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area[time_level] );
	    }
	}
    }

    // work on ifk face as a BOTTOM face
    // t1 in the i-ordinate direction
    // t2 in the j-ordinate direction
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax + 1; ++k ) {
		iface = get_ifk(i,j,k);
		p0 = &(get_vtx(i,j,k)->pos[time_level]);
		p1 = &(get_vtx(i+1,j,k)->pos[time_level]);
		p2 = &(get_vtx(i+1,j+1,k)->pos[time_level]);
		p3 = &(get_vtx(i,j+1,k)->pos[time_level]);
		quad_properties( *p0, *p1, *p2, *p3,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area[time_level] );
	    }
	}
    }

    /* Propagate cross-cell lengths into the ghost cells. */
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    i = imin;
	    cell = get_cell(i,j,k);
	    get_cell(i-1,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	    get_cell(i-2,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	    i = imax;
	    cell = get_cell(i,j,k);
	    get_cell(i+1,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	    get_cell(i+2,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    j = jmin;
	    cell = get_cell(i,j,k);
	    get_cell(i,j-1,k)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	    get_cell(i,j-2,k)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	    j = jmax;
	    cell = get_cell(i,j,k);
	    get_cell(i,j+1,k)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	    get_cell(i,j+2,k)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    k = kmin;
	    cell = get_cell(i,j,k);
	    get_cell(i,j,k-1)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	    get_cell(i,j,k-2)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	    k = kmax;
	    cell = get_cell(i,j,k);
	    get_cell(i,j,k+1)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	    get_cell(i,j,k+2)->copy_values_from(*cell, COPY_CELL_LENGTHS, time_level);
	}
    }

    /* Extrapolate (with first-order) cell positions and volumes to ghost cells. */
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    i = imin;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i+1,j,k);
	    ghost_cell = get_cell(i-1,j,k);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i-2,j,k);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_2->volume[time_level] - cell_2->volume[time_level];
	    i = imax;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i-1,j,k);
	    ghost_cell = get_cell(i+1,j,k);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i+2,j,k);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    j = jmin;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j+1,k);
	    ghost_cell = get_cell(i,j-1,k);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j-2,k);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	    j = jmax;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j-1,k);
	    ghost_cell = get_cell(i,j+1,k);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j+2,k);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    k = kmin;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j,k+1);
	    ghost_cell = get_cell(i,j,k-1);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j,k-2);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	    k = kmax;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j,k-1);
	    ghost_cell = get_cell(i,j,k+1);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j,k+2);
	    ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	    ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	}
    }
    
    return SUCCESS;
} // end compute_primary_cell_geometric_data()

int Block::compute_distance_to_nearest_wall_for_all_cells(size_t dimensions, size_t time_level)
{
    FV_Cell *cell_at_wall[6];
    double dist[6], half_width[6];
    FV_Interface *face_at_wall;

    for ( FV_Cell *cellp : active_cells ) {
	std::vector<size_t> ijk = to_ijk_indices(cellp->id);
	size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2];
	// Step 1: get distances to all boundaries along all index directions.
	// If the block is not too distorted, these directions should take us
	// straight to the bounding walls.
	// North
	face_at_wall = get_ifj(i,jmax+1,k);
	dist[NORTH] = vabs(cellp->pos[time_level] - face_at_wall->pos);
	cell_at_wall[NORTH] = get_cell(i,jmax,k);
	half_width[NORTH] = vabs(cell_at_wall[NORTH]->pos[time_level] - face_at_wall->pos);
	// East
	face_at_wall = get_ifi(imax+1,j,k);
	dist[EAST] = vabs(cellp->pos[time_level] - face_at_wall->pos);
	cell_at_wall[EAST] = get_cell(imax,j,k);
	half_width[EAST] = vabs(cell_at_wall[EAST]->pos[time_level] - face_at_wall->pos);
	// South
	face_at_wall = get_ifj(i,jmin,k);
	dist[SOUTH] = vabs(cellp->pos[time_level] - face_at_wall->pos);
	cell_at_wall[SOUTH] = get_cell(i,jmin,k);
	half_width[SOUTH] = vabs(cell_at_wall[SOUTH]->pos[time_level] - face_at_wall->pos);
	// West
	face_at_wall = get_ifi(imin,j,k);
	dist[WEST] = vabs(cellp->pos[time_level] - face_at_wall->pos);
	cell_at_wall[WEST] = get_cell(imin,j,k);
	half_width[WEST] = vabs(cell_at_wall[WEST]->pos[time_level] - face_at_wall->pos);
	if ( dimensions == 3 ) {
	    // Top
	    face_at_wall = get_ifk(i,j,kmax+1);
	    dist[TOP] = vabs(cellp->pos[time_level] - face_at_wall->pos);
	    cell_at_wall[TOP] = get_cell(i,j,kmax);
	    half_width[TOP] = vabs(cell_at_wall[TOP]->pos[time_level] - face_at_wall->pos);
	    // Bottom
	    face_at_wall = get_ifk(i,j,kmin);
	    dist[BOTTOM] = vabs(cellp->pos[time_level] - face_at_wall->pos);
	    cell_at_wall[BOTTOM] = get_cell(i,j,kmin);
	    half_width[BOTTOM] = vabs(cell_at_wall[BOTTOM]->pos[time_level] - face_at_wall->pos);
	} else {
	    cell_at_wall[TOP] = NULL;
	    dist[TOP] = 0.0;
	    half_width[TOP] = 0.0;
	    cell_at_wall[BOTTOM] = NULL;
	    dist[BOTTOM] = 0.0;
	    half_width[BOTTOM] = 0.0;
	}

	// Step 2: Just in case there are no real walls for this block...
	//
	// We'll start by storing the largest distance and 
	// corresponding wall-cell half-width so that we have 
	// a relatively large distance in case there are no walls
	// on the boundary of the block.
	cellp->distance_to_nearest_wall = dist[0];
	cellp->half_cell_width_at_wall = half_width[0];
	for ( size_t iface = 1; iface < 6; ++iface ) {
	    if ( cell_at_wall[iface] != NULL &&
		 dist[iface] > cellp->distance_to_nearest_wall ) {
		cellp->distance_to_nearest_wall = dist[iface];
		cellp->half_cell_width_at_wall = half_width[iface];
		cellp->cell_at_nearest_wall = cell_at_wall[iface];
	    }
	}
		
	// Step 3: find the closest real wall.
	for ( size_t iface = 0; iface < 6; ++iface ) {
	    if ( dimensions == 2 && iface >= 4 ) break; // Top,Bottom are 4,5
	    if ( bcp[iface]->is_wall() && 
		 bcp[iface]->type_code != SLIP_WALL &&
		 dist[iface] < cellp->distance_to_nearest_wall && 
		 cell_at_wall[iface] != NULL ) {
		cellp->distance_to_nearest_wall = dist[iface];
		cellp->half_cell_width_at_wall = half_width[iface];
		cellp->cell_at_nearest_wall = cell_at_wall[iface];
	    }
	}
    } // end for *cellp
    
    return SUCCESS;
} // end compute_distance_to_nearest_wall_for_all_cells()

int Block::compute_secondary_cell_geometric_data(size_t dimensions, size_t time_level)
// Compute secondary-cell and interface geometric properties.
// Will be used for computing gradients for viscous terms.
//
// To do: The code for the 3D cells has been ported from eilmer2 without
//        taking advantage of the eilmer3 structure where ghost-cells have
//        been filled in with useful geometric information.
//        We should make use of this information.
{
    size_t i, j, k;
    Vector3 dummy;
    double iLen, jLen, kLen;
    FV_Vertex *vertex;
    FV_Interface *iface;
    Vector3 *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7;
    Vector3 ds;

    if ( dimensions == 2 ) {
	secondary_areas_2D(time_level);
	return SUCCESS;
    }

    /*
     * Internal secondary cell geometry information
     */
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    for ( k = kmin; k <= kmax-1; ++k ) {
		vertex = get_vtx(i+1,j+1,k+1);
		p0 = &(get_cell(i,j,k)->pos[time_level]);
		p1 = &(get_cell(i+1,j,k)->pos[time_level]);
		p2 = &(get_cell(i+1,j+1,k)->pos[time_level]);
		p3 = &(get_cell(i,j+1,k)->pos[time_level]);
		p4 = &(get_cell(i,j,k+1)->pos[time_level]);
		p5 = &(get_cell(i+1,j,k+1)->pos[time_level]);
		p6 = &(get_cell(i+1,j+1,k+1)->pos[time_level]);
		p7 = &(get_cell(i,j+1,k+1)->pos[time_level]);
		hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				     dummy, vertex->volume, iLen, jLen, kLen );
	    }
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    for ( k = kmin; k <= kmax-1; ++k ) {
		iface = get_sifi(i,j,k);
		p0 = &(get_cell(i,j,k)->pos[time_level]);
		p3 = &(get_cell(i,j+1,k)->pos[time_level]);
		p7 = &(get_cell(i,j+1,k+1)->pos[time_level]);
		p4 = &(get_cell(i,j,k+1)->pos[time_level]);
		quad_properties( *p0, *p3, *p7, *p4,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area[time_level] );
	    }
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax-1; ++k ) {
		iface = get_sifj(i,j,k);
		p0 = &(get_cell(i,j,k)->pos[time_level]);
		p4 = &(get_cell(i,j,k+1)->pos[time_level]);
		p5 = &(get_cell(i+1,j,k+1)->pos[time_level]);
		p1 = &(get_cell(i+1,j,k)->pos[time_level]);
		quad_properties( *p0, *p4, *p5, *p1,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area[time_level] );
	    }
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		iface = get_sifk(i,j,k);
		p0 = &(get_cell(i,j,k)->pos[time_level]);
		p1 = &(get_cell(i+1,j,k)->pos[time_level]);
		p2 = &(get_cell(i+1,j+1,k)->pos[time_level]);
		p3 = &(get_cell(i,j+1,k)->pos[time_level]);
		quad_properties( *p0, *p1, *p2, *p3,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area[time_level] );
	    }
	}
    }

    /*
     * East boundary secondary cell geometry information
     */
    i = imax;
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_cell(i,j,k)->pos[time_level]);
	    p1 = &(get_ifi(i+1,j,k)->pos);
	    p2 = &(get_ifi(i+1,j+1,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos[time_level]);
	    p4 = &(get_cell(i,j,k+1)->pos[time_level]);
	    p5 = &(get_ifi(i+1,j,k+1)->pos);
	    p6 = &(get_ifi(i+1,j+1,k+1)->pos);
	    p7 = &(get_cell(i,j+1,k+1)->pos[time_level]);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifi(i+1,j,k);
	    p1 = &(get_ifi(i+1,j,k)->pos);
	    p2 = &(get_ifi(i+1,j+1,k)->pos);
	    p6 = &(get_ifi(i+1,j+1,k+1)->pos);
	    p5 = &(get_ifi(i+1,j,k+1)->pos);
	    quad_properties( *p1, *p2, *p6, *p5,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos[time_level]);
	    p4 = &(get_cell(i,j,k+1)->pos[time_level]);
	    p5 = &(get_ifi(i+1,j,k+1)->pos);
	    p1 = &(get_ifi(i+1,j,k)->pos);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos[time_level]);
	    p1 = &(get_ifi(i+1,j,k)->pos);
	    p2 = &(get_ifi(i+1,j+1,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos[time_level]);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }

    /*
     * West boundary secondary cell geometry information
     */
    i = imin - 1;
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_ifi(i+1,j,k)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos[time_level]);
	    p2 = &(get_cell(i+1,j+1,k)->pos[time_level]);
	    p3 = &(get_ifi(i+1,j+1,k)->pos);
	    p4 = &(get_ifi(i+1,j,k+1)->pos);
	    p5 = &(get_cell(i+1,j,k+1)->pos[time_level]);
	    p6 = &(get_cell(i+1,j+1,k+1)->pos[time_level]);
	    p7 = &(get_ifi(i+1,j+1,k+1)->pos);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_ifi(i+1,j,k)->pos);
	    p3 = &(get_ifi(i+1,j+1,k)->pos);
	    p7 = &(get_ifi(i+1,j+1,k+1)->pos);
	    p4 = &(get_ifi(i+1,j,k+1)->pos);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_ifi(i+1,j,k)->pos);
	    p4 = &(get_ifi(i+1,j,k+1)->pos);
	    p5 = &(get_cell(i+1,j,k+1)->pos[time_level]);
	    p1 = &(get_cell(i+1,j,k)->pos[time_level]);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_ifi(i+1,j,k)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos[time_level]);
	    p2 = &(get_cell(i+1,j+1,k)->pos[time_level]);
	    p3 = &(get_ifi(i+1,j+1,k)->pos);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }

    /*
     * North boundary secondary cell geometry information
     */
    j = jmax;
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_cell(i,j,k)->pos[time_level]);
	    p1 = &(get_cell(i+1,j,k)->pos[time_level]);
	    p2 = &(get_ifj(i+1,j+1,k)->pos);
	    p3 = &(get_ifj(i,j+1,k)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos[time_level]);
	    p5 = &(get_cell(i+1,j,k+1)->pos[time_level]);
	    p6 = &(get_ifj(i+1,j+1,k+1)->pos);
	    p7 = &(get_ifj(i,j+1,k+1)->pos);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for (k = kmin; k <= kmax-1; ++k) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos[time_level]);
	    p3 = &(get_ifj(i,j+1,k)->pos);
	    p7 = &(get_ifj(i,j+1,k+1)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos[time_level]);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for (k = kmin; k <= kmax-1; ++k) {
	    iface = get_sifj(i,j+1,k);
	    p3 = &(get_ifj(i,j+1,k)->pos);
	    p7 = &(get_ifj(i,j+1,k+1)->pos);
	    p6 = &(get_ifj(i+1,j+1,k+1)->pos);
	    p2 = &(get_ifj(i+1,j+1,k)->pos);
	    quad_properties( *p3, *p7, *p6, *p2,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos[time_level]);
	    p1 = &(get_cell(i+1,j,k)->pos[time_level]);
	    p2 = &(get_ifj(i+1,j+1,k)->pos);
	    p3 = &(get_ifj(i,j+1,k)->pos);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }

    /*
     * South boundary secondary cell geometry information
     */
    j = jmin - 1;
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_ifj(i,j+1,k)->pos);
	    p1 = &(get_ifj(i+1,j+1,k)->pos);
	    p2 = &(get_cell(i+1,j+1,k)->pos[time_level]);
	    p3 = &(get_cell(i,j+1,k)->pos[time_level]);
	    p4 = &(get_ifj(i,j+1,k+1)->pos);
	    p5 = &(get_ifj(i+1,j+1,k+1)->pos);
	    p6 = &(get_cell(i+1,j+1,k+1)->pos[time_level]);
	    p7 = &(get_cell(i,j+1,k+1)->pos[time_level]);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_ifj(i,j+1,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos[time_level]);
	    p7 = &(get_cell(i,j+1,k+1)->pos[time_level]);
	    p4 = &(get_ifj(i,j+1,k+1)->pos);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_ifj(i,j+1,k)->pos);
	    p4 = &(get_ifj(i,j+1,k+1)->pos);
	    p5 = &(get_ifj(i+1,j+1,k+1)->pos);
	    p1 = &(get_ifj(i+1,j+1,k)->pos);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_ifj(i,j+1,k)->pos);
	    p1 = &(get_ifj(i+1,j+1,k)->pos);
	    p2 = &(get_cell(i+1,j+1,k)->pos[time_level]);
	    p3 = &(get_cell(i,j+1,k)->pos[time_level]);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }

    /*
     * Top boundary secondary cell geometry information
     */
    k = kmax;
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_cell(i,j,k)->pos[time_level]);
	    p1 = &(get_cell(i+1,j,k)->pos[time_level]);
	    p2 = &(get_cell(i+1,j+1,k)->pos[time_level]);
	    p3 = &(get_cell(i,j+1,k)->pos[time_level]);
	    p4 = &(get_ifk(i,j,k+1)->pos);
	    p5 = &(get_ifk(i+1,j,k+1)->pos);
	    p6 = &(get_ifk(i+1,j+1,k+1)->pos);
	    p7 = &(get_ifk(i,j+1,k+1)->pos);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos[time_level]);
	    p3 = &(get_cell(i,j+1,k)->pos[time_level]);
	    p7 = &(get_ifk(i,j+1,k+1)->pos);
	    p4 = &(get_ifk(i,j,k+1)->pos);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos[time_level]);
	    p4 = &(get_ifk(i,j,k+1)->pos);
	    p5 = &(get_ifk(i+1,j,k+1)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos[time_level]);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    iface = get_sifk(i,j,k+1);
	    p4 = &(get_ifk(i,j,k+1)->pos);
	    p5 = &(get_ifk(i+1,j,k+1)->pos);
	    p6 = &(get_ifk(i+1,j+1,k+1)->pos);
	    p7 = &(get_ifk(i,j+1,k+1)->pos);
	    quad_properties( *p4, *p5, *p6, *p7,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }

    /*
     * Bottom boundary secondary cell geometry information
     */
    k = kmin - 1;
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_ifk(i,j,k+1)->pos);
	    p1 = &(get_ifk(i+1,j,k+1)->pos);
	    p2 = &(get_ifk(i+1,j+1,k+1)->pos);
	    p3 = &(get_ifk(i,j+1,k+1)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos[time_level]);
	    p5 = &(get_cell(i+1,j,k+1)->pos[time_level]);
	    p6 = &(get_cell(i+1,j+1,k+1)->pos[time_level]);
	    p7 = &(get_cell(i,j+1,k+1)->pos[time_level]);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( i = imin; i <= imax; ++i) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_ifk(i,j,k+1)->pos);
	    p3 = &(get_ifk(i,j+1,k+1)->pos);
	    p7 = &(get_cell(i,j+1,k+1)->pos[time_level]);
	    p4 = &(get_cell(i,j,k+1)->pos[time_level]);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_ifk(i,j,k+1)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos[time_level]);
	    p5 = &(get_cell(i+1,j,k+1)->pos[time_level]);
	    p1 = &(get_ifk(i+1,j,k+1)->pos);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_ifk(i,j,k+1)->pos);
	    p1 = &(get_ifk(i+1,j,k+1)->pos);
	    p2 = &(get_ifk(i+1,j+1,k+1)->pos);
	    p3 = &(get_ifk(i,j+1,k+1)->pos);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area[time_level] );
	}
    }

    return SUCCESS;
} // end compute_secondary_cell_geometric_data_3D()

int Block::calc_volumes_2D(size_t time_level)
/// \brief Compute the PRIMARY cell volumes, areas, and centers 
///        from the vertex positions.
///
/// For 2D planar, assume unit length in the Z-direction.
/// For axisymmetry, compute a volume per radian.
///
/// Determine minimum length and aspect ratio, also.
{
    size_t i, j;
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double xN, yN, xS, yS, xE, yE, xW, yW;
    double vol, max_vol, min_vol, xyarea;
    double dx, dy, dxN, dyN, dxE, dyE;
    double lengthN, lengthE, length_max, length_min, length_cross;
    double max_aspect, aspect_ratio;
    FV_Cell *cell, *source_cell, *target_cell;

    // Cell layout
    // C-----B     3-----2
    // |     |     |     |
    // |  c  |     |  c  |
    // |     |     |     |
    // D-----A     0-----1
    
    max_vol = 0.0;
    min_vol = 1.0e6;    /* arbitrarily large */
    max_aspect = 0.0;
    for ( i = imin; i <= imax; ++i ) {
        for ( j = jmin; j <= jmax; ++j ) {
            cell = get_cell(i,j);
	    // These are the corners.
	    xA = cell->vtx[1]->pos[time_level].x;
	    yA = cell->vtx[1]->pos[time_level].y;
	    xB = cell->vtx[2]->pos[time_level].x;
	    yB = cell->vtx[2]->pos[time_level].y;
	    xC = cell->vtx[3]->pos[time_level].x;
	    yC = cell->vtx[3]->pos[time_level].y;
	    xD = cell->vtx[0]->pos[time_level].x;
	    yD = cell->vtx[0]->pos[time_level].y;
	    // Cell area in the (x,y)-plane.
            xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                            (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
	    // Cell Centroid.
            cell->pos[time_level].x = 1.0 / (xyarea * 6.0) * 
                      ((yB - yA) * (xA * xA + xA * xB + xB * xB) + 
                       (yC - yB) * (xB * xB + xB * xC + xC * xC) +
                       (yD - yC) * (xC * xC + xC * xD + xD * xD) + 
                       (yA - yD) * (xD * xD + xD * xA + xA * xA));
            cell->pos[time_level].y = -1.0 / (xyarea * 6.0) * 
                     ((xB - xA) * (yA * yA + yA * yB + yB * yB) + 
                      (xC - xB) * (yB * yB + yB * yC + yC * yC) +
                      (xD - xC) * (yC * yC + yC * yD + yD * yD) + 
                      (xA - xD) * (yD * yD + yD * yA + yA * yA));
	    cell->pos[time_level].z = 0.0;
	    // Cell Volume.
            if (get_axisymmetric_flag() == 1) {
                // Volume per radian = centroid y-ordinate * cell area
                vol = xyarea * cell->pos[time_level].y;
            } else {
                // Assume unit depth in the z-direction.
                vol = xyarea;
            }
            if (vol < 0.0) {
                printf("Negative cell volume: vol[%d][%d] = %e\n", 
		       static_cast<int>(i), static_cast<int>(j), vol);
                return 1;
            }
            if (vol > max_vol) max_vol = vol;
            if (vol < min_vol) min_vol = vol;
            cell->volume[time_level] = vol;
            cell->area[time_level] = xyarea;
	    // Check cell length scale using North and East boundaries.
	    // Also, save the minimum length for later use in the CFL
	    // checking routine.  Record max aspect ratio over all cells.
            dxN = xC - xB;
            dyN = yC - yB;
            dxE = xA - xB;
            dyE = yA - yB;
            lengthN = sqrt(dxN * dxN + dyN * dyN);
            lengthE = sqrt(dxE * dxE + dyE * dyE);

            length_max = lengthN;
            if (lengthE > length_max) length_max = lengthE;
            length_cross = xyarea / length_max; 
            // estimate of minimum width of cell
            length_min = lengthN;
            if (lengthE < length_min) length_min = lengthE;
            if (length_cross < length_min) length_min = length_cross;
            cell->L_min = length_min;
            aspect_ratio = length_max / length_min;
            if (aspect_ratio > max_aspect) max_aspect = aspect_ratio;

	    // Record the cell widths in the i- and j-index directions.
	    // The widths are measured between corresponding midpoints of
            // the bounding interfaces.
            // This data is later used by the high-order reconstruction.
            xN = 0.5 * (xC + xB);
            yN = 0.5 * (yC + yB);
            xS = 0.5 * (xD + xA);
            yS = 0.5 * (yD + yA);
            xE = 0.5 * (xA + xB);
            yE = 0.5 * (yA + yB);
            xW = 0.5 * (xD + xC);
            yW = 0.5 * (yD + yC);
            dx = xN - xS;
            dy = yN - yS;
            cell->jLength = sqrt(dx * dx + dy * dy);
            dx = xE - xW;
            dy = yE - yW;
            cell->iLength = sqrt(dx * dx + dy * dy);
	    cell->kLength = 0.0;
        } // j loop
    } // i loop

    // We now need to mirror the cell iLength and jLength
    // around the boundaries.
    // Those boundaries that are adjacent to another block
    // will be updated later with the other-block's cell lengths.
    for ( i = imin; i <= imax; ++i ) {
        // North boundary
        j = jmax;
        source_cell = get_cell(i,j);
	target_cell = get_cell(i,j+1);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
	target_cell = get_cell(i,j+2);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        // South boundary
        j = jmin;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i,j-1);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i,j-2);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
    } // end for i

    for ( j = jmin; j <= jmax; ++j ) {
        // East boundary
        i = imax;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i+1,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i+2,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        // West boundary
        i = imin;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i-1,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i-2,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
    } // end for j
	
    return SUCCESS;
} // end calc_volumes_2D()


/// \brief Compute the secondary cell cell areas in the (x,y)-plane.
///
/// The secondary cells are centred on the vertices of the 
/// primary cells and have primary cell centres as their corners.
/// For this particular secondary cell, centred on a vertex v (i,j),
/// the near-by primary cell i,j is centred on B'.
///
///          +-----+
///          |     |
///       C'-+--B' |
///       |  |  |  |
///       |  v--+--+
///       |     |
///       D'----A'
///
int Block::secondary_areas_2D(size_t time_level)
{
    size_t i, j;
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double xyarea, max_area, min_area;

    max_area = 0.0;
    min_area = 1.0e6;   // arbitrarily large
    // First, do all of the internal secondary cells.
    // i.e. The ones centred on primary vertices which 
    // are not on a boundary.
    for (i = imin+1; i <= imax; ++i) {
        for (j = jmin+1; j <= jmax; ++j) {
	    // These are the corners.
            xA = get_cell(i,j-1)->pos[time_level].x;
            yA = get_cell(i,j-1)->pos[time_level].y;
            xB = get_cell(i,j)->pos[time_level].x;
            yB = get_cell(i,j)->pos[time_level].y;
            xC = get_cell(i-1,j)->pos[time_level].x;
            yC = get_cell(i-1,j)->pos[time_level].y;
            xD = get_cell(i-1,j-1)->pos[time_level].x;
            yD = get_cell(i-1,j-1)->pos[time_level].y;
	    // Cell area in the (x,y)-plane.
            xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                            (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
            if (xyarea < 0.0) {
                printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		       static_cast<int>(id), static_cast<int>(i), static_cast<int>(j),
		       xyarea);
                return BAD_CELLS_ERROR;
            }
            if (xyarea > max_area) max_area = xyarea;
            if (xyarea < min_area) min_area = xyarea;
            get_vtx(i,j)->area = xyarea;
        } // j loop
    } // i loop

    // Note that the secondary cells along block boundaries are HALF cells.
    //
    // East boundary.
    i = imax+1;
    for (j = jmin+1; j <= jmax; ++j) {
	xA = get_ifi(i,j-1)->pos.x;
	yA = get_ifi(i,j-1)->pos.y;
	xB = get_ifi(i,j)->pos.x;
	yB = get_ifi(i,j)->pos.y;
        xC = get_cell(i-1,j)->pos[time_level].x;
        yC = get_cell(i-1,j)->pos[time_level].y;
        xD = get_cell(i-1,j-1)->pos[time_level].x;
        yD = get_cell(i-1,j-1)->pos[time_level].y;
	// Cell area in the (x,y)-plane.
        xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                        (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
        if (xyarea < 0.0) {
            printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		   static_cast<int>(id), static_cast<int>(i), static_cast<int>(j),
		   xyarea);
            return BAD_CELLS_ERROR;
        }
        if (xyarea > max_area) max_area = xyarea;
        if (xyarea < min_area) min_area = xyarea;
        get_vtx(i,j)->area = xyarea;
    } // j loop 

    // Fudge corners -- not expecting to use this data.
    get_vtx(i,jmin)->area = 0.5 * get_vtx(i,jmin+1)->area;
    get_vtx(i,jmax+1)->area = 0.5 * get_vtx(i,jmax)->area;
    
    // West boundary.
    i = imin;
    for (j = jmin+1; j <= jmax; ++j) {
        xA = get_cell(i,j-1)->pos[time_level].x;
        yA = get_cell(i,j-1)->pos[time_level].y;
        xB = get_cell(i,j)->pos[time_level].x;
        yB = get_cell(i,j)->pos[time_level].y;
	xC = get_ifi(i,j)->pos.x;
	yC = get_ifi(i,j)->pos.y;
	xD = get_ifi(i,j-1)->pos.x;
	yD = get_ifi(i,j-1)->pos.y;
	// Cell area in the (x,y)-plane.
        xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                        (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
        if (xyarea < 0.0) {
	    printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		   static_cast<int>(id), static_cast<int>(i), static_cast<int>(j),
		   xyarea);
            return BAD_CELLS_ERROR;
        }
        if (xyarea > max_area) max_area = xyarea;
        if (xyarea < min_area) min_area = xyarea;
        get_vtx(i,j)->area = xyarea;
    } // j loop 

    // Fudge corners.
    get_vtx(i,jmin)->area = 0.5 * get_vtx(i,jmin+1)->area;
    get_vtx(i,jmax+1)->area = 0.5 * get_vtx(i,jmax)->area;

    // North boundary.
    j = jmax+1;
    for (i = imin+1; i <= imax; ++i) {
	// These are the corners.
        xA = get_cell(i,j-1)->pos[time_level].x;
        yA = get_cell(i,j-1)->pos[time_level].y;
	xB = get_ifj(i,j)->pos.x;
	yB = get_ifj(i,j)->pos.y;
	xC = get_ifj(i-1,j)->pos.x;
	yC = get_ifj(i-1,j)->pos.y;
        xD = get_cell(i-1,j-1)->pos[time_level].x;
        yD = get_cell(i-1,j-1)->pos[time_level].y;
	// Cell area in the (x,y)-plane.
        xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                        (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
        if (xyarea < 0.0) {
            printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		   static_cast<int>(id), static_cast<int>(i), static_cast<int>(j),
		   xyarea);
            return BAD_CELLS_ERROR;
        }
        if (xyarea > max_area) max_area = xyarea;
        if (xyarea < min_area) min_area = xyarea;
        get_vtx(i,j)->area = xyarea;
    } // i loop 

    // Fudge corners.
    get_vtx(imin,j)->area = 0.5 * get_vtx(imin+1,j)->area;
    get_vtx(imax+1,j)->area = 0.5 * get_vtx(imax,j)->area;

    // South boundary.
    j = jmin;
    for (i = imin+1; i <= imax; ++i) {
	xA = get_ifj(i,j)->pos.x;
	yA = get_ifj(i,j)->pos.y;
        xB = get_cell(i,j)->pos[time_level].x;
        yB = get_cell(i,j)->pos[time_level].y;
        xC = get_cell(i-1,j)->pos[time_level].x;
        yC = get_cell(i-1,j)->pos[time_level].y;
	xD = get_ifj(i-1,j)->pos.x;
	yD = get_ifj(i-1,j)->pos.y;
	// Cell area in the (x,y)-plane.
        xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                        (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
        if (xyarea < 0.0) {
            printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		   static_cast<int>(id), static_cast<int>(i), static_cast<int>(j),
		   xyarea);
            return BAD_CELLS_ERROR;
        }
        if (xyarea > max_area) max_area = xyarea;
        if (xyarea < min_area) min_area = xyarea;
        get_vtx(i,j)->area = xyarea;
    } // i loop

    // Fudge corners.
    get_vtx(imin,j)->area = 0.5 * get_vtx(imin+1,j)->area;
    get_vtx(imax+1,j)->area = 0.5 * get_vtx(imax,j)->area;

    //printf("Max area = %e, Min area = %e\n", max_area, min_area);

    return SUCCESS;
} // end secondary_areas_2D()

/// \brief Compute the interface lengths and direction cosines for interfaces.
///
int Block::calc_faces_2D(size_t time_level)
{
    FV_Interface *iface;
    size_t i, j;
    double xA, xB, yA, yB, xC, yC;
    double LAB, LBC;

    // East-facing interfaces.
    for (i = imin; i <= imax+1; ++i) {
        for (j = jmin; j <= jmax; ++j) {
            iface = get_ifi(i,j);
	    // These are the corners.
            xA = get_vtx(i,j)->pos[time_level].x; 
	    yA = get_vtx(i,j)->pos[time_level].y;
            xB = get_vtx(i,j+1)->pos[time_level].x; 
	    yB = get_vtx(i,j+1)->pos[time_level].y;
            LAB = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA));
            if (LAB < 1.0e-9) {
                printf("Zero length ifi[%d,%d]: %e\n",
		       static_cast<int>(i), static_cast<int>(j), LAB);
            }
            // Direction cosines for the unit normal.
            iface->n.x = (yB - yA) / LAB;
            iface->n.y = -(xB - xA) / LAB;
	    iface->n.z = 0.0;  // 2D plane
	    iface->t2 = Vector3(0.0, 0.0, 1.0);
	    iface->t1 = cross(iface->n, iface->t2);
            // Length in the XY-plane.
            iface->length = LAB;
            // Mid-point and area.
            iface->Ybar = 0.5 * (yA + yB);
            if (get_axisymmetric_flag() == 1) {
                // Interface area per radian.
                iface->area[time_level] = LAB * iface->Ybar;
            } else {
                // Assume unit depth in the Z-direction.
                iface->area[time_level] = LAB;
            }
	    iface->pos = (get_vtx(i,j)->pos[time_level] + get_vtx(i,j+1)->pos[time_level])/2.0;
	    
        } // j loop
    } // i loop
    
    // North-facing interfaces.
    for (i = imin; i <= imax; ++i) {
        for (j = jmin; j <= jmax+1; ++j) {
            iface = get_ifj(i,j);
	    // These are the corners.
            xB = get_vtx(i+1,j)->pos[time_level].x;
            yB = get_vtx(i+1,j)->pos[time_level].y;
            xC = get_vtx(i,j)->pos[time_level].x;
            yC = get_vtx(i,j)->pos[time_level].y;
            LBC = sqrt((xC - xB) * (xC - xB) + (yC - yB) * (yC - yB));
            if (LBC < 1.0e-9) {
                printf("Zero length ifj[%d,%d]: %e\n",
		       static_cast<int>(i), static_cast<int>(j), LBC);
            }
            // Direction cosines for the unit normal.
            iface->n.x = (yC - yB) / LBC;
            iface->n.y = -(xC - xB) / LBC;
	    iface->n.z = 0.0;  // 2D plane
	    iface->t2 = Vector3(0.0, 0.0, 1.0);
	    iface->t1 = cross(iface->n, iface->t2);
            // Length in the XY-plane.
            iface->length = LBC;
            // Mid-point and area.
            iface->Ybar = 0.5 * (yC + yB);
            if (get_axisymmetric_flag() == 1) {
                // Interface area per radian.
                iface->area[time_level] = LBC * iface->Ybar;
            } else {
                // Assume unit depth in the Z-direction.
                iface->area[time_level] = LBC;
            }
	    iface->pos = (get_vtx(i+1,j)->pos[time_level] + get_vtx(i,j)->pos[time_level])/2.0;
        } // j loop
    } // i loop
    return SUCCESS;
} // end calc_faces_2D()

int Block::calc_ghost_cell_geom_2D(size_t time_level)
/// \brief Compute the ghost cell positions and volumes.
///
/// 'Compute' is a bit too strong to describe what we do here.
///  Rather this is a first-order extrapolation
/// from interior cells to estimate the position
/// and volume of the ghost cells.
{
    size_t i, j;
    FV_Cell *cell_1, *cell_2, *ghost_cell;
    // East boundary
    i = imax;
    for ( j = jmin; j <= jmax; ++j ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i-1,j);
	ghost_cell = get_cell(i+1,j);
	ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i+2,j);
	ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
    }
    // West boundary
    i = imin;
    for ( j = jmin; j <= jmax; ++j ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i+1,j);
	ghost_cell = get_cell(i-1,j);
	ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i-2,j);
	ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
    }
    // North boundary
    j = jmax;
    for ( i = imin; i <= imax; ++i ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i,j-1);
	ghost_cell = get_cell(i,j+1);
	ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i,j+2);
	ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
    }
    // South boundary
    j = jmin;
    for ( i = imin; i <= imax; ++i ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i,j+1);
	ghost_cell = get_cell(i,j-1);
	ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i,j-2);
	ghost_cell->pos[time_level] = 2.0*cell_1->pos[time_level] - cell_2->pos[time_level];
	ghost_cell->volume[time_level] = 2.0*cell_1->volume[time_level] - cell_2->volume[time_level];
    }
    return SUCCESS;
} // end Block::calc_ghost_cell_geom_2D()

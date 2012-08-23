/** \file exch_mpi.cxx
 * \ingroup eilmer3
 * \brief Functions to exchange boundary data via MPI.
 *
 * \author PA Jacobs and RJ Goozee
 * \version 05-Mar-08 Eilmer3 port
 */

// Intel MPI requires mpi.h included BEFORE stdio.h
#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
extern "C" {
#   include "../../../lib/util/source/useful.h"
}
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "exch2d.hh"
#include "exch3d.hh"
#include "exch_mpi.hh"

double *send_buffer[6];
double *receive_buffer[6];


/** \brief Returns the number of double values that will be sent
 *   in the buffer for MPI communication.
 */
int number_of_double_values( int face, int ni, int nj, int nk, int nv )
{
    int ne;
    if ( face == NORTH || face == SOUTH ) {
	ne = 2 * ni * nk * nv;
    } else if ( face == EAST || face == WEST ) {
	ne = 2 * nj * nk * nv;
    } else if ( face == TOP || face == BOTTOM ) {
	ne = 2 * ni * nj * nv;
    } else {
	printf("Wrong boundary index! = %d\n", face );
	ne = 0;
    }
    return ne;
}


/** \brief Allocate memory for the MPI send and receive buffers.
 *
 * Allocate space to two times the linear dimension of each boundary.
 * This provides space for ghost cells two deep along each edge.
 * \returns 0 if successful, 1 otherwise.
 */
int allocate_send_and_receive_buffers( Block *bd ) 
{
    int flag = 0;
    int total_bytes = 0;
    if (bd == NULL) {
        printf("allocate_send_and_receive_buffers: NULL pointer\n");
        return VALUE_ERROR;
    }
    int nv = number_of_values_in_cell_copy(COPY_ALL_CELL_DATA);
    for ( int face = 0; face < 6; ++face ) {
	int ne = number_of_double_values(face, bd->nni, bd->nnj, bd->nnk, nv);
	send_buffer[face] = (double *) calloc(ne, sizeof(double));
	if (send_buffer[face] == NULL)
	    ++flag;
	else
	    total_bytes += ne * sizeof(double);
	receive_buffer[face] = (double *) calloc(ne, sizeof(double));
	if (receive_buffer[face] == NULL)
	    ++flag;
	else
	    total_bytes += ne * sizeof(double);
    }
    if ( flag > 0 ) {
        printf("Block %d: buffer allocation failed.\n", bd->id);
	return MEMORY_ERROR;
    }
    return SUCCESS;
} // end allocate_send_and_receive_buffers()


/** \brief A tag should be unique for each block and boundary combination. 
 */
int make_tag( int block_id, int face ) 
{
    return 10 * block_id + face;
}


/** \brief Ensure that all boundary data is exchanged between
 *         connected boundaries on adjacent blocks.
 *
 */
int mpi_exchange_boundary_data( int jb, int type_of_copy )
{
    int other_block, other_face;
    MPI_Status status[6];
    MPI_Request request[6];
    int tag, ne, nfaces;

    global_data *G = get_global_data_ptr();
    Block *bdp = get_block_data_ptr(jb);
    int nv  = number_of_values_in_cell_copy( type_of_copy );
    if ( G->dimensions == 3 ) {
	nfaces = 6;
    } else {
	nfaces = 4;
    }

    // Post non-blocking receives.
    for ( int face = 0; face < nfaces; ++face ) {
	if ( bdp->bcp[face]->type_code == ADJACENT ||
	     bdp->bcp[face]->type_code == ADJACENT_PLUS_UDF ) {
	    other_block = bdp->bcp[face]->neighbour_block;
	    other_face = bdp->bcp[face]->neighbour_face;
	    ne = number_of_double_values(face, bdp->nni, bdp->nnj, bdp->nnk, nv);
	    tag = make_tag( other_block, other_face );
	    MPI_Irecv( receive_buffer[face], ne, MPI_DOUBLE, other_block, 
		       tag, MPI_COMM_WORLD, &(request[face]) );
#           if 0
	    printf("Post receive: block[%d] face[%d] from block[%d] face[%d]\n",
		   jb, face, other_block, other_face );
#           endif
	}
    } // end for face...

    // Blocking sends (to corresponding receives on other processes).
    for ( int face = 0; face < nfaces; ++face ) {
	if ( bdp->bcp[face]->type_code == ADJACENT ||
	     bdp->bcp[face]->type_code == ADJACENT_PLUS_UDF ) {
	    other_block = bdp->bcp[face]->neighbour_block;
	    other_face = bdp->bcp[face]->neighbour_face;
	    ne = number_of_double_values(face, bdp->nni, bdp->nnj, bdp->nnk, nv);
	    tag = make_tag( jb, face );
	    if ( G->dimensions == 2 ) {
		copy_into_send_buffer_2D(bdp, face, type_of_copy, send_buffer[face]);
	    } else {
		copy_into_send_buffer_3D(bdp, face, type_of_copy, send_buffer[face]);
	    }
#           if 0
	    printf("Send: block[%d] face[%d] to block[%d] face[%d]\n",
		   jb, face, other_block, other_face );
#           endif
	    MPI_Send( send_buffer[face], ne, MPI_DOUBLE, other_block, tag, MPI_COMM_WORLD );
	}
    } // end for face...

    // Wait for receives to complete.
    // Once they complete, copy the data back into the ghost cells.
    for ( int face = 0; face < nfaces; ++face ) {
	if ( bdp->bcp[face]->type_code == ADJACENT ||
	     bdp->bcp[face]->type_code == ADJACENT_PLUS_UDF ) {
	    other_block = bdp->bcp[face]->neighbour_block;
	    other_face = bdp->bcp[face]->neighbour_face;
	    MPI_Wait( &(request[face]), &(status[face]) );
	    if ( G->dimensions == 2 ) {
		copy_from_receive_buffer_2D(bdp, face, type_of_copy, receive_buffer[face]);
	    } else {
		copy_from_receive_buffer_3D(bdp, face, type_of_copy, receive_buffer[face]);
	    }
#           if 0
	    printf("Received OK: block[%d] face[%d] from block[%d] face[%d]\n",
		   jb, face, other_block, other_face );
#           endif
	}
    } // end for ( face...

    return SUCCESS;
} /* end mpi_exchange_boundary_data() */


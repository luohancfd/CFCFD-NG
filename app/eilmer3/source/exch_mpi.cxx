/// \file exch_mpi.cxx
/// \ingroup eilmer3
/// \brief Functions to exchange boundary data via MPI.
///
/// \author PA Jacobs and RJ Goozee
/// \version 05-Mar-08 Eilmer3 port
/// \version 26-Aug-2012 Multiple blocks per MPI process.

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

std::vector <double *> send_buffer;
std::vector <double *> receive_buffer;
std::vector<MPI_Status> status;
std::vector<MPI_Request> request;


/// \brief Returns the number of double values that will be sent
///  in the buffer for MPI communication.
int number_of_double_values(int face, int ni, int nj, int nk, int nv)
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


/// \brief Allocate memory for the MPI send and receive buffers.
///
/// Allocate space to two times the linear dimension of each boundary.
/// This provides space for ghost cells two deep along each edge.
/// \returns 0 if successful, 1 otherwise.
int allocate_send_and_receive_buffers(void) 
{
    global_data &G = *get_global_data_ptr();
    Block *bdp;
    int flag = 0;
    int total_bytes = 0;
    int n_local_blocks = G.my_blocks.size(); 
    status.resize(n_local_blocks*6);
    request.resize(n_local_blocks*6);
    send_buffer.resize(n_local_blocks*6);
    receive_buffer.resize(n_local_blocks*6);

    if ( get_verbose_flag() ) printf("Begin allocate_send_and_receive_buffers()");
    int nv = number_of_values_in_cell_copy(COPY_ALL_CELL_DATA);
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	for ( int face = 0; face < 6; ++face ) {
	    int ne = number_of_double_values(face, bdp->nni, bdp->nnj, bdp->nnk, nv);
	    send_buffer[jb*6+face] = (double *) calloc(ne, sizeof(double));
	    if (send_buffer[jb*6+face] == NULL)
		++flag;
	    else
		total_bytes += ne * sizeof(double);
	    receive_buffer[jb*6+face] = (double *) calloc(ne, sizeof(double));
	    if (receive_buffer[jb*6+face] == NULL)
		++flag;
	    else
		total_bytes += ne * sizeof(double);
	}
    } // end for jb...
    if ( get_verbose_flag() ) printf("    allocated %d bytes for data buffers.", total_bytes);
    if ( flag > 0 ) {
        printf("Block %d: buffer allocation failed.\n", bdp->id);
	return MEMORY_ERROR;
    }
    return SUCCESS;
} // end allocate_send_and_receive_buffers()


/// \brief Clean up memory for the MPI send and receive buffers.
int delete_send_and_receive_buffers(void) 
{
    global_data &G = *get_global_data_ptr();
    if ( get_verbose_flag() ) printf("Begin delete_send_and_receive_buffers()");
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	for ( int face = 0; face < 6; ++face ) {
	    free(send_buffer[jb*6+face]);
	    send_buffer[jb*6+face] = 0;
	    free(receive_buffer[jb*6+face]);
	    receive_buffer[jb*6+face] = 0;
	}
    } // end for jb...
    status.clear();
    request.clear();
    send_buffer.clear();
    receive_buffer.clear();
    if ( get_verbose_flag() ) printf("    done deleting buffers.");
    return SUCCESS;
} // end delete_send_and_receive_buffers()


/// \brief A tag should be unique for each block and boundary combination. 
int make_tag(int block_id, int face) 
{
    return 10 * block_id + face;
}


/// \brief Ensure that all boundary data is exchanged between
///        connected boundaries on adjacent blocks.
int mpi_exchange_boundary_data(int type_of_copy, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    Block *bdp;
    int other_block, other_face;
    int tag, ne, nfaces;

    int nv  = number_of_values_in_cell_copy(type_of_copy);
    if ( G.dimensions == 3 ) {
	nfaces = 6;
    } else {
	nfaces = 4;
    }

    // Post non-blocking receives.
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	for ( int face = 0; face < nfaces; ++face ) {
	    if ( bdp->bcp[face]->type_code == ADJACENT ||
		 bdp->bcp[face]->type_code == ADJACENT_PLUS_UDF ) {
		other_block = bdp->bcp[face]->neighbour_block;
		other_face = bdp->bcp[face]->neighbour_face;
		ne = number_of_double_values(face, bdp->nni, bdp->nnj, bdp->nnk, nv);
		tag = make_tag(other_block, other_face);
		MPI_Irecv(receive_buffer[jb*6+face], ne, MPI_DOUBLE, 
			  G.mpi_rank_for_block[other_block], 
			  tag, MPI_COMM_WORLD, &(request[jb*6+face]));
#               if 0
		printf("Post receive: block[%d] face[%d] from block[%d] face[%d]\n",
		       bdp->id, face, other_block, other_face );
#               endif
	    }
	} // end for face...
    } // end for jb...

    // Blocking sends (to corresponding receives on other processes).
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	for ( int face = 0; face < nfaces; ++face ) {
	    if ( bdp->bcp[face]->type_code == ADJACENT ||
		 bdp->bcp[face]->type_code == ADJACENT_PLUS_UDF ) {
		other_block = bdp->bcp[face]->neighbour_block;
		other_face = bdp->bcp[face]->neighbour_face;
		ne = number_of_double_values(face, bdp->nni, bdp->nnj, bdp->nnk, nv);
		tag = make_tag(bdp->id, face);
		if ( G.dimensions == 2 ) {
		    copy_into_send_buffer_2D(bdp, face, type_of_copy, send_buffer[jb*6+face], gtl);
		} else {
		    copy_into_send_buffer_3D(bdp, face, type_of_copy, send_buffer[jb*6+face], gtl);
		}
#               if 0
		printf("Send: block[%d] face[%d] to block[%d] face[%d]\n",
		       bdp->id, face, other_block, other_face );
#               endif
		MPI_Send(send_buffer[jb*6+face], ne, MPI_DOUBLE, 
			 G.mpi_rank_for_block[other_block],
			 tag, MPI_COMM_WORLD);
	    }
	} // end for face...
    } // end for jb...

    // Wait for receives to complete.
    // Once they complete, copy the data back into the ghost cells.
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	for ( int face = 0; face < nfaces; ++face ) {
	    if ( bdp->bcp[face]->type_code == ADJACENT ||
		 bdp->bcp[face]->type_code == ADJACENT_PLUS_UDF ) {
		other_block = bdp->bcp[face]->neighbour_block;
		other_face = bdp->bcp[face]->neighbour_face;
		MPI_Wait(&(request[jb*6+face]), &(status[jb*6+face]));
		if ( G.dimensions == 2 ) {
		    copy_from_receive_buffer_2D(bdp, face, type_of_copy, receive_buffer[jb*6+face], gtl);
		} else {
		    copy_from_receive_buffer_3D(bdp, face, type_of_copy, receive_buffer[jb*6+face], gtl);
		}
#               if 0
		printf("Received OK: block[%d] face[%d] from block[%d] face[%d]\n",
		       bdp->id, face, other_block, other_face);
#               endif
	    }
	} // end for ( face...
    } // end for jb...

    return SUCCESS;
} // end mpi_exchange_boundary_data()


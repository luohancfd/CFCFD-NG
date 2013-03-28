/** \file exch_mpi.h
 * \ingroup eilmer3
 * \brief Header for MPI specific items.
 */

#ifndef EXCH_MPI_HH
#define EXCH_MPI_HH

#include "block.hh"

int number_of_double_values(int face, int ni, int nj, int nk, int nv);
int allocate_send_and_receive_buffers(void);
int delete_send_and_receive_buffers(void);
int make_tag(int block_id, int face);
int mpi_exchange_boundary_data(int type_of_copy, size_t time_level);

#endif


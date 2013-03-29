/** \file exch2d.hh
 * \ingroup eilmer3
 * \brief Prototypes for the boundary exchange functions. 
 */

#ifndef EXCH2D_HH
#define EXCH2D_HH

int exchange_shared_boundary_data(int jb, int type_of_copy, size_t gtl);

int copy_boundary_data_2D(int jb, int type_of_copy, size_t gtl);
int copy_to_north_boundary_2D(Block *A, Block *B, int B_bndry, int type_of_copy, size_t gtl);
int copy_to_east_boundary_2D(Block *A, Block *B, int B_bndry, int type_of_copy, size_t gtl);
int copy_to_south_boundary_2D(Block *A, Block *B, int B_bndry, int type_of_copy, size_t gtl);
int copy_to_west_boundary_2D(Block *A, Block *B, int B_bndry, int type_of_copy, size_t gtl);
int copy_to_east_boundary_diaphragm_2D(Block *A, Block *B, 
				       int B_bndry, int type_of_copy,
				       double diaphragm_time_fraction,
				       double diaphragm_rupture_diameter, 
				       double sim_time, size_t gtl);
int copy_to_west_boundary_diaphragm_2D(Block *A, Block *B, 
				       int B_bndry, int type_of_copy,
				       double diaphragm_time_fraction, 
				       double diaphragm_rupture_diameter, 
				       double sim_time, size_t gtl);

int copy_into_send_buffer_2D(Block *bd, int bndry, int type_of_copy,
			     double *send_buffer, size_t gtl);
int copy_from_receive_buffer_2D(Block *bd, int bndry, int type_of_copy,
				double *receive_buffer, size_t gtl);

#endif

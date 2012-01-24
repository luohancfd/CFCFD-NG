/** \file exch2d.hh
 * \ingroup eilmer3
 * \brief Prototypes for the boundary exchange functions. 
 */

#ifndef EXCH2D_HH
#define EXCH2D_HH

int exchange_shared_boundary_data(int jb, int type_of_copy);

int copy_boundary_data_2D(int jb, int type_of_copy);
int copy_to_north_boundary_2D(struct Block *A, struct Block *B,
			      int B_bndry, int type_of_copy);
int copy_to_east_boundary_2D(struct Block *A, struct Block *B,
			     int B_bndry, int type_of_copy);
int copy_to_south_boundary_2D(struct Block *A, struct Block *B,
			      int B_bndry, int type_of_copy);
int copy_to_west_boundary_2D(struct Block *A, struct Block *B,
			     int B_bndry, int type_of_copy);
int copy_to_east_boundary_diaphragm_2D(struct Block *A, 
				       struct Block *B, 
				       int B_bndry, int type_of_copy,
				       double diaphragm_time_fraction,
				       double diaphragm_rupture_diameter, 
				       double sim_time);
int copy_to_west_boundary_diaphragm_2D(struct Block *A, 
				       struct Block *B, 
				       int B_bndry, int type_of_copy,
				       double diaphragm_time_fraction, 
				       double diaphragm_rupture_diameter, 
				       double sim_time);

int copy_into_send_buffer_2D( Block *bd, int bndry, int type_of_copy, double *send_buffer );
int copy_from_receive_buffer_2D( Block *bd, int bndry, int type_of_copy, double *receive_buffer );

#endif

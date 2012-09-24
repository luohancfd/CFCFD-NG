// l_adapt.hh

#ifndef L_ADAPT_HH
#define L_ADAPT_HH

int L_adapt_cells(struct slug_data *A);
int L_fuse_cell_data(struct L_cell *source1,
                     struct L_cell *source2, struct L_cell *destination);
int L_split_cell_data_roughly( struct L_cell *source, double xL, 
			       struct L_cell *dest1, struct L_cell *dest2);
int L_split_cell_data_smoothly( struct L_cell *ci,   /* the cell to be split */
				struct L_cell *cim1, /* the cell to left     */
				struct L_cell *cip1, /* the cell to right    */
				struct L_cell *ca,   /* left new cell        */ 
				struct L_cell *cb ); /* right new cell       */

#endif

// l_adapt.hh

#ifndef L_ADAPT_HH
#define L_ADAPT_HH

#include "l_slug.hh"
#include "l1d.hh"
#include "l_cell.hh"

int L_adapt_cells(GasSlug *A);
int L_fuse_cell_data(LCell *source1, LCell *source2, LCell *destination);
int L_split_cell_data_roughly(LCell *source, double xL, LCell *dest1, LCell *dest2);
int L_split_cell_data_smoothly(LCell *ci,   /* the cell to be split */
			       LCell *cim1, /* the cell to left     */
			       LCell *cip1, /* the cell to right    */
			       LCell *ca,   /* left new cell        */ 
			       LCell *cb ); /* right new cell       */

#endif

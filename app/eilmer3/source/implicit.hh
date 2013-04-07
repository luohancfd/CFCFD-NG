/** \file implicit.hh
 * \ingroup eilmer3
 * \brief Prototypes for functions to compute point implicit viscous updates.
 *
 * \author OJ
 * \version Oct 2009
 */

#include <stdio.h>
#include <stdlib.h>
/*---------------------------------------------------------------------*/
int gasdynamic_point_implicit_inviscid_increment(double dt);
int gasdynamic_fully_implicit_inviscid_increment(double dt);
int inviscid_point_implicit_update_for_cell(FV_Cell *cell);
int calculate_M_inviscid(FV_Cell *cell, int dimensions);
int calculate_inviscid_jacobian(FV_Cell *cell, FV_Interface *iface);
int calculate_h_inviscid(FV_Cell *cell, int dimensions);
int gasdynamic_point_implicit_viscous_increment(void);
int gasdynamic_fully_implicit_viscous_increment(void);
int point_implicit_update_for_cell(FV_Cell *cell);
int calculate_M(FV_Cell *cell, int dimensions);
int calculate_viscous_jacobian(FV_Cell *cell, FV_Interface *iface);
int calculate_h(FV_Cell *cell, int dimensions);
void gaussj(FV_Cell *cell, int n, int m);

/*---------------------------------------------------------------------*/

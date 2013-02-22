/* td_br_filter.i 
 * SWIG interface header file for td_br_filter.c
 */

%module td_br_filter 
%{
#include "td_br_filter.h"
%}

extern int cfilter( char *vname, int nhalf );
extern int c_heat_transfer_emf( char *vname, double dt,
				double rhock, double E, double a );
extern int c_heat_transfer_temp( char *vname, double dt,
				 double rhock, double E, double a );
extern int c_emf_temp( char *vname, double dt,
		       double rhock, double E, double a);
extern int c_surface_temperature( char *vname, double dt,
				  double rhock, double E, double a );

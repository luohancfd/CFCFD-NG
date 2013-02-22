/** \file td_br_filter.h
 * \ingroup tds
 * \brief Fast filter functions for td_browser.tcl.
 */

int cfilter( Tcl_Interp *interp, char *vname, int nhalf );
int c_heat_transfer_emf( Tcl_Interp *interp, char *vname, double dt,
		     double rhock, double E, double a );
int c_heat_transfer_temp( Tcl_Interp *interp, char *vname, double dt,
		     double rhock, double E, double a );
int c_emf_temp( Tcl_Interp *interp, char *vname, double dt,
		     double rhock, double E, double a );
int c_surface_temperature( Tcl_Interp *interp, char *vname, double dt,
			   double rhock, double E, double a );

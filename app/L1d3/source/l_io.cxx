/** \file l_io.cxx
 * \ingroup l1d3
 * \brief I/O functions for l1d.c.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <sstream>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/kinetics/reaction-update.hh"
#include "l1d.hh"
#include "l_kernel.hh"
#include "l_diaph.hh"
#include "l_piston.hh"
#include "l_cell.hh"
#include "l_slug.hh"
using namespace std;


int print_simulation_status(FILE *strm, const char* efname, int step, SimulationData& SD,
			    vector<GasSlug>& A, vector<DiaphragmData>& Diaph,
			    vector<PistonData>& Pist, double cfl_max, 
			    double cfl_tiny, double time_tiny) {
    /*
     * Print the simulation status to the specified stream.
     * If the specified stream is NULL, then send the data to the 
     * event log file.
     */

    char msg_string[256];
    int close_stream, js, jd, jp;
    double p_max, x_max, E_tot;

    UNUSED_VARIABLE(cfl_tiny);
    UNUSED_VARIABLE(time_tiny);
    
    close_stream = 0;

    if ( strm == NULL ) {
        /* Presume that we want to open the events file. */
        strm = fopen( efname, "a+" );
        if ( strm != NULL ) { close_stream = 1; }
    }

    if ( strm != NULL ) {
	sprintf( msg_string,
		 "Step= %7d, time= %e, dt= %10.3e, CFL= %10.3e\n",
                 step, SD.sim_time, SD.dt_global, cfl_max);
	fputs( msg_string, strm );
	for (js = 0; js < SD.nslug; ++js) {
	    A[js].maximum_p(&p_max, &x_max);
	    E_tot = A[js].total_energy();
	    sprintf( msg_string,
		     "Slug %d: p_max=%10.3e, x=%10.3e, Etot=%10.3e, dt=%10.3e, nnx=%d\n",
		     js, p_max, x_max, E_tot, A[js].dt_allow, A[js].nnx);
	    fputs( msg_string, strm );
	}
	for (jd = 0; jd < SD.ndiaphragm; ++jd) {
	    sprintf( msg_string, "Diaph[%d].is_burst = %d, trigger_time = %e\n",
		     jd, Diaph[jd].is_burst, Diaph[jd].trigger_time);
	    fputs( msg_string, strm );
        }
	for (jp = 0; jp < SD.npiston; ++jp) {
	    sprintf( msg_string, "Piston %d: flags=%d%d%d",
		     jp, Pist[jp].is_restrain, Pist[jp].on_buffer,
		     Pist[jp].brakes_on);
	    fputs( msg_string, strm );
	    sprintf( msg_string, " x=%10.3e V=%10.3e a=%10.3e KE=%10.3e mass=%10.3e\n",
		     Pist[jp].x, Pist[jp].V, Pist[jp].DVDt[0],
		     0.5 * Pist[jp].mass * Pist[jp].V * Pist[jp].V, Pist[jp].mass);
	    fputs( msg_string, strm );
	    sprintf( msg_string, 
		     "          hit_count= %d, speed at last strike= %e\n",
		     Pist[jp].hit_buffer_count, Pist[jp].V_buffer);
	    fputs( msg_string, strm );
	} // end for
    }

    if ( close_stream ) {
        fclose( strm );
    }

    return SUCCESS;
}

int log_event(const char *efname, const char* event_message ) {
    /*
     * Write a message to the events log file.
     * This file is opened and closed each time so that
     * Windows users can see the data.
     */
    FILE *efp;
    efp = fopen(efname, "a+");
    if ( efp != NULL ) {
        fputs(event_message, efp);
        /* print_status( efp ); */
        fclose( efp );
    }
    return SUCCESS;
}


int L_write_cell_history(GasSlug& A, FILE* hisfile)
// Write out the flow solution in a (small) subset of cells
// at a different (often smaller) time interval to the full
// flow solution.
// As of 30-Sep-2012, we assume indexing 0..nn-1 in the slug.
{
    // The output format for this function needs to be kept the 
    // same as that for L_write_x_history().
    for ( int i = 0; i < A.hncell; ++i) {
        int ix = A.hxcell[i] + A.ixmin;
	fprintf(hisfile, "%s\n", A.Cell.at(ix).write_cell_values_to_string().c_str());
    }
    fflush(hisfile);
    return SUCCESS;
} // end function L_write_cell_history


int L_write_x_history(double xloc, std::vector<GasSlug>& A, FILE* hisfile)
// Write out the flow solution at a specified x-location
// at a different (often smaller) time interval to the full
// flow solution.
// To get values at the end of a slug, say at a reflecting
// wall, it may be necessary to specify an x-location which
// will always fall within the last cell and not at the very edge.
// Input...
// xloc    : x-location
// A       : reference to the vector of GasSlugs
// hisfile : file to which data is to be written
{
    // The output format for this function needs to be kept the 
    // same as that for L_write_cell_history().
    Gas_model *gmodel = get_gas_model_ptr();
    LCell icell = LCell(gmodel);
    // Find the gas slug containing the x-location
    int found=0;
    for ( size_t js = 0; js < A.size(); ++js ) {
        found = A[js].interpolate_cell_data(xloc, icell);
        if (found == 1) break;
    }
    fprintf(hisfile, "%s\n", icell.write_cell_values_to_string().c_str() );
    fflush(hisfile);
    return SUCCESS;
} // end function L_write_x_history

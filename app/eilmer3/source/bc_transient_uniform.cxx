// bc_transient_uniform.cxx

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <valarray>

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_transient_uniform.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

TransientUniformBC::TransientUniformBC( Block *bdp, int which_boundary,
					std::string filename )
    : BoundaryCondition(bdp, which_boundary, TRANSIENT_UNI, "TransientUniformBC",
			0, false, false, -1, -1, 0),
      filename(filename)
{
    // Reads the time-record flow state data from a previously written file. 
    //
    // The format expected of this file is space-separated values for 
    //
    //   t, p, u, v, w, T, tke, omega, f0, f1 ... 
    //
    // for each cell with an extra line at the top to specify the variable names.  
    // If there is only one species, it is sufficient to supply values for p, u, v, T only.
    // If vibrational or electronic energies are requested, they are
    // assumed to be in equilibrium with the transrotational temp T.
    char line[512], token[512];
    FILE *fp;
    std::vector<double> massf_line;
    gmodel = get_gas_model_ptr();
    nsp = gmodel->get_number_of_species();
    nmodes = gmodel->get_number_of_modes();
#   if ECHO_ALL
    cout << "TransientUniformBC(): Read data file." << endl;
#   endif
    fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
	cerr << "TransientUniformBC(): cannot open file " << filename << endl;
        exit(FILE_ERROR);
    }
    if ( fgets(line, sizeof(line), fp) == NULL ) {
	cerr << "TransientUniformBC(): failure of fgets()" << endl;
	cerr << "    Couldn't read first line of file." << endl;
	exit(FILE_ERROR);
    }
    nsample = 0;
    // For each line in the file, store the flow state data.
    while ( 1 ) {
	if ( fgets(line, sizeof(line), fp) == NULL ) {
	    break;
	}
	++nsample;
	// Pull the line apart with the string tokenizer.
	strcpy( token, strtok(line, " \t") ); tta.push_back(atof(token));
	strcpy( token, strtok(NULL, " \t") ); pa.push_back(atof(token));
	strcpy( token, strtok(NULL, " \t") ); ua.push_back(atof(token));
	strcpy( token, strtok(NULL, " \t") ); va.push_back(atof(token));
	strcpy( token, strtok(NULL, " \t") ); wa.push_back(atof(token));
	strcpy( token, strtok(NULL, " \t") ); Ta.push_back(atof(token));
	strcpy( token, strtok(NULL, " \t") ); tkea.push_back(atof(token));
	strcpy( token, strtok(NULL, " \t") ); omegaa.push_back(atof(token));
#       if ECHO_ALL
        printf( "Sample[%d]: time=%e p=%e u=%e v=%e w=%e T=%e\n", 
		nsample, tta[nsample-1], pa[nsample-1], ua[nsample-1],
		va[nsample-1], wa[nsample-1], Ta[nsample-1] );
	printf( "            tke=%e omega=%e\n", tkea[nsample-1], omegaa[nsample-1] );
        fflush(stdout);
#       endif

	massf_line.resize(nsp, 0.0);
        if ( nsp == 1 ) {
            massf_line[0] = 1.0;
        } else {
	    for ( int isp = 0; isp < nsp; ++isp ) {
		strcpy( token, strtok(NULL, " \t") );
		massf_line[isp] = atof(token);
	    }
        }
	massfa.push_back(massf_line);
#       if ECHO_ALL
	printf( "            " );
	for ( int isp = 0; isp < nsp; ++isp ) {
            printf( "massf[%d]=%e, ", isp, massfa[nsample-1][isp] );
	}
        printf( "\n" ); fflush(stdout);
#       endif
    } // end for

#   if ECHO_ALL
    cout << "TransientUniformBC(): finished reading " << nsample << " time samples." << endl;
#   endif
} // end TransientUniformBC constructor

TransientUniformBC::TransientUniformBC( const TransientUniformBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation), 
      filename(bc.filename)
{
    nsp = bc.nsp;
    nmodes = bc.nmodes;
    cerr << "TransientUniformBC copy constructor isn't complete!!!!!!!" << endl;
    exit( NOT_IMPLEMENTED_ERROR );
}

TransientUniformBC::TransientUniformBC()
    : BoundaryCondition(0, 0, TRANSIENT_UNI, "TransientUniformBC",
			0, false, false, -1, -1, 0),
      filename("")
{ /* Cannot do much without more information. */ }

TransientUniformBC & 
TransientUniformBC::operator=(const TransientUniformBC &bc)
{
    BoundaryCondition::operator=(bc);
    if ( this != &bc ) {
	filename = bc.filename;
	tta = bc.tta; pa = bc.pa; ua = bc.ua;
	va = bc.va; wa = bc.wa; Ta = bc.Ta;
	massfa = bc.massfa;
	gmodel = bc.gmodel;
	nsample = bc.nsample; nsp = bc.nsp; nmodes=bc.nmodes;
    }
    return *this;
}

TransientUniformBC::~TransientUniformBC() 
{}

int TransientUniformBC::apply_inviscid( double t )
{
    int ii;
    double p, u, v, w, mu_t, k_t, tke, omega;
    std::vector<double> T;
    std::vector<double> mf;
    double alpha;
    Block & bd = *bdp;

    // Check the ends of the time range.
    if ( t <= tta[0] || t >= tta[nsample - 1] ) {
        if ( t <= tta[0] ) ii = 0; else ii = nsample - 1;
        p = pa[ii];
        for ( int isp = 0; isp < nsp; ++isp ) mf.push_back(massfa[ii][isp]);
        for ( int imode = 0; imode < nmodes; ++imode ) T.push_back(Ta[ii]);
        u = ua[ii];
        v = va[ii];
        w = wa[ii];
	mu_t = 0.0;
	k_t = 0.0;
        tke = tkea[ii];
        omega = omegaa[ii];
    } else {
        // Do a linear search to find the bracketing time segment.
        ii = 1;
        while ( t > tta[ii] ) ++ii;
        alpha = (t - tta[ii-1]) / (tta[ii] - tta[ii-1]);
        if ( alpha > 1.0 ) alpha = 1.0;
        if ( alpha < 0.0 ) alpha = 0.0;
        p = (1.0 - alpha) * pa[ii-1] + alpha * pa[ii];
        for ( int isp = 0; isp < nsp; ++isp ) {
            mf.push_back( (1.0 - alpha) * massfa[ii-1][isp] + alpha * massfa[ii][isp] );
        }
        for ( int imode = 0; imode < nmodes; ++imode ) {
            T.push_back( (1.0 - alpha) * Ta[ii-1] + alpha * Ta[ii] );
        }
        u = (1.0 - alpha) * ua[ii-1] + alpha * ua[ii];
        v = (1.0 - alpha) * va[ii-1] + alpha * va[ii];
        w = (1.0 - alpha) * wa[ii-1] + alpha * wa[ii];
	mu_t = 0.0;
	k_t = 0.0;
        tke = (1.0 - alpha) * tkea[ii-1] + alpha * tkea[ii];
        omega = (1.0 - alpha) * omegaa[ii-1] + alpha * omegaa[ii];
    } 
    int S = 0;
    CFlowCondition *gsp = new CFlowCondition(gmodel, p, u, v, w, T, mf, "", 
					     tke, omega, mu_t, k_t, S);
    mf.clear();
    T.clear();
    // Now, fill in the ghost cells, assuming that the flow is
    // essentially like SupersonicIn.
    int i, j, k;
    FV_Cell *dest_cell;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for k
 	break;
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    } // end switch

    delete gsp;
    return SUCCESS;
}

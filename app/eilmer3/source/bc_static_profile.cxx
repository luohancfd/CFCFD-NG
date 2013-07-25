// bc_static_profile.cxx

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
#include "bc_static_profile.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

StaticProfileBC::StaticProfileBC( Block *bdp, int which_boundary,
				  const std::string filename, size_t n_profile )
    : BoundaryCondition(bdp, which_boundary, STATIC_PROF, "StaticProfileBC",
			0, false, false, -1, -1, 0),
      filename(filename), n_profile(n_profile)
{
    // Reads the flow state data from a previously written profile file.
    //
    // The default is to take in one profile slice and to copy this profile
    // into the two ghost cell slices. However, if the user wishes to use
    // individual profile slices for each slice of ghost cell, the user can
    // input two profile slices by setting the n_profile option to 2.
    //
    // The format expected of this file is that written by the Python
    // code found in e3_flow.py, as used by the postprocessing program e3_post.py.
    // Look for the functions variable_list_for_cell() and write_cell_data().
    // Also look at variable_list_for_cell() in cell.cxx.
    //
    // The first line in the file specifies the variable names on the 
    // remaining lines.  The elements on all lines are space separated.
    //
    // For the input of two profile slices, the first profile listed in the
    // file must be for the inner ghost cells, and the second profile for
    // the outer ghost cells.
    //
    char line[512], token[512];
    double x, y, z, volume, rho, u, v, w, p, a, mu, k,mu_t, k_t;
    int S;
    double Q_rad_org, f_rad_org, Q_rE_rad, tke, omega, dt_chem, dt_therm;
    std::vector<double> massf, e, T;
    size_t ncell, ncell_for_profile;
    size_t ncell_read_from_file = 0;
    FILE *fp;
    global_data *G = get_global_data_ptr();
    //---------------------------------------------------------------------------
    // FIX-ME: Update this BC for 3D, using C++ tokenizing stream and something 
    // like the usual block indexing for the storage order.
    if ( G->dimensions == 3 ) {
	cerr << "StaticProfileBC is not implemented for 3D." << endl;
	exit(NOT_IMPLEMENTED_ERROR);
    }
    //--------------------------------------------------------------------------
    Gas_model *gmodel = get_gas_model_ptr();
    nsp = gmodel->get_number_of_species();
    massf.resize(nsp);
    nmodes = gmodel->get_number_of_modes();
    e.resize(nmodes);
    T.resize(nmodes);
        
    fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        cerr << "StaticProfileBC() constructor:"
	     << " cannot open file " << filename;
	exit(FILE_ERROR);
    }
    // Read and ignore the comment line containing the variable names.
    if ( fgets(line, sizeof(line), fp) == NULL ) {
	cerr << "StaticProfileBC(): failure of fgets()" << endl;
	cerr << "    Had expected to read a line of variable names." << endl;
	cerr << "    Quitting program." << endl;
	exit(FILE_ERROR);
    }
    // For data each line in the file, store the flow state data.
    while ( 1 ) {
	if ( fgets(line, sizeof(line), fp) == NULL ) break;
	/* Pull the line apart with the string tokenizer. */
	strcpy(token, strtok(line, " ")); sscanf(token, "%lf", &x);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &y);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &z);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &volume);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &rho);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &u);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &v);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &w);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &p);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &a);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &mu);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &k); // labelled as k[0], only one present
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &mu_t);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &k_t);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%d", &S);
	if ( get_radiation_flag() == 1 ) {
	    strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &Q_rad_org);
	    strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &f_rad_org);
	    strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &Q_rE_rad);
	}
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &tke);
	strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &omega);
        if ( nsp == 1 ) {
	    strcpy( token, strtok(NULL, " ") ); massf[0] = 1.0;
	    // ignore the mass-fraction value in the file.
        } else {
	    for ( size_t isp = 0; isp < nsp; ++isp ) {
		strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &(massf[isp]));
	    }
        }
	if ( nsp > 1 ) {
	    strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &dt_chem );
	} else {
	    dt_chem = -1.0;
	}
	for ( size_t imode = 0; imode < nmodes; ++imode ) {
	    strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &(e[imode]));
	    strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &(T[imode]));
	}
	if ( nmodes > 1 ) {
	    strcpy(token, strtok(NULL, " ")); sscanf(token, "%lf", &dt_therm );
	} else {
	    dt_therm = -1.0;
	}
#       if 1
	cout << "x=" << x << " y=" << y << " z=" << z 
	     << " volume=" << volume << " rho=" << rho 
	     << " u=" << u << " v=" << v << " w=" << w 
	     << " p=" << p << " a=" << a << " mu=" << mu << " k=" << k 
	     << " mu_t=" << mu_t << " k_t=" << k_t << " S=" << S 
	     << " tke=" << tke << " omega=" << omega 
	     << " dt_chem=" << dt_chem << " e[0]=" << e[0] << " T[0]=" << T[0] << endl;
#       endif
	flow_profile.push_back(new CFlowCondition(gmodel, p, u, v, w, T, massf, "", tke, omega, mu_t, k_t, S));
	++ncell_read_from_file;
    } // end while
#   if ECHO_ALL
    cout << "StaticProfileBC() constructor: read " << ncell_read_from_file << " cells." << endl; 
#   endif
    // For the case with two input profiles, check that the number of cells read is an even number
    if ( ncell_read_from_file % 2 != 0 ) {
        cerr << "StaticProfileBC() constructor:" << endl
             << "    For 2 input profiles, the number of cells read should be an even number: " << endl
             << ", ncell_read_from_file=" << ncell_read_from_file << endl;
        exit(BAD_INPUT_ERROR);
    }
    // Check for that the number of cells is appropriate for this boundary
    if ( which_boundary == NORTH || which_boundary == SOUTH ) {
	ncell = bdp->nni;
    } else {
	ncell = bdp->nnj;
    }
    ncell_for_profile = ncell_read_from_file / n_profile;
    if ( ncell != ncell_for_profile ) {
        cerr << "StaticProfileBC() constructor:" << endl 
	     << "    Inconsistent numbers of cells: ncell=" << ncell
	     << ", ncell_for_profile=" << ncell_for_profile << endl;
        exit(BAD_INPUT_ERROR);
    }
    massf.clear();
    e.clear();
    T.clear();
#   if ECHO_ALL
    cout << "StaticProfileBC() constructor: done." << endl;
#   endif
}

StaticProfileBC::StaticProfileBC( const StaticProfileBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face),
      filename(bc.filename), n_profile(bc.n_profile)
{
    cerr << "StaticProfileBC() copy constructor is not implemented." << endl;
    exit( NOT_IMPLEMENTED_ERROR );
}

StaticProfileBC::StaticProfileBC()
    : BoundaryCondition(0, 0, STATIC_PROF, "StaticProfileBC",
			0, false, false, -1, -1, 0),
      filename(""), n_profile(0)
{}

StaticProfileBC & StaticProfileBC::operator=(const StaticProfileBC &bc)
{
    if ( this != &bc ) {
	BoundaryCondition::operator=(bc);
	filename = bc.filename;
	n_profile = bc.n_profile;
	nsp = bc.nsp;
	ncell_for_profile = bc.ncell_for_profile;
	for ( size_t i = 0; i < bc.flow_profile.size(); ++i ) {
	    flow_profile.push_back(new CFlowCondition(*(bc.flow_profile[i])));
	}
    }
    return *this;
}

StaticProfileBC::~StaticProfileBC() 
{
    for ( size_t i = 0; i < flow_profile.size(); ++i ) {
	delete flow_profile[i];
	flow_profile[i] = 0;
    }
    flow_profile.clear();
}

int StaticProfileBC::apply_convective( double t )
{
    size_t i, ifirst, ilast, j, jfirst, jlast, ncell_for_profile;
    FV_Cell *dest_cell;
    CFlowCondition *gsp;
    Block & bd = *bdp;

    ncell_for_profile = flow_profile.size() / n_profile;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	ifirst = bd.imin;
	ilast = bd.imax;
        for (i = ifirst; i <= ilast; ++i) {
            gsp = flow_profile[i - ifirst];
            dest_cell = bd.get_cell(i,j+1);
            dest_cell->copy_values_from(*gsp);
            if (n_profile == 2) {
                gsp = flow_profile[(i - ifirst) + ncell_for_profile];
            }
            dest_cell = bd.get_cell(i,j+2);
            dest_cell->copy_values_from(*gsp);
        } // end i loop
	break;
    case EAST:
	i = bd.imax;
	jfirst = bd.jmin;
	jlast = bd.jmax;
        for (j = jfirst; j <= jlast; ++j) {
            gsp = flow_profile[j - jfirst];
            dest_cell = bd.get_cell(i+1,j);
            dest_cell->copy_values_from(*gsp);
            if (n_profile == 2) {
                gsp = flow_profile[(j - jfirst) + ncell_for_profile];
            }
            dest_cell = bd.get_cell(i+2,j);
            dest_cell->copy_values_from(*gsp);
	} // end j loop
	break;
    case SOUTH:
	j = bd.jmin;
	ifirst = bd.imin;
	ilast = bd.imax;
        for (i = ifirst; i <= ilast; ++i) {
            gsp = flow_profile[i - ifirst];
            dest_cell = bd.get_cell(i,j-1);
            dest_cell->copy_values_from(*gsp);
            if (n_profile == 2) {
                gsp = flow_profile[(i - ifirst) + ncell_for_profile];
            }
            dest_cell = bd.get_cell(i,j-2);
            dest_cell->copy_values_from(*gsp);
        } // end i loop
	break;
    case WEST:
	i = bd.imin;
	jfirst = bd.jmin;
	jlast = bd.jmax;
        for (j = jfirst; j <= jlast; ++j) {
            gsp = flow_profile[j - jfirst];
            dest_cell = bd.get_cell(i-1,j);
            dest_cell->copy_values_from(*gsp);
            if (n_profile == 2) {
                gsp = flow_profile[(j - jfirst) + ncell_for_profile];
            }
            dest_cell = bd.get_cell(i-2,j);
            dest_cell->copy_values_from(*gsp);
        } // end j loop
 	break;
    // TODO: the TOP and BOTTOM boundaries.
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    } // end switch

    return SUCCESS;
}

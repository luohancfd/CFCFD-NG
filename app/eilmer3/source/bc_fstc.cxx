// bc_fstc.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "bc_fstc.hh"
#include "bc_catalytic.hh"
#include "bc_menter_correction.hh"

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

fstcBC::fstcBC( Block *bdp, int which_boundary, const std::string filename )
    : BoundaryCondition(bdp, which_boundary, FSTC, "fstcBC",
			0, true, false, -1, -1, 0),
      filename(filename)
{
    // Reads the temperature profile from the solid solver output.

    char line[512], token[512];
    double T;
    //std::vector<double> fstcT;
    FILE *fp;
    size_t ncell, nread;

    fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        cerr << "fstcBC() constructor:"
	     << " cannot open file " << filename << endl;
	exit(FILE_ERROR);
    }

    if ( fgets(line, sizeof(line), fp) == NULL ) {
	cerr << "fstcBC(): failure of fgets()" << endl;
	cerr << "Quitting program." << endl;
	exit(FILE_ERROR);
    }

    nread = sscanf(line, "%zu", &ncell_for_profile);
    if ( nread != 1 ) {
        cerr << "fstcBC() constructor:"
	     << "Could not read ncell_for_profile from line:" << endl
	     << line;
        exit(BAD_INPUT_ERROR);
    }

    // Check for that the number of cells is appropriate for this boundary
    if ( which_boundary == NORTH || which_boundary == SOUTH ) {
	ncell = bdp->nni;
    } else {
	ncell = bdp->nnj;
    }

    if ( ncell != ncell_for_profile ) {
        cerr << "fstcBC() constructor:" << endl
	     << "    Inconsistent numbers of cells: ncell=" << ncell
	     << ", ncell_for_profile=" << ncell_for_profile << endl;
        exit(BAD_INPUT_ERROR);
    }

    /* For each line in the file, store the temperature data. */
    for ( size_t ii = 0; ii < ncell_for_profile; ++ii ) {
	if ( fgets(line, sizeof(line), fp) == NULL ) {
	    cerr << "fstcBC(): failure of fgets()" << endl;
	    cerr << "Quitting program." << endl;
	    exit(FILE_ERROR);
	}
    /* Pull the line apart with the string tokenizer. */
    strcpy( token, strtok(line, " ") );
    sscanf( token, "%lf", &T );
    fstc_TProfile.push_back(T);
    } // end for
} //end of constructor

fstcBC::fstcBC( const fstcBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation),
      filename(bc.filename)
{}

fstcBC::fstcBC()
    : BoundaryCondition(0, 0, FSTC, "fstcBC",
			0, true, false, -1, -1, 0),
      filename("")
{ /* Cannot do anything useful here. */ }

fstcBC & fstcBC::operator=(const fstcBC &bc)
{
    BoundaryCondition::operator=(bc);
    filename = bc.filename;
    cerr << "fstcBC() assignment operator with file: " << filename
	 << "Not implemented. " << endl;
    exit(NOT_IMPLEMENTED_ERROR);
    return *this;
}

fstcBC::~fstcBC() {}

int fstcBC::apply_viscous( double t )
{
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    Block & bd = *bdp;
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) {
            fs.gas->T[imode] = fstc_TProfile[i-bd.imin];
		}
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[NORTH]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) {
            fs.gas->T[imode] = fstc_TProfile[j-bd.jmin];
		}
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[EAST]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) {
            fs.gas->T[imode] = fstc_TProfile[i-bd.imin];
		}
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[SOUTH]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) {
            fs.gas->T[imode] = fstc_TProfile[j-bd.jmin];
		}
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[WEST]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[TOP];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) {
            fs.gas->T[imode] = fstc_TProfile[j-bd.jmin];
		}
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[TOP]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[BOTTOM];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) {
            fs.gas->T[imode] = fstc_TProfile[j-bd.jmin];
		}
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[BOTTOM]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_viscous not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
} // end fstcBC::apply_viscous()

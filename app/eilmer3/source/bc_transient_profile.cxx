// bc_transient_profile.cxx
// PJ 2014-nov-17 adapted from bc_static_profile as a starting point.

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>

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
#include "bc_transient_profile.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

TransientProfileBC::TransientProfileBC(Block *bdp, int which_boundary,
				 const std::string _filename, size_t _n_profile)
    : BoundaryCondition(bdp, which_boundary, TRANSIENT_PROF),
      filename(_filename), n_profile(_n_profile)
{
    // Reads the flow state data from a previously written transient profile file.
    //
    // The default is to take in a double-cell profile slice and to copy this profile
    // into the two ghost cells at each corresponding location on the inflow face. 
    //
    // The format expected of this file is that written by the Block::write_profile
    // code found in block_io.cxx.
    // The first line in the file specifies the variable names for the data
    // that appears on the remaining lines.  The data on all lines are space separated.
    // For the input of two profile slices, two cells (ghost cell 1 followed by ghost cell 2)
    // of data are read at each point across the surface.  
    //
    double x, y, z, volume, rho, u, v, w, p, a, mu, mu_t, k_t, Sfloat;
    int S;
    double Q_rad_org, f_rad_org, Q_rE_rad, tke, omega, dt_chem, dt_therm;
    std::vector<double> massf, e, T, k;
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    nsp = gmodel->get_number_of_species();
    massf.resize(nsp);
    nmodes = gmodel->get_number_of_modes();
    e.resize(nmodes);
    T.resize(nmodes);
    k.resize(nmodes);

    if ( G.verbosity_level >= 2 ) {
	cout << "TransientProfileBC() constructor: filename= " << filename << endl;
    }
    std::ifstream fstrm(filename.c_str());
    std::string text;
    getline(fstrm, text); // Read and ignore the comment line saying where the data has come from.
    getline(fstrm, text); // Read and ignore the comment line containing the variable names.
    getline(fstrm, text); // Read and ignore the comment line containing the numbers of cells across the profile.
    // [TODO] would be good to cross-check these numbers.
    switch ( which_boundary ) {
    case NORTH:
    case SOUTH:
	ncell_for_profile = bdp->nni * bdp->nnk;
	break;
    case EAST:
    case WEST:
	ncell_for_profile = bdp->nnj * bdp->nnk;
	break;
    case TOP:
    case BOTTOM:
	ncell_for_profile = bdp->nni * bdp->nnj;
	break;
    default:
	throw std::runtime_error("Should not have arrived here -- invalid boundary.");
    }

    // The remaining lines in the file are in sets for particular time instances.
    // [TODO] This first cut of the code will read only the first time instance.
    // We need to decide the mechanics of reading instances as needed and how to interpolate them.

    getline(fstrm, text); // Read the line containing the time stamp.
    double time_stamp;
    sscanf(text.c_str(), "time= %lg", &time_stamp);
    if ( G.verbosity_level >= 2 ) {
	cout << "TransientProfileBC: reading profile for time_stamp= " << time_stamp << endl;
    }
    // For data each line in the file, store the flow state data for later use in the ghost cells.
    for ( size_t ncell = 0; ncell < ncell_for_profile; ++ncell ) {
	for ( size_t i_profile = 0; i_profile < n_profile; ++i_profile ) {
	    // Suck in a line of data and save it as a flow condition.
	    fstrm >> x >> y >> z >> volume >> rho >> u >> v >> w >> p >> a >> mu;
	    for ( size_t imode = 0; imode < nmodes; ++imode ) { fstrm >> k[imode]; }
	    fstrm >> mu_t >> k_t;
	    fstrm >> Sfloat; S = int(Sfloat); // beware of reading int from float string!!
	    if ( G.radiation ) { fstrm >> Q_rad_org >> f_rad_org >> Q_rE_rad; }
	    fstrm >> tke >> omega;
	    if ( nsp == 1 ) {
		fstrm >> massf[0]; massf[0] = 1.0; // ignore the mass-fraction value in the file.
	    } else {
		for ( size_t isp = 0; isp < nsp; ++isp ) { fstrm >> massf[isp]; }
	    }
	    if ( nsp > 1 ) { fstrm >> dt_chem; } else { dt_chem = -1.0; }
	    for ( size_t imode = 0; imode < nmodes; ++imode ) {
		fstrm >> e[imode] >> T[imode];
	    }
	    if ( nmodes > 1 ) { fstrm >> dt_therm; } else { dt_therm = -1.0; }
	    if ( G.verbosity_level >= 3 ) {
		cout << "x=" << x << " y=" << y << " z=" << z 
		     << " volume=" << volume << " rho=" << rho 
		     << " u=" << u << " v=" << v << " w=" << w 
		     << " p=" << p << " a=" << a << " mu=" << mu << " k[0]=" << k[0] 
		     << " mu_t=" << mu_t << " k_t=" << k_t << " S=" << S 
		     << " tke=" << tke << " omega=" << omega 
		     << " dt_chem=" << dt_chem << " e[0]=" << e[0] << " T[0]=" << T[0] << endl;
	    }
	    flow_profile.push_back(new CFlowCondition(gmodel,p,u,v,w,T,massf,"",tke,omega,mu_t,k_t,S));
	} // end for i_profile
    } // end for ncell

    massf.clear(); e.clear(); T.clear();
    if ( G.verbosity_level >= 2 ) {
	cout << "TransientProfileBC() constructor is done." << endl;
    }
} // end TransientProfileBC constructor

TransientProfileBC::TransientProfileBC(const TransientProfileBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      filename(bc.filename), n_profile(bc.n_profile)
{
    cerr << "TransientProfileBC() copy constructor is not implemented." << endl;
    exit( NOT_IMPLEMENTED_ERROR );
}

TransientProfileBC::TransientProfileBC()
    : BoundaryCondition(0, 0, STATIC_PROF),
      filename(""), n_profile(0)
{}

TransientProfileBC & TransientProfileBC::operator=(const TransientProfileBC &bc)
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

TransientProfileBC::~TransientProfileBC() 
{
    for ( size_t i = 0; i < flow_profile.size(); ++i ) {
	delete flow_profile[i];
	flow_profile[i] = 0;
    }
    flow_profile.clear();
}

void TransientProfileBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "filename= " << filename << endl;
    cout << lead_in << "n_profile= " << n_profile << endl;
    return;
}

int TransientProfileBC::apply_convective( double t )
{
    size_t i, j, k, ncell;
    FV_Cell *dest_cell;
    CFlowCondition *gsp;
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		ncell = ((k - bd.kmin)*bd.nni + (i - bd.imin)) * n_profile;
		gsp = flow_profile[ncell];
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*gsp);
		if (n_profile == 2) gsp = flow_profile[ncell+1];
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // end k loop
	break;
    case EAST:
	i = bd.imax;
        for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	    for ( j = bd.jmin; j <= bd.jmax; ++j ) {
		ncell = ((k - bd.kmin)*bd.nnj + (j - bd.jmin)) * n_profile;
		gsp = flow_profile[ncell];
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*gsp);
		if (n_profile == 2) gsp = flow_profile[ncell+1];
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // end k loop
	break;
    case SOUTH:
	j = bd.jmin;
        for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		ncell = ((k - bd.kmin)*bd.nni + (i - bd.imin)) * n_profile;
		gsp = flow_profile[ncell];
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*gsp);
		if (n_profile == 2) gsp = flow_profile[ncell+1];
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // end k loop
	break;
    case WEST:
	i = bd.imin;
        for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	    for ( j = bd.jmin; j <= bd.jmax; ++j ) {
		ncell = ((k - bd.kmin)*bd.nnj + (j - bd.jmin)) * n_profile;
		gsp = flow_profile[ncell];
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*gsp);
		if (n_profile == 2) gsp = flow_profile[ncell+1];
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // end k loop
	break;
    case TOP:
	k = bd.kmax;
	for ( j = bd.jmin; j <= bd.jmax; ++j ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		ncell = ((j - bd.jmin)*bd.nni + (i - bd.imin)) * n_profile;
		gsp = flow_profile[ncell];
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*gsp);
		if (n_profile == 2) gsp = flow_profile[ncell+1];
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // end j loop
	break;
    case BOTTOM:
	k = bd.kmin;
	for ( j = bd.jmin; j <= bd.jmax; ++j ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		ncell = ((j - bd.jmin)*bd.nni + (i - bd.imin)) * n_profile;
		gsp = flow_profile[ncell];
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*gsp);
		if (n_profile == 2) gsp = flow_profile[ncell+1];
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // end j loop
    } // end switch which_boundary

    return SUCCESS;
} // end TransientProfileBC::apply_convective

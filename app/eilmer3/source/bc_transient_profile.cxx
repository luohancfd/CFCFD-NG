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
    global_data &G = *get_global_data_ptr();
    if ( G.verbosity_level >= 2 ) {
	cout << "TransientProfileBC() constructor: filename= " << filename << endl;
    }
    fstrm.open(filename.c_str());
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
    // We will read in two timeslices and store them as t0 and t1 with t0 < t1.

    fstrm >> text;
    fstrm >> t0; 

    if ( G.verbosity_level >= 2 ) {
	cout << "TransientProfileBC: reading profile for time_stamp= " << t0 << endl;
    }
    read_profile_into_t0();
    
    fstrm >> text;
    fstrm >> t1; 

    if ( G.verbosity_level >= 2 ) {
	cout << "TransientProfileBC: reading profile for time_stamp= " << t1 << endl;
    }
    read_profile_into_t1();

    // Just copy anything to get started.
    cf_interp = new CFlowCondition(*(flow_profile_t0[0]));

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
	for ( size_t i = 0; i < bc.flow_profile_t0.size(); ++i ) {
	    flow_profile_t0.push_back(new CFlowCondition(*(bc.flow_profile_t0[i])));
	    flow_profile_t1.push_back(new CFlowCondition(*(bc.flow_profile_t1[i])));
	}
    }
    return *this;
}

TransientProfileBC::~TransientProfileBC() 
{
    for ( size_t i = 0; i < flow_profile_t0.size(); ++i ) {
	delete flow_profile_t0[i];
	flow_profile_t0[i] = 0;
	delete flow_profile_t1[i];
	flow_profile_t1[i] = 0;
    }
    flow_profile_t0.clear();
    flow_profile_t1.clear();
    delete cf_interp;
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
    global_data &G = *get_global_data_ptr();
    size_t i, j, k, ncell;
    FV_Cell *dest_cell;
    CFlowCondition *cf0, *cf1;
    Block & bd = *bdp;

    // We'll need where the simulation time 't' is in relation to the
    // time slices in t0 and t1.
    double w; // this is weighting to apply to linearly blend t0 and t1 profiles.
    if ( t < t0 ) {
	// If the simulation time is earlier than the values we have for
	// t0, then we'll just assume that the profile is constant in
	// time up until t0.
	w = 0.0; // full weight to profile at t0
    }
    
    if ( t >= t0 && t < t1 ) {
	w = (t - t0)/(t1 - t0);
    }

    if ( t >= t1 ) {
	if ( fstrm.eof() ) {
	    // We've run out of profiles to read. So we just give full weighting to
	    // profile at t1.
	    w = 1.0;
	}
	else {
	    // We'll try to keep reading the file until we find a t1 > t.
	    while ( !fstrm.eof() ) {
		fstrm >> text;
		fstrm >> t1; 
		if ( G.verbosity_level >= 2 ) {
		    cout << "TransientProfileBC: reading profile for time_stamp= " << t1 << endl;
		}
		if ( fstrm.eof() ) {
		    // Get out of this loop now if we did read past the end of the file.
		    w = 1.0;
		    break;
		}
		update_profile_t0_with_t1();
		read_profile_into_t1();
		if ( t1 > t ) {
		    // We found a new profile to use.
		    w = (t - t0)/(t1 - t0);
		    break;
		}
	    }
	    // if we reached this far, we got to the end of the file.
	    // so just give full weighting to profile at t1.
	    w = 1.0;
	}
    }

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		ncell = ((k - bd.kmin)*bd.nni + (i - bd.imin)) * n_profile;
		cf0 = flow_profile_t0[ncell];
		cf1 = flow_profile_t1[ncell];
		interpolate_condition(*cf0, *cf1, w, *cf_interp);
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*cf_interp);
		if (n_profile == 2) {
		    cf0 = flow_profile_t0[ncell+1];
		    cf1 = flow_profile_t1[ncell+1];
		    interpolate_condition(*cf0, *cf1, w, *cf_interp);
		}
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*cf_interp);
	    } // end i loop
	} // end k loop
	break;
    case EAST:
	i = bd.imax;
        for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	    for ( j = bd.jmin; j <= bd.jmax; ++j ) {
		ncell = ((k - bd.kmin)*bd.nnj + (j - bd.jmin)) * n_profile;
		cf0 = flow_profile_t0[ncell];
		cf1 = flow_profile_t1[ncell];
		interpolate_condition(*cf0, *cf1, w, *cf_interp);
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*cf_interp);
		if (n_profile == 2) {
		    cf0 = flow_profile_t0[ncell+1];
		    cf1 = flow_profile_t1[ncell+1];
		    interpolate_condition(*cf0, *cf1, w, *cf_interp);
		}
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*cf_interp);
	    } // end j loop
	} // end k loop
	break;
    case SOUTH:
	j = bd.jmin;
        for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		ncell = ((k - bd.kmin)*bd.nni + (i - bd.imin)) * n_profile;
		cf0 = flow_profile_t0[ncell];
		cf1 = flow_profile_t1[ncell];
		interpolate_condition(*cf0, *cf1, w, *cf_interp);
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*cf_interp);
		if (n_profile == 2) {
		    cf0 = flow_profile_t0[ncell+1];
		    cf1 = flow_profile_t1[ncell+1];
		    interpolate_condition(*cf0, *cf1, w, *cf_interp);
		}
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*cf_interp);
	    } // end i loop
	} // end k loop
	break;
    case WEST:
	i = bd.imin;
        for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	    for ( j = bd.jmin; j <= bd.jmax; ++j ) {
		ncell = ((k - bd.kmin)*bd.nnj + (j - bd.jmin)) * n_profile;
		cf0 = flow_profile_t0[ncell];
		cf1 = flow_profile_t1[ncell];
		interpolate_condition(*cf0, *cf1, w, *cf_interp);
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*cf_interp);
		if (n_profile == 2) {
		    cf0 = flow_profile_t0[ncell+1];
		    cf1 = flow_profile_t1[ncell+1];
		    interpolate_condition(*cf0, *cf1, w, *cf_interp);
		}
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*cf_interp);
	    } // end j loop
	} // end k loop
	break;
    case TOP:
	k = bd.kmax;
	for ( j = bd.jmin; j <= bd.jmax; ++j ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		ncell = ((j - bd.jmin)*bd.nni + (i - bd.imin)) * n_profile;
		cf0 = flow_profile_t0[ncell];
		cf1 = flow_profile_t1[ncell];
		interpolate_condition(*cf0, *cf1, w, *cf_interp);
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*cf_interp);
		if (n_profile == 2) {
		    cf0 = flow_profile_t0[ncell+1];
		    cf1 = flow_profile_t1[ncell+1];
		    interpolate_condition(*cf0, *cf1, w, *cf_interp);
		}
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*cf_interp);
	    } // end i loop
	} // end j loop
	break;
    case BOTTOM:
	k = bd.kmin;
	for ( j = bd.jmin; j <= bd.jmax; ++j ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		ncell = ((j - bd.jmin)*bd.nni + (i - bd.imin)) * n_profile;
		cf0 = flow_profile_t0[ncell];
		cf1 = flow_profile_t1[ncell];
		interpolate_condition(*cf0, *cf1, w, *cf_interp);
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*cf_interp);
		if (n_profile == 2) {
		    cf0 = flow_profile_t0[ncell+1];
		    cf1 = flow_profile_t1[ncell+1];
		    interpolate_condition(*cf0, *cf1, w, *cf_interp);
		}
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*cf_interp);
	    } // end i loop
	} // end j loop
    } // end switch which_boundary

    return SUCCESS;
} // end TransientProfileBC::apply_convective

void
TransientProfileBC::
read_profile_into_t0()
{
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
	    flow_profile_t0.push_back(new CFlowCondition(gmodel,p,u,v,w,T,massf,"",tke,omega,mu_t,k_t,S));
	} // end for i_profile
    } // end for ncell
}

void
TransientProfileBC::
read_profile_into_t1()
{
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

    // First clear out t1 flow profiles in case we've previously used the vector.
    flow_profile_t1.clear();
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
	    flow_profile_t1.push_back(new CFlowCondition(gmodel,p,u,v,w,T,massf,"",tke,omega,mu_t,k_t,S));
	} // end for i_profile
    } // end for ncell
}

void
TransientProfileBC::
update_profile_t0_with_t1()
{
    t0 = t1;
    for ( size_t i = 0; i < flow_profile_t0.size(); ++i ) {
	    flow_profile_t0[i] = flow_profile_t1[i];
    }
}

void
TransientProfileBC::
interpolate_condition(CFlowCondition &cf0, CFlowCondition &cf1, double w, CFlowCondition &cf_interp)
{
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();

    cf_interp.gas->p = (1.0 - w) * cf0.gas->p + w * cf1.gas->p;
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	cf_interp.gas->massf[isp] = (1.0 - w) * cf0.gas->massf[isp] + w * cf1.gas->massf[isp];
    }
    // We'll need to scale after this interpolation step
    scale_mass_fractions(cf_interp.gas->massf);
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	cf_interp.gas->T[imode] = (1.0 - w) * cf0.gas->T[imode] + w * cf1.gas->T[imode];
    }
    gmodel->eval_thermo_state_pT(*(cf_interp.gas));
    gmodel->eval_transport_coefficients(*(cf_interp.gas));
    cf_interp.u = (1.0 - w) * cf0.u + w * cf1.u;
    cf_interp.v = (1.0 - w) * cf0.v + w * cf1.v;
    cf_interp.w = (1.0 - w) * cf0.w + w * cf1.w;
    cf_interp.tke = (1.0 - w) * cf0.tke + w * cf1.tke;
    cf_interp.omega = (1.0 - w) * cf0.omega + w * cf1.omega;
    cf_interp.mu_t = (1.0 - w) * cf0.mu_t + w * cf1.mu_t;
    cf_interp.k_t = (1.0 - w) * cf0.k_t + w * cf1.k_t;
    // Can't interpolate on S.
    // If either t0 or t1 is 1, then set intermediate as 1.
    if ( cf0.S == 1 || cf1.S == 1 )
	cf_interp.S = 1;
    else
	cf_interp.S = 0;
    
}

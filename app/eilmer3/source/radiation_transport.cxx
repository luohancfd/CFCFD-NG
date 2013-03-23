/** \file radiation_transport.cxx
 *  \ingroup eilmer3
 *  \brief Definitions for the radiation transport model class and functions.
 **/

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../../lib/radiation/source/spectral_model.hh"
#include "../../../lib/radiation/source/LOS_pieces.hh"
#include "../../../lib/util/source/useful.h"
#include "../../../lib/nm/source/exponential_integrals.hh"

#include "radiation_transport.hh"
#include "cell.hh"
#include "kernel.hh"
#include "bc.hh"

const int VERBOSE_RADIATION_TRANSPORT = 1;
const int NO_CLUSTERING = 0;
const int CLUSTERING_BY_VOLUME = 1;
const int CLUSTERING_BY_AREA = 2;
const int STANDARD_ABSORPTION = 0;
const int PARTITIONED_ENERGY_ABSORPTION = 1;


using namespace std;

/* --------- Model: "RadiationTransportModel: ----------- */
RadiationTransportModel::RadiationTransportModel( lua_State *L ) {}

RadiationTransportModel::~RadiationTransportModel()
{
    for ( size_t ithread=0; ithread<rsm_.size(); ++ithread )
    	delete rsm_[ithread];
}

int RadiationTransportModel::set_radiation_spectral_model( string file_name, size_t nthreads )
{
    // Create one spectral model per thread
    int flag = SUCCESS;
    for ( size_t ithread=0; ithread<nthreads; ++ithread ) {
    	rsm_.push_back( create_radiation_spectral_model( file_name ) );
    	if ( rsm_.back()==0 ) flag = FAILURE;
    }
    
    return flag;
}

/* --------- Model: "NoRadiationTransportModel" --------- */

NoRadiationTransportModel::NoRadiationTransportModel( void ) {}

NoRadiationTransportModel::~NoRadiationTransportModel() {}

string
NoRadiationTransportModel::str() const
{
    return "NoRadiationTransportModel";
}

int
NoRadiationTransportModel::initialise()
{
    // Do nothing
    
    return SUCCESS;
}

void NoRadiationTransportModel::compute_Q_rad_for_flowfield()
{
    // Do nothing
    
    return;
}

/* --------- Model: "OpticallyThin" --------- */

OpticallyThin::OpticallyThin( lua_State *L )
: RadiationTransportModel(L)
{
    // Why aren't we using get_int()? Booleans makes the input file more neat
    spectrally_resolved_ = static_cast<int>(get_boolean(L,-1,"spectrally_resolved"));
}

OpticallyThin::~OpticallyThin() {}

string
OpticallyThin::str() const
{
    return "OpticallyThin";
}

int
OpticallyThin::initialise()
{
    // Do nothing
    
    return SUCCESS;
}

void OpticallyThin::compute_Q_rad_for_flowfield()
{
    Block * bdp;
    global_data &G = *get_global_data_ptr();  // set up a reference

    size_t jb;
#   ifdef _OPENMP
#   pragma omp barrier
#   pragma omp parallel for private(jb) schedule(runtime)
#   endif
    for ( jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	if ( bdp->active != 1 ) continue;
	this->compute_Q_rad_for_block(bdp);
    }
    
    return;
}

void OpticallyThin::compute_Q_rad_for_block( Block * A )
{
    // No flowfield reabsorption, 100% emission
    FV_Cell *cellp;
    int ithread = omp_get_thread_num();
    
#   if VERBOSE_RADIATION_TRANSPORT
    size_t total_cells = A->nnk * A->nnj * A->nni;
    size_t count = 0;
#   endif

    for ( size_t k = A->kmin; k <= A->kmax; ++k ) {
	for ( size_t i = A->imin; i <= A->imax; ++i ) {
	    for ( size_t j = A->jmin; j <= A->jmax; ++j ) {
#		if VERBOSE_RADIATION_TRANSPORT
		if ( A->id==0 ) {
		    cout << "Computing spectra for cell: " << count << " of: " << total_cells << " in block: " << A->id << endl;
		    ++count;
		}
#               endif
		cellp = A->get_cell(i,j,k);
		cellp->Q_rE_rad = ( - 4.0 * M_PI ) * \
		rsm_[ithread]->radiative_integrated_emission_for_gas_state(*cellp->fs->gas, spectrally_resolved_);
	    }
	}
    }
    
    cout << "Finished block: " << A->id << endl;
    
    return;
}

/* --------- Model: "TangentSlab" --------- */

TangentSlab::TangentSlab( lua_State *L )
 : RadiationTransportModel(L)
{
    exact_ = get_boolean( L, -1, "exact_formulation" );
}

TangentSlab::~TangentSlab()
{}

string
TangentSlab::str() const
{
    return "TangentSlab";
}

int
TangentSlab::initialise()
{
    // Do nothing
    
    return SUCCESS;
}

void TangentSlab::compute_Q_rad_for_flowfield()
{
    // FIXME: eventually want to be able to do calculation across blocks
    
    Block * bdp;
    global_data &G = *get_global_data_ptr();  // set up a reference
    
    size_t jb;
#   ifdef _OPENMP
#   pragma omp barrier
#   pragma omp parallel for private(jb) schedule(runtime)
#   endif
    for ( jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	if ( bdp->active != 1 ) continue;
	this->compute_Q_rad_for_block(bdp);
    }
    
    return;
}

void TangentSlab::compute_Q_rad_for_block( Block * A )
{
#   if VERBOSE_RADIATION_TRANSPORT
    cout << "TangentSlab::calculate_source_term_for_block()" << endl;
#   endif
    FV_Cell * cellp;
    Vector3 * coords_im1;
    Vector3 * coords_i;
    double s = 0.0;
    size_t index;
    
    // Perform TS calc for lines of cells in i space
    // NOTE: assuming left to right flow
    for ( size_t k = A->kmin; k <= A->kmax; ++k ) {
    	for ( size_t j = A->jmin; j <= A->jmax; ++j ) {
#           if VERBOSE_RADIATION_TRANSPORT
    	    cout << "performing TS calc for j = " << j << endl;
#           endif
            // Initialise TS_data with the number of points to be used
    	    TS_data TS = TS_data( rsm_[omp_get_thread_num()], A->nni );
    	    // west IFace of first cell is s=0.0
    	    coords_im1 = &(A->get_cell(A->imin,j,k)->iface[WEST]->pos);
    	    TS.T_i_ = A->get_cell(A->imin,j,k)->iface[WEST]->fs->gas->T[0];
    	    for ( size_t i = A->imin; i <= A->imax; ++i ) {
    	    	cellp = A->get_cell(i,j,k);
    	    	coords_i = &(cellp->pos);
    	    	s += vabs( *coords_i - *coords_im1 );
    	    	TS.set_rad_point( i - A->imin, cellp->fs->gas, &(cellp->Q_rE_rad), s, cellp->iLength );
    	    	coords_im1 = coords_i;
    	    }
    	    // use east IFace of last cell for s_f and T_f
    	    TS.T_f_ = A->get_cell(A->imax,j,k)->iface[EAST]->fs->gas->T[0];
	    // index into 1D heat-flux vectors (i parts are omitted as we assume this is always the easterly face)
	    index = (A->jmax-A->jmin+1)*(k-A->kmin) + (j-A->jmin);
	    // perform tangent-slab integration
            if ( exact_ ) {
	        // Exact but slow version
	        A->bcp[EAST]->q_rad[index] = TS.exact_solve_for_divq();
            }
            else {
	        // Approximate but speedy version
    	        A->bcp[EAST]->q_rad[index] = TS.quick_solve_for_divq();
            }
        }
    }
    
    return;
}

/* --------- Model: "DiscreteTransfer" --------- */

DiscreteTransfer::DiscreteTransfer( lua_State *L )
 : RadiationTransportModel(L)
{
    nrays_ = get_int(L,-1,"nrays");
    
    dl_lmin_ratio_ = 0.9;	// ratio of RT step to min cell length scale
    
    dl_min_ = 1.0e-10;		// Minimum step length in metres
    
    E_min_ = 1.0e-30;		// min energy per photon packet to be considered
    
    string clustering = get_string(L,-1,"clustering");
    if ( clustering == "by volume" )
    	clustering_ = CLUSTERING_BY_VOLUME;
    else if ( clustering == "by area" )
    	clustering_ = CLUSTERING_BY_AREA;
    else
    	clustering_ = NO_CLUSTERING;
    
    string binning = get_string(L,-1,"binning");
    if ( binning == "frequency" )
	binning_ = FREQUENCY_BINNING;
    else if ( binning == "opacity" )
	binning_ = OPACITY_BINNING;
    else
	binning_ = NO_BINNING;
    if ( binning_ )
	N_bins_ = get_int(L,-1,"N_bins");

    // Done.
}

DiscreteTransfer::~DiscreteTransfer()
{
    for ( size_t ib=0; ib<cells_.size(); ++ib ) {
    	for ( size_t ic=0; ic<cells_[ib].size(); ++ic ) {
    	    delete cells_[ib][ic];
    	}
    }
    
    for ( size_t ib=0; ib<interfaces_.size(); ++ib ) {
    	for ( size_t iface=0; iface<interfaces_[ib].size(); ++iface ) {
    	    delete interfaces_[ib][iface];
    	}
    }
    
    delete cf_;
}

string
DiscreteTransfer::str() const
{
    return "DiscreteTransfer";
}

int
DiscreteTransfer::initialise()
{
    // initialise the cell finder
    global_data &G = *get_global_data_ptr();
    if ( G.dimensions == 2 ) {
    	cf_ = new CellFinder2D();	// use default nvertices
    }
    else {
    	cf_ = new CellFinder3D();	// use default nvertices
    }
    
    // initialise cells and interfaces from block data
    // NOTE: not doing in parallel as this should be very quick
    cells_.resize( G.nblock );
    interfaces_.resize( G.nblock );
    size_t nthreads = omp_get_max_threads();
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	Block * bdp = G.my_blocks[jb];
#       if VERBOSE_RADIATION_TRANSPORT
    	cout << "Thread " << omp_get_thread_num() << ": Initialising cells in block: " << jb << endl;
#       endif
	// if ( bdp->active != 1 ) continue;
	for ( size_t k = bdp->kmin; k <= bdp->kmax; ++k ) {
	    for ( size_t j = bdp->jmin; j <= bdp->jmax; ++j ) {
	    	for ( size_t i = bdp->imin; i <= bdp->imax; ++i ) {
	    	    FV_Cell * cell = bdp->get_cell(i,j,k);
	    	    cells_[jb].push_back( new DiscreteTransferCell( cell->fs->gas, &(cell->Q_rE_rad), cell->pos, cell->volume ) );
	    	    cells_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    cells_[jb].back()->Q_rE_rad_temp_.resize( nthreads );
	    	    // test for and create wall interfaces
	    	    if ( j==bdp->jmin && bdp->bcp[SOUTH]->is_wall_flag ) {
	    	    	interfaces_[jb].push_back( new DiscreteTransferInterface( cell->iface[SOUTH]->fs->gas, cell->iface[SOUTH]->pos, cell->iface[SOUTH]->area, cell->iface[SOUTH]->length ) );
	    	    	interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    }
	    	    else if ( j==bdp->jmax && bdp->bcp[NORTH]->is_wall_flag ) {
	    	    	interfaces_[jb].push_back( new DiscreteTransferInterface( cell->iface[NORTH]->fs->gas, cell->iface[NORTH]->pos, cell->iface[NORTH]->area, cell->iface[NORTH]->length ) );
	    	    	interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    }
	    	    if ( i==bdp->imin && bdp->bcp[WEST]->is_wall_flag ) {
	    	    	interfaces_[jb].push_back( new DiscreteTransferInterface( cell->iface[WEST]->fs->gas, cell->iface[WEST]->pos, cell->iface[WEST]->area, cell->iface[WEST]->length ) );
	    	    	interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    }
	    	    else if ( i==bdp->imax && bdp->bcp[EAST]->is_wall_flag ) {
	    	    	interfaces_[jb].push_back( new DiscreteTransferInterface( cell->iface[EAST]->fs->gas, cell->iface[EAST]->pos, cell->iface[EAST]->area, cell->iface[EAST]->length ) );
	    	    	interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    }
	    	    if ( G.dimensions == 3 ) {
	    	    	if ( k==bdp->kmin && bdp->bcp[BOTTOM]->is_wall_flag ) {
	    	    	    interfaces_[jb].push_back( new DiscreteTransferInterface( cell->iface[BOTTOM]->fs->gas, cell->iface[BOTTOM]->pos, cell->iface[BOTTOM]->area, cell->iface[BOTTOM]->length ) );
	    	    	    interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	    interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    	}
	    	    	else if ( k==bdp->kmax && bdp->bcp[TOP]->is_wall_flag ) {
	    	    	    interfaces_[jb].push_back( new DiscreteTransferInterface( cell->iface[TOP]->fs->gas, cell->iface[TOP]->pos, cell->iface[TOP]->area, cell->iface[TOP]->length ) );
	    	    	    interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	    interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    	}
	    	    }
	    	}
	    }
	}
	
#       if VERBOSE_RADIATION_TRANSPORT
	cout << "- created " << cells_[jb].size() << " DiscreteTransferCell's" << endl;
	cout << "- created " << interfaces_[jb].size() << " DiscreteTransferInterface's" << endl;
#       endif
    }
    
    return SUCCESS;
}

void DiscreteTransfer::compute_Q_rad_for_flowfield()
{
    // Initialise geometry data
    global_data &G = *get_global_data_ptr();
    bool planar_flag = false;
    if ( G.dimensions == 2 && !get_axisymmetric_flag() ) planar_flag = true;
    
    // 0. Loop over spectral blocks (CHECKME)
    for ( size_t isb=0; isb < static_cast<size_t>(rsm_[0]->get_spectral_blocks()); ++isb ) {
	for ( size_t irsm=0; irsm<rsm_.size(); ++irsm ) {
	    rsm_[irsm]->reset_spectral_params( isb );
	}
    
	// 1. Set all source terms to zero and store all spectra
	for ( size_t ib=0; ib<cells_.size(); ++ib ) {
	    size_t ic;
#	    ifdef _OPENMP
#	    pragma omp parallel for private(ic) schedule(runtime)
#	    endif
	    for ( ic=0; ic<cells_[ib].size(); ++ic ) {
	    	RayTracingCell * cell = cells_[ib][ic];
		*(cell->Q_rE_rad_) = 0.0;
		// Also make sure thread vector is zero
		for ( size_t iQ=0; iQ<cell->Q_rE_rad_temp_.size(); ++ iQ ) {
		    cell->Q_rE_rad_temp_[iQ] = 0.0;
		}
		cout << "Thread " << omp_get_thread_num() << ": Recomputing spectra for cell: " << ic << " in block: " << ib;
		cells_[ib][ic]->recompute_spectra( rsm_[omp_get_thread_num()] );
		cout << " - j_total = " << cells_[ib][ic]->X_->integrate_emission_spectra() << endl;
	    }
	    size_t iface;
#	    ifdef _OPENMP
#	    pragma omp barrier
#	    pragma omp parallel for private(iface) schedule(runtime)
#	    endif
	    for ( iface=0; iface<interfaces_[ib].size(); ++iface ) {
	    	RayTracingInterface * interface = interfaces_[ib][iface];
		cout << "Thread " << omp_get_thread_num() << ": Recomputing spectra for interface: " << iface << " in block: " << ib;
		interface->recompute_spectra( rsm_[omp_get_thread_num()] );
		cout << " - I_total = " << interface->S_->integrate_intensity_spectra() << endl;
	    }
	}
	
	// 2. Compute cell_E_rad_total_max and interface_E_rad_total_max
	// NOTE: easiest to not do this in parallel
	double cell_E_rad_total_max = 0.0;
	double interface_E_rad_total_max = 0.0;
	for ( size_t ib=0; ib<cells_.size(); ++ib ) {
	    for ( size_t ic=0; ic<cells_[ib].size(); ++ic ) {
	    	RayTracingCell * cell = cells_[ib][ic];
		if ( clustering_==CLUSTERING_BY_VOLUME )
		    cell->E_rad_ = cell->X_->j_int.back() * cell->vol_;
		else if ( clustering_==CLUSTERING_BY_AREA )
		    cell->E_rad_ = cell->X_->j_int.back() * cell->area_;
		if ( cell->E_rad_ > cell_E_rad_total_max ) cell_E_rad_total_max = cell->E_rad_;
	    }
	    for ( size_t iface=0; iface<interfaces_[ib].size(); ++iface ) {
	    	RayTracingInterface * interface = interfaces_[ib][iface];
		if ( clustering_==CLUSTERING_BY_VOLUME )
		    interface->E_rad_ = interface->S_->I_int.back() * interface->area_;
		else if ( clustering_==CLUSTERING_BY_AREA )
		    interface->E_rad_ = interface->S_->I_int.back() * interface->length_;
		if ( interface->E_rad_ > interface_E_rad_total_max ) interface_E_rad_total_max = interface->E_rad_;
	    }
	}
	
	// 3. Create spectral bins
	if ( binning_==FREQUENCY_BINNING ) {
	    N_bins_ = create_spectral_bin_vector( cells_[0][0]->X_->nu, binning_, N_bins_, B_ );
	}
	else if ( binning_==OPACITY_BINNING ) {
	    // We need to solve for the spatially independent mean opacity/absorption (see Eq 2.2 of Wray, Ripoll and Prabhu)
	    vector<double> kappa_mean(rsm_[0]->get_spectral_points());
	    size_t inu;
#	    ifdef _OPENMP
#   	    pragma omp barrier
#	    pragma omp parallel for private(inu) schedule(runtime)
#	    endif
	    for ( inu=0; inu < static_cast<size_t>(rsm_[omp_get_thread_num()]->get_spectral_points()); inu++ ) {
		double j_nu_dV = 0.0, S_nu_dV = 0.0;
		for ( size_t ib=0; ib<cells_.size(); ++ib ) {
		    for ( size_t ic=0; ic<cells_[ib].size(); ++ic ) {
			RayTracingCell * cell = cells_[ib][ic];
			double dj_nu_dV = cell->X_->j_nu[inu] * cell->vol_;
			j_nu_dV += dj_nu_dV;
			if ( cell->X_->kappa_nu[inu] > 0.0 )
			    S_nu_dV += dj_nu_dV / cell->X_->kappa_nu[inu];
		    }
		}
		// If the source function is zero then kappa should also be zero
		kappa_mean[inu] = ( S_nu_dV==0.0 ) ? 0.0 : j_nu_dV / S_nu_dV;
	    }
	    N_bins_ = create_spectral_bin_vector( kappa_mean, binning_, N_bins_, B_ );
	}

	// 4. Create binned spectra
	if ( binning_ ) {
	    for ( size_t ib=0; ib<cells_.size(); ++ib ) {
		size_t ic;
#	        ifdef _OPENMP
#	        pragma omp parallel for private(ic) schedule(runtime)
#  	        endif
		for ( ic=0; ic<cells_[ib].size(); ++ic ) {
		    RayTracingCell * cell = cells_[ib][ic];
		    cout << "Thread " << omp_get_thread_num() << ": Creating binned spectra for cell: " << ic << " in block: " << ib;
		    cell->Y_ = new BinnedCoeffSpectra( cell->X_, B_ );
		    cout << " - j_total from spectra = " << cell->X_->integrate_emission_spectra() << endl;
		    cout << " - j_total from binning = " << cell->Y_->sum_emission() << endl;
		}
		size_t iface;
#	        ifdef _OPENMP
#	        pragma omp barrier
#	        pragma omp parallel for private(iface) schedule(runtime)
#	        endif
		for ( iface=0; iface<interfaces_[ib].size(); ++iface ) {
		    RayTracingInterface * interface = interfaces_[ib][iface];
		    cout << "Thread " << omp_get_thread_num() << ": Creating binned spectra for interface: " << iface << " in block: " << ib;
		    interface->U_ = new BinnedSpectralIntensity( interface->S_, B_ );
		    cout << " - I_total from spectra = " <<  interface->S_->integrate_intensity_spectra() << endl;
		    cout << " - I_total from binning = " << interface->U_->sum_intensity() << endl;
		}
	    }
	}

	// 3. create rays
	for ( size_t ib=0; ib<cells_.size(); ++ib ) {
	    size_t ii, jj, kk;
	    for ( size_t ic=0; ic<cells_[ib].size(); ++ic ) {
		RayTracingCell * cell = cells_[ib][ic];
		size_t nrays = nrays_;
		if ( clustering_ ) {
		    nrays = int( double(nrays_) * cell->E_rad_ / cell_E_rad_total_max );
		}
		if ( nrays==0 ) nrays = 1;
		this->initialise_rays_for_cell( cell, nrays, G.dimensions, planar_flag );
#               if VERBOSE_RADIATION_TRANSPORT
		cout << "Thread " << omp_get_thread_num() << ": Tracing rays for cell: " << ic << " in block: " << ib << endl;
#               endif
		cell->get_CFD_cell_indices( ii, jj, kk );
		// Cycle over all rays
		size_t ir;
#   	        ifdef _OPENMP
#               pragma omp parallel for private(ir) schedule(runtime)
#               endif
		for ( ir=0; ir<cell->rays_.size(); ++ir ) {
		    RayTracingRay * ray = cell->rays_[ir];
		    // use ii and jj as initial cell guess
		    this->trace_ray( ray, ib, ii, jj, kk );
		}
	    }
	    for ( size_t iface=0; iface<interfaces_[ib].size(); ++iface ) {
		RayTracingInterface * interface = interfaces_[ib][iface];
	        size_t nrays = nrays_;
	        if ( clustering_ ) {
		    nrays = int( double(nrays_) * interface->E_rad_ / interface_E_rad_total_max );
	        }
		if ( nrays==0 ) nrays = 1;
		this->initialise_rays_for_interface( interface, nrays_, G.dimensions, planar_flag );
#               if VERBOSE_RADIATION_TRANSPORT
		cout << "Thread " << omp_get_thread_num() << ": Tracing rays for interface: " << iface << " in block: " << ib << endl;
#               endif
		interface->get_CFD_cell_indices( ii, jj, kk );
		// Cycle over all rays
		size_t ir;
#   	        ifdef _OPENMP
#               pragma omp parallel for private(ir) schedule(runtime)
#               endif
		for ( ir=0; ir<interface->rays_.size(); ++ir ) {
		    RayTracingRay * ray = interface->rays_[ir];
		    // use ii and jj as initial cell guess
		    this->trace_ray( ray, ib, ii, jj, kk );
		    ++ray;
		}
	    }
	}
	
#       if 0
	// write a random cells rays to file 
	srand((unsigned)time(0));
	ib = rand()%int(cells_.size());
	size_t ic = rand()%int(cells_[ib].size());
	cout << "writing rays to file from block: " << ib << ", cell: " << ic << endl;
	cells_[ib][ic]->write_rays_to_file("random_cell_rays.txt");
#       endif
	
	// 4. perform transport of radiative energy throughout the grid
	for ( size_t ib=0; ib<cells_.size(); ++ib ) {
#	    if VERBOSE_RADIATION_TRANSPORT
	    cout << "Thread " << omp_get_thread_num() << ": Integrating along rays for block: " << ib << endl;
#           endif
	    for ( size_t ic=0; ic<cells_[ib].size(); ++ic ) {
#		if VERBOSE_RADIATION_TRANSPORT
	    	cout << "Thread " << omp_get_thread_num() << ": Integrating along rays for cell: " << ic << " in block: " << ib << endl;
#               endif
		RayTracingCell * cell = cells_[ib][ic];
		size_t ir;                                              
#		ifdef _OPENMP
#   	    	pragma omp barrier
#		pragma omp parallel for private(ir) schedule(runtime)
#		endif
		for ( ir=0; ir<cell->rays_.size(); ++ir ) {
		    RayTracingRay * ray = cell->rays_[ir];
		    if ( binning_ ) {
			// perform radiation transport with the binned spectra
			for ( size_t iB=0; iB<N_bins_; ++iB ) {
			    // calculate frequency interval for this frequency point
			    // NOTE: we are taking care here to remain consistent with a trapezoidal discretisation of the spectra
			    //       where the trapezoid heights are taken as the average of the two bounding points
			    // calculate energy of this photon packet
			    double E = cell->Y_->j_bin[iB] * cell->vol_ * ray->domega_;
			    if ( E < E_min_ ) continue;
			    // subtract emitted energy from origin cell
			    cell->Q_rE_rad_temp_[omp_get_thread_num()] -= E / cell->vol_;
			    // step along the ray path
			    double L = 0.0, dL;
			    for ( size_t ip=0; ip<ray->points_.size(); ++ip ) {
				RayTracingPoint * point = ray->points_[ip];
				// calculate absorbed energy
				double kappa_bin = point->Y_->kappa_bin[iB];
				dL = point->s_ - L; L += dL;
				double dE = ( 1.0 - exp( - kappa_bin * dL ) ) * E; E -= dE;
				// add absorbed energy to traversed cell
				(*point->Q_rE_rad_temp_)[omp_get_thread_num()] += dE / point->vol_;
			    }
			    // store the remaining energy to be dumped onto the wall element
			    ray->E_exit_ += E;
			}
		    }
		    else {
			size_t nnu = rsm_[omp_get_thread_num()]->get_spectral_points();
			for ( size_t inu=0; inu<nnu; ++inu ) {
			    // calculate frequency interval for this frequency point
			    // NOTE: we are taking care here to remain consistent with a trapezoidal discretisation of the spectra
			    //       where the trapezoid heights are taken as the average of the two bounding points
			    double dnu=0.0;
			    if      ( inu==0 )     dnu = 0.5 * fabs( cell->X_->nu[inu+1] - cell->X_->nu[inu] );
			    else if ( inu==nnu-1 ) dnu = 0.5 * fabs( cell->X_->nu[inu] - cell->X_->nu[inu-1] );
			    else                   dnu = 0.5 * fabs( cell->X_->nu[inu] - cell->X_->nu[inu-1] ) + 0.5 * fabs( cell->X_->nu[inu+1] - cell->X_->nu[inu] );
			    // calculate energy of this photon packet
			    double E = cell->X_->j_nu[inu] * cell->vol_ * ray->domega_ * dnu;
			    if ( E < E_min_ ) continue;
			    // subtract emitted energy from origin cell
			    cell->Q_rE_rad_temp_[omp_get_thread_num()] -= E / cell->vol_;
			    // step along the ray path
			    double L = 0.0, dL;
			    for ( size_t ip=0; ip<ray->points_.size(); ++ip ) {
				RayTracingPoint * point = ray->points_[ip];
				// calculate absorbed energy
				double kappa_nu = point->X_->kappa_nu[inu];
				dL = point->s_ - L; L += dL;
				double dE = ( 1.0 - exp( - kappa_nu * dL ) ) * E; E -= dE;
				// add absorbed energy to traversed cell
				(*point->Q_rE_rad_temp_)[omp_get_thread_num()] += dE / point->vol_;
			    }
			    // store the remaining energy to be dumped onto the wall element
			    ray->E_exit_ += E;
			}
		    }
		}
	    }
	    for ( size_t iface=0; iface<interfaces_[ib].size(); ++iface ) {
#		if VERBOSE_RADIATION_TRANSPORT
	    	cout << "Thread " << omp_get_thread_num() << ": Integrating along rays for interface: " << iface << " in block: " << ib << endl;
#               endif
		RayTracingInterface * interface = interfaces_[ib][iface];
		size_t ir;
#		ifdef _OPENMP
#   	    	pragma omp barrier
#		pragma omp parallel for private(ir) schedule(runtime)
#		endif
		for ( ir=0; ir<interface->rays_.size(); ++ir ) {
		    RayTracingRay * ray = interface->rays_[ir];
		    if ( binning_ ) {
			for ( size_t iB=0; iB<N_bins_; ++iB ) {
			    // calculate energy of this photon packet
			    double E = interface->U_->I_bin[iB] * interface->area_ * ray->domega_;
			    if ( E < E_min_ ) continue;
			    // step along the ray path
			    double L = 0.0, dL;
			    for ( size_t ip=0; ip<ray->points_.size(); ++ip ) {
				RayTracingPoint * point = ray->points_[ip];
				// calculate absorbed energy
				double kappa_bin = point->Y_->kappa_bin[iB];
				dL = point->s_ - L; L += dL;
				double dE = ( 1.0 - exp( - kappa_bin * dL ) ) * E; E -= dE;
				// add absorbed energy to traversed cell
				(*point->Q_rE_rad_temp_)[omp_get_thread_num()] += dE / point->vol_;
			    }
			    // store the remaining energy to be dumped onto the wall element
			    // NOTE: only if the number of points is non-zero to avoid those
			    //       rays that are directed into the wall
			    if ( ray->points_.size() > 0 ) ray->E_exit_ += E;
			}
		    }
		    else {
			size_t nnu = rsm_[omp_get_thread_num()]->get_spectral_points();
			for ( size_t inu=0; inu<nnu; ++inu ) {
			    // calculate frequency interval for this frequency point
			    // NOTE: we are taking care here to remain consistent with a trapezoidal discretisation of the spectra
			    //       where the trapezoid heights are taken as the average of the two bounding points
			    double dnu=0.0;
			    if      ( inu==0 )     dnu = 0.5 * fabs( interface->S_->nu[inu+1] - interface->S_->nu[inu] );
			    else if ( inu==nnu-1 ) dnu = 0.5 * fabs( interface->S_->nu[inu] - interface->S_->nu[inu-1] );
			    else                   dnu = 0.5 * fabs( interface->S_->nu[inu] - interface->S_->nu[inu-1] ) + 0.5 * fabs( interface->S_->nu[inu+1] - interface->S_->nu[inu] );
			    // calculate energy of this photon packet
			    double E = interface->S_->I_nu[inu] * interface->area_ * ray->domega_ * dnu;
			    if ( E < E_min_ ) continue;
			    // step along the ray path
			    double L = 0.0, dL;
			    for ( size_t ip=0; ip<ray->points_.size(); ++ip ) {
				RayTracingPoint * point = ray->points_[ip];
				// calculate absorbed energy
				double kappa_nu = point->X_->kappa_nu[inu];
				dL = point->s_ - L; L += dL;
				double dE = ( 1.0 - exp( - kappa_nu * dL ) ) * E; E -= dE;
				// add absorbed energy to traversed cell
				(*point->Q_rE_rad_temp_)[omp_get_thread_num()] += dE / point->vol_;
			    }
			    // store the remaining energy to be dumped onto the wall element
			    // NOTE: only if the number of points is non-zero to avoid those
			    //       rays that are directed into the wall
			    if ( ray->points_.size() > 0 ) ray->E_exit_ += E;
			}
		    }
		}
	    }
	}
	
    } // end spectral block loop
    
    
    // 5. Sum Q_rE_rad_temp values and set Q_rE_rad in CFD cells
    for ( size_t ib=0; ib<cells_.size(); ++ib ) {
    	size_t ic;
#   	ifdef _OPENMP
#   	pragma omp barrier
#   	pragma omp parallel for private(ic) schedule(runtime)
#  	endif
    	for ( ic=0; ic<cells_[ib].size(); ++ic ) {
    	    RayTracingCell * cell = cells_[ib][ic];
    	    for ( size_t iQ=0; iQ<cell->Q_rE_rad_temp_.size(); ++iQ ) {
    	    	*(cell->Q_rE_rad_) += cell->Q_rE_rad_temp_[iQ];
    	    }
    	    // cout << "*(cell->Q_rE_rad_) [ ib = " << ib << ", ic = " << ic << " ] = " << *(cell->Q_rE_rad_) << endl;
    	}
    }
    
    // 6. dump all exiting energy onto appropriate wall elements
    // NOTE: cannot do this with multiple threads due to possible race condition
    for ( size_t ib=0; ib<cells_.size(); ++ib ) {
    	for ( size_t ic=0; ic<cells_[ib].size(); ++ic ) {
    	    for ( size_t iray=0; iray<cells_[ib][ic]->rays_.size(); ++iray ) {
    	    	RayTracingRay * ray = cells_[ib][ic]->rays_[iray];
    	    	*(ray->q_rad_) += ray->E_exit_ / ray->exit_area_;
    	    }
    	}
    }
    
    for ( size_t ib=0; ib<interfaces_.size(); ++ib ) {
    	for ( size_t iface=0; iface<interfaces_[ib].size(); ++iface ) {
    	    for ( size_t iray=0; iray<interfaces_[ib][iface]->rays_.size(); ++iray ) {
    	    	RayTracingRay * ray = interfaces_[ib][iface]->rays_[iray];
    	    	*(ray->q_rad_) += ray->E_exit_ / ray->exit_area_;
    	    }
    	}
    }
    
    // clean up
    // all done in deconstructor for now
    
    return;
}


void DiscreteTransfer::initialise_rays_for_cell( RayTracingCell * cell, size_t nrays, size_t ndim, bool planar )
{
    // clear any existing rays
    cell->rays_.resize(0);
    
    if ( planar == true ) {
	// uniform planar ray distribution
	double alpha = 0.0, theta = 0.0, phi = 0.0, domega = 0.0;
	for ( size_t iray=0; iray<nrays; ++iray ) {
	    alpha = 2.0 * M_PI * iray / nrays;
	    if ( alpha>=0.0 && alpha<M_PI/2.0 ) {
		theta = 0.0; phi = alpha;
	    } else if ( alpha>=M_PI/2.0 && alpha<1.5*M_PI ) {
		theta = M_PI; phi = M_PI - alpha;
	    } else if ( alpha>=1.5*M_PI && alpha<2.0*M_PI ) {
		theta = 0.0; phi = alpha - M_PI*2.0;
	    }
	    domega = fabs( cos( phi ) ) * 2.0 * M_PI / double( nrays ) * M_PI ;
	    cell->rays_.push_back( new RayTracingRay2D( theta, phi, domega, cell->origin_ ) );
	}
    } else {
	// uniform 3D ray distribution - create rays via golden spiral method
	vector<Vector3> pts;
	double inc = M_PI * ( 3.0 - sqrt(5.0) );
	double off = 2.0 / double(nrays);
	double phi, theta = 0.0;
	for ( size_t iray=0; iray<nrays; ++iray ) {
	    double y = double(iray)*off - 1.0 + ( off/2.0 );
	    double r = sqrt( 1.0 - y*y );
	    phi = double(iray)*inc;
	    pts.push_back( Vector3(cos(phi)*r, y, sin(phi)*r) );
	}
	
	// Extract angles in our coordinate system
	double domega = 4.0 * M_PI / double ( nrays );
	for ( size_t iray=0; iray<nrays; ++iray ) {
	    double x=pts[iray].x, y=pts[iray].y, z=pts[iray].z;
	    double L = vabs( pts[iray] );
	    phi = asin( y / L );
	    // calculate theta
	    if      ( x >0.0 && z >0.0 ) theta = atan(z/x);
	    else if ( x==0.0 && z >0.0 ) theta = 0.5*M_PI;
	    else if ( x <0.0 && z >0.0 ) theta = M_PI - atan(z/-x);
	    else if ( x <0.0 && z==0.0 ) theta = 1.0*M_PI;
	    else if ( x <0.0 && z <0.0 ) theta = M_PI + atan(-z/x);
	    else if ( x==0.0 && z <0.0 ) theta = 1.5*M_PI;
	    else if ( x >0.0 && z <0.0 ) theta = 2.0*M_PI - atan(-z/x);
	    else if ( x >0.0 && z==0.0 ) theta = 0.0;
	    // else: impossible!
	    
	    if ( ndim==2 ) {
		cell->rays_.push_back( new RayTracingRay2D( theta, phi, domega, cell->origin_ ) );
	    }
	    else if ( ndim==3 ) {
		cell->rays_.push_back( new RayTracingRay3D( theta, phi, domega, cell->origin_ ) );
	    }
	}
    }
    
    return;
}

void DiscreteTransfer::initialise_rays_for_interface( RayTracingInterface * interface, size_t nrays, size_t ndim, bool planar )
{
    // clear any existing rays
    interface->rays_.resize(0);
    
    if ( planar == true ) {
	// uniform planar ray distribution
	double alpha = 0.0, theta = 0.0, phi = 0.0, domega = 0.0;
	for ( size_t iray=0; iray<nrays; ++iray ) {
	    alpha = 2.0 * M_PI * iray / nrays;
	    if ( alpha>=0.0 && alpha<M_PI/2.0 ) {
		theta = 0.0; phi = alpha;
	    } else if ( alpha>=M_PI/2.0 && alpha<1.5*M_PI ) {
		theta = M_PI; phi = M_PI - alpha;
	    } else if ( alpha>=1.5*M_PI && alpha<2.0*M_PI ) {
		theta = 0.0; phi = alpha - M_PI*2.0;
	    }
	    domega = fabs( cos( phi ) ) * 2.0 * M_PI / double( nrays ) * M_PI ;
	    interface->rays_.push_back( new RayTracingRay2D( theta, phi, domega, interface->origin_ ) );
	}
    } else {
	// uniform 3D ray distribution - create rays via golden spiral method
	vector<Vector3> pts;
	double inc = M_PI * ( 3.0 - sqrt(5.0) );
	double off = 2.0 / double(nrays);
	double phi, theta = 0.0;
	for ( size_t iray=0; iray<nrays; ++iray ) {
	    double y = double(iray)*off - 1.0 + ( off/2.0 );
	    double r = sqrt( 1.0 - y*y );
	    phi = double(iray)*inc;
	    pts.push_back( Vector3(cos(phi)*r, y, sin(phi)*r) );
	}
	
	// Extract angles in our coordinate system
	double domega = 4.0 * M_PI / double ( nrays );
	for ( size_t iray=0; iray<nrays; ++iray ) {
	    double x=pts[iray].x, y=pts[iray].y, z=pts[iray].z;
	    double L = vabs( pts[iray] );
	    phi = asin( y / L );
	    // calculate theta
	    if      ( x >0.0 && z >0.0 ) theta = atan(z/x);
	    else if ( x==0.0 && z >0.0 ) theta = 0.5*M_PI;
	    else if ( x <0.0 && z >0.0 ) theta = M_PI - atan(z/-x);
	    else if ( x <0.0 && z==0.0 ) theta = 1.0*M_PI;
	    else if ( x <0.0 && z <0.0 ) theta = M_PI + atan(-z/x);
	    else if ( x==0.0 && z <0.0 ) theta = 1.5*M_PI;
	    else if ( x >0.0 && z <0.0 ) theta = 2.0*M_PI - atan(-z/x);
	    else if ( x >0.0 && z==0.0 ) theta = 0.0;
	    // else: impossible!
	    
	    if ( ndim==2 ) {
		interface->rays_.push_back( new RayTracingRay2D( theta, phi, domega, interface->origin_ ) );
	    }
	    else if ( ndim==3 ) {
		interface->rays_.push_back( new RayTracingRay3D( theta, phi, domega, interface->origin_ ) );
	    }
	}
    }
    
    return;
}

int DiscreteTransfer::trace_ray( RayTracingRay * ray, size_t ib, size_t ic, size_t jc, size_t kc )
{
    Block * A = get_block_data_ptr( ib );
    FV_Cell * cell = A->get_cell(ic,jc,kc);		// start at origin cell
    
    double L = 0.5 * dl_lmin_ratio_ * cell->L_min;	// small initial step length
    Vector3 p = ray->get_point_on_line( L );
    
    RayTracingCell * RTcell;
    
    // step along the ray, creating the points as we go
    while ( ( ray->status_ = cf_->find_cell( p, ib, ic, jc, kc ) ) == INSIDE_GRID ) {
    	// Get pointers to the new block and cell
    	A = get_block_data_ptr( ib );
    	cell = A->get_cell(ic,jc,kc);
    	// create the RayTracingPoint
    	RTcell = cells_[ib][get_cell_index(A,ic,jc,kc)];
    	ray->points_.push_back( new RayTracingPoint( cell->fs->gas, &(cell->Q_rE_rad), L, cell->volume, RTcell->X_, RTcell->Y_, &(RTcell->Q_rE_rad_temp_) ) );
	// calculate next position on ray
	L += dl_lmin_ratio_ * cell->L_min;
	p = ray->get_point_on_line( L );
    }

    // Make sure block and cell pointers are up to date
    A = get_block_data_ptr( ib );
    cell = A->get_cell(ic,jc,kc);

    if ( ray->status_ == ERROR ) {
    	cout << "DiscreteTransfer::trace_ray()" << endl
    	     << "CellFinder failed at location: " << p.str() << endl
    	     << "Exiting program" << endl;
    	exit( FAILURE);
    }
    else {
    	// set q_rad pointer and exit area
    	size_t index = A->bcp[ ray->status_ ]->get_heat_flux_index( ic, jc, kc );
    	ray->q_rad_ = &(A->bcp[ ray->status_ ]->q_rad[ index ]);
    	ray->exit_area_ = cell->iface[ ray->status_ ]->area;
    }
    
    return SUCCESS;
}

/* --------- Model: "MonteCarlo" --------- */

MonteCarlo::MonteCarlo( lua_State *L )
 : RadiationTransportModel(L)
{
    nrays_ = get_int(L,-1,"nrays");	// average number of rays/photons per cell/interface
    
    dl_lmin_ratio_ = 0.3;	// ratio of RT step to min cell length scale
    
    dl_min_ = 1.0e-15;		// Minimum step length in metres
    
    E_min_ = 1.0e-30;		// min energy per photon packet to be considered
    
    string absorption = get_string(L,-1,"absorption");
    if ( absorption == "standard" )
	absorption_ = STANDARD_ABSORPTION;
    else if ( absorption == "partitioned energy" )
	absorption_ = PARTITIONED_ENERGY_ABSORPTION;
    else {
	cout << "MonteCarlo::MonteCarlo()" << endl
	     << "Absorption method: " << absorption << " not implemented" << endl
	     << "Exiting program." << endl;
	exit( NOT_IMPLEMENTED_ERROR );
    }

    string clustering = get_string(L,-1,"clustering");
    if ( clustering == "by volume" )
    	clustering_ = CLUSTERING_BY_VOLUME;
    else if ( clustering == "by area" )
    	clustering_ = CLUSTERING_BY_AREA;
    else
    	clustering_ = NO_CLUSTERING;
    
    // Done.
}

MonteCarlo::~MonteCarlo()
{
    for ( size_t ib=0; ib<cells_.size(); ++ib ) {
    	for ( size_t ic=0; ic<cells_[ib].size(); ++ic ) {
    	    delete cells_[ib][ic];
    	}
    }
    
    for ( size_t ib=0; ib<interfaces_.size(); ++ib ) {
    	for ( size_t iface=0; iface<interfaces_[ib].size(); ++iface ) {
    	    delete interfaces_[ib][iface];
    	}
    }
    
    delete cf_;
    
    delete rg_;
}

string
MonteCarlo::str() const
{
    return "MonteCarlo";
}

int
MonteCarlo::initialise()
{
    // initialise geometry data
    global_data &G = *get_global_data_ptr();
    size_t nfaces = 0;
    if ( G.dimensions == 2 ) {
    	cf_ = new CellFinder2D();	// use default nvertices
    	nfaces = 4;
    }
    else {
    	cf_ = new CellFinder3D();	// use default nvertices
    	nfaces = 6;
    }
    
    // Initialise random number generator
    int32 seed = (int32)time(0);
    rg_ = new CRandomMersenne(seed);
    
    // initialise cells and interfaces from block data
    // NOTE: not doing in parallel as this should be very quick
    cells_.resize( G.nblock );
    interfaces_.resize( G.nblock );
    size_t nthreads = omp_get_max_threads();
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	Block * bdp = G.my_blocks[jb];
#       if VERBOSE_RADIATION_TRANSPORT
    	cout << "Thread " << omp_get_thread_num() << ": Initialising cells in block: " << jb << endl;
#       endif
	// if ( bdp->active != 1 ) continue;
	for ( size_t k = bdp->kmin; k <= bdp->kmax; ++k ) {
	    for ( size_t j = bdp->jmin; j <= bdp->jmax; ++j ) {
	    	for ( size_t i = bdp->imin; i <= bdp->imax; ++i ) {
	    	    FV_Cell * cell = bdp->get_cell(i,j,k);
	    	    cells_[jb].push_back( new MonteCarloCell( cell->fs->gas, &(cell->Q_rE_rad), cell->pos, cell->volume, cell->area ) );
	    	    cells_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    cells_[jb].back()->Q_rE_rad_temp_.resize( nthreads );
	    	    cells_[jb].back()->interfaces_.resize( nfaces );
	    	    // test for and create wall interfaces
	    	    if ( j==bdp->jmin && bdp->bcp[SOUTH]->is_wall_flag ) {
	    	    	interfaces_[jb].push_back( new MonteCarloInterface( cell->iface[SOUTH]->fs->gas, cell->iface[SOUTH]->pos, cell->iface[SOUTH]->area, cell->iface[SOUTH]->length ) );
	    	    	interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    	size_t hf_index_ =  bdp->bcp[SOUTH]->get_heat_flux_index( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_ = &(bdp->bcp[SOUTH]->q_rad[hf_index_]);
	    	    	cells_[jb].back()->interfaces_[SOUTH] = interfaces_[jb].back();
	    	    }
	    	    else if ( j==bdp->jmax && bdp->bcp[NORTH]->is_wall_flag ) {
	    	    	interfaces_[jb].push_back( new MonteCarloInterface( cell->iface[NORTH]->fs->gas, cell->iface[NORTH]->pos, cell->iface[NORTH]->area, cell->iface[NORTH]->length ) );
	    	    	interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    	size_t hf_index_ =  bdp->bcp[NORTH]->get_heat_flux_index( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_ = &(bdp->bcp[NORTH]->q_rad[hf_index_]);
	    	    	cells_[jb].back()->interfaces_[NORTH] = interfaces_[jb].back();
	    	    }
	    	    if ( i==bdp->imin && bdp->bcp[WEST]->is_wall_flag ) {
	    	    	interfaces_[jb].push_back( new MonteCarloInterface( cell->iface[WEST]->fs->gas, cell->iface[WEST]->pos, cell->iface[WEST]->area, cell->iface[WEST]->length ) );
	    	    	interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    	size_t hf_index_ =  bdp->bcp[WEST]->get_heat_flux_index( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_ = &(bdp->bcp[WEST]->q_rad[hf_index_]);
	    	    	cells_[jb].back()->interfaces_[WEST] = interfaces_[jb].back();
	    	    }
	    	    else if ( i==bdp->imax && bdp->bcp[EAST]->is_wall_flag ) {
	    	    	interfaces_[jb].push_back( new MonteCarloInterface( cell->iface[EAST]->fs->gas, cell->iface[EAST]->pos, cell->iface[EAST]->area, cell->iface[EAST]->length) );
	    	    	interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    	size_t hf_index_ =  bdp->bcp[EAST]->get_heat_flux_index( i, j, k );
	    	    	interfaces_[jb].back()->q_rad_ = &(bdp->bcp[EAST]->q_rad[hf_index_]);
	    	    	cells_[jb].back()->interfaces_[EAST] = interfaces_[jb].back();
	    	    }
	    	    if ( G.dimensions == 3 ) {
	    	    	if ( k==bdp->kmin && bdp->bcp[BOTTOM]->is_wall_flag ) {
	    	    	    interfaces_[jb].push_back( new MonteCarloInterface( cell->iface[BOTTOM]->fs->gas, cell->iface[BOTTOM]->pos, cell->iface[BOTTOM]->area, cell->iface[BOTTOM]->length ) );
	    	    	    interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	    interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    	    size_t hf_index_ =  bdp->bcp[BOTTOM]->get_heat_flux_index( i, j, k );
	    	    	    interfaces_[jb].back()->q_rad_ = &(bdp->bcp[BOTTOM]->q_rad[hf_index_]);
	    	    	    cells_[jb].back()->interfaces_[BOTTOM] = interfaces_[jb].back();
	    	    	}
	    	    	else if ( k==bdp->kmax && bdp->bcp[TOP]->is_wall_flag ) {
	    	    	    interfaces_[jb].push_back( new MonteCarloInterface( cell->iface[TOP]->fs->gas, cell->iface[TOP]->pos, cell->iface[TOP]->area, cell->iface[TOP]->length ) );
	    	    	    interfaces_[jb].back()->set_CFD_cell_indices( i, j, k );
	    	    	    interfaces_[jb].back()->q_rad_temp_.resize( nthreads );
	    	    	    size_t hf_index_ =  bdp->bcp[TOP]->get_heat_flux_index( i, j, k );
	    	    	    interfaces_[jb].back()->q_rad_ = &(bdp->bcp[TOP]->q_rad[hf_index_]);
	    	    	    cells_[jb].back()->interfaces_[TOP] = interfaces_[jb].back();
	    	    	}
	    	    }
	    	}
	    }
	}
	
#       if VERBOSE_RADIATION_TRANSPORT
	cout << "- created " << cells_[jb].size() << " MonteCarloCell's" << endl;
	cout << "- created " << interfaces_[jb].size() << " MonteCarloInterface's" << endl;
#       endif
    }
    
    return SUCCESS;
}

void MonteCarlo::compute_Q_rad_for_flowfield()
{
    // Initialise geometry data
    global_data &G = *get_global_data_ptr();
    bool planar_flag = false;
    if ( G.dimensions == 2 && !get_axisymmetric_flag() ) planar_flag = true;
    
    // 0. Loop over spectral blocks (CHECKME)
    for ( size_t isb=0; isb < static_cast<size_t>(rsm_[0]->get_spectral_blocks()); ++isb ) {
	for ( size_t irsm=0; irsm<rsm_.size(); ++irsm ) {
	    rsm_[irsm]->reset_spectral_params( isb );
	}

	// 1. Set all source terms to zero, store all spectra
	for ( size_t ib=0; ib<cells_.size(); ++ib ) {
	    size_t ic;
#	    ifdef _OPENMP
#	    pragma omp parallel for private(ic) schedule(runtime)
#	    endif
	    for ( ic=0; ic<cells_[ib].size(); ++ic ) {
	    	RayTracingCell * cell = cells_[ib][ic];
		*(cell->Q_rE_rad_) = 0.0;
		// Also make sure thread vector is zero
		for ( size_t iQ=0; iQ<cell->Q_rE_rad_temp_.size(); ++ iQ ) {
		    cell->Q_rE_rad_temp_[iQ] = 0.0;
		}
		cout << "Thread " << omp_get_thread_num() << ": Recomputing spectra for cell: " << ic << " in block: " << ib;
		cell->recompute_spectra( rsm_[omp_get_thread_num()] );
		if ( rsm_[omp_get_thread_num()]->get_spectral_points()==1 )
		    cout << " - j_total = " << cell->X_->j_int[0] << endl;
		else
		    cout << " - j_total = " << cell->X_->integrate_emission_spectra() << endl;
	    }
	    size_t iface;
#	    ifdef _OPENMP
#	    pragma omp parallel for private(iface) schedule(runtime)
#	    endif
	    for ( iface=0; iface<interfaces_[ib].size(); ++iface ) {
	    	RayTracingInterface * interface = interfaces_[ib][iface];
		cout << "Thread " << omp_get_thread_num() << ": Recomputing spectra for interface: " << iface << " in block: " << ib;
		interface->recompute_spectra( rsm_[omp_get_thread_num()] );
		if ( rsm_[omp_get_thread_num()]->get_spectral_points()==1 )
		    cout << " - I_total = " << interface->S_->I_int[0] << endl;
		else
		    cout << " - I_total = " << interface->S_->integrate_intensity_spectra() << endl;
	    }
	}
	
	// 2. Compute cell_E_rad_total_max and interface_E_rad_total_max
	// NOTE: easiest to not do this in parallel
	double cell_E_rad_total_max = 0.0;
	double interface_E_rad_total_max = 0.0;
	for ( size_t ib=0; ib<cells_.size(); ++ib ) {
	    for ( size_t ic=0; ic<cells_[ib].size(); ++ic ) {
	    	RayTracingCell * cell = cells_[ib][ic];
		if ( clustering_==CLUSTERING_BY_VOLUME )
		    cell->E_rad_ = cell->X_->j_int.back() * cell->vol_;
		else if ( clustering_==CLUSTERING_BY_AREA )
		    cell->E_rad_ = cell->X_->j_int.back() * cell->area_;
		if ( cell->E_rad_ > cell_E_rad_total_max ) cell_E_rad_total_max = cell->E_rad_;
	    }
	    for ( size_t iface=0; iface<interfaces_[ib].size(); ++iface ) {
	    	RayTracingInterface * interface = interfaces_[ib][iface];
		if ( clustering_==CLUSTERING_BY_VOLUME )
		    interface->E_rad_ = interface->S_->I_int.back() * interface->area_;
		else if ( clustering_==CLUSTERING_BY_AREA )
		    interface->E_rad_ = interface->S_->I_int.back() * interface->length_;
		if ( interface->E_rad_ > interface_E_rad_total_max ) interface_E_rad_total_max = interface->E_rad_;
	    }
	}
	
	// 3. create ray one-by-one and transport radiation along it
	for ( size_t ib=0; ib<cells_.size(); ++ib ) {
	    size_t ii, jj, kk;
	    for ( size_t ic=0; ic<cells_[ib].size(); ++ic ) {
#           	if VERBOSE_RADIATION_TRANSPORT
    	    	cout << "Thread " << omp_get_thread_num() << ": Tracing rays for cell: " << ic << " in block: " << ib;
#           	endif
		RayTracingCell * cell = cells_[ib][ic];
		cell->get_CFD_cell_indices( ii, jj, kk );
		size_t nrays = nrays_;
		if ( clustering_ ) {
		    nrays = int( double(nrays_) * cell->E_rad_ / cell_E_rad_total_max );
		}
		if ( nrays==0 ) nrays = 1;
		cout << " - nrays: " << nrays << endl;
		// Cycle over all rays
		size_t iray;
#   	   	ifdef _OPENMP
#          	pragma omp parallel for private(iray) schedule(runtime)
#          	endif
		for ( iray=0; iray<nrays; ++iray ) {
		    RayTracingRay * ray = this->create_new_ray_for_cell( cell, nrays, G.dimensions, planar_flag );
		    double nu = cell->X_->random_frequency( rg_->Random() );
		    // Each photon for a given cell has the same energy (domega is constant)
		    double E = cell->X_->j_int.back() * cell->vol_ * ray->domega_;
		    if ( E < E_min_ ) continue;
		    // subtract emitted energy from origin cell
		    cell->Q_rE_rad_temp_[omp_get_thread_num()] -= E / cell->vol_;
		    if ( absorption_ == STANDARD_ABSORPTION )
		        this->trace_ray_standard( ray, ib, ii, jj, kk, nu, E );
		    else
			this->trace_ray_partitioned_energy( ray, ib, ii, jj, kk, nu, E );
		    delete ray;
		}
	    }
	    for ( size_t iface=0; iface<interfaces_[ib].size(); ++iface ) {
#           	if VERBOSE_RADIATION_TRANSPORT
    	    	cout << "Thread " << omp_get_thread_num() << ": Tracing rays for interface: " << iface << " in block: " << ib;
#           	endif
	    	RayTracingInterface * interface = interfaces_[ib][iface];
	    	interface->get_CFD_cell_indices( ii, jj, kk );
	        size_t nrays = nrays_;
	        if ( clustering_ ) {
		    nrays = int( double(nrays_) * interface->E_rad_ / interface_E_rad_total_max );

	        }
		if ( nrays<=0 ) nrays = 1;
		cout << " - nrays: " << nrays << endl;
		// Cycle over all rays
		size_t iray;
#   	   	ifdef _OPENMP
#          	pragma omp parallel for private(iray) schedule(runtime)
#          	endif
		for ( iray=0; iray<nrays; ++iray ) {
		    RayTracingRay * ray = this->create_new_ray_for_interface( interface, nrays, G.dimensions, planar_flag );
                    double nu = interface->S_->random_frequency( rg_->Random() );
		    double E = interface->S_->I_int.back() * interface->area_ * ray->domega_;
		    if ( E < E_min_ ) continue;
		    if ( absorption_ == STANDARD_ABSORPTION )
		        this->trace_ray_standard( ray, ib, ii, jj, kk, nu, E );
		    else
			this->trace_ray_partitioned_energy( ray, ib, ii, jj, kk, nu, E );
		    delete ray;
		}
	    }
	}
    }
	
    // 5. Sum Q_rE_rad_temp values and set Q_rE_rad in CFD cells
    for ( size_t ib=0; ib<cells_.size(); ++ib ) {
    	size_t ic;
#   	ifdef _OPENMP
#   	pragma omp barrier
#   	pragma omp parallel for private(ic) schedule(runtime)
#  	endif
    	for ( ic=0; ic<cells_[ib].size(); ++ic ) {
    	    RayTracingCell * cell = cells_[ib][ic];
    	    for ( size_t iQ=0; iQ<cell->Q_rE_rad_temp_.size(); ++iQ ) {
    	    	*(cell->Q_rE_rad_) += cell->Q_rE_rad_temp_[iQ];
    	    }
    	    // cout << "*(cell->Q_rE_rad_) [ ib = " << ib << ", ic = " << ic << " ] = " << *(cell->Q_rE_rad_) << endl;
    	}
    }
    
    // 6. Sum q_rad_temp values and set q_rad in CFD interfaces
    for ( size_t ib=0; ib<cells_.size(); ++ib ) {
    	size_t iface;
#   	ifdef _OPENMP
#   	pragma omp barrier
#   	pragma omp parallel for private(iface) schedule(runtime)
#  	endif
    	for ( iface=0; iface<interfaces_[ib].size(); ++iface ) {
    	    RayTracingInterface * interface = interfaces_[ib][iface];
    	    for ( size_t iq=0; iq<interface->q_rad_temp_.size(); ++iq ) {
    	    	*(interface->q_rad_) += interface->q_rad_temp_[iq];
    	    }
    	}
    }

    // clean up
    // all done in deconstructor for now
    
    
    return;
}

RayTracingRay * MonteCarlo::create_new_ray_for_cell( RayTracingCell * cell, size_t nrays, size_t ndim, bool planar )
{
    // Initialise uniformly random rays
    double domega = 4.0 * M_PI / double( nrays );
    double theta = 2.0 * M_PI * rg_->Random();
    double phi = acos( 1.0 - 2.0 * rg_->Random() ) - M_PI / 2.0;
    // theta = 0 or 180 degrees for planar geometry
    if ( planar ) theta = ( theta < M_PI ) ? 0.0 : M_PI;
    // Randomize origin location (for the moment just use the origin)
#   if 0
    Vector3 tmp;
#   else
    Vector3 origin( cell->origin_ );
#   endif
    RayTracingRay * ray = 0;
    if ( ndim==2 ) {
	ray = new RayTracingRay2D( theta, phi, domega, origin );
    }
    else {
	ray = new RayTracingRay3D( theta, phi, domega, origin );
    }
    
    return ray;
}

RayTracingRay * MonteCarlo::create_new_ray_for_interface( RayTracingInterface * interface, size_t nrays, size_t ndim, bool planar )
{
    // Initialise uniformly random rays
    double domega = 4.0 * M_PI / double( nrays );
    double theta = 2.0 * M_PI * rg_->Random();
    double phi = acos( 1.0 - 2.0 * rg_->Random() ) - M_PI / 2.0;
    // theta = 0 or 180 degrees for planar geometry
    if ( planar ) theta = ( theta < M_PI ) ? 0.0 : M_PI;
    // Randomize origin location (for the moment just use the origin)
#   if 0
    Vector3 tmp;
#   else
    Vector3 origin( interface->origin_ );
#   endif
    RayTracingRay * ray = 0;
    if ( ndim==2 ) {
	ray = new RayTracingRay2D( theta, phi, domega, origin );
    }
    else {
	ray = new RayTracingRay3D( theta, phi, domega, origin );
    }
    
    return ray;
}


int MonteCarlo::trace_ray_standard( RayTracingRay * ray, size_t ib, size_t ic, size_t jc, size_t kc, double nu, double E )
{
    // 1. Initialise pointers
    Block * A = get_block_data_ptr( ib );
    FV_Cell * cell = A->get_cell(ic,jc,kc);		// start at origin cell
    
    double dL = 0.5 * dl_lmin_ratio_ * cell->L_min;	// small initial step length
    double L = dL;
    Vector3 p = ray->get_point_on_line( L );
    
    RayTracingCell * RTcell = 0;
    
    size_t count = 0;

    // limiting optical thickness
    double tau_lim = 1.0 / exp( rg_->Random() );

    // accumulative optical thickness
    double tau_acc = 0.0;

    // step along the ray, dumping energy as we go
    while ( ( ray->status_ = cf_->find_cell( p, ib, ic, jc, kc ) ) == INSIDE_GRID ) {
    	// Get pointers to the new block and cell
    	A = get_block_data_ptr( ib );
    	cell = A->get_cell(ic,jc,kc);
    	RTcell = cells_[ib][get_cell_index(A,ic,jc,kc)];
	// calculate optical thickness
	double kappa = RTcell->X_->kappa_from_nu(nu);
	tau_acc += kappa * dL;
	// add all energy to traversed cell if limiting optical thickness exceeded
	if ( tau_acc >= tau_lim ) {
	    RTcell->Q_rE_rad_temp_[omp_get_thread_num()] += E / RTcell->vol_;
	    return SUCCESS;
	}
    	// calculate next position on ray
    	dL = dl_lmin_ratio_ * cell->L_min;
    	L += dL;
    	++count;
	p = ray->get_point_on_line( L );
    }

    // Make sure block and cell pointers are up to date
    A = get_block_data_ptr( ib );
    cell = A->get_cell(ic,jc,kc);
    
    if ( ray->status_ == ERROR ) {
    	cout << "MonteCarlo::trace_ray()" << endl
    	     << "CellFinder failed at location: " << p.str() << endl
    	     << "Exiting program" << endl;
    	exit( FAILURE);
    }
    else if ( count>0 && A->bcp[ray->status_]->is_wall_flag ) {
    	// dump energy onto interface (only if more than one point and this is a wall)
    	// size_t index = A->bcp[ ray->status_ ]->get_heat_flux_index( ic, jc, kc );
	// FIXME: need a way of finding the nearest interface to an exited photons
	//        location when the photon exits on a corner
    	RayTracingInterface * RTinterface =  RTcell->interfaces_[ray->status_];
    	bool skip_interface = false;
    	if ( RTinterface==0 ) {
    	    cout << "MonteCarlo::trace_ray" << endl
    	         << "bad RayTracingInterface pointer on RayTracingCell" << endl;
    	    skip_interface = true;
#           if EXIT_ON_RT_FAILURE
    	    cout << "Bailing out!" << endl;
    	    cout << "ray->status_ = " << ray->status_ << endl;
    	    cout << "RTcell->interfaces_[NORTH] = " << RTcell->interfaces_[NORTH] << endl;
    	    cout << "RTcell->interfaces_[SOUTH] = " << RTcell->interfaces_[SOUTH] << endl;
    	    cout << "RTcell->interfaces_[EAST] = " << RTcell->interfaces_[EAST] << endl;
    	    cout << "RTcell->interfaces_[WEST] = " << RTcell->interfaces_[WEST] << endl;
    	    cout << "cell->pos = " << RTcell->origin_.str() << endl;
    	    cout << "point = " << p.str() << endl;
    	    exit( BAD_INPUT_ERROR );
#           endif
    	}
    	if ( ! skip_interface )
    	    RTinterface->q_rad_temp_[omp_get_thread_num()] += E / RTinterface->area_;
    }
    
    return SUCCESS;
}

int MonteCarlo::trace_ray_partitioned_energy( RayTracingRay * ray, size_t ib, size_t ic, size_t jc, size_t kc, double nu, double E )
{
    // 1. Initialise pointers
    Block * A = get_block_data_ptr( ib );
    FV_Cell * cell = A->get_cell(ic,jc,kc);             // start at origin cell

    double dL = 0.5 * dl_lmin_ratio_ * cell->L_min;     // small initial step length
    double L = dL;
    Vector3 p = ray->get_point_on_line( L );

    RayTracingCell * RTcell = 0;

    size_t count = 0;

    // step along the ray, dumping energy as we go
    while ( ( ray->status_ = cf_->find_cell( p, ib, ic, jc, kc ) ) == INSIDE_GRID ) {
        // Get pointers to the new block and cell
        A = get_block_data_ptr( ib );
        cell = A->get_cell(ic,jc,kc);
        RTcell = cells_[ib][get_cell_index(A,ic,jc,kc)];
        // calculate absorbed energy
        double kappa = RTcell->X_->kappa_from_nu(nu);
        double dE = ( 1.0 - exp( - kappa * dL ) ) * E; E -= dE;
        // add absorbed energy to traversed cell
        RTcell->Q_rE_rad_temp_[omp_get_thread_num()] += dE / RTcell->vol_;
        // calculate next position on ray
        dL = dl_lmin_ratio_ * cell->L_min;
        L += dL;
        ++count;
        p = ray->get_point_on_line( L );
    }

    // Make sure block and cell pointers are up to date
    A = get_block_data_ptr( ib );
    cell = A->get_cell(ic,jc,kc);

    if ( ray->status_ == ERROR ) {
        cout << "MonteCarlo::trace_ray()" << endl
             << "CellFinder failed at location: " << p.str() << endl
             << "Exiting program" << endl;
        exit( FAILURE);
    }
    else if ( count>0 && A->bcp[ray->status_]->is_wall_flag ) {
        // dump energy onto interface (only if more than one point and this is a wall)
        // size_t index = A->bcp[ ray->status_ ]->get_heat_flux_index( ic, jc, kc );
	// FIXME: need a way of finding the nearest interface to an exited photons
	//        location when the photon exits on a corner
        RayTracingInterface * RTinterface =  RTcell->interfaces_[ray->status_];
        bool skip_interface = false;
        if ( RTinterface==0 ) {
            cout << "MonteCarlo::trace_ray" << endl
                 << "bad RayTracingInterface pointer on RayTracingCell" << endl;
            skip_interface = true;
#           if EXIT_ON_RT_FAILURE
            cout << "Bailing out!" << endl;
            cout << "ray->status_ = " << ray->status_ << endl;
            cout << "RTcell->interfaces_[NORTH] = " << RTcell->interfaces_[NORTH] << endl;
            cout << "RTcell->interfaces_[SOUTH] = " << RTcell->interfaces_[SOUTH] << endl;
            cout << "RTcell->interfaces_[EAST] = " << RTcell->interfaces_[EAST] << endl;
            cout << "RTcell->interfaces_[WEST] = " << RTcell->interfaces_[WEST] << endl;
            cout << "cell->pos = " << RTcell->origin_.str() << endl;
            cout << "point = " << p.str() << endl;
            exit( BAD_INPUT_ERROR );
#           endif
        }
        if ( ! skip_interface )
            RTinterface->q_rad_temp_[omp_get_thread_num()] += E / RTinterface->area_;
    }

    return SUCCESS;
}

RadiationTransportModel * create_radiation_transport_model( const string file_name )
{
    // 0. Initialise a Radiation_transport_model pointer
    RadiationTransportModel * rtm;
    
    // 1. Get transport model name from lua file
    lua_State *L = initialise_radiation_lua_State();

    if( luaL_dofile(L, file_name.c_str()) != 0 ) {
	ostringstream ost;
	ost << "set_radiation_transport_model():\n";
	ost << "Error in input file: " << file_name << endl;
	input_error(ost);
    }
    
    lua_getglobal(L,"transport_data");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "set_radiation_transport_model():\n";
	ost << "Error in the declaration of transport_data - a table is expected.\n";
	input_error(ost);
    }
    
    string transport_model = get_string(L, -1, "transport_model");
    
    // 2. Create the transport model
    if ( ECHO_RAD_INPUT > 0 ) 
	cout << "--- Creating a new " << transport_model << " transport model"
   	     << " from file: " << file_name << " ---" << endl;
        
    if ( transport_model == "none" ) {
	rtm = new NoRadiationTransportModel();
    }
    else if( transport_model == "optically thin" ) {
	rtm = new OpticallyThin(L);
    }
    else if( transport_model == "tangent slab" ) {
	rtm = new TangentSlab(L);
    }
    else if( transport_model == "discrete transfer" ) {
	rtm = new DiscreteTransfer(L);
    }
    else if( transport_model == "monte carlo" ) {
	rtm = new MonteCarlo(L);
    }
    else {
	cout << "The specified radiation transport model: " << transport_model << endl;
	cout << "is not available of no yet implemented.\n";
	cout << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    lua_pop(L,1);	// pop transport_data
    lua_close(L);
    
    // 3. Create the spectral model(s)
    rtm->set_radiation_spectral_model( file_name, omp_get_max_threads() );
    
    return rtm;
}


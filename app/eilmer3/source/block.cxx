/// \file block.cxx
/// \ingroup eilmer3
/// \brief Functions that apply to the block as a whole.
///
/// \author PJ
/// \version 23-Jun-2006 Abstracted from cns_tstp.cxx
///

#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <unistd.h>
#include <stdexcept>
#include <math.h>
extern "C" {
#include <zlib.h>
}
#include "cell.hh"
#include "kernel.hh"
#include "block.hh"
#include "bc.hh"

//-----------------------------------------------------------------------------

Block::Block()
    : bcp(N_INTERFACE,NULL) // Let everything else default initialize.
{}

Block::Block(const Block &b)
    : id(b.id), active(b.active), omegaz(b.omegaz),
      dt_allow(b.dt_allow),
      cfl_min(b.cfl_min), cfl_max(b.cfl_max),
      mass_residual(b.mass_residual),
      energy_residual(b.energy_residual),
      mass_residual_loc(b.mass_residual_loc),
      energy_residual_loc(b.energy_residual_loc),
      hncell(b.hncell),
      hicell(b.hicell), hjcell(b.hjcell), hkcell(b.hkcell),
      nidim(b.nidim), njdim(b.njdim), nkdim(b.nkdim),
      nni(b.nni), nnj(b.nnj), nnk(b.nnk),
      imin(b.imin), imax(b.imax),
      jmin(b.jmin), jmax(b.jmax),
      kmin(b.kmin), kmax(b.kmax),
      active_cells(b.active_cells),
      baldwin_lomax_iturb(b.baldwin_lomax_iturb),
      bcp(b.bcp),
      ctr_(b.ctr_),
      ifi_(b.ifi_), ifj_(b.ifj_), ifk_(b.ifk_),
      vtx_(b.vtx_),
      sifi_(b.sifi_), sifj_(b.sifj_), sifk_(b.sifk_)
{}

Block & Block::operator=(const Block &b)
{
    if ( this != &b ) { // Avoid aliasing
	id = b.id; active = b.active; omegaz = b.omegaz;
	dt_allow = b.dt_allow;
	cfl_min = b.cfl_min; cfl_max = b.cfl_max;
	mass_residual = b.mass_residual;
	energy_residual = b.energy_residual;
	mass_residual_loc = b.mass_residual_loc;
	energy_residual_loc = b.energy_residual_loc;
	hncell = b.hncell;
	hicell = b.hicell; hjcell = b.hjcell; hkcell = b.hkcell;
	nidim = b.nidim; njdim = b.njdim; nkdim = b.nkdim;
	nni = b.nni; nnj = b.nnj; nnk = b.nnk;
	imin = b.imin; imax = b.imax;
	jmin = b.jmin; jmax = b.jmax;
	kmin = b.kmin; kmax = b.kmax;
	active_cells = b.active_cells;
	baldwin_lomax_iturb = b.baldwin_lomax_iturb;
	bcp = b.bcp;
	ctr_ = b.ctr_;
	ifi_ = b.ifi_; ifj_ = b.ifj_; ifk_ = b.ifk_;
	vtx_ = b.vtx_;
	sifi_ = b.sifi_; sifj_ = b.sifj_; sifk_ = b.sifk_;
    }
    return *this;
}

Block::~Block() 
{
    for ( size_t i = 0; i < bcp.size(); ++i )
	delete bcp[i];
    // for ( size_t i = 0; i < shock_iface_pos.size(); ++i )
    //	delete shock_iface_pos[i];
}

/// \brief Allocate memory for the internal arrays of the block.
/// Returns 0 if successful, 1 otherwise.
int Block::array_alloc(size_t dimensions)
{
    global_data &G = *get_global_data_ptr();
    if ( G.verbose_init_messages ) cout << "array_alloc(): Begin for block " <<  id << endl;
    // Check for obvious errors.
    if ( nidim <= 0 || njdim <= 0 || nkdim <= 0 ) {
        cerr << "array_alloc(): Declared dimensions are zero or negative: " 
	     << nidim << ", " << njdim << ", " << nkdim << endl;
	exit( VALUE_ERROR );
    }
    Gas_model *gm = get_gas_model_ptr();

    // Allocate vectors of pointers for the entire block.
    size_t ntot = nidim * njdim * nkdim;
    ctr_.resize(ntot); 
    ifi_.resize(ntot);
    ifj_.resize(ntot);
    if ( dimensions == 3 ) ifk_.resize(ntot);
    vtx_.resize(ntot);
    sifi_.resize(ntot);
    sifj_.resize(ntot);
    if ( dimensions == 3 ) sifk_.resize(ntot);
    // Now, create the actual objects.
    for (size_t gid = 0; gid < ntot; ++gid) {
        ctr_[gid] = new FV_Cell(gm);
	ctr_[gid]->id = gid;
	std::vector<size_t> ijk = to_ijk_indices(gid);
	size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2];
	if ( i >= imin && i <= imax && j >= jmin && j <= jmax && k >= kmin && k <= kmax ) {
	    active_cells.push_back(ctr_[gid]);
	}
        ifi_[gid] = new FV_Interface(gm);
	ifi_[gid]->id = gid;
        ifj_[gid] = new FV_Interface(gm);
	ifj_[gid]->id = gid;
        if ( dimensions == 3 ) {
	    ifk_[gid] = new FV_Interface(gm);
	    ifk_[gid]->id = gid;
	}
        vtx_[gid] = new FV_Vertex(gm);
	vtx_[gid]->id = gid;
        sifi_[gid] = new FV_Interface(gm);
	sifi_[gid]->id = gid;
        sifj_[gid] = new FV_Interface(gm);
	sifj_[gid]->id = gid;
        if ( dimensions == 3 ) {
	    sifk_[gid] = new FV_Interface(gm);
	    sifk_[gid]->id = gid;
	}
    } // gid loop

    if ( G.verbose_init_messages || id == 0 ) {
	cout << "Block " << id << ": finished creating " << ntot << " cells." << endl;
    }
    return SUCCESS;
} // end of array_alloc()


int Block::array_cleanup(size_t dimensions)
{
    global_data &G = *get_global_data_ptr();
    if ( G.verbose_init_messages || id == 0 ) {
	cout << "array_cleanup(): Begin for block " <<  id << endl;
    }
    // Need to clean up allocated memory.
    size_t ntot = nidim * njdim * nkdim;
    for (size_t gid = 0; gid < ntot; ++gid) {
	delete ctr_[gid];
	delete ifi_[gid];
	delete ifj_[gid];
	if ( dimensions == 3 ) delete ifk_[gid];
	delete vtx_[gid];
	delete sifi_[gid];
	delete sifj_[gid];
	if ( dimensions == 3 ) delete sifk_[gid];
    } // gid loop
    ctr_.clear(); 
    ifi_.clear();
    ifj_.clear();
    ifk_.clear();
    vtx_.clear();
    sifi_.clear();
    sifj_.clear();
    sifk_.clear();
    active_cells.clear();
    return SUCCESS;
} // end of array_cleanup()

//-----------------------------------------------------------------------------

int Block::bind_interfaces_to_cells( size_t dimensions )
{
    FV_Cell *cellp;
    size_t kstart, kend;

    if ( dimensions == 3 ) {
	kstart = kmin-1;
	kend = kmax+1;
    } else {
	kstart = 0;
	kend = 0;
    }
    for ( size_t k = kstart; k <= kend; ++k ) {
	for ( size_t j = jmin-1; j <= jmax+1; ++j ) {
	    for ( size_t i = imin-1; i <= imax+1; ++i ) {
		cellp = get_cell(i,j,k);
		cellp->iface[NORTH] = get_ifj(i,j+1,k);
		cellp->iface[EAST] = get_ifi(i+1,j,k);
		cellp->iface[SOUTH] = get_ifj(i,j,k);
		cellp->iface[WEST] = get_ifi(i,j,k);
		cellp->vtx[0] = get_vtx(i,j,k);
		cellp->vtx[1] = get_vtx(i+1,j,k);
		cellp->vtx[2] = get_vtx(i+1,j+1,k);
		cellp->vtx[3] = get_vtx(i,j+1,k);
		if ( dimensions == 3 ) {
		    cellp->iface[TOP] = get_ifk(i,j,k+1);
		    cellp->iface[BOTTOM] = get_ifk(i,j,k);
		    cellp->vtx[4] = get_vtx(i,j,k+1);
		    cellp->vtx[5] = get_vtx(i+1,j,k+1);
		    cellp->vtx[6] = get_vtx(i+1,j+1,k+1);
		    cellp->vtx[7] = get_vtx(i,j+1,k+1);
		} else {
		    cellp->iface[TOP] = NULL;
		    cellp->iface[BOTTOM] = NULL;
		    cellp->vtx[4] = NULL;
		    cellp->vtx[5] = NULL;
		    cellp->vtx[6] = NULL;
		    cellp->vtx[7] = NULL;
		} // end if
	    }
	}
    }
    return SUCCESS;
} // end bind_interfaces_to_cells()


/// \brief Set the base heat source values for this block.
int Block::set_base_qdot(global_data &gd, size_t gtl)
{
    double total_qdot_for_block = 0.0;
    total_qdot_for_block = 0.0;
    for ( FV_Cell *cellp: active_cells ) {
	cellp->base_qdot = 0.0;
	for ( CHeatZone &hz : gd.heat_zone ) {
	    if ( cellp->pos[gtl].x >= hz.x0 && cellp->pos[gtl].x <= hz.x1 &&
		 cellp->pos[gtl].y >= hz.y0 && cellp->pos[gtl].y <= hz.y1 &&
		 (gd.dimensions == 2 || (cellp->pos[gtl].z >= hz.z0 && cellp->pos[gtl].z <= hz.z1)) ) {
		cellp->base_qdot += hz.qdot;
	    }
	} // for hz
	total_qdot_for_block += cellp->base_qdot * cellp->volume[gtl];
    } // for cellp
    if ( total_qdot_for_block > 1.0e-10 ) {
	cout << "set_base_qdot(): block " << id
	     << " base_qdot= " << total_qdot_for_block << " Watts" << endl;
    }
    return SUCCESS;
} // end set_base_qdot()


/// \brief Set the reactions-allowed flag for cells in this block.
int Block::identify_reaction_zones(global_data &gd, size_t gtl)
{
    size_t total_cells_in_reaction_zones = 0;
    size_t total_cells = 0;
    for ( FV_Cell *cellp: active_cells ) {
	if ( gd.n_reaction_zone > 0 ) {
	    // User-specified reaction zones; mask off reacting/nonreacting zones.
	    cellp->fr_reactions_allowed = false;
	    for ( CReactionZone &rz : gd.reaction_zone ) {
		if ( cellp->pos[gtl].x >= rz.x0 && cellp->pos[gtl].x <= rz.x1 &&
		     cellp->pos[gtl].y >= rz.y0 && cellp->pos[gtl].y <= rz.y1 &&
		     (gd.dimensions == 2 || 
		      (cellp->pos[gtl].z >= rz.z0 && cellp->pos[gtl].z <= rz.z1)) ) {
		    cellp->fr_reactions_allowed = true;
		}
	    } // end for( &rz
	} else {
	    // No user-specified zones; always allow reactions.
	    cellp->fr_reactions_allowed = true;
	}
	total_cells_in_reaction_zones += (cellp->fr_reactions_allowed ? 1: 0);
	total_cells += 1;
    } // for cellp
    if ( gd.reacting ) {
	cout << "identify_reaction_zones(): block " << id
	     << " cells inside zones = " << total_cells_in_reaction_zones 
	     << " out of " << total_cells << endl;
	if ( gd.n_reaction_zone == 0 ) {
	    cout << "Note that for no user-specified zones,"
		 << " the whole domain is allowed to be reacting." << endl;
	}
    }
    return SUCCESS;
} // end identify_reaction_zones()


/// \brief Set the in-turbulent-zone flag for cells in this block.
int Block::identify_turbulent_zones(global_data &gd, size_t gtl)
{
    size_t total_cells_in_turbulent_zones = 0;
    size_t total_cells = 0;
    for ( FV_Cell *cellp: active_cells ) {
	if ( gd.n_turbulent_zone > 0 ) {
	    cellp->in_turbulent_zone = false;
	    for ( CTurbulentZone &tz : gd.turbulent_zone ) {
		if ( cellp->pos[gtl].x >= tz.x0 && cellp->pos[gtl].x <= tz.x1 &&
		     cellp->pos[gtl].y >= tz.y0 && cellp->pos[gtl].y <= tz.y1 &&
		     (gd.dimensions == 2 || 
		      (cellp->pos[gtl].z >= tz.z0 && cellp->pos[gtl].z <= tz.z1)) ) {
		    cellp->in_turbulent_zone = true;
		}
	    } // for tz
	} else {
	    cellp->in_turbulent_zone = true;
	}
	total_cells_in_turbulent_zones += (cellp->in_turbulent_zone ? 1: 0);
	total_cells += 1;
    } // for cellp
    if ( gd.turbulence_model != TM_NONE ) {
	cout << "identify_turbulent_zones(): block " << id
	     << " cells inside zones = " << total_cells_in_turbulent_zones 
	     << " out of " << total_cells << endl;
	if ( gd.n_turbulent_zone == 0 ) {
	    cout << "Note that for no user-specified zones,"
		 << " the whole domain is allowed to be turbulent." << endl;
	}
    }
    return SUCCESS;
} // end identify_turbulent_zones()


int Block::clear_fluxes_of_conserved_quantities(size_t dimensions)
{
    FV_Interface *IFace;

    for ( size_t k = kmin; k <= kmax; ++k ) {
	for (size_t j = jmin; j <= jmax; ++j) {
	    for (size_t i = imin; i <= imax+1; ++i) {
		IFace = get_ifi(i,j,k);
		IFace->F->clear_values();
	    } // for i
	} // for j
    } // for k
    for ( size_t k = kmin; k <= kmax; ++k ) {
	for (size_t j = jmin; j <= jmax+1; ++j) {
	    for (size_t i = imin; i <= imax; ++i) {
		IFace = get_ifj(i,j,k);
		IFace->F->clear_values();
	    } // for i
	} // for j
    } // for k
    if ( dimensions == 3 ) {
	for ( size_t k = kmin; k <= kmax+1; ++k ) {
	    for (size_t j = jmin; j <= jmax; ++j) {
		for (size_t i = imin; i <= imax; ++i) {
		    IFace = get_ifk(i,j,k);
		    IFace->F->clear_values();
		} // for i
	    } // for j
	} // for k
    } // end if G.dimensions == 3
    return SUCCESS;
}

int Block::propagate_data_west_to_east(size_t dimensions)
// Propagate data from the west ghost cell, right across the block.
// This is a useful starting state for the block-sequenced calculation
// where the final flow is expected to be steady-state.
{
    FV_Cell *src, *dest;
    Gas_model *gm = get_gas_model_ptr();
    for ( size_t k = kmin; k <= kmax; ++k ) {
	for ( size_t j = jmin; j <= jmax; ++j) {
	    src = get_cell(imin-1,j);
	    for ( size_t i = imin; i <= imax; ++i ) {
		dest = get_cell(i,j,k);
		dest->copy_values_from(*src, COPY_FLOW_STATE, 0);
		Gas_data *gas = dest->fs->gas;
		if ( gm->eval_thermo_state_pT(*gas) != SUCCESS ||
		     gm->eval_transport_coefficients(*gas) != SUCCESS ) {
		    printf( "propagate_data_west_to_east(): Duff call to thermo model.\n" );
		    printf( "   i=%d, j=%d, k=%d\n", static_cast<int>(i),
			    static_cast<int>(j), static_cast<int>(k) );
		    gas->print_values();
		    throw std::runtime_error("Block::propagate_data_west_to_east(): "
					     "Duff call to thermo model.");
		}
	    } // for i
	} // for j
    } // for k
    return SUCCESS;
} // end propagate_data_west_to_east()


int Block::count_invalid_cells(size_t dimensions, size_t gtl)
/// \brief Returns the number of cells that contain invalid data.
///
/// This data can be identified by the density of internal energy 
/// being on the minimum limit or the velocity being very large.
//
// To do: We should probably make this function more 3D friendly, however,
// it should not be invoked (ever) if the code is working well!
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    size_t number_of_invalid_cells = 0;
    for ( FV_Cell *cellp: active_cells ) {
	if ( cellp->check_flow_data() == false ) {
	    ++number_of_invalid_cells;
	    std::vector<size_t> ijk = to_ijk_indices(cellp->id);
	    size_t i = ijk[0]; size_t j = ijk[1]; size_t k = ijk[2];
	    printf("count_invalid_cells: block_id = %d, cell[%d,%d,%d]\n", 
		   static_cast<int>(id), static_cast<int>(i), 
		   static_cast<int>(j), static_cast<int>(k));
	    cellp->print();
	    if ( adjust_invalid_cell_data ) {
		// We shall set the cell data to something that
		// is valid (and self consistent).
		FV_Cell *other_cellp;
		std::vector<FV_Cell *> neighbours;
		if ( prefer_copy_from_left ) {
		    printf( "Adjusting cell data by copying data from left.\n" );
		    // Presently, only searches around in the i,j plane
		    other_cellp = get_cell(i-1,j,k);
		    if ( other_cellp->check_flow_data() ) neighbours.push_back(other_cellp);
		} else {
		    printf( "Adjusting cell data to a local average.\n" );
		    other_cellp = get_cell(i+1,j,k);
		    if ( other_cellp->check_flow_data() ) neighbours.push_back(other_cellp);
		    other_cellp = get_cell(i,j-1,k);
		    if ( other_cellp->check_flow_data() ) neighbours.push_back(other_cellp);
		    other_cellp = get_cell(i,j+1,k);
		    if ( other_cellp->check_flow_data() ) neighbours.push_back(other_cellp);
		}
		if ( neighbours.size() == 0 ) {
		    throw std::runtime_error("Block::count_invalid_cells(): "
					     "There were no valid neighbours to replace cell data.\n");
		}
		cellp->replace_flow_data_with_average(neighbours);
		cellp->encode_conserved(gtl, 0, omegaz, with_k_omega);
		cellp->decode_conserved(gtl, 0, omegaz, with_k_omega);
		printf("after flow-data replacement: block_id = %d, cell[%d,%d,%d]\n", 
		       static_cast<int>(id), static_cast<int>(i),
		       static_cast<int>(j), static_cast<int>(k));
		cellp->print();
	    } // end adjust_invalid_cell_data 
	} // end of if invalid data...
    } // for cellp
    return number_of_invalid_cells;
} // end of count_invalid_cells()


int Block::init_residuals(size_t dimensions)
/// \brief Initialization of data for later computing residuals.
{
    mass_residual = 0.0;
    mass_residual_loc = Vector3(0.0, 0.0, 0.0);
    energy_residual = 0.0;
    energy_residual_loc = Vector3(0.0, 0.0, 0.0);
    for ( FV_Cell *cellp: active_cells ) {
	cellp->rho_at_start_of_step = cellp->fs->gas->rho;
	cellp->rE_at_start_of_step = cellp->U[0]->total_energy;
    }
    return SUCCESS;
} // end of init_residuals()


int Block::compute_residuals(size_t dimensions, size_t gtl)
/// \brief Compute the residuals using previously stored data.
///
/// The largest residual of density for all cells was the traditional way
/// mbcns/Elmer estimated the approach to steady state.
/// However, with the splitting up of the increments for different physical
/// processes, this residual calculation code needed a bit of an update.
/// Noting that the viscous-stress, chemical and radiation increments
/// do not affect the mass within a cell, we now compute the residuals 
/// for both mass and (total) energy for all cells, the record the largest
/// with their location. 
{
    mass_residual = 0.0;
    mass_residual_loc = Vector3(0.0, 0.0, 0.0);
    energy_residual = 0.0;
    energy_residual_loc = Vector3(0.0, 0.0, 0.0);
    for ( FV_Cell *cellp: active_cells ) {
	double local_residual = (cellp->fs->gas->rho - cellp->rho_at_start_of_step) 
	    / cellp->fs->gas->rho;
	local_residual = fabs(local_residual);
	if ( local_residual > mass_residual ) {
	    mass_residual = local_residual;
	    mass_residual_loc.x = cellp->pos[gtl].x;
	    mass_residual_loc.y = cellp->pos[gtl].y;
	    mass_residual_loc.z = cellp->pos[gtl].z;
	}
	// In the following line, the zero index is used because,
	// at the end of the gas-dynamic update, that index holds
	// the updated data.
	local_residual = (cellp->U[0]->total_energy - cellp->rE_at_start_of_step) 
	    / cellp->U[0]->total_energy;
	local_residual = fabs(local_residual);
	if ( local_residual > energy_residual ) {
	    energy_residual = local_residual;
	    energy_residual_loc.x = cellp->pos[gtl].x;
	    energy_residual_loc.y = cellp->pos[gtl].y;
	    energy_residual_loc.z = cellp->pos[gtl].z;
	}
    } // for cellp
    return SUCCESS;
} // end of compute_residuals()


int Block::determine_time_step_size()
/// \brief Compute the local time step limit for all cells in the block.
///
/// The overall time step is limited by the worst-case cell.
/// \returns 0 on success, DT_SEARCH_FAILED otherwise.
///
/// \verbatim
/// Some Definitions...
/// ----------------
/// dt_global  : global time step for the block
/// G.cfl_target : desired CFL number
/// cfl_min  : approximate minimum CFL number in the block
/// cfl_max  : approximate maximum CFL number in the block
/// dt_allow : allowable time step (i.e. the maximum dt that
///            satisfies both the CFL target and the viscous
///            time step limit)
/// cfl_allow : allowable CFL number, t_order dependent
/// \endverbatim
{
    global_data *gdp = get_global_data_ptr();
    bool with_k_omega = (gdp->turbulence_model == TM_K_OMEGA && 
			 !gdp->separate_update_for_k_omega_source);
    bool first;
    double dt_local, cfl_local, signal, cfl_allow;
    // These limits allow the simulation of the sod shock tube
    // to get just a little wobbly around the shock.
    // Lower values of cfl should be used for a smooth solution.
    switch ( number_of_stages_for_update_scheme(get_gasdynamic_update_scheme()) ) {
    case 1: cfl_allow = 0.9; break;
    case 2: cfl_allow = 1.2; break;
    case 3: cfl_allow = 1.6; break;
    default: cfl_allow = 0.9;
    }

    first = true;
    for ( FV_Cell *cp: active_cells ) {
	signal = cp->signal_frequency(gdp->dimensions, with_k_omega);
	cfl_local = gdp->dt_global * signal; // Current (Local) CFL number
	dt_local = gdp->cfl_target / signal; // Recommend a time step size.
	if ( first ) {
	    cfl_min = cfl_local;
	    cfl_max = cfl_local;
	    dt_allow = dt_local;
	    first = false;
	} else {
	    if (cfl_local < cfl_min) cfl_min = cfl_local;
	    if (cfl_local > cfl_max) cfl_max = cfl_local;
	    if (dt_local < dt_allow) dt_allow = dt_local;
	}
    } // for cp
    if ( cfl_max > 0.0 && cfl_max < cfl_allow ) {
	return SUCCESS;
    } else {
	printf( "determine_time_step_size(): bad CFL number was encountered\n" );
	printf( "    cfl_max = %e for Block %d\n", cfl_max, static_cast<int>(id) );
	printf( "    If this cfl_max value is not much larger than 1.0,\n" );
	printf( "    your simulation could probably be restarted successfully\n" );
	printf( "    with some minor tweaking." );
	printf( "    That tweaking should probably include a reduction\n");
	printf( "    in the size of the initial time-step, dt\n");
	printf( "    If this job is a restart/continuation of an old job, look in\n");
	printf( "    the old-job.finish file for the value of dt at termination.\n");
	return DT_SEARCH_FAILED;
    }
} // end of determine_time_step_size()


int Block::detect_shock_points(size_t dimensions)
/// \brief Detects shocks by looking for compression between adjacent cells.
///
/// The velocity component normal to the cell interfaces
/// is used as the indicating variable.
{
    FV_Cell *cL, *cR;
    FV_Interface *IFace;
    double uL, uR, aL, aR, a_min;

    // Change in normalised velocity to indicate a shock.
    // A value of -0.05 has been found suitable to detect the levels of
    // shock compression observed in the "sod" and "cone20" test cases.
    // It may need to be tuned for other situations, especially when
    // viscous effects are important.
    global_data &G = *get_global_data_ptr();
    double tol = G.compression_tolerance;

    // First, work across North interfaces and
    // locate shocks using the (local) normal velocity.
    for ( size_t k = kmin; k <= kmax; ++k ) {
	for ( size_t i = imin; i <= imax; ++i ) {
	    for ( size_t j = jmin-1; j <= jmax; ++j ) {
		cL = get_cell(i,j,k);
		cR = get_cell(i,j+1,k);
		IFace = cL->iface[NORTH];
		uL = cL->fs->vel.x * IFace->n.x + cL->fs->vel.y * IFace->n.y + cL->fs->vel.z * IFace->n.z;
		uR = cR->fs->vel.x * IFace->n.x + cR->fs->vel.y * IFace->n.y + cR->fs->vel.z * IFace->n.z;
		aL = cL->fs->gas->a;
		aR = cR->fs->gas->a;
		if (aL < aR)
		    a_min = aL;
		else
		    a_min = aR;
		IFace->fs->S = ((uR - uL) / a_min < tol);
	    } // j loop
	} // i loop
    } // for k
    
    // Second, work across East interfaces and
    // locate shocks using the (local) normal velocity.
    for ( size_t k = kmin; k <= kmax; ++k ) {
	for ( size_t i = imin-1; i <= imax; ++i ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		cL = get_cell(i,j,k);
		cR = get_cell(i+1,j,k);
		IFace = cL->iface[EAST];
		uL = cL->fs->vel.x * IFace->n.x + cL->fs->vel.y * IFace->n.y + cL->fs->vel.z * IFace->n.z;
		uR = cR->fs->vel.x * IFace->n.x + cR->fs->vel.y * IFace->n.y + cR->fs->vel.z * IFace->n.z;
		aL = cL->fs->gas->a;
		aR = cR->fs->gas->a;
		if (aL < aR)
		    a_min = aL;
		else
		    a_min = aR;
		IFace->fs->S = ((uR - uL) / a_min < tol);
	    } // j loop
	} // i loop
    } // for k
    
    if ( dimensions == 3 ) {
	// Third, work across Top interfaces.
	for ( size_t i = imin; i <= imax; ++i ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t k = kmin-1; k <= kmax; ++k ) {
		    cL = get_cell(i,j,k);
		    cR = get_cell(i,j,k+1);
		    IFace = cL->iface[TOP];
		    uL = cL->fs->vel.x * IFace->n.x + cL->fs->vel.y * IFace->n.y + cL->fs->vel.z * IFace->n.z;
		    uR = cR->fs->vel.x * IFace->n.x + cR->fs->vel.y * IFace->n.y + cR->fs->vel.z * IFace->n.z;
		    aL = cL->fs->gas->a;
		    aR = cR->fs->gas->a;
		    if (aL < aR)
			a_min = aL;
		    else
			a_min = aR;
		    IFace->fs->S = ((uR - uL) / a_min < tol);
		} // for k
	    } // j loop
	} // i loop
    } // if ( dimensions == 3 )
    
    // Finally, mark cells as shock points if any of their
    // interfaces are shock points.
    for ( size_t k = kmin; k <= kmax; ++k ) {
	for ( size_t i = imin; i <= imax; ++i ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		FV_Cell *cellp = get_cell(i,j,k);
		cellp->fs->S = cellp->iface[EAST]->fs->S || cellp->iface[WEST]->fs->S ||
		    cellp->iface[NORTH]->fs->S || cellp->iface[SOUTH]->fs->S ||
		    ( dimensions == 3 && (cellp->iface[BOTTOM]->fs->S || cellp->iface[TOP]->fs->S) );
	    } // j loop
	} // i loop
    } // for k
    
    return SUCCESS;
} // end of detect_shock_points()

//-----------------------------------------------------------------------------

int find_nearest_cell(double x, double y, double z, 
		      size_t *jb_near,
		      size_t *i_near, size_t *j_near, size_t *k_near,
		      size_t gtl)
/// \brief Given an (x,y,z) position, locate the nearest cell centre.  
///
/// @param x, y, z : coordinates of the desired point
/// @param i_near, j_near, k_near : pointers to indices of the cell centre are stored here
/// @returns 1 for a close match, 0 if there were no close-enough cells.
{
    global_data *gd = get_global_data_ptr();
    Block *bdp;
    size_t ig, jg, kg, jbg;
    double dx, dy, dz, nearest, distance;

    jbg = 0;
    bdp = get_block_data_ptr(jbg);
    ig = bdp->imin; jg = bdp->jmin; kg = bdp->kmin;
    FV_Cell *mycp = bdp->get_cell(ig,jg,kg);
    dx = x - mycp->pos[gtl].x; dy = y - mycp->pos[gtl].y; dz = z - mycp->pos[gtl].z;
    nearest = sqrt(dx*dx + dy*dy + dz*dz);

    for ( size_t jb = 0; jb < gd->nblock; jb++ ) {
	bdp = get_block_data_ptr(jb);
	for ( FV_Cell *cp: bdp->active_cells ) {
	    dx = x - cp->pos[gtl].x; dy = y - cp->pos[gtl].y; dz = z - cp->pos[gtl].z;
	    distance = sqrt(dx*dx + dy*dy + dz*dz);
	    if (distance < nearest) {
		std::vector<size_t> ijk = bdp->to_ijk_indices(cp->id);
		ig = ijk[0]; jg = ijk[1]; kg = ijk[2];
		nearest = distance; jbg = jb;
	    }
	} // for cp
    } // for jb
    bdp = get_block_data_ptr(jbg);
    *jb_near = jbg; *i_near = ig; *j_near = jg; *k_near = kg;
    if ( nearest > bdp->get_ifi(ig,jg,kg)->length || 
	 nearest > bdp->get_ifj(ig,jg,kg)->length || 
	 (gd->dimensions == 3 && (nearest > bdp->get_ifk(ig,jg,kg)->length)) ) {
        printf("We really did not get close to (%e, %e, %e)\n", x, y, z);
        printf("Nearest %e lengths %e %e %e \n", nearest, 
	       bdp->get_ifi(ig,jg,kg)->length, bdp->get_ifj(ig,jg,kg)->length, 
	       bdp->get_ifk(ig,jg,kg)->length);
	return 0;
    } else {
        return 1;
    }
} // end find_nearest_cell()


// Some global data for locate_cell().
// To shorten the search, starting_block is the first block searched
// when requested locate_cell() is requested to locate the cell
// containing a specified point.
size_t starting_block = 0;

int locate_cell(double x, double y, double z,
	        size_t *jb_found, size_t *i_found, size_t *j_found, size_t *k_found,
		size_t gtl)
// Returns 1 if a cell containing the sample point (x,y,z) is found, else 0.
// The indices of the containing cell are recorded, if found.
//
// To consider: maybe we should use *jb_found as the starting_block value.
{
    global_data *gd = get_global_data_ptr();
    Block *bdp;
    Vector3 p = Vector3(x,y,z);
    *i_found = 0; *j_found = 0; *k_found = 0; *jb_found = 0;
    // Search the blocks, starting from the block in which the last point was found.
    for ( size_t jb = starting_block; jb < gd->nblock; jb++ ) {
	bdp = get_block_data_ptr(jb);
	for ( FV_Cell *cp: bdp->active_cells ) {
	    if ( cp->point_is_inside(p, gd->dimensions, gtl) ) {
		std::vector<size_t> ijk = bdp->to_ijk_indices(cp->id);
		*i_found = ijk[0]; *j_found = ijk[1]; *k_found = ijk[2]; 
		*jb_found = jb;
		starting_block = jb; // remember for next time
		return 1;
	    }
	} // for cp
    } // jb-loop
    // If we reach this point, then the point may be in one of the other blocks.
    for ( size_t jb = 0; jb < starting_block; jb++ ) {
	bdp = get_block_data_ptr(jb);
	for ( FV_Cell *cp: bdp->active_cells ) {
	    if ( cp->point_is_inside(p, gd->dimensions, gtl) ) {
		std::vector<size_t> ijk = bdp->to_ijk_indices(cp->id);
		*i_found = ijk[0]; *j_found = ijk[1]; *k_found = ijk[2]; 
		*jb_found = jb;
		starting_block = jb; // remember for next time
		return 1;
	    }
	} // for cp
    } // for jb
    // If we arrive here, we have not located the containing cell.
    return 0;
} // end locate_cell()

/** \file bc.cxx
 * \ingroup eilmer3
 * \brief Core Boundary Conditions for eilmer3.
 *
 * \author PAJ, Rowan and Dan.
 *
 */

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
#include <unistd.h>
#include <vector>

extern "C" {
#include <zlib.h>
}

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_catalytic.hh"
#include "bc_adjacent.hh"
#include "bc_supersonic_in.hh"
#include "bc_extrapolate_out.hh"
#include "bc_shock_fitting_in.hh"
#include "bc_slip_wall.hh"
#include "bc_adiabatic.hh"
#include "bc_fixed_t.hh"
#include "bc_sliding_t.hh"
#include "bc_subsonic_in.hh"
#include "bc_transient_uniform.hh"
#include "bc_static_profile.hh"
#include "bc_fixed_p_out.hh"
#include "bc_ablating.hh"
#include "bc_surface_energy_balance.hh"
#include "bc_user_defined.hh"
#include "bc_fstc.hh"
#include "kernel.hh"
#include "diffusion.hh"

using namespace std;

const int VERBOSE_BCS = 0; // Set to 1 to print mem usage, writing of data etc

//-----------------------------------------------------------------
// Class-based boundary conditions in which all information about 
// each boundary condition is contained within a single class 
// rather than scattered across separate functions.
// The trade-off is that we have to have boundary selection within 
// the new functions.

/// \brief Set up ghost-cell values for the inviscid flux calculations.
///
/// This is a coordinating function for use in the main time-stepping loop.
int apply_inviscid_bc( Block &bd, double t, int dimensions )
{
    bd.bcp[NORTH]->apply_inviscid(t);
    bd.bcp[EAST]->apply_inviscid(t);
    bd.bcp[SOUTH]->apply_inviscid(t);
    bd.bcp[WEST]->apply_inviscid(t);
    if ( dimensions == 3 ) {
	bd.bcp[TOP]->apply_inviscid(t);
	bd.bcp[BOTTOM]->apply_inviscid(t);
    }
    return SUCCESS;
}

int apply_viscous_bc( Block &bd, double t, int dimensions )
{
    bd.bcp[NORTH]->apply_viscous(t);
    bd.bcp[EAST]->apply_viscous(t);
    bd.bcp[SOUTH]->apply_viscous(t);
    bd.bcp[WEST]->apply_viscous(t);
    if ( dimensions == 3 ) {
	bd.bcp[TOP]->apply_viscous(t);
	bd.bcp[BOTTOM]->apply_viscous(t);
    }
    return SUCCESS;
}

//-----------------------------------------------------------------------

BoundaryCondition::
BoundaryCondition( Block *bdp, int which_boundary, int type_code,
		   std::string name_of_BC, int x_order,
		   bool is_wall, bool use_udf_flux,
		   int neighbour_block, int neighbour_face,
		   int neighbour_orientation,
		   int wc_bc, int sponge_flag, int xforce_flag )
    : bdp(bdp), which_boundary(which_boundary), type_code(type_code),
      name_of_BC(name_of_BC), x_order(x_order),
      is_wall_flag(is_wall), use_udf_flux_flag(use_udf_flux),
      neighbour_block(neighbour_block), neighbour_face(neighbour_face),
      neighbour_orientation(neighbour_orientation),
      wc_bc(wc_bc), sponge_flag(sponge_flag), xforce_flag(xforce_flag)
{
    Block & bd = *bdp;
    // 1. Determine size heat flux vectors
    int dim = 1;
    switch ( which_boundary ) {
    case NORTH:
	// j = constant at jmax
	dim = bd.nni * bd.nnk;
	break;
    case EAST:
	// i = constant at imax
	dim = bd.nnj * bd.nnk;
	break;
    case SOUTH:
	// j = constant at jmin
	dim = bd.nni * bd.nnk;
	break;
    case WEST:
	// i = constant at imin
	dim = bd.nnj * bd.nnk;
	break;
    case TOP:
	// k = constant at kmax
	dim = bd.nni * bd.nnj;
	break;
    case BOTTOM:
	// k = constant at kmin
	dim = bd.nni * bd.nnj;
	break;
    default:
	printf( "Error: which_boundary %d was not understood\n", 
		which_boundary );
	exit(NOT_IMPLEMENTED_ERROR);
    }
    
    // 2. Size heat flux vectors
    q_cond.resize( dim );
    q_diff.resize( dim );
    q_rad.resize( dim );
    
    // 3. Determine boundary cell indice limits
    switch( which_boundary ) {
    case NORTH:
	kmin = bd.kmin; kmax = bd.kmax;
	jmin = bd.jmax; jmax = bd.jmax;
	imin = bd.imin; imax = bd.imax;
	break;
    case EAST:
	kmin = bd.kmin; kmax = bd.kmax;
	jmin = bd.jmin; jmax = bd.jmax;
	imin = bd.imax; imax = bd.imax;
	break;
    case SOUTH:
	kmin = bd.kmin; kmax = bd.kmax;
	jmin = bd.jmin; jmax = bd.jmin;
	imin = bd.imin; imax = bd.imax;
	break;
    case WEST:
	kmin = bd.kmin; kmax = bd.kmax;
	jmin = bd.jmin; jmax = bd.jmax;
	imin = bd.imin; imax = bd.imin;
	break;
    case TOP:
	kmin = bd.kmax; kmax = bd.kmax;
	jmin = bd.jmin; jmax = bd.jmax;
	imin = bd.imin; imax = bd.imax;
	break;
    case BOTTOM:
	kmin = bd.kmin; kmax = bd.kmin;
	jmin = bd.jmin; jmax = bd.jmax;
	imin = bd.imin; imax = bd.imax;
	break;
    default:
	printf( "Error: Boundary code %d not understood.\n", which_boundary );
	exit(NOT_IMPLEMENTED_ERROR);
    }
    
    // 4. Set the catalytic wall BC pointer to null
    cw = 0;
    
#   if VERBOSE_BCS
    int total_bytes = 3 * dim * sizeof(double);
    
    cout << "Block " << bd.id << ", boundary " << which_boundary
	 << ": Have allocated " << total_bytes << " bytes of memory." << endl;
#   endif
}

// Shouldn't really have a boundary condition object created without
// reference to a particular block but, just in case the compiler wants it...
BoundaryCondition::
BoundaryCondition()
    : bdp(0), which_boundary(0), type_code(0),
      name_of_BC("Unspecified"), x_order(0),
      is_wall_flag(false), use_udf_flux_flag(false),
      neighbour_block(-1), neighbour_face(-1),
      neighbour_orientation(0),
      wc_bc(NON_CATALYTIC), sponge_flag(0), xforce_flag(0),
      cw(0)
{}

BoundaryCondition::
BoundaryCondition( const BoundaryCondition &bc )
    : bdp(bc.bdp), // Still bound to the original block.
      which_boundary(bc.which_boundary), type_code(bc.type_code),
      name_of_BC(bc.name_of_BC), x_order(bc.x_order),
      is_wall_flag(bc.is_wall_flag), use_udf_flux_flag(bc.use_udf_flux_flag),
      neighbour_block(bc.neighbour_block), neighbour_face(bc.neighbour_face),
      neighbour_orientation(bc.neighbour_orientation),
      wc_bc(bc.wc_bc), sponge_flag(bc.sponge_flag), xforce_flag(bc.xforce_flag),
      q_cond(bc.q_cond), q_diff(bc.q_diff), q_rad(bc.q_rad),
      imin(bc.imin), imax(bc.imax), jmin(bc.jmin), jmax(bc.jmax)
{
    // if ( bc.cw ) cw = new CatalyticWallBC(*bc.cw);
    cw = 0;
}

BoundaryCondition & BoundaryCondition::operator=( const BoundaryCondition &bc )
{
    if ( this != &bc ) {
	bdp = bc.bdp; // This new BC is still bound to the original block.
	which_boundary = bc.which_boundary; 
	type_code = bc.type_code;
	name_of_BC = bc.name_of_BC;
	x_order = bc.x_order;
	is_wall_flag = bc.is_wall_flag;
	use_udf_flux_flag = bc.use_udf_flux_flag;
	neighbour_block = bc.neighbour_block;
	neighbour_face = bc.neighbour_face;
	neighbour_orientation = bc.neighbour_orientation;
	wc_bc = bc.wc_bc;
	sponge_flag = bc.sponge_flag;
	xforce_flag = bc.xforce_flag;
	q_cond = bc.q_cond;
	q_diff = bc.q_diff;
	q_rad = bc.q_rad;
	imin = bc.imin; imax = bc.imax; jmin = bc.jmin; jmax = bc.jmax;
	// if ( bc.cw ) cw = new CatalyticWallBC(*bc.cw);
	cw = 0;
    }
    return *this;
}

BoundaryCondition::~BoundaryCondition()
{
    // Delete the catalytic BC object
    if ( cw ) delete cw;
}

void BoundaryCondition::print_info( std::string lead_in )
{
    cout << lead_in << "block_id" << bdp->id << endl;
    cout << lead_in << "which_boundary=" << which_boundary 
	 << "(" << get_face_name(which_boundary) << ")" << endl;
    cout << lead_in << "type_code=" << type_code 
	 << "(" << name_of_BC << ")" << endl;
    cout << lead_in << "x_order=" << x_order << endl;
    cout << lead_in << "is_wall_flag=" << is_wall_flag << endl;
    cout << lead_in << "use_udf_flux=" << use_udf_flux_flag << endl;
    cout << lead_in << "neighbour_block=" << neighbour_block << endl;
    cout << lead_in << "neighbour_face=" << neighbour_face 
	 << "(" << get_face_name(neighbour_face) << ")" << endl;
    cout << lead_in << "neighbour_orientation=" << neighbour_orientation << endl;
    cout << lead_in << "wc_bc=" << wc_bc << endl;
    cout << lead_in << "sponge_flag=" << sponge_flag << endl;
    cout << lead_in << "xforce_flag=" << xforce_flag << endl;
    return;
}

int BoundaryCondition::apply_inviscid( double t )
{
    // The default inviscid boundary condition is to reflect
    // the normal component of the velocity at the ghost-cell
    // centres -- slip-wall. 
    int i, j, k;
    FV_Cell *src_cell, *dest_cell;
    FV_Interface *IFace;
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[NORTH];
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i,j-1,k);
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[EAST];
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i-1,j,k);
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[SOUTH];
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i,j+1,k);
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[WEST];
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i+1,j,k);
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[TOP];
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i,j,k-1);
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[BOTTOM];
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		// ghost cell 2.
		src_cell = bd.get_cell(i,j,k+1);
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		reflect_normal_velocity(dest_cell, IFace);
		if (get_mhd_flag() == 1) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
}

int BoundaryCondition::apply_viscous( double t )
{
    // The default behaviour for viscous terms
    // is to do nothing at the boundary.
    return SUCCESS;
}

bool BoundaryCondition::is_wall()
{
    return is_wall_flag;
}

bool BoundaryCondition::use_udf_flux()
{
    return use_udf_flux_flag;
}

int BoundaryCondition::compute_surface_heat_flux( void ) 
/// \brief Calculate the heat flux values at all cell interfaces bounding walls
{
    FV_Interface * IFace;
    FV_Cell * cell_one;
    Vector3 cc, ic, i0c;
    int i, j, k;
    int index;
    Block & bd = *bdp;
    
    // 1. Check if this BC represents a wall
    if ( is_wall_flag ) {
	// 2. Loop over all bounding cells
	// NOTE: - restricting derivatives to first-order for now
	//       - using conductivities at wall (from interface)
	for ( k=kmin; k<=kmax; ++k ) {
	    for ( j=jmin; j<=jmax; ++j ) {
		for ( i=imin; i<=imax; ++i ) {
		    // calc. index into 1D heat-flux vectors
		    index = get_heat_flux_index( i, j, k );
		    cell_one = bd.get_cell(i,j,k);
		    IFace = cell_one->iface[which_boundary];
		    if ( compute_cell_interface_surface_heat_flux( IFace, cell_one, index ) ) {
			cerr << "BoundaryCondition::compute_surface_heat_flux()" << endl
			     << "compute_cell_interface_surface_heat_flux() failed for index: "
			     << index << endl;
			return FAILURE;
		    }
		} // end i loop
	    } // end j loop
	} // end k loop
    } // end if ( iswall )

    return SUCCESS;
}

int BoundaryCondition::
compute_cell_interface_surface_heat_flux(FV_Interface * IFace, 
					 FV_Cell * cell_one, int index) 
/// \brief Calculate the heat flux values for a single cell interface
{
    Vector3 cc, ic, i0c;
    double d1;
    int iT, isp;
    Gas_model * gm = get_gas_model_ptr();
    int nTs = gm->get_number_of_modes();
    int nsp = gm->get_number_of_species();
    double dTds;
    vector<double> dfds(nsp), dfd0(nsp), dfd00(nsp), js(nsp), j0(nsp), j00(nsp);
    double viscous_factor = get_viscous_factor();
    
    cc = cell_one->pos;
    ic = IFace->pos;
    i0c = cc - ic;
    d1 = - dot(i0c,IFace->n);
    // 1. Calculate conductive heat flux
    q_cond[index] = 0.0;
    for ( iT=0; iT<nTs; ++iT ){
	dTds = ( cell_one->fs->gas->T[iT] - IFace->fs->gas->T[iT] ) / d1;
	q_cond[index] += dTds * IFace->fs->gas->k[iT];
    }
    // 2. Calculate diffusive heat flux
    q_diff[index] = 0.0;
    if ( get_diffusion_flag() ) {
	for ( isp=0; isp<nsp; ++isp ){
	    dfds[isp] = ( cell_one->fs->gas->massf[isp] - IFace->fs->gas->massf[isp] ) / d1;
	}
	// Apply a diffusion model
	double D_t = 0.0;
	if ( get_k_omega_flag() == 1 ) {
	    double Sc_t = get_turbulence_schmidt_number();
	    D_t = IFace->fs->mu_t / (IFace->fs->gas->rho * Sc_t);
	}
	calculate_diffusion_fluxes(*(IFace->fs->gas), D_t, dfds, dfd0, dfd00, js, j0, j00 );
	for ( isp=0; isp<nsp; ++isp ){
	    q_diff[index] -= viscous_factor * js[isp] * gm->enthalpy(*(IFace->fs->gas), isp);
	}
    }

	// Need to do some modifying to the heat flux values for the purpose of outputting heat fluxes in the post processing
	// In the rest of the code q is used based on the normal vectors of the cells however here we want to specify that
	// q is always going into the boundary. The reasoning behind this is that the extraction of heat flux is usually for
	// the purpose of thermal structural modelling. This code simply switches the sign of the output so that it is going
	// out of the cell and into the structure for all faces
	if ( which_boundary == SOUTH or which_boundary == WEST or which_boundary == BOTTOM ) {
		q_cond[index] = q_cond[index] * (-1.0);
	}

    // 3. Calculate radiative heat flux - should be already stored in q_rad[]

    return SUCCESS;
}

int BoundaryCondition::write_surface_heat_flux( string filename, double sim_time ) 
/// \brief Write the heat flux values at all cell interfaces bounding walls
///
/// Basic formatting borrowed from Block::write_solution
{
    FV_Cell * cell;
    FV_Interface * IFace;
    int i, j, k;
    int index;
    double Re_wall;
    Block & bd = *bdp;
    
    FILE *fp;
#   if VERBOSE_BCS
    if ( bd.id == 0 && which_boundary == 0 ) {
	printf( "write_surface_heat_flux(): At t = %e, start block = %d, boundary = %d.\n",
	    sim_time, bd.id, which_boundary );
    }
#   endif
    if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	cerr << "write_solution(): Could not open " << filename << "; BAILING OUT" << endl;
	exit( FILE_ERROR );
    }
    fprintf(fp, "%20.12e\n", sim_time);
    string var_list = "";
    var_list += "\"index.i\" \"index.j\" \"index.k\" ";
    var_list += "\"pos.x\" \"pos.y\" \"pos.z\" ";
    var_list += "\"q_cond\" \"q_diff\" \"q_rad\" \"T_wall\"";
    var_list += "\"T_cell\" \"rho_cell\" \"un_cell\" \"Re_wall\"";
    fprintf(fp, "%s\n", var_list.c_str());
    
    // 1. Check if this BC represents a wall
    if ( is_wall_flag ) {
	// 2. Write dimensions of surface data to file
	fprintf(fp, "%d %d %d\n", (imax-imin+1), (jmax-jmin+1), (kmax-kmin+1));
	// 3. Loop over all interfaces
	// NOTE - assuming viscous BC's have been applied, as Twall is taken 
	//        to be equal to IFace->T[0]
	for ( k=kmin; k<=kmax; ++k ) {
	    for ( j=jmin; j<=jmax; ++j ) {
		for ( i=imin; i<=imax; ++i ) {
		    // calc. index into 1D heat-flux vectors
		    index = (jmax-jmin+1)*(imax-imin+1)*(k-kmin) + 
			    (imax-imin+1)*(j-jmin) + (i-imin);
		    cell = bd.get_cell(i,j,k);
		    IFace = cell->iface[which_boundary];
		    cell->fs->vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
		    Re_wall = cell->calculate_wall_Reynolds_number(which_boundary);
		    // 5. Write heat-flux data for interface to file
		    fprintf(fp, "%d %d %d ", i, j, k);
		    fprintf(fp, "%20.12e %20.12e %20.12e ", 
			    IFace->pos.x, IFace->pos.y, IFace->pos.z);
		    fprintf(fp, "%20.12e %20.12e %20.12e ", 
			    q_cond[index], q_diff[index], q_rad[index]);
		    fprintf(fp, "%20.12e ", IFace->fs->gas->T[0]);
		    fprintf(fp, "%20.12e %20.12e %20.12e %20.12e \n", 
		    	    cell->fs->gas->T[0], cell->fs->gas->rho, cell->fs->vel.x, Re_wall );
		    cell->fs->vel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
		} // end i loop
	    } // end j loop
	} // end k loop
    } // end if ( iswall )
    else {
	fprintf(fp, "0 0 0\n");
    } // end else
    fclose(fp);
    
    return SUCCESS;
} // end of Block::write_heat_flux_data()

int BoundaryCondition::write_fstc_heat_flux( string filename, double sim_time )
/// \brief Write the heat flux values at all cell interfaces bounding walls
///        Used for fluid-structure thermal coupling
/// Basic formatting borrowed from Block::write_solution
{
    FV_Cell * cell;
    FV_Interface * IFace;
    int i, j, k;
    int index;
    Block & bd = *bdp;

    FILE *fp;
#   if VERBOSE_BCS
    if ( bd.id == 0 && which_boundary == 0 ) {
	printf( "write_surface_heat_flux(): At t = %e, start block = %d, boundary = %d.\n",
	    sim_time, bd.id, which_boundary );
    }
#   endif
    if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	cerr << "write_solution(): Could not open " << filename << "; BAILING OUT" << endl;
	exit( FILE_ERROR );
    }
    fprintf(fp, "%20.12e\n", sim_time);
    string var_list = "";
    var_list += "\"index.i\" \"index.j\" \"index.k\" ";
    var_list += "\"pos.x\" \"pos.y\" \"pos.z\" ";
    var_list += "\"q_cond\" \"q_diff\" \"q_rad\" \"q_tot\" ";
    var_list += "\"T_wall\" \"T_cell\" ";
    fprintf(fp, "%s\n", var_list.c_str());

    // 1. Check if this BC represents a wall
    if ( is_wall_flag ) {
	// 2. Write dimensions of surface data to file
	fprintf(fp, "%d %d %d\n", (imax-imin+1), (jmax-jmin+1), (kmax-kmin+1));
	// 3. Loop over all interfaces
	// NOTE - assuming viscous BC's have been applied, as Twall is taken
	//        to be equal to IFace->T[0]
	for ( k=kmin; k<=kmax; ++k ) {
	    for ( j=jmin; j<=jmax; ++j ) {
		for ( i=imin; i<=imax; ++i ) {
		    // calc. index into 1D heat-flux vectors
		    index = (jmax-jmin+1)*(imax-imin+1)*(k-kmin) +
			    (imax-imin+1)*(j-jmin) + (i-imin);
		    cell = bd.get_cell(i,j,k);
		    IFace = cell->iface[which_boundary];
		    // 5. Write heat-flux data for interface to file
		    fprintf(fp, "%d %d %d ", i, j, k);
		    fprintf(fp, "%12.4e %12.4e %12.4e ",
			    IFace->pos.x, IFace->pos.y, IFace->pos.z);
		    fprintf(fp, "%12.4e %12.4e %12.4e ",
			    q_cond[index], q_diff[index], q_rad[index]);
                    fprintf(fp, "%12.4e ", q_cond[index]+q_diff[index]+q_rad[index]);
		    fprintf(fp, "%12.4e ", IFace->fs->gas->T[0]);
		    fprintf(fp, "%12.4e \n", cell->fs->gas->T[0] );
		} // end i loop
	    } // end j loop
	} // end k loop
    } // end if ( iswall )
    else {
	fprintf(fp, "0 0 0\n");
    } // end else
    fclose(fp);

    return SUCCESS;
} // end of BoundaryCondition::write_fstc_heat_flux()

double BoundaryCondition::
read_surface_heat_flux( string filename, int dimensions, int zip_files )
{
#   define NCHAR 4000
    char line[NCHAR];
    FILE *fp;
    gzFile zfp;
    char *gets_result;
    int i, j, k;
    int index;
    double sim_time;

    if ( get_verbose_flag() && which_boundary == 0 ) 
	printf("read_surface_heat_flux(): Start surface %d.\n", which_boundary);
    if (zip_files) {
	fp = NULL;
	filename += ".gz";
	if ( access(filename.c_str(), F_OK) != 0 ) return 0;
	if ((zfp = gzopen(filename.c_str(), "r")) == NULL) {
	    cerr << "read_surface_heat_flux(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit(FILE_ERROR);
	}
    } else {
	zfp = NULL;
	if ( access(filename.c_str(), F_OK) != 0 ) return 0;
	if ((fp = fopen(filename.c_str(), "r")) == NULL) {
	    cerr << "read_surface_heat_flux(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit(FILE_ERROR);
	}
    }
    if (zip_files) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if (gets_result == NULL) {
	printf("read_surface_heat_flux(): Empty surface heat flux file file while looking for sim_time value.\n");
	exit(BAD_INPUT_ERROR);
    }
    sscanf(line, "%lf", &sim_time);
    if ( get_verbose_flag() && which_boundary == 0 ) 
	printf("read_surface_heat_flux(): Time = %e\n", sim_time);
    if (zip_files) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if (gets_result == NULL) {
	printf("read_surface_heat_flux(): Empty flow field file while looking for line of variable names.\n");
	exit(BAD_INPUT_ERROR);
    }
    // The line just read should be the list of variable names, double-quoted.
    if (zip_files) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if ( gets_result == NULL ) {
	printf("read_surface_heat_flux(): Empty flow field file while looking for surface data dimensions.\n");
	exit(BAD_INPUT_ERROR);
    }
    sscanf(line, "%d %d %d", &i, &j, &k);
    if ( i==0 && j==0 && k==0 ) return 0.0;
    if ( i != (imax-imin+1) || j != (jmax-jmin+1) || k != ((dimensions == 3) ? (kmax-kmin+1) : 1) ) {
	printf("read_surface_heat_flux(): surface %d, mismatch in surface data dimensions\n", which_boundary);
	printf("    This misalignment could be caused by a having a different number\n");
	printf("    of fields for each cell's entry.\n");
	exit(BAD_INPUT_ERROR);
    }
    for ( k = kmin; k <= kmax; ++k ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( i = imin; i <= imax; ++i ) {
		// calc. index into 1D heat-flux vectors
		index = (jmax-jmin+1)*(imax-imin+1)*(k-kmin) + 
			(imax-imin+1)*(j-jmin) + (i-imin);
		// All surface element data is on one line.
		if (zip_files) {
		    gets_result = gzgets(zfp, line, NCHAR);
		} else {
		    gets_result = fgets(line, NCHAR, fp);
		}
		if (gets_result == NULL) {
		    printf("read_surface_heat_flux(): Empty flow field file while reading surface element data.\n");
		    exit(BAD_INPUT_ERROR);
		}
		scan_string_for_surface_heat_flux( q_cond[index], q_diff[index], q_rad[index], line );
	    }
	}
    }
    if (zip_files) {
	gzclose(zfp);
    } else {
	fclose(fp);
    }
    return sim_time;
#   undef NCHAR
}

int BoundaryCondition::write_vertex_velocities( std::string filename, double sim_time, int dimensions )
{
    int i, j, k, irangemax, jrangemax, krangemax;
    FV_Vertex *vtx;
    Block & bd = *bdp;
    
    FILE *fp;
#   if VERBOSE_BCS
    if ( bd.id == 0 && which_boundary == 0 ) {
	printf( "write_vertex_velocities(): At t = %e, start block = %d, boundary = %d.\n",
	    sim_time, bd.id, which_boundary );
    }
#   endif
    if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	cerr << "write_solution(): Could not open " << filename << "; BAILING OUT" << endl;
	exit( FILE_ERROR );
    }
    fprintf(fp, "%20.12e\n", sim_time);
    string var_list = "";
    var_list += "\"index.i\" \"index.j\" \"index.k\" ";
    var_list += "\"pos.x\" \"pos.y\" \"pos.z\" ";
    var_list += "\"vel.x\" \"vel.y\" \"vel.z\" ";
    fprintf(fp, "%s\n", var_list.c_str());
    
    // 2. Write dimensions of surface data to file
    if ( ( which_boundary == NORTH ) || ( which_boundary == SOUTH ) ) {
	fprintf(fp, "%d %d %d\n", (imax+1-imin+1), (jmax-jmin+1), (kmax-kmin+1));
	irangemax = imax+1;
	jrangemax = jmax;
	krangemax = kmax+1;
    } else if ( ( which_boundary == EAST ) || ( which_boundary == WEST ) ) {
	fprintf(fp, "%d %d %d\n", (imax-imin+1), (jmax+1-jmin+1), (kmax-kmin+1));
	irangemax = imax;
	jrangemax = jmax+1;
	krangemax = kmax+1;
    } else if ( ( which_boundary == TOP ) || ( which_boundary == BOTTOM ) ) {
	fprintf(fp, "%d %d %d\n", (imax-imin+1), (jmax-jmin+1), (kmax+1-kmin+1));
	irangemax = imax+1;
	jrangemax = jmax+1;
	krangemax = kmax;
    } else {
	printf( "Error: Boundary %d\n not implemented.", 
		which_boundary );
	exit(NOT_IMPLEMENTED_ERROR);
    }
// 3. Loop over all vertices
    if ( dimensions == 2 ) {
    	krangemax = kmax;
    }
    for ( k = kmin; k <= krangemax; ++k ) {
	for ( j = jmin; j <= jrangemax; ++j ) {
	    for ( i = imin; i <= irangemax; ++i ) {
		vtx = bd.get_vtx(i,j,k);
		fprintf(fp, "%d %d %d ", i, j, k);
		fprintf(fp, "%20.12e %20.12e %20.12e ", 
			vtx->pos.x, vtx->pos.y, vtx->pos.z);
		fprintf(fp, "%20.12e %20.12e %20.12e \n", 
			vtx->vel.x, vtx->vel.y, vtx->vel.z);
	    } // end i loop
	} // end j loop
    } // end k loop
    fclose(fp);
    return SUCCESS;
}

//------------------------------------------------------------------------

BoundaryCondition *create_BC( Block *bdp, int which_boundary, int type_of_BC, 
			      int inflow_condition_id, std::string filename, int n_profile,
			      double Twall, double Pout, int x_order, int is_wall, int use_udf_flux,
			      int other_block, int other_face, int neighbour_orientation,
			      int sponge_flag, int xforce_flag, 
			      vector<double> &mdot, double epsilon,
			      int wc_bc, string wcbc_fname, vector<double> f_wall,
			      double Twall_i, double Twall_f, double t_i, double t_f,
			      int assume_ideal)
// Return the appropriate BoundaryCondition object based on the
// integer identity of the boundary condition.
// Eventually, we would like to be able to use more descriptive strings
// as the identities.
//
{
    BoundaryCondition *newBC;

    switch ( type_of_BC ) {
    case ADJACENT:
	newBC = new AdjacentBC( bdp, which_boundary, other_block,
				other_face, neighbour_orientation );
	break;
    case SUP_IN:
	newBC = new SupersonicInBC( bdp, which_boundary, inflow_condition_id );
	break;
    case EXTRAPOLATE_OUT:
	newBC = new ExtrapolateOutBC( bdp, which_boundary, x_order, sponge_flag );
	break;
    case SHOCK_FITTING_IN:
	newBC = new ShockFittingInBC( bdp, which_boundary, inflow_condition_id );
	break;
    case SLIP_WALL:
	newBC = new SlipWallBC( bdp, which_boundary );
	break;
    case ADIABATIC:
	newBC = new AdiabaticBC( bdp, which_boundary );
	break;
    case FIXED_T:
	newBC = new FixedTBC( bdp, which_boundary, Twall );
	break;
    case SLIDING_T:
	newBC = new SlidingTBC( bdp, which_boundary, Twall_i, Twall_f, t_i, t_f );
	break;
    case SUBSONIC_IN:
	newBC = new SubsonicInBC( bdp, which_boundary, inflow_condition_id, assume_ideal );
	break;
    case TRANSIENT_UNI:
	newBC = new TransientUniformBC( bdp, which_boundary, filename );
	break;
    case STATIC_PROF:
	newBC = new StaticProfileBC( bdp, which_boundary, filename, n_profile );
	break;
    case FIXED_P_OUT:
	newBC = new FixedPOutBC( bdp, which_boundary, Pout, x_order );
	break;
    case USER_DEFINED:
	newBC = new UserDefinedBC( bdp, which_boundary, filename, 
				   is_wall==1, use_udf_flux==1 );
	break;
    case ADJACENT_PLUS_UDF:
	newBC = new AdjacentPlusUDFBC( bdp, which_boundary, other_block, 
				       other_face, neighbour_orientation,
				       filename, is_wall==1, use_udf_flux==1);
	break;
    case SEB:
    	newBC = new SurfaceEnergyBalanceBC( bdp, which_boundary, epsilon );
    	break;
    case ABLATING:
    	newBC = new AblatingBC( bdp, which_boundary, Twall, mdot, filename );
    	break;
    case FSTC:
	newBC = new fstcBC( bdp, which_boundary, filename );
        break;
    default:
	cerr << "create_BC() error: boundary condition \"" << type_of_BC 
	     << "\" is not available." << endl;
	exit( FAILURE );
	newBC = new SlipWallBC( bdp, which_boundary );
	break;
    } // end switch
    
    if ( wc_bc != NON_CATALYTIC && type_of_BC != FIXED_T 
    	                        && type_of_BC != ADIABATIC
    	                        && type_of_BC != SEB
    	                        && type_of_BC != ABLATING
    	                        && type_of_BC != SLIDING_T
                                && type_of_BC != FSTC) {
        cerr << "create_BC() error: cannot use wc_bc " << wc_bc 
             << " with type_of_BC " << type_of_BC << endl
    	     << "Bailing out!" << endl;
        exit( FAILURE );
    }
    else if ( wc_bc == NON_CATALYTIC ) {
    	// do nothing (newBC->cw should already be a null pointer)
    }
    else if ( wc_bc == EQUIL_CATALYTIC ) {
    	newBC->cw = new EquilibriumCatalyticWallBC( wcbc_fname );
    }
    else if ( wc_bc == SUPER_CATALYTIC ) {
    	newBC->cw = new SuperCatalyticWallBC( f_wall );
    }
    else {
    	cerr << "create_BC() error: wc_bc " << wc_bc << " not available." << endl
    	     << "Bailing out!" << endl;
    	exit( FAILURE );
    }
    
    newBC->wc_bc = wc_bc;
    newBC->xforce_flag = xforce_flag;
    return newBC;
}

//------------------------------------------------------------------------

/** \brief Check that the blocks are reasonably connected.
 *
 * Each boundary should be connected to a maximum of one
 * other boundary.
 * Each connection should have two entries in the table.
 * (One each way.) This routine actually steps through each
 * connection and makes sure that the converse connection
 * is present and correct.
 * The number of cells along connected boundaries should match.
 *
 */
int check_connectivity()
{
    global_data &G = *get_global_data_ptr();
    Block *bdp;
    Block *other_bdp;
    int jb, face, other_block, other_face;
    int nnA, nnB;
    int fail = 0;

    // Check both forward and reverse connections.
    for ( jb = 0; jb < G.nblock; ++jb ) {
	bdp = get_block_data_ptr(jb);
	for ( face = NORTH; face <= ((G.dimensions == 2) ? WEST : BOTTOM); ++face ) {
	    other_block = bdp->bcp[face]->neighbour_block;
	    other_face = bdp->bcp[face]->neighbour_face;
	    if ( other_block >= 0 ) {
		other_bdp = get_block_data_ptr(other_block);
		if ( other_bdp->bcp[other_face]->neighbour_block != jb ||
		     other_bdp->bcp[other_face]->neighbour_face != face ) {
		    cerr << "blocks " << jb << " and " << other_block 
			 << " incorrectly connected" << endl;
		    fail = 1;
		}
	    }
	} // end for face
    } // end for jb

    if (fail == 0) {
        if ( get_verbose_flag() ) cout << "Forward and Backward connections are OK." << endl;
    } else {
        cout << "Block connections fail." << endl;
        exit( VALUE_ERROR );
    }

    if ( G.dimensions == 2 ) {
	// Check numbers of cells along connected boundaries.
	for ( jb = 0; jb < G.nblock; ++jb ) {
	    bdp = get_block_data_ptr(jb);
	    for ( face = NORTH; face <= WEST; ++face ) {
		other_block = bdp->bcp[face]->neighbour_block;
		other_face = bdp->bcp[face]->neighbour_face;
		if ( other_block >= 0 ) {
		    other_bdp = get_block_data_ptr(other_block);
		    if ( face == NORTH || face == SOUTH ) 
			nnA = bdp->nni; 
		    else 
			nnA = bdp->nnj;
		    if ( other_face == NORTH || other_face == SOUTH ) 
			nnB = other_bdp->nni; 
		    else 
			nnB = other_bdp->nnj;
		    if ( nnA != nnB ) {
			cout << "blocks " << jb << " and "<< other_block 
			     <<" have mismatched cells " << nnA << " " << nnB << endl;
			fail = 1;
		    }
		} // end if other_block
	    } // end for face
	} // end for jb
	if (fail == 0) {
	    if ( get_verbose_flag() ) cout << "Numbers of cells for adjacent boundaries are OK." << endl;
	} else {
	    cerr << "Numbers of cells for adjacent boundaries FAIL." << endl;
	    exit( VALUE_ERROR );
	}
    } else {
	// FIX-ME bring the code over from e3prep
	if ( get_verbose_flag() ) cout << "Numbers of cells on joined faces are not checked for 3D." << endl;
    }

    return SUCCESS;
} // end check_connectivity()

//------------------------------------------------------------------------

/// \brief Read a line from a boundaries ".heat" file
///
/// DFP, June 2010
int
scan_string_for_surface_heat_flux( double &q_cond, double &q_diff, double &q_rad, char *bufptr )
// There isn't any checking of the file content.
// If anything gets out of place, the result is wrong data.
{
    // Look for a new-line character and truncate the string there.
    char *cptr = strchr(bufptr, '\n');
    if ( cptr != NULL ) cptr = '\0'; 
    // Now, we should have a string with only numbers separated by spaces.
    int i = atoi(strtok( bufptr, " " )); // tokenize on space characters
    int j = atoi(strtok( NULL, " " )); 
    int k = atoi(strtok( NULL, " " )); 
    double x = atof(strtok( NULL, " " )); 
    double y = atof(strtok( NULL, " " )); 
    double z = atof(strtok( NULL, " " ));
    q_cond = atof(strtok( NULL, " " ));
    q_diff = atof(strtok( NULL, " " ));
    q_rad = atof(strtok( NULL, " " ));

    UNUSED_VARIABLE(i); UNUSED_VARIABLE(j); UNUSED_VARIABLE(k);
    UNUSED_VARIABLE(x); UNUSED_VARIABLE(y); UNUSED_VARIABLE(z);
    
    return SUCCESS;
} 

#undef VERBOSE_BCS


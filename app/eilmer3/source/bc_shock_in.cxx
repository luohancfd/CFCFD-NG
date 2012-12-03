// bc_shock_in.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_shock_in.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

ShockInBC::ShockInBC( Block &bdp, int which_boundary, 
				int inflow_condition_id )
    : BoundaryCondition(bdp, which_boundary, SHOCK_IN, "ShockIn", 
			0, false, false, -1, -1, 0), 
      inflow_condition_id(inflow_condition_id) {}

ShockInBC::ShockInBC( const ShockInBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation),
      inflow_condition_id(bc.inflow_condition_id) {}

ShockInBC::~ShockInBC() {}

int ShockInBC::apply_inviscid( double t )
{
    // Set up ghost cells with inflow state. 
    int i, j, k;
    FV_Cell *dest_cell;
    FV_Interface *dest_face;
    global_data &gdp = *get_global_data_ptr();
    CFlowCondition *gsp = gdp.gas_state[inflow_condition_id];

    switch ( which_boundary ) {
    case NORTH:
	j = bdp.jmax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		dest_cell = bdp.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifj(i,j+1,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*gsp);
		shock_inflow_fluxes(dest_face, bdp.omegaz);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bdp.imax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		dest_cell = bdp.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifi(i+1,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*gsp);
		shock_inflow_fluxes(dest_face, bdp.omegaz);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bdp.jmin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		dest_cell = bdp.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifj(i,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*gsp);
		shock_inflow_fluxes(dest_face, bdp.omegaz);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bdp.imin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		dest_cell = bdp.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifi(i,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*gsp);
		shock_inflow_fluxes(dest_face, bdp.omegaz);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bdp.kmax;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		dest_cell = bdp.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifk(i,j,k+1);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*gsp);
		shock_inflow_fluxes(dest_face, bdp.omegaz);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bdp.kmin;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		dest_cell = bdp.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifk(i,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*gsp);
		shock_inflow_fluxes(dest_face, bdp.omegaz);
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    } // end switch

    return SUCCESS;
}

/// \brief Calculate inviscid fluxes based only on freestream conditions. 
///
/// It is assumed that the shock is stationary.

int ShockInBC::shock_inflow_fluxes(FV_Interface *IFace, double omegaz)
{
    double rho, e, p, ke;
    double un, vt1, vt2;
    ConservedQuantities &F = *(IFace->F);
    FlowState *IFace_flow_state = IFace->fs;
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = get_gas_model_ptr()->get_number_of_modes();
    //
	rho = IFace_flow_state->gas->rho;
	un = IFace_flow_state->vel.x;
	vt1 = IFace_flow_state->vel.y;
	vt2 = IFace_flow_state->vel.z;
	p = IFace_flow_state->gas->p;
	e = IFace_flow_state->gas->e[0];
	ke = 0.5 * (un*un + vt1*vt1 + vt2*vt2);
	/* Kinetic energy per unit volume. */

	/* Mass flux (mass / unit time / unit area) */
	F.mass = rho * un;
	/* Flux of normal momentum */
	F.momentum.x = rho * un * un + p;
	/* Flux of tangential momentum */
	F.momentum.y = rho * un * vt1;
	F.momentum.z = rho * un * vt2;
	/* Flux of Total Energy */
	F.total_energy = rho * un * (e + ke) + p * un;
	F.tke = rho * un * IFace_flow_state->tke;  // turbulence kinetic energy
	F.omega = rho * un * IFace_flow_state->omega;  // pseudo vorticity
	/* Species mass flux */
	for ( int isp = 0; isp < nsp; ++isp ) {
	    F.massf[isp] = F.mass * IFace_flow_state->gas->massf[isp];
	}
	/* Individual energies. */
	// NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = F.mass * IFace_flow_state->gas->e[imode];
	}
	
    if ( omegaz != 0.0 ) {
	    // Rotating frame.
	    double x = IFace->pos.x;
	    double y = IFace->pos.y;
	    double rsq = x*x + y*y;
	    // The conserved quantity is rothalpy,
	    // so we need to take -(u**2)/2 off the total energy flux.
	    // Note that rotating frame velocity u = omegaz * r.
	    F.total_energy -= F.mass * 0.5*omegaz*omegaz*rsq;
    }
  
	return SUCCESS;
} // end ShockInBC::shock_inflow_fluxes()

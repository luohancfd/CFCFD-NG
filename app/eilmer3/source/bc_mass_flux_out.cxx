// bc_mass_flux_out.cxx
// PJ, 13-Jan-2014

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_mass_flux_out.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

MassFluxOutBC::MassFluxOutBC(Block *bdp, int which_boundary, 
			     double _mass_flux, double p_init, double _relax_factor)
    : BoundaryCondition(bdp, which_boundary, MASS_FLUX_OUT), 
      mass_flux(_mass_flux), external_pressure(p_init), relax_factor(_relax_factor) 
{}

MassFluxOutBC::MassFluxOutBC(const MassFluxOutBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      mass_flux(bc.mass_flux),
      external_pressure(bc.external_pressure),
      relax_factor(bc.relax_factor) 
{}

MassFluxOutBC::MassFluxOutBC()
    : BoundaryCondition(0, 0, MASS_FLUX_OUT),
      mass_flux(0.0), relax_factor(0.05) 
{}

MassFluxOutBC & MassFluxOutBC::operator=(const MassFluxOutBC &bc)
{
    BoundaryCondition::operator=(bc);
    mass_flux = bc.mass_flux; // Ok for self-assignment.
    external_pressure = bc.external_pressure;
    relax_factor = bc.relax_factor;
    return *this;
}

MassFluxOutBC::~MassFluxOutBC() {}

void MassFluxOutBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "mass_flux= " << mass_flux << " kg/s/m**2" << endl;
    cout << lead_in << "(initial)external_pressure= " << external_pressure << " Pa" <<endl;
    cout << lead_in << "relax_factor= " << relax_factor << endl;
    return;
}

int MassFluxOutBC::apply_convective(double t)
{
    // Fill ghost cells with data from just inside the boundary
    // using zero-order extrapolation (i.e. just copy the data)
    // and then impose a pressure that has been incremented to
    // eliminate the error between the specified mass_flux and
    // the current estimated average mass_flux across the boundary.
    //
    // We assume that this boundary is an outflow boundary.
    //
    // See notes on page 82,83 of PJ's Eilmer3 workbook, Jan 2014.
    // Assuming approximately constant density, the pressure has
    // been linked to the velocity with the Bernoulli equation.
    size_t i, j, k;
    FV_Cell *src_cell, *dest_cell;
    Gas_model *gmodel = get_gas_model_ptr();
    Block & bd = *bdp;
    double area = 0.0;
    double rhoUA = 0.0; // current mass_flux through boundary
    double rhoA = 0.0;
    FV_Interface *face;
    double dp = 0.0;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	// First, estimate current mass flux across boundary.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[NORTH];
		area += face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA += src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end i loop
	} // for k
	dp = relax_factor * 0.5 * rhoA/area * (rhoUA*rhoUA/(area*area) - mass_flux*mass_flux);
	external_pressure += dp;
	// Apply ghost-cell conditions.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
	// First, estimate current mass flux across boundary.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[EAST];
		area += face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA += src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end j loop
	} // for k
	dp = relax_factor * 0.5 * rhoA/area * (rhoUA*rhoUA/(area*area) - mass_flux*mass_flux);
	external_pressure += dp;
	// Apply ghost-cell conditions.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
	// First, estimate current mass flux across boundary.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[SOUTH];
		area += face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA -= src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end i loop
	} // for k
	dp = relax_factor * 0.5 * rhoA/area * (rhoUA*rhoUA/(area*area) - mass_flux*mass_flux);
	external_pressure += dp;
	// Apply ghost-cell conditions.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
	// First, estimate current mass flux across boundary.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[WEST];
		area += face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA -= src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end j loop
	} // for k
	dp = relax_factor * 0.5 * rhoA/area * (rhoUA*rhoUA/(area*area) - mass_flux*mass_flux);
	external_pressure += dp;
	// Apply ghost-cell conditions.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
	// First, estimate current mass flux across boundary.
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[TOP];
		area += face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA += src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end j loop
	} // for i
	dp = relax_factor * 0.5 * rhoA/area * (rhoUA*rhoUA/(area*area) - mass_flux*mass_flux);
	external_pressure += dp;
	// Apply ghost-cell conditions.
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
	// First, estimate current mass flux across boundary.
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[BOTTOM];
		area += face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA -= src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end j loop
	} // for i
	dp = relax_factor * 0.5 * rhoA/area * (rhoUA*rhoUA/(area*area) - mass_flux*mass_flux);
	external_pressure += dp;
	// Apply ghost-cell conditions.
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = external_pressure;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for i
 	break;
    default:
	throw std::runtime_error("apply_inviscid not implemented for this boundary");
    } // end switch

    return SUCCESS;
}

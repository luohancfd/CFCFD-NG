// bc_subsonic_in.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_subsonic_in.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

SubsonicInBC::SubsonicInBC( Block &bdp, int which_boundary, 
			    int inflow_condition_id, int assume_ideal )
    : BoundaryCondition(bdp, which_boundary, SUBSONIC_IN, "SubsonicInBC",
			false, false, -1, -1, 0), 
      inflow_condition_id(inflow_condition_id),
      use_ideal_gas_relations(assume_ideal)
{}

SubsonicInBC::SubsonicInBC( const SubsonicInBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation), 
      inflow_condition_id(bc.inflow_condition_id) 
{}

SubsonicInBC::~SubsonicInBC() {}

int SubsonicInBC::apply_inviscid( double t )
{
    int i, j, k;
    FV_Cell *src_cell, *dest_cell;
    double u, v, w, velocity;
    global_data &gdp = *get_global_data_ptr();
    CFlowCondition *gstagp = gdp.gas_state[inflow_condition_id];
    CFlowCondition *gsp = new CFlowCondition(*gstagp);

    switch ( which_boundary ) {
    case NORTH:
	j = bdp.jmax;
	for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		src_cell = bdp.get_cell(i,j,k);
		u = src_cell->fs->vel.x;
		v = src_cell->fs->vel.y;
		w = src_cell->fs->vel.z;
		velocity = sqrt(u*u + v*v + w*w);
		subsonic_inflow_properties(gstagp, gsp, velocity);
		dest_cell = bdp.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bdp.imax;
	for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		u = src_cell->fs->vel.x;
		v = src_cell->fs->vel.y;
		w = src_cell->fs->vel.z;
		velocity = sqrt(u*u + v*v + w*w);
		subsonic_inflow_properties(gstagp, gsp, velocity);
		dest_cell = bdp.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bdp.jmin;
	for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		src_cell = bdp.get_cell(i,j,k);
		u = src_cell->fs->vel.x;
		v = src_cell->fs->vel.y;
		w = src_cell->fs->vel.z;
		velocity = sqrt(u*u + v*v + w*w);
		subsonic_inflow_properties(gstagp, gsp, velocity);
		dest_cell = bdp.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bdp.imin;
	for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		u = src_cell->fs->vel.x;
		v = src_cell->fs->vel.y;
		w = src_cell->fs->vel.z;
		velocity = sqrt(u*u + v*v + w*w);
		subsonic_inflow_properties(gstagp, gsp, velocity);
		dest_cell = bdp.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bdp.kmax;
	for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		u = src_cell->fs->vel.x;
		v = src_cell->fs->vel.y;
		w = src_cell->fs->vel.z;
		velocity = sqrt(u*u + v*v + w*w);
		subsonic_inflow_properties(gstagp, gsp, velocity);
		dest_cell = bdp.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bdp.kmin;
	for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		u = src_cell->fs->vel.x;
		v = src_cell->fs->vel.y;
		w = src_cell->fs->vel.z;
		velocity = sqrt(u*u + v*v + w*w);
		subsonic_inflow_properties(gstagp, gsp, velocity);
		dest_cell = bdp.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    } // end switch

    delete gsp;
    return SUCCESS;
}

/// \brief Estimate the free-stream conditions at the boundary cell. 
///
/// It is assumed that the velocity passed to this function 
/// is the inflow velocity (magnitude) of the cell just inside the boundary.
///
/// Boundary condition developed by Rowan Gollan (2001) for a perfect gas, only.
/// Extended by PJ to work with arbitrary (equilibrium) gases, July-2011.
/// Still doesn't work for nonequilibrium chemistry and vibrational-nonequilibrium
/// is likewise ignored.
///
/// TODO: make this like Hannes' subsonic inflow UDF for use in turbomachinery calcs.
int SubsonicInBC::subsonic_inflow_properties(CFlowCondition *stagnation, 
					     CFlowCondition *inflow_state, 
					     double inflow_velocity)
{
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int status_flag;
    //
    // Carry over some of the properties without further calculation.
    inflow_state->tke = stagnation->tke;
    inflow_state->omega = stagnation->omega;
    for ( int isp = 0; isp < nsp; ++isp ) {
	inflow_state->gas->massf[isp] = stagnation->gas->massf[isp];
    }
    //
    // Compute free-stream properties from stagnation conditions
    // and free-stream velocity.
    double a_0 = stagnation->gas->a;
    double p_0 = stagnation->gas->p;
    double rho_0 = stagnation->gas->rho;
    if ( use_ideal_gas_relations ) {
	double R = gmodel->R(*(stagnation->gas), status_flag);
	double C_v = gmodel->Cv(*(stagnation->gas), status_flag);
	double C_p = gmodel->Cp(*(stagnation->gas), status_flag);
	double GAMMA = C_p / C_v;
	double gm1 = GAMMA - 1.0;
	//
	// TODO: check these formulae; this approach does not match the stepping
	// approach used below.
	inflow_state->u = inflow_velocity;
	inflow_state->gas->a = sqrt(pow(a_0, 2) - gm1 * pow(inflow_velocity, 2) / 2);
	inflow_state->gas->T[0] = pow(inflow_state->gas->a, 2) / (GAMMA * R);
	//
	double M = inflow_state->u / inflow_state->gas->a;
	double T_ratio = 1 + gm1 / 2 * pow(M, 2);
	//
	inflow_state->gas->p = p_0 / pow(T_ratio, (GAMMA / gm1));
	inflow_state->gas->rho = rho_0 / pow(T_ratio, (1 / gm1));
	inflow_state->gas->e[0] = C_v * inflow_state->gas->T[0];
    } else {
	// Assume an internally-reversible process and use Gibbs equation
	// to take small isentropic steps from the stagnation conditions
	// down to the inflow condition that has to have the same total enthalpy.
	double e_0 = stagnation->gas->e[0];
	double H_0 = e_0 + p_0/rho_0; // Total enthalpy, at stagnation.
	// Now expand the free-stream gas to achieve the same total enthalpy.
	double kinetic_energy = 0.5 * inflow_velocity * inflow_velocity;
	double p = p_0;
	double rho = rho_0;
	double e = e_0;
	double H = e + p/rho + kinetic_energy;
	double delta = -0.005; // relative increment in density, small.
	int count = 0;
	int count_limit = 500;
	// If the density has to drop too far, we probably have 
	// a supersonic inflow condition and shouldn't be using this one.
	// TODO: revisit this stepping and follow up with an interpolation
	// stage to get very close to the correct inflow condition.
	while ( H > H_0 and count < count_limit ) {
	    double drho = rho * delta;
	    double de = p*drho/(rho*rho);
	    rho += drho;
	    e += de;
	    inflow_state->gas->rho = rho;
	    inflow_state->gas->e[0] = e;
	    gmodel->eval_thermo_state_rhoe(*(inflow_state->gas));
	    p = inflow_state->gas->p;
	    H = e + p/rho + kinetic_energy;
	    ++count;
	}
	if ( count >= count_limit ) {
	    cout << "SubsonicInBC: count_limit exceeded for generalized stepping." << endl;
	    cout << "    rho_0=" << rho_0 << ", p=" << p_0 << ", e_0=" << e_0 
		 << ", T_0=" << stagnation->gas->T[0] << ", H_0=" << H_0 << endl;
	    cout << "    rho=" << rho << ", p=" << p << ", e=" << e 
		 << ", v=" << inflow_velocity << ", T=" << inflow_state->gas->T[0]
		 << ", H=" << H << endl;
	}
    }
    if ( get_viscous_flag() ) gmodel->eval_transport_coefficients(*(inflow_state->gas));
    if ( get_diffusion_flag() ) gmodel->eval_diffusion_coefficients(*(inflow_state->gas));
    return SUCCESS;
} // end SubsonicInBC::subsonic_inflow_properties()

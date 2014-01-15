// bc_subsonic_in.cxx

#include <numeric>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_subsonic_in.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

SubsonicInBC::SubsonicInBC(Block *bdp, int which_boundary, 
			   int inflow_condition_id, int assume_ideal)
    : BoundaryCondition(bdp, which_boundary, SUBSONIC_IN), 
      inflow_condition_id(inflow_condition_id),
      use_ideal_gas_relations(assume_ideal)
{
    Gas_model *gmodel = get_gas_model_ptr();
    global_data &gd = *get_global_data_ptr();
    CFlowCondition *gstagp = gd.gas_state[inflow_condition_id];
    // Set all temperatures at T[0]
    for ( int imode = 1; imode < gmodel->get_number_of_modes(); ++imode ) {
	gstagp->gas->T[imode] = gstagp->gas->T[0];
    }
    gmodel->eval_thermo_state_pT(*(gstagp->gas));
    s0 = gmodel->mixture_entropy(*(gstagp->gas));
    h0 = gmodel->mixture_enthalpy(*(gstagp->gas));
    cout << "SubsonicInBC: set up boundary condition." << endl;
    cout << "    p0= " << gstagp->gas->p << " Pa" << endl;
    cout << "    T0= " << gstagp->gas->T[0] << " degrees K" << endl;
    cout << "    s0= " << s0 << " J/kg/K" << endl;
    cout << "    h0= " << h0 << " J/kg" << endl;
}

SubsonicInBC::SubsonicInBC(const SubsonicInBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code), 
      inflow_condition_id(bc.inflow_condition_id),
      use_ideal_gas_relations(bc.use_ideal_gas_relations),
      s0(bc.s0), h0(bc.h0)
{}

SubsonicInBC::SubsonicInBC()
    : BoundaryCondition(0, 0, SUBSONIC_IN), 
      inflow_condition_id(0), use_ideal_gas_relations(0), s0(0.0), h0(0.0)
{}

SubsonicInBC & SubsonicInBC::operator=(const SubsonicInBC &bc)
{
    // Benign for self-assignment.
    BoundaryCondition::operator=(bc);
    inflow_condition_id = bc.inflow_condition_id;
    use_ideal_gas_relations = bc.use_ideal_gas_relations;
    s0 = bc.s0;
    h0 = bc.h0;
    return *this;
}

SubsonicInBC::~SubsonicInBC() {}

void SubsonicInBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "inflow_condition_id= " << inflow_condition_id << endl;
    cout << lead_in << "use_ideal_gas_relations= " << use_ideal_gas_relations << endl;
    cout << lead_in << "entropy= " << s0 << endl;
    return;
}

int SubsonicInBC::apply_convective(double t)
{
    Block & bd = *bdp;
    size_t i, j, k;
    FV_Cell *src_cell, *dest_cell;
    FV_Interface *face;
    double p;
    global_data &gd = *get_global_data_ptr();
    CFlowCondition &gstagp = *(gd.gas_state[inflow_condition_id]);
    CFlowCondition gsp(gstagp);

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[NORTH];
		p = src_cell->fs->gas->p;
		// Choose a flow direction that is locally-inward at this location on the boundary.
		subsonic_inflow_properties(gstagp, -(face->n.x), -(face->n.y), -(face->n.z), gsp, p);
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(gsp);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[EAST];
		p = src_cell->fs->gas->p;
		subsonic_inflow_properties(gstagp, -(face->n.x), -(face->n.y), -(face->n.z), gsp, p);
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(gsp);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[SOUTH];
		p = src_cell->fs->gas->p;
		subsonic_inflow_properties(gstagp, face->n.x, face->n.y, face->n.z, gsp, p);
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(gsp);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[WEST];
		p = src_cell->fs->gas->p;
		subsonic_inflow_properties(gstagp, face->n.x, face->n.y, face->n.z, gsp, p);
#               if 0
		// Some debug... PJ 15-Jan-2014
		cout << "subsonic conditions:" << endl;
		cout << "    p= " << p 
		     << " gsp.p= " << gsp.gas->p
		     << " T= " << gsp.gas->T[0] << endl;
		cout << "    gsp.u= " << gsp.u 
		     << " gsp.v= " << gsp.v
		     << " gsp.w= " << gsp.w << endl; 
#               endif
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(gsp);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[TOP];
		p = src_cell->fs->gas->p;
		subsonic_inflow_properties(gstagp, -(face->n.x), -(face->n.y), -(face->n.z), gsp, p);
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(gsp);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[BOTTOM];
		p = src_cell->fs->gas->p;
		subsonic_inflow_properties(gstagp, face->n.x, face->n.y, face->n.z, gsp, p);
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(gsp);
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

/// \brief Estimate the free-stream conditions at the boundary cell. 
///
/// We use the pressure from the cell just inside the boundary
/// and compute the isentropic expansion from stagnation conditions
/// to the cell pressure to compute the gas state.
///
/// Boundary condition developed by Rowan Gollan (2001) for a perfect gas, only.
/// Extended by PJ to work with arbitrary (equilibrium) gases, July-2011.
/// Still doesn't work for nonequilibrium chemistry 
/// (and vibrational-nonequilibrium is likewise ignored).
/// Generalise and revised by RJG for all gases, August 2013.
/// 14-Jan-2014 PJ make direction an input so that we can specify
///                the normal to the boundary.
///
/// TODO: make this like Hannes' subsonic inflow UDF for use in turbomachinery calcs.
int SubsonicInBC::subsonic_inflow_properties(const CFlowCondition &stagnation,
					     double dir_x, double dir_y, double dir_z,
					     CFlowCondition &inflow_state, 
					     double inflow_pressure)
{
    global_data &G = *get_global_data_ptr();
    Gas_model &gmodel = *get_gas_model_ptr();
    size_t nsp = gmodel.get_number_of_species();
    size_t nmodes = gmodel.get_number_of_modes();

    // Carry over some of the properties without further calculation.
    inflow_state.tke = stagnation.tke;
    inflow_state.omega = stagnation.omega;
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	inflow_state.gas->massf[isp] = stagnation.gas->massf[isp];
    }
    // Compute inflow state by assuming isentropic expansion from stagnation conditions
    // Hence, we pass the stagnation entropy
    inflow_state.gas->p  = inflow_pressure;
    gmodel.eval_thermo_state_ps(*(inflow_state.gas), inflow_pressure, s0);
    // Compute gas speed from energy conservation
    double h = gmodel.mixture_enthalpy(*(inflow_state.gas));
    // If the internal energy is greater than stagnation conditions,
    // it probably means there is a wave trying to push back into
    // our reservior.  We'll stop that by just setting stagnated conditions.
    // then just return stagnated conditions
    if ( h >= h0 ) {
	inflow_state.gas->p = stagnation.gas->p;
	for ( size_t imode = 0; imode < nmodes; ++imode ) {
	    inflow_state.gas->T[imode] = stagnation.gas->T[0];
	}
	gmodel.eval_thermo_state_pT(*(inflow_state.gas));
	inflow_state.u = 0.0;
	inflow_state.v = 0.0;
	inflow_state.w = 0.0;
    } else {
	double speed = sqrt(2.0*(h0 - h));
	inflow_state.u = speed * dir_x;
	inflow_state.v = speed * dir_y;
	inflow_state.w = speed * dir_z;
    }

    if ( G.viscous ) gmodel.eval_transport_coefficients(*(inflow_state.gas));
    if ( G.diffusion ) gmodel.eval_diffusion_coefficients(*(inflow_state.gas));
    return SUCCESS;
} // end SubsonicInBC::subsonic_inflow_properties()

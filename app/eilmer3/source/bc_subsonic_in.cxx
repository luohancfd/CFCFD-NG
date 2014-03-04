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

void SubsonicInBC::setup_stagnation_condition()
{
    Gas_model &gmodel = *get_gas_model_ptr();
    // Assume equilibrium, set all temperatures at T[0]
    for ( int imode = 1; imode < gmodel.get_number_of_modes(); ++imode ) {
	gstagp.gas->T[imode] = gstagp.gas->T[0];
    }
    gmodel.eval_thermo_state_pT(*(gstagp.gas));
    s0 = gmodel.mixture_entropy(*(gstagp.gas));
    h0 = gmodel.mixture_enthalpy(*(gstagp.gas));
    return;
}

SubsonicInBC::SubsonicInBC(Block *bdp, int which_boundary, int inflow_condition_id_,
			   double mass_flux_, double relax_factor_,
			   subsonic_in_direction_t direction_type_,
			   std::vector<double> direction_vector_,
			   double direction_alpha_, double direction_beta_,
			   bool assume_ideal)
    : BoundaryCondition(bdp, which_boundary, SUBSONIC_IN), 
      inflow_condition_id(inflow_condition_id_),
      mass_flux(mass_flux_),
      relax_factor(relax_factor_),
      direction_type(direction_type_),
      direction_vector(direction_vector_),
      direction_alpha(direction_alpha_),
      direction_beta(direction_beta_),
      use_ideal_gas_relations(assume_ideal)
{
    global_data &gd = *get_global_data_ptr();
    gstagp = *(gd.gas_state[inflow_condition_id]);
    setup_stagnation_condition();
    double mag = 0.0;
    for ( double elem: direction_vector ) mag += elem * elem;
    if ( mag > 1.0e-6 ) {
	// Ensure that we have a unit vector.
	mag = sqrt(mag);
	for ( double &elem: direction_vector ) elem /= mag;
    } else {
	// An arbitrary, but well-defined, unit vector.
	direction_vector.resize(3);
	direction_vector[0] = 1.0;
	direction_vector[1] = 0.0;
	direction_vector[2] = 0.0;
    }	
    cout << "SubsonicInBC: set up boundary condition." << endl;
    cout << "    p0= " << gstagp.gas->p << " Pa" << endl;
    cout << "    T0= " << gstagp.gas->T[0] << " degrees K" << endl;
    cout << "    s0= " << s0 << " J/kg/K" << endl;
    cout << "    h0= " << h0 << " J/kg" << endl;
    cout << "    mass_flux= " << mass_flux << " kg/s/m**2" << endl;
    cout << "    relax_factor= " << relax_factor << endl;
    cout << "    direction_type= " << get_subsonic_in_direction_name(direction_type) << endl;
    cout << "    direction_vector= (" << direction_vector[0] << ", " << direction_vector[1]
	 << ", " << direction_vector[2] << ")" << endl;
    cout << "    direction_alpha= " << direction_alpha << endl;
    cout << "    direction_beta= " << direction_beta << endl;
    cout << "    use_ideal_gas_relations= " << use_ideal_gas_relations << endl;
}

SubsonicInBC::SubsonicInBC(const SubsonicInBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code), 
      inflow_condition_id(bc.inflow_condition_id),
      mass_flux(bc.mass_flux),
      relax_factor(bc.relax_factor),
      direction_type(bc.direction_type),
      direction_vector(bc.direction_vector),
      direction_alpha(bc.direction_alpha),
      direction_beta(bc.direction_beta),
      use_ideal_gas_relations(bc.use_ideal_gas_relations),
      gstagp(bc.gstagp)
{
    setup_stagnation_condition();
}

SubsonicInBC::SubsonicInBC()
    : BoundaryCondition(0, 0, SUBSONIC_IN), 
      inflow_condition_id(0), mass_flux(0.0), relax_factor(0.05), 
      direction_type(SUBSONIC_IN_NORMAL),
      direction_vector(std::vector<double>(3, 0.0)),
      direction_alpha(0.0),
      direction_beta(0.0),
      use_ideal_gas_relations(false)
{
    gstagp = CFlowCondition();
    s0 = 0.0;
    h0 = 0.0;
}

SubsonicInBC & SubsonicInBC::operator=(const SubsonicInBC &bc)
{
    if ( this != &bc ) { // Avoid self-assignment.
	BoundaryCondition::operator=(bc);
	inflow_condition_id = bc.inflow_condition_id;
	mass_flux = bc.mass_flux;
	relax_factor = bc.relax_factor;
	use_ideal_gas_relations = bc.use_ideal_gas_relations;
	gstagp = bc.gstagp;
	direction_type = bc.direction_type;
	direction_vector = bc.direction_vector;
	direction_alpha = bc.direction_alpha;
	direction_beta = bc.direction_beta;
    }
    setup_stagnation_condition();
    return *this;
}

SubsonicInBC::~SubsonicInBC() {}

void SubsonicInBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "inflow_condition_id= " << inflow_condition_id << endl;
    cout << lead_in << "mass_flux= " << mass_flux << endl;
    cout << lead_in << "relax_factor= " << relax_factor << endl;
    cout << lead_in << "direction_type= " 
	 << get_subsonic_in_direction_name(direction_type) << endl;
    cout << lead_in << "direction_vector= [" << direction_vector[0]
	 << ", " << direction_vector[1] << ", " << direction_vector[2] << "]" << endl;
    cout << lead_in << "direction_alpha= " << direction_alpha << endl;
    cout << lead_in << "direction_beta= " << direction_beta << endl;
    cout << lead_in << "use_ideal_gas_relations= " << use_ideal_gas_relations << endl;
    cout << lead_in << "entropy= " << s0 << endl;
    cout << lead_in << "enthalpy= " << h0 << endl;
    return;
}

int SubsonicInBC::apply_convective(double t)
{
    Block & bd = *bdp;
    size_t i, j, k;
    FV_Cell *src_cell, *dest_cell;
    FV_Interface *face;
    // global_data &gd = *get_global_data_ptr();
    CFlowCondition gsp(gstagp);
    double area = 0.0;
    double rhoUA = 0.0; // current mass_flux through boundary
    double rhoA = 0.0;
    double pA = 0.0;
    double p;
    double dp_over_p = 0.0;
    double speed, vt, vr, vz;
    double x, y, rxy;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	// First, estimate current inflow condition and mass flux into the block, across boundary.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[NORTH];
		area += face->area[0];
		pA += src_cell->fs->gas->p * face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA -= src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end i loop
	} // for k
	p = pA / area; // Average pressure across boundary.
	if ( mass_flux > 0.0 ) {
	    // Adjust the pressure to better achieve the specified mass flux.
	    dp_over_p = relax_factor * 0.5 * rhoA/area * (mass_flux*mass_flux - rhoUA*rhoUA/(area*area)) / p;
	    gstagp.gas->p *= 1.0 + dp_over_p;
	    setup_stagnation_condition();
	}
	speed = subsonic_inflow_properties(gstagp, gsp, p);
	// For turbo inflow, beta sets the axial flow.
	vz = speed * sin(direction_beta);
	// ...and alpha sets the angle of the flow in the plane of rotation.
	vt = speed * cos(direction_beta) * sin(direction_alpha);
	vr = -speed * cos(direction_beta) * cos(direction_alpha);
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[NORTH];
		switch ( direction_type ) {
		case SUBSONIC_IN_UNIFORM:
		    // We're given the flow direction.
		    gsp.u = speed * direction_vector[0];
		    gsp.v = speed * direction_vector[1];
		    gsp.w = speed * direction_vector[2];
		    break;
		case SUBSONIC_IN_AXIAL:
		    // Axial-flow through a presumably circular surface.
		    // [TODO] 27-Feb-2014 PJ: check that this fall-through is OK.
		case SUBSONIC_IN_RADIAL:
		    // Radial-in through a presumably cylindrical surface.
		    // We also presume that the grid is orthogonal to the boundary
		    // so that we don't have to recompute angles for each of the ghost cells.
		    x = face->pos.x; y = face->pos.y; rxy = sqrt(x*x + y*y);
		    gsp.u = vr * x/rxy - vt * y/rxy;
		    gsp.v = vt * x/rxy + vr * y/rxy;
		    gsp.w = vz;
		    break;
		case SUBSONIC_IN_NORMAL:
		default:
		    // Choose a flow direction that is locally-inward at this location on the boundary.
		    gsp.u = speed * -(face->n.x);
		    gsp.v = speed * -(face->n.y);
		    gsp.w = speed * -(face->n.z);
		}
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(gsp);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
	// First, estimate current inflow conditions and mass flux across boundary.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[EAST];
		area += face->area[0];
		pA += src_cell->fs->gas->p * face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA -= src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end j loop
	} // for k
	p = pA / area; // Average pressure across boundary.
	if ( mass_flux > 0.0 ) {
	    // Adjust the pressure to better achieve the specified mass flux.
	    dp_over_p = relax_factor * 0.5 * rhoA/area * (mass_flux*mass_flux - rhoUA*rhoUA/(area*area)) / p;
	    gstagp.gas->p *= 1.0 + dp_over_p;
	    setup_stagnation_condition();
	}
	speed = subsonic_inflow_properties(gstagp, gsp, p);
	// For turbo inflow, beta sets the axial flow.
	vz = speed * sin(direction_beta);
	// ...and alpha sets the angle of the flow in the plane of rotation.
	vt = speed * cos(direction_beta) * sin(direction_alpha);
	vr = -speed * cos(direction_beta) * cos(direction_alpha);
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[EAST];
		switch ( direction_type ) {
		case SUBSONIC_IN_UNIFORM:
		    // We're given the flow direction.
		    gsp.u = speed * direction_vector[0];
		    gsp.v = speed * direction_vector[1];
		    gsp.w = speed * direction_vector[2];
		    break;
		case SUBSONIC_IN_AXIAL:
		    // Axial-flow through a presumably circular surface.
		    // [TODO] 27-Feb-2014 PJ: check that this fall-through is OK.
		case SUBSONIC_IN_RADIAL:
		    // Radial-in through a presumably cylindrical surface.
		    // We also presume that the grid is orthogonal to the boundary
		    // so that we don't have to recompute angles for each of the ghost cells.
		    x = face->pos.x; y = face->pos.y; rxy = sqrt(x*x + y*y);
		    gsp.u = vr * x/rxy - vt * y/rxy;
		    gsp.v = vt * x/rxy + vr * y/rxy;
		    gsp.w = vz;
		    break;
		case SUBSONIC_IN_NORMAL:
		default:
		    // Choose a flow direction that is locally-inward at this location on the boundary.
		    gsp.u = speed * -(face->n.x);
		    gsp.v = speed * -(face->n.y);
		    gsp.w = speed * -(face->n.z);
		}
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(gsp);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
	// First, estimate current inflow conditions and mass flux across boundary.
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[SOUTH];
		area += face->area[0];
		pA += src_cell->fs->gas->p * face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA += src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end i loop
	} // for k
	p = pA / area; // Average pressure across boundary.
	if ( mass_flux > 0.0 ) {
	    // Adjust the pressure to better achieve the specified mass flux.
	    dp_over_p = relax_factor * 0.5 * rhoA/area * (mass_flux*mass_flux - rhoUA*rhoUA/(area*area)) / p;
	    gstagp.gas->p *= 1.0 + dp_over_p;
	    setup_stagnation_condition();
	}
	speed = subsonic_inflow_properties(gstagp, gsp, p);
	// For turbo inflow, beta sets the axial flow.
	vz = speed * sin(direction_beta);
	// ...and alpha sets the angle of the flow in the plane of rotation.
	vt = speed * cos(direction_beta) * sin(direction_alpha);
	vr = -speed * cos(direction_beta) * cos(direction_alpha);
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[SOUTH];
		switch ( direction_type ) {
		case SUBSONIC_IN_UNIFORM:
		    // We're given the flow direction.
		    gsp.u = speed * direction_vector[0];
		    gsp.v = speed * direction_vector[1];
		    gsp.w = speed * direction_vector[2];
		    break;
		case SUBSONIC_IN_AXIAL:
		    // Axial-flow through a presumably circular surface.
		    // [TODO] 27-Feb-2014 PJ: check that this fall-through is OK.
		case SUBSONIC_IN_RADIAL:
		    // Radial-in through a presumably cylindrical surface.
		    // We also presume that the grid is orthogonal to the boundary
		    // so that we don't have to recompute angles for each of the ghost cells.
		    x = face->pos.x; y = face->pos.y; rxy = sqrt(x*x + y*y);
		    gsp.u = vr * x/rxy - vt * y/rxy;
		    gsp.v = vt * x/rxy + vr * y/rxy;
		    gsp.w = vz;
		    break;
		case SUBSONIC_IN_NORMAL:
		default:
		    // Choose a flow direction that is locally-inward at this location on the boundary.
		    gsp.u = speed * face->n.x;
		    gsp.v = speed * face->n.y;
		    gsp.w = speed * face->n.z;
		}
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(gsp);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
	// First, estimate current inflow condition and mass flux across boundary.
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[WEST];
		area += face->area[0];
		pA += src_cell->fs->gas->p * face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA += src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end j loop
	} // for k
	p = pA / area; // Average pressure across boundary.
	if ( mass_flux > 0.0 ) {
	    // Adjust the pressure to better achieve the specified mass flux.
	    dp_over_p = relax_factor * 0.5 * rhoA/area * (mass_flux*mass_flux - rhoUA*rhoUA/(area*area)) / p;
	    gstagp.gas->p *= 1.0 + dp_over_p;
	    setup_stagnation_condition();
#           if 0
	    // Some debug... PJ 25-Feb-2014
	    cout << "Adjusting for mass flux:" << endl;
	    cout << "    area= " << area 
		 << " p= " << p
		 << " rhoA= " << rhoA
		 << " rhoUA= " << rhoUA << endl;
	    cout << "    dp_over_p= " << dp_over_p
		 << " gstagp_p= " << gstagp.gas->p << endl;
#           endif
	}
	speed = subsonic_inflow_properties(gstagp, gsp, p);
	// For turbo inflow, beta sets the axial flow.
	vz = speed * sin(direction_beta);
	// ...and alpha sets the angle of the flow in the plane of rotation.
	vt = speed * cos(direction_beta) * sin(direction_alpha);
	vr = -speed * cos(direction_beta) * cos(direction_alpha);
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[WEST];
		switch ( direction_type ) {
		case SUBSONIC_IN_UNIFORM:
		    // We're given the flow direction.
		    gsp.u = speed * direction_vector[0];
		    gsp.v = speed * direction_vector[1];
		    gsp.w = speed * direction_vector[2];
		    break;
		case SUBSONIC_IN_AXIAL:
		    // Axial-flow through a presumably circular surface.
		    // [TODO] 27-Feb-2014 PJ: check that this fall-through is OK.
		case SUBSONIC_IN_RADIAL:
		    // Radial-in through a presumably cylindrical surface.
		    // We also presume that the grid is orthogonal to the boundary
		    // so that we don't have to recompute angles for each of the ghost cells.
		    x = face->pos.x; y = face->pos.y; rxy = sqrt(x*x + y*y);
		    gsp.u = vr * x/rxy - vt * y/rxy;
		    gsp.v = vt * x/rxy + vr * y/rxy;
		    gsp.w = vz;
		    break;
		case SUBSONIC_IN_NORMAL:
		default:
		    // Choose a flow direction that is locally-inward at this location on the boundary.
		    gsp.u = speed * face->n.x;
		    gsp.v = speed * face->n.y;
		    gsp.w = speed * face->n.z;
		}
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
	// First, estimate current inflow conditions and mass flux across boundary.
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[TOP];
		area += face->area[0];
		pA += src_cell->fs->gas->p * face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA -= src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end j loop
	} // for i
	p = pA / area; // Average pressure across boundary.
	if ( mass_flux > 0.0 ) {
	    // Adjust the pressure to better achieve the specified mass flux.
	    dp_over_p = relax_factor * 0.5 * rhoA/area * (mass_flux*mass_flux - rhoUA*rhoUA/(area*area)) / p;
	    gstagp.gas->p *= 1.0 + dp_over_p;
	    setup_stagnation_condition();
	}
	speed = subsonic_inflow_properties(gstagp, gsp, p);
	// For turbo inflow, beta sets the axial flow.
	vz = speed * sin(direction_beta);
	// ...and alpha sets the angle of the flow in the plane of rotation.
	vt = speed * cos(direction_beta) * sin(direction_alpha);
	vr = -speed * cos(direction_beta) * cos(direction_alpha);
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[TOP];
		switch ( direction_type ) {
		case SUBSONIC_IN_UNIFORM:
		    // We're given the flow direction.
		    gsp.u = speed * direction_vector[0];
		    gsp.v = speed * direction_vector[1];
		    gsp.w = speed * direction_vector[2];
		    break;
		case SUBSONIC_IN_AXIAL:
		    // Axial-flow through a presumably circular surface.
		    // [TODO] 27-Feb-2014 PJ: check that this fall-through is OK.
		case SUBSONIC_IN_RADIAL:
		    // Radial-in through a presumably cylindrical surface.
		    // We also presume that the grid is orthogonal to the boundary
		    // so that we don't have to recompute angles for each of the ghost cells.
		    x = face->pos.x; y = face->pos.y; rxy = sqrt(x*x + y*y);
		    gsp.u = vr * x/rxy - vt * y/rxy;
		    gsp.v = vt * x/rxy + vr * y/rxy;
		    gsp.w = vz;
		    break;
		case SUBSONIC_IN_NORMAL:
		default:
		    // Choose a flow direction that is locally-inward at this location on the boundary.
		    gsp.u = speed * -(face->n.x);
		    gsp.v = speed * -(face->n.y);
		    gsp.w = speed * -(face->n.z);
		}
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(gsp);
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(gsp);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
	// First, estimate current inflow conditions and mass flux across boundary.
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[BOTTOM];
		area += face->area[0];
		pA += src_cell->fs->gas->p * face->area[0];
		rhoA += src_cell->fs->gas->rho * face->area[0];
		rhoUA += src_cell->fs->gas->rho * dot(src_cell->fs->vel, face->n) * face->area[0];
	    } // end j loop
	} // for i
	p = pA / area; // Average pressure across boundary.
	if ( mass_flux > 0.0 ) {
	    // Adjust the pressure to better achieve the specified mass flux.
	    dp_over_p = relax_factor * 0.5 * rhoA/area * (mass_flux*mass_flux - rhoUA*rhoUA/(area*area)) / p;
	    gstagp.gas->p *= 1.0 + dp_over_p;
	    setup_stagnation_condition();
	}
	speed = subsonic_inflow_properties(gstagp, gsp, p);
	// For turbo inflow, beta sets the axial flow.
	vz = speed * sin(direction_beta);
	// ...and alpha sets the angle of the flow in the plane of rotation.
	vt = speed * cos(direction_beta) * sin(direction_alpha);
	vr = -speed * cos(direction_beta) * cos(direction_alpha);
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		face = src_cell->iface[BOTTOM];
		switch ( direction_type ) {
		case SUBSONIC_IN_UNIFORM:
		    // We're given the flow direction.
		    gsp.u = speed * direction_vector[0];
		    gsp.v = speed * direction_vector[1];
		    gsp.w = speed * direction_vector[2];
		    break;
		case SUBSONIC_IN_AXIAL:
		    // Axial-flow through a presumably circular surface.
		    // [TODO] 27-Feb-2014 PJ: check that this fall-through is OK.
		case SUBSONIC_IN_RADIAL:
		    // Radial-in through a presumably cylindrical surface.
		    // We also presume that the grid is orthogonal to the boundary
		    // so that we don't have to recompute angles for each of the ghost cells.
		    x = face->pos.x; y = face->pos.y; rxy = sqrt(x*x + y*y);
		    gsp.u = vr * x/rxy - vt * y/rxy;
		    gsp.v = vt * x/rxy + vr * y/rxy;
		    gsp.w = vz;
		    break;
		case SUBSONIC_IN_NORMAL:
		default:
		    // Choose a flow direction that is locally-inward at this location on the boundary.
		    gsp.u = speed * face->n.x;
		    gsp.v = speed * face->n.y;
		    gsp.w = speed * face->n.z;
		}
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
/// 27-Feb-2014: PJ: return the scalar speed rather than set the velocity components.
///
double SubsonicInBC::subsonic_inflow_properties(const CFlowCondition &stagnation,
						CFlowCondition &inflow_state, 
						double inflow_pressure)
{
    global_data &G = *get_global_data_ptr();
    Gas_model &gmodel = *get_gas_model_ptr();
    size_t nsp = gmodel.get_number_of_species();
    size_t nmodes = gmodel.get_number_of_modes();
    double speed = 0.0;

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
	speed = sqrt(2.0*(h0 - h));
	inflow_state.u = speed;
	inflow_state.v = 0.0;
	inflow_state.w = 0.0;
    }

    if ( G.viscous ) gmodel.eval_transport_coefficients(*(inflow_state.gas));
    if ( G.diffusion ) gmodel.eval_diffusion_coefficients(*(inflow_state.gas));
    return speed;
} // end SubsonicInBC::subsonic_inflow_properties()

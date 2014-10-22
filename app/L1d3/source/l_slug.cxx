// l_slug.cxx
// Refactored from l1d code 25-Sep-2012

#include <vector>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <numeric>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "l1d.hh"
#include "l_cell.hh"
#include "l_slug.hh"
#include "l_io.hh"
#include "l_rivp.hh"

GasSlug::GasSlug(int indx, SimulationData& SD, 
		 std::string config_file_name, int echo_input)
{
    cout << "GasSlug constructor using config file." << endl;
    ConfigParser dict = ConfigParser(config_file_name);
    std::stringstream tag;
    tag << indx;
    std::string section = "slug-" + tag.str();
    if (echo_input == 1) cout << "Reading slug " << indx << " parameters..." << endl;
    Gas_model *gmodel = get_gas_model_ptr();

    dict.parse_int(section, "nn", nnx, 10);
    dict.parse_int(section, "cluster_to_end_L", cluster_to_end_1, 0);
    dict.parse_int(section, "cluster_to_end_R", cluster_to_end_2, 0);
    dict.parse_double(section, "cluster_strength", cluster_strength, 0.0);
    if (echo_input == 1) {
	cout << "    nn = " << nnx << endl;
	cout << "    cluster_to_end_L = " << cluster_to_end_1 << endl;
	cout << "    cluster_to_end_R = " << cluster_to_end_2 << endl;
	cout << "    cluster_strength = " << cluster_strength << endl;
    }
    // Adaptivity parameters for this gas slug. 
    dict.parse_int(section, "nnmax", nxmax, nnx);
    dict.parse_int(section, "adaptive", adaptive, 0);
    dict.parse_double(section, "dxmin", dxmin, 0.0);
    dict.parse_double(section, "dxmax", dxmax, 0.0);
    if (echo_input == 1) {
	cout << "    nnmax = " << nxmax << endl;
	cout << "    adaptive = " << adaptive << endl;
	cout << "    dxmin = " << dxmin << endl;
	cout << "    dxmax = " << dxmax << endl;
    }
    // Allow for ghost cells 
    nghost = 2;
    nxdim = nxmax + 2 * nghost;
    if ( nxdim > NDIM ) {
	cout << "**** Problem with Slug[" << indx << "]" << endl;
	cout << "     NDIM=" << NDIM 
	     << " is not large enough for nxdim=" << nxdim << endl;
        cout << "     nnx=" << nnx << " nxmax=" << nxmax << endl;
	cout << "Quitting." << endl;
	exit(1);
    }
    // An array of cells with internal structures that need to be allocated.
    for ( int i = 0; i < nxdim; ++i ) {
	Cell.push_back(LCell(gmodel));
    }
    set_index_range();

    // Flag for viscous effects:
    // =1: include them
    // =0: do not include them
    // Flag for wall temperature
    // =0: use specified wall temperature
    // =1: use adiabatic wall temperature
    dict.parse_int(section, "viscous_effects", viscous_effects, 0);
    dict.parse_int(section, "adiabatic_flag", adiabatic, 0);
    if (echo_input == 1) {
	cout << "    viscous_effects = " << viscous_effects << endl;
	cout << "    adiabatic_flag = " << adiabatic << endl;
    }

    // Boundary conditions.
    // By default, the gas slug is unconnected.
    left_slug_id = -1;
    left_slug_end_id = -1;
    right_slug_id = -1;
    right_slug_end_id = -1;
    //
    left_piston_id = -1;
    right_piston_id = -1;
    //
    left_diaphragm_id = -1;
    right_diaphragm_id = -1;
    //
    set_left_end_ustar = 0;
    set_right_end_ustar = 0;
    left_ustar = 0.0;
    right_ustar = 0.0;

    // Process BC data for left boundary. 
    std::string line, control_string, label;
    dict.parse_string(section, "BC_L", line, "");
    stringstream ss_L(line);
    ss_L >> control_string;
    if (control_string == "S") {
        left_bc_type = SLUG;
	ss_L >> left_slug_id;
	ss_L >> label;
        if (label[0] == 'L' || label[0] == 'l' || label[0] == '0')
            left_slug_end_id = LEFT;
        if (label[0] == 'R' || label[0] == 'r' || label[0] == '1')
            left_slug_end_id = RIGHT;
        if (echo_input == 1) {
            cout << "   left_boundary: neighbour slug = " 
		 << left_slug_id << " end = " << left_slug_end_id << endl;
        }
    } else if (control_string == "SD") {
        left_bc_type = SLUG_DIAPHRAGM;
	ss_L >> left_slug_id;
	ss_L >> label;
        if (label[0] == 'L' || label[0] == 'l' || label[0] == '0')
            left_slug_end_id = LEFT;
        if (label[0] == 'R' || label[0] == 'r' || label[0] == '1')
            left_slug_end_id = RIGHT;
	ss_L >> left_diaphragm_id;
        if (echo_input == 1) {
            cout << "    left_boundary: neighbour slug = "
		 << left_slug_id << " end = " <<  left_slug_end_id
		 << " diaphragm = " << left_diaphragm_id << endl;
        }
    } else if (control_string == "P") {
        left_bc_type = PISTON;
        ss_L >> left_piston_id;
        if (echo_input == 1) {
            cout << "    left_boundary: neighbour piston = "
		 << left_piston_id << endl;
        }
        set_left_end_ustar = 1;
    } else if (control_string == "V") {
        left_bc_type = SOLID_BOUNDARY;
        ss_L >> left_ustar;
        if (echo_input == 1) {
            cout << "    left_boundary: imposed velocity = "
		 << left_ustar << endl;
        }
        set_left_end_ustar = 1;
    } else if (control_string == "F") {
        left_bc_type = FREE_END;
        if (echo_input == 1) {
            cout << "    left_boundary: free-end" << endl;
        }
    } else {
        cout << "    Invalid control string: " << control_string << endl;
        exit(-1);
    } // end if control_string...

    // Process BC data for right boundary. 
    dict.parse_string(section, "BC_R", line, "");
    stringstream ss_R(line);
    ss_R >> control_string;
    if (control_string == "S") {
        right_bc_type = SLUG;
	ss_R >> right_slug_id;
	ss_R >> label;
        if (label[0] == 'L' || label[0] == 'l' || label[0] == '0')
            right_slug_end_id = LEFT;
        if (label[0] == 'R' || label[0] == 'r' || label[0] == '1')
            right_slug_end_id = RIGHT;
        if (echo_input == 1) {
            cout << "   right_boundary: neighbour slug = " 
		 << right_slug_id << " end = " << right_slug_end_id << endl;
        }
    } else if (control_string == "SD") {
        right_bc_type = SLUG_DIAPHRAGM;
	ss_R >> right_slug_id;
	ss_R >> label;
        if (label[0] == 'L' || label[0] == 'l' || label[0] == '0')
            right_slug_end_id = LEFT;
        if (label[0] == 'R' || label[0] == 'r' || label[0] == '1')
            right_slug_end_id = RIGHT;
	ss_R >> right_diaphragm_id;
        if (echo_input == 1) {
            cout << "    right_boundary: neighbour slug = "
		 << right_slug_id << " end = " <<  right_slug_end_id
		 << " diaphragm = " << right_diaphragm_id << endl;
        }
    } else if (control_string == "P") {
        right_bc_type = PISTON;
        ss_R >> right_piston_id;
        if (echo_input == 1) {
            cout << "    right_boundary: neighbour piston = "
		 << right_piston_id << endl;
        }
        set_right_end_ustar = 1;
    } else if (control_string == "V") {
        right_bc_type = SOLID_BOUNDARY;
        ss_R >> right_ustar;
        if (echo_input == 1) {
            cout << "    right_boundary: imposed velocity = "
		 << right_ustar << endl;
        }
        set_right_end_ustar = 1;
    } else if (control_string == "F") {
        right_bc_type = FREE_END;
        if (echo_input == 1) {
            cout << "    right_boundary: free-end" << endl;
        }
    } else {
        cout << "    Invalid control string: " << control_string << endl;
        exit(-1);
    } // end if control_string...

    // Time stepping and order of reconstruction.
    dt = SD.dt_init;
    cfl_target = SD.CFL;
    dt0 = dt;
    Torder = SD.Torder;
    Xorder = SD.Xorder;
    // Plotting and history output events.
    dict.parse_int(section, "hncell", hncell, 0);
    std::vector<int> vint_default;
    vint_default.resize(hncell);
    for ( size_t i = 0; i < vint_default.size(); ++i ) vint_default[i] = 0;
    dict.parse_vector_of_ints(section, "hxcell", hxcell, vint_default);
    if (echo_input == 1) {
	cout << "    hncell = " << hncell << endl;
	cout << "    hxcell =";
	for ( int i = 0; i < hncell; ++i ) cout << " " << hxcell[i];
	cout << endl;
    }
    // Initial slug state. 
    int nsp = gmodel->get_number_of_species();
    initial_flow_state = new LFlowState(gmodel);
    dict.parse_double(section, "initial_xL", xbegin, 0.0);
    dict.parse_double(section, "initial_xR", xend, 0.0);
    dict.parse_double(section, "initial_p", initial_flow_state->gas->p, 100.0e3);
    dict.parse_double(section, "initial_u", initial_flow_state->u, 0.0);
    dict.parse_double(section, "initial_T", initial_flow_state->gas->T[0], 300.0);
    std::vector<double> vdbl_default;
    vdbl_default.resize(nsp);
    for ( size_t i = 0; i < vdbl_default.size(); ++i ) vdbl_default[i] = 0.0;
    dict.parse_vector_of_doubles(section, "massf", initial_flow_state->gas->massf, vdbl_default);
    if (echo_input == 1) {
	cout << "    initial_xL = " << xbegin << endl;
	cout << "    initial_xR = " << xend << endl;
	cout << "    initial_p = " << initial_flow_state->gas->p << endl;
	cout << "    initial_u = " << initial_flow_state->u << endl;
	cout << "    initial_T = " << initial_flow_state->gas->T[0] << endl;
	cout << "    massf =";
	for ( int i = 0; i < nsp; ++i ) cout << " " << initial_flow_state->gas->massf[i];
	cout << endl;
    }
    double f_sum = initial_flow_state->gas->massf[0];
    for ( int isp = 1; isp < nsp; ++isp ) {
	f_sum += initial_flow_state->gas->massf[isp];
    }
    if ( fabs(f_sum - 1.0) > 1.0e-4 ) {
	printf( "Species mass fractions do not sum to 1.0: %e\n", f_sum );
	exit(-1);
    }
    // Density, Internal energy, Speed of Sound, and 
    // molecular transport coefficients. 
    gmodel->eval_thermo_state_pT(*(initial_flow_state->gas));
    gmodel->eval_transport_coefficients(*(initial_flow_state->gas), gmodel);
    double e = accumulate(initial_flow_state->gas->e.begin(),
			  initial_flow_state->gas->e.end(), 0.0);
    if (echo_input == 1) {
	cout << "    rho = " << initial_flow_state->gas->rho
	     << " e = " << e
	     << " a = " << initial_flow_state->gas->a << endl;
	cout << "    R = " << gmodel->R(*(initial_flow_state->gas))
	     << " Cv = " << gmodel->Cv(*(initial_flow_state->gas)) << endl;
	cout << "    mu = " << initial_flow_state->gas->mu
	     << " k = " << initial_flow_state->gas->k[0] << endl;
    }
} // end GasSlug constructor


GasSlug::GasSlug(const GasSlug& gs) // copy constructor
{
    dt = gs.dt;
    dt0 = gs.dt0;
    dt_allow = gs.dt_allow;
    cfl_target = gs.cfl_target;
    sim_time = gs.sim_time;
    cfl_min = gs.cfl_min;
    cfl_max = gs.cfl_max;
    cfl_tiny = gs.cfl_tiny;
    time_tiny = gs.time_tiny;
    residual = gs.residual;
    max_steps = gs.max_steps;
    Xorder = gs.Xorder;
    Torder = gs.Torder;
    hncell = gs.hncell;
    for ( size_t i = 0; i < gs.hxcell.size(); ++i ) hxcell.push_back(gs.hxcell[i]);
    viscous_effects = gs.viscous_effects;
    adiabatic = gs.adiabatic;
    initial_flow_state = new LFlowState(*(gs.initial_flow_state));
    adaptive = gs.adaptive;
    nxdim = gs.nxdim;
    nxmax = gs.nxmax;
    nnx = gs.nnx;
    nghost = gs.nghost;
    ixmin = gs.ixmin;
    ixmax = gs.ixmax;
    xbegin = gs.xbegin;
    xend = gs.xend;
    dxmin = gs.dxmin;
    dxmax = gs.dxmax;
    cluster_to_end_1 = gs.cluster_to_end_1;
    cluster_to_end_2 = gs.cluster_to_end_2;
    cluster_strength = gs.cluster_strength;
    left_bc_type = gs.left_bc_type;
    right_bc_type = gs.right_bc_type;
    left_slug_id = gs.left_slug_id;
    right_slug_id = gs.right_slug_id;
    left_slug_end_id = gs.left_slug_end_id;
    right_slug_end_id = gs.right_slug_end_id;
    left_piston_id = gs.left_piston_id;
    right_piston_id = gs.right_piston_id;
    left_diaphragm_id = gs.left_diaphragm_id;
    right_diaphragm_id = gs.right_diaphragm_id;
    set_left_end_ustar = gs.set_left_end_ustar;
    set_right_end_ustar = gs.set_right_end_ustar;
    left_ustar = gs.left_ustar;
    left_pstar = gs.left_pstar;
    right_pstar = gs.right_pstar;
    for ( int i = 0; i < nxdim; ++i ) {
	Cell.push_back(LCell(gs.Cell[i]));
    }
} // end GasSlug copy constructor


GasSlug::~GasSlug(void)
{
    Cell.clear();
    delete initial_flow_state;
    hxcell.clear();
}


int GasSlug::set_index_range(void)
// Set up min and max indices for convenience in later work.
// Active cells should then be addressible as 
// cell[ix], ixmin <= ix <= ixmax.
{
    ixmin = nghost;
    ixmax = ixmin + nnx - 1;
    return SUCCESS;
}


int GasSlug::read_state(FILE* infile)
// Read the flow solution for all faces and cells in a slug. 
{
#   define NCHAR 3200
    char line[NCHAR];
    int ix, isp, imode;
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();
    if (fgets(line, NCHAR, infile) == NULL) {
        printf("Empty flow field file.\n");
        return FAILURE;
    }
    sscanf(line, "%lf", &(sim_time));
    if (fgets(line, NCHAR, infile) == NULL) {
        printf("Empty flow field file.\n");
        return FAILURE;
    }
    sscanf(line, "%d %d %d", &ix, &isp, &imode);
    if (ix <= nxmax) {
        nnx = ix;
    } else {
        printf("Trying to read too many cells into this slug.\n");
        return FAILURE;
    }
    set_index_range();
    if ( isp != nsp ) {
        printf("Inconsistent number of species: expected %d, read %d\n", nsp, isp);
	return FAILURE;
    }
    if ( imode != nmodes ) {
        printf("Inconsistent number of energy modes: expected %d, read %d\n", nmodes, imode);
	return FAILURE;
    }
    // From here on, we hope that the solution file is ok.
    // Interfaces between cells.
    for ( ix = ixmin - 1; ix <= ixmax; ++ix ) {
	LCell* c = &( Cell[ix] );
	if (fgets(line, NCHAR, infile) == NULL) {
	    printf("Problem reading flow field file.\n");
	    return FAILURE;
	}
	if ( c->scan_iface_values_from_string(line) ) {
	    printf("IFace for Cell[%d] failed to read from line:\n", ix);
	    printf("%s\n", line);
	    return FAILURE;
	}
    } // end for ix
    // The cells with flow data.
    for ( ix = ixmin; ix <= ixmax; ++ix ) {
	LCell* c = &( Cell[ix] );
	if (fgets(line, NCHAR, infile) == NULL) {
	    printf("Problem reading flow field file.\n");
	    return FAILURE;
	}
	if ( c->scan_cell_values_from_string(line) ) {
	    printf("Cell[%d] failed to read from line:\n", ix);
	    printf("%s\n", line);
	    return FAILURE;
	}
	double f_sum = c->gas->massf[0];
	for ( isp = 1; isp < nsp; ++isp ) {
	    f_sum += c->gas->massf[isp];
	}
	if ( fabs(f_sum - 1.0) > 0.001 ) {
	    printf("Species don't sum correctly %g\n", f_sum );
	    for ( isp = 1; isp < nsp; ++isp ) {
		printf("    %d %e\n", isp, c->gas->massf[isp]);
	    }
	    return FAILURE;
	}
    } // end for ix
    return SUCCESS;
#   undef NCHAR
} // end read_state


int GasSlug::write_state(FILE* outfile)
// Write the flow solution (i.e. the interface positions and
// the primary variables at the cell centers) to a file.
{
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();
    fprintf(outfile, "%e  # begin slug data: sim_time\n", sim_time);
    fprintf(outfile, "%d %d %d # nnx, nsp\n", nnx, nsp, nmodes);
    // Interfaces between cells.
    for ( int ix = ixmin - 1; ix <= ixmax; ++ix ) {
        fprintf(outfile, "%s\n", Cell[ix].write_iface_values_to_string().c_str());
    }
    // The actual cells.
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        fprintf(outfile, "%s\n", Cell[ix].write_cell_values_to_string().c_str());
    }
    fflush(outfile);
    return SUCCESS;
} // end write_state()

int GasSlug::compute_areas(TubeModel *tube)
// Compute the tube cross-section area, local loss coefficient
// and specified wall temperature as a function of position.
// Also, compute the cell mid-points.
{
    for ( int ix = ixmin - 1; ix <= ixmax; ++ix ) {
        // Compute appropriate interval.
        int im1 = (int) ((Cell[ix].x - tube->x1) / tube->dx);
        int i = im1 + 1;
        // Now interpolate the area.
        if ( i <= 0 ) {
            Cell[ix].area = tube->area[0];
        } else if ( i >= tube->n ) {
            Cell[ix].area = tube->area[tube->n - 1];
        } else {
            double alpha = (Cell[ix].x - tube->x1 - im1 * tube->dx) / tube->dx;
            Cell[ix].area = tube->area[im1] * (1.0 - alpha) + tube->area[i] * alpha;
        }
    }
    // Compute cell volumes and cell midpoints 
    // from the current interface areas and positions.
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        Cell[ix].volume = 0.5*(Cell[ix].area + Cell[ix-1].area)
            * (Cell[ix].x - Cell[ix-1].x);
        Cell[ix].xmid = 0.5*(Cell[ix].x + Cell[ix-1].x);
    }
    // Interpolate the Loss coefficients and specified wall temperatures
    // at the cell midpoints.
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
	// Compute appropriate interval.
        int im1 = (int) ((Cell[ix].xmid - tube->x1) / tube->dx);
        int i = im1 + 1;
	// Now interpolate the loss coefficient.
        if ( i <= 0 ) {
            Cell[ix].K_over_L = tube->K_over_L[0];
            Cell[ix].T_Wall = tube->T_Wall[0];
        } else if ( i >= tube->n ) {
            Cell[ix].K_over_L = tube->K_over_L[tube->n - 1];
            Cell[ix].T_Wall = tube->T_Wall[tube->n - 1];
        } else {
            double alpha = (Cell[ix].xmid - tube->x1 - im1*tube->dx) / tube->dx;
            Cell[ix].K_over_L = tube->K_over_L[im1] * (1.0 - alpha)
                + tube->K_over_L[i] * alpha;
            Cell[ix].T_Wall = tube->T_Wall[im1] * (1.0 - alpha)
                + tube->T_Wall[i] * alpha;
        }
    }
    return SUCCESS;
} // end compute_areas()


double GasSlug::total_energy()
{
    double te = 0.0;
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
	double e = accumulate(Cell[ix].gas->e.begin(), Cell[ix].gas->e.end(), 0.0);
        te += Cell[ix].mass * (e + 0.5*Cell[ix].u*Cell[ix].u);
    }
    return te;
} // end total_energy()


int GasSlug::maximum_p(double *p_max, double *x_max)
// For the given slug of gas, find the maximum pressure and its location.
{
    double pp = 0.0;
    double xx = 0.0;
    int something_done = 0;
    for ( int ix = ixmin+1; ix <= ixmax-1; ++ix ) {
        double p = (Cell[ix-1].gas->p + Cell[ix].gas->p + Cell[ix+1].gas->p)/3.0;
        if ( p > pp ) {
            pp = p;
            xx = Cell[ix].xmid;
        }
        something_done = 1;
    }
    if ( !something_done ) {
        // Arbitrary choice as there are not enough cells, anyway.
        pp = Cell[ixmin].gas->p;
        xx = Cell[ixmin].xmid;
    }
    *p_max = pp;
    *x_max = xx;
    return SUCCESS;
} // end maximum_p()


int GasSlug::fill_data()
// Fill the flow field for this block with uniform (initial) conditions.
{
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        LCell* target = &(Cell[ix]);
        target->gas->copy_values_from(*(initial_flow_state->gas));
        target->u = initial_flow_state->u;
    }
    return SUCCESS;
} // end fill_data()

int GasSlug::encode_conserved()
// Compute the conserved quantities in each cell from the primary variables.
{
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
	Cell[ix].encode_conserved();
    }
    return SUCCESS;
} // end encode_conserved


int GasSlug::decode_conserved()
// Compute the primary variables from the conserved quantities in each cell.
{
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
	Cell[ix].decode_conserved();
    }
    return SUCCESS;
} // end function L_decode_conserved


int GasSlug::set_chemistry_timestep(double dt)
{
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        LCell* c = &( Cell[ix] );
	c->dt_chem = dt;
    }
    return SUCCESS;
}


int GasSlug::set_thermal_timestep(double dt)
{
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        LCell* c = &( Cell[ix] );
	c->dt_therm = dt;
    }
    return SUCCESS;
}


int GasSlug::chemical_increment(double dt)
// Use Rowan's chemistry module to update the species mass fractions.
{
    Gas_model *gmodel = get_gas_model_ptr();
    Reaction_update *rupdate = get_reaction_update_ptr();
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        LCell* c = &( Cell[ix] );
	int flag = rupdate->update_state(*(c->gas), dt, c->dt_chem, gmodel);
	if ( flag != 0 ) {
	    cout << "chemical update failed" << endl
		 << "    for cell " << ix << " at x=" << c->xmid << endl;
	}
	// The update only changes mass fractions, we need to impose
	// a thermodynamic constraint based on a call to the equation
	// of state.
	gmodel->eval_thermo_state_rhoe(*(c->gas));
	// Ensure viscous properties are up-to-date.
	gmodel->eval_transport_coefficients(*(c->gas), gmodel);
    } // end ix loop
    return SUCCESS;
} // end chemical_increment

/*
 * -----------------
 * Production Terms.
 * -----------------
 */

double f_darcy_weisbach(double Re, double D)
{
    // Darcy-Weisbach friction factor for fully-developed pipe flow.
    // Re is based on pipe diameter.
    double f, lgtmp;
    // size of wall roughness elements in metres.
    double eps = 0.025e-3;
    if (Re < 10.0) {
	// A reasonable limit for very low speeds.
	f = 6.4;
    } else if (Re < 2000.0) {
	// Laminar regime.
	f = 64.0 / Re;
    } else if (Re < 4000.0) {
	// Transition regime.
	f = 0.032 / pow(Re/2000.0, 0.3187);
    } else {
	// Fully turbulent
#       define  SMOOTH  1
#       if (SMOOTH == 1)
	lgtmp = 1.327359 - 0.9 * log10(Re);
	UNUSED_VARIABLE(D);
	UNUSED_VARIABLE(eps);
#       else
	lgtmp = log10(21.25 * pow(Re, -0.9) + eps / D);
#       endif
	f = 1.14 - 2.0 * lgtmp;
	f = 1.0 / (f * f);
    }
    return f;
}

double f_flat_plate( double Re )
{
    // Flat-plate friction factor with Re based on length along plate.
    // See Holman (1986) Heat Transfer, Eq 5.125 to 5.127
    double f;
    if (Re < 1.0e4) { 
	// somewhat arbitrary lower Re limit to stop 
	// possibility of near-infinite friction and heat flux
	f = 8 * 0.332 / sqrt(10000.0);
    } else if (Re < 5.0e5) { // Laminar regime
	f = 8 * 0.332 / sqrt(Re); 
    } else if (Re < 1.0e7) { // lower Re turbulent regime.
	f = 8 * 0.0296 * pow(Re, -0.2);
    } else { // high Re tubulent regime.
	f = 8 * 0.185 * pow(log10(Re), -2.584);
    }
    return f;
}


int GasSlug::source_vector()
// Compute the components of the source vector, Q.
// This vector is used to include viscous losses and
// heat transfer and mass sources/sinks.
{
    double lambda, f, Re_D, Re_L, D, abs_u, area, M;
    double T_wall, T_aw, F_LOSS, T_ref, St;
    double omega, tau0, h, length;
    double w_dot, q, gam;
    double sigma, alpha, Tgas, q_rad;
    double r2r1, m_loss, mom_loss, E_loss, kt, kl;

    Gas_model *gmodel = get_gas_model_ptr();
    double mypi = 4.0*atan(1.0);

    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
	LCell* cell = &(Cell[ix]);
	// Initialize the source terms to zero.
	// Each physical effect modelled will add or subtract 
	// something from these quantities.
	cell->Q_m = 0.0;
	cell->Q_mom = 0.0;
	cell->Q_E = 0.0;
	// The following variables are used to record the magnitude
	// of the viscous effects -- clear them in case they don't get
	// set on this pass..
	cell->shear_stress = 0.0;
	cell->heat_flux = 0.0;
	if ( viscous_effects == INVISCID ) {
	    // Nothing further to be done for this cell.
	    continue;
	}
	if ( L_get_case_id() == DRUMMOND_WITH_M4_NOZZLE ) {
	    // Treat nozzle flow as inviscid because 
	    // turbulent pipe flow model or other shocked gas models
	    // are unrealistic for losses in the (assumed steady) nozzle flow.
	    if (cell->x >= 0.049) { // start of nozzle contraction
		continue; // Nothing further to be done for this cell.
	    }
	}

	// Prepare to compute the viscous effects / losses.
	gam = gmodel->Cp(*(cell->gas)) / gmodel->Cv(*(cell->gas));
	// Local tube geometry
	area = 0.5 * (cell->area + Cell[ix-1].area);
	D = 2.0 * sqrt(area / mypi); // effective diameter
	length = cell->x - Cell[ix-1].x; // between interfaces
	// Flow state
	abs_u = fabs(cell->u); // Local gas speed
	double Prandtl = 0.75; // a constant value for the moment -- FIX-ME
	omega = pow(Prandtl, 0.333); // Recovery factor -- assume turbulent
	M = abs_u / cell->gas->a; // Local Mach number
#       define APPLY_COMPRESSIBILITY_FACTOR 1
	if ( APPLY_COMPRESSIBILITY_FACTOR == 1 ) {
	    lambda = 1.0 + (gam - 1.0) * 0.5 * omega * M * M;
	} else {
	    lambda = 1.0;
	}
	T_aw = lambda * cell->gas->T[0]; // Adiabatic wall temperature
	// Local wall temperature 
	if (adiabatic == 1) {
	    T_wall = T_aw;
	} else if (adiabatic == 2) {
	    // John Hunter's suggestion.  Heat transfer will tend to 
	    // increase the wall temperature to something close to the
	    // core temperature of the gas.
	    T_wall = cell->gas->T[0] - 400.0;
	    if (T_wall < cell->T_Wall) T_wall = cell->T_Wall;
	} else {
	    T_wall = cell->T_Wall;
	}
	// Transport properties based on Eckert reference conditions.
	T_ref = cell->gas->T[0] + 0.5 * (T_wall - cell->gas->T[0]) +
	    0.22 * (T_aw - cell->gas->T[0]);
	cell->ref->copy_values_from(*(cell->gas));
	cell->ref->T[0] = T_ref;
	cell->ref->rho = cell->gas->rho * cell->gas->T[0] / T_ref;
	cell->ref->p = cell->ref->rho * gmodel->R(*(cell->ref)) * cell->ref->T[0];
	gmodel->eval_thermo_state_pT(*(cell->ref));
	gmodel->eval_transport_coefficients(*(cell->ref), gmodel);
	// Local Reynolds number based on diameter and reference conditions.
	Re_D = cell->ref->rho * D * abs_u / cell->ref->mu;
	// use distance moved as reference distance for Re
	Re_L = cell->ref->rho * cell->L_bar * abs_u / cell->ref->mu;

	if ( ( viscous_effects == VISCOUS_LOSS_FLAT_PLATE_F || 
	       viscous_effects == DRB_MASS_LOSS_FLAT_PLATE_F )
	     && cell->L_bar > 0.0 ) {
	    // Friction and heat flux based on a flat plate calculation 
	    // No compressibility correction apart from
	    // property evaluation at Eckert reference conditions
	    // and adiabatic wall temperature based on recovery factor.
	    f = f_flat_plate( Re_L );
	} else {
	    // Default: friction factor determined from 
	    // fully-developed pipe flow correlation.
	    f = f_darcy_weisbach( Re_D, D ) / lambda;
	}

	if ( L_get_case_id() == DRUMMOND_WITH_M4_NOZZLE && sim_time < 0.02 ) {
	    // After shock-reflection, revert to the pipe-flow friction factor.
	    f = f_darcy_weisbach( Re_D, D ) / lambda;
	}

	// Local shear stress, in Pa, computed from the friction factor. 
	tau0 = -(0.125 * f) * cell->gas->rho * cell->u * abs_u;
	cell->shear_stress = tau0;

	// Rate of energy transfer into the cell.
	// First, shaft work (shear stress).
	// I think that this term should be zero not the following...
	// w_dot = tau0 * cell->u * length * mypi * D;
	w_dot = 0.0;
  
	// Convective heat transfer coefficient.
	St = (f * 0.125) * pow(Prandtl, -0.667);
	h = cell->ref->rho * gmodel->Cp(*(cell->gas)) * abs_u * St;

	// DJM special
	if ( viscous_effects == VISCOUS_LOSS_PIPE_F_HALF_H ) h *= 0.5;

	// Convective heat transfer from the wall into the gas cell.
	if ( adiabatic == 1 ) {
	    q = 0.0;
	} else {
	    q = h * mypi * D * length * (T_wall - T_aw);
	}

	// Compute the energy radiated by the gas cell to the tube wall.
#       define  RADIATION  0
	if (RADIATION == 1) {
	    Tgas = cell->gas->T[0];
	    alpha = 10.0;   /* 1/metres */
	    sigma = 5.6696e-8;  /* W/(m**2.K**4) */
	    q_rad = 2.0 * length * D * D * mypi *
		alpha * exp(-alpha * 0.5 * D) *
		sigma * Tgas * Tgas * Tgas * Tgas;
	} else {
	    q_rad = 0.0;
	}

	// Record the heat flux for the output file.
	cell->heat_flux = (q - q_rad) / (mypi * D * length);

	// Now, we actually apply the viscous effects to the 
	// conserved quantities in the cell.
	//
	// Decide how to remove the mass, momentum and energy
	// from the cell -- with, or without associated mass loss.
	//
	// DEFAULT: Momentum and Energy loss without associated mass loss.
	// Rate of change of Momentum by shear stress and 
	// energy loss by convective (and possibly radiative) heat transfer.
	m_loss   = 0.0;
	mom_loss = -(tau0 * length * mypi * D);
	E_loss   = -(w_dot + q - q_rad);

	if ( viscous_effects == CJD_MASS_LOSS_LAMINAR ||
	     viscous_effects == CJD_MASS_LOSS_TURBULENT ) {
	    // Con Doolan's mass-loss model for shocked gas.
	    // m_loss = mass loss per unit time
	    // This is based on the incompressible flat plate boundary layer
	    // equations and transformed to compresssible via the Howarth 
	    // transformation.
	    r2r1 = cell->gas->rho / initial_flow_state->gas->rho;
	    if (cell->L_bar > 0.0 && r2r1 > 1.001) {
		// Only formulated for strong shocks going forward.
		if ( viscous_effects == CJD_MASS_LOSS_TURBULENT ) {
		    kt = 0.375;
		    m_loss = 0.70 * mypi * D * length * cell->ref->rho * abs_u *
			kt / pow(Re_L, 1.0 / 5.0);
		    mom_loss = 0.889 * m_loss * cell->u;
		    E_loss = gmodel->Cp(*(cell->gas)) * cell->ref->T[0] * m_loss;
		} else { // Laminar
		    kl = 5.0;
		    m_loss = 0.333 * mypi * D * length * cell->ref->rho * 
			abs_u * kl / sqrt(Re_L);
		    mom_loss = m_loss * cell->u;
		    double e = accumulate(cell->gas->e.begin(),
					  cell->gas->e.end(), 0.0);
		    E_loss = m_loss * (e + abs_u * abs_u * 0.5);
		}
		if (m_loss < MINIMUM_MASS) { // Can't remember why.
		    E_loss   = 0.0;
		    mom_loss = 0.0;
		}
		// *** FIX ME ***
		// DRB had the following modificaton on the CJD mass loss: 
		// cell->Q_m -= m_loss * (-tau0 * length * mypi * D / mom_loss);
		// Is it important?
	    }
	} else if ( viscous_effects == DRB_MASS_LOSS_PIPE_F ||
		    viscous_effects == DRB_MASS_LOSS_FLAT_PLATE_F ) {
	    // Use default momentum and energy loss (via friction factor).
	    // Obtain mass loss from the momentum loss, noting that
	    // cell velocity is the momentum per unit mass of the core flow. 
	    m_loss = mom_loss / cell->u;
	}
	cell->Q_m   -= m_loss;
	cell->Q_mom -= mom_loss;
	cell->Q_E   -= E_loss;

	// Pipe fitting loss is a momentum loss on top of other viscous effects. 
	F_LOSS = -(cell->K_over_L) * 0.5 * cell->gas->rho * cell->u * abs_u;
	cell->Q_mom += area * length * F_LOSS;
        
    } // end for ix
    return SUCCESS;
} // end function L_source_vector


int GasSlug::axial_heat_flux(double k)
// Compute an "axial heat flux" to kill off glitches
// which seem to appear at the contact surfaces between gas slugs.
// Note that the form of the flux is quadratic in dT/dx rather than linear.
{
    // The thermal conductivity, k, is not intended to be accurate. 
    // For air at nominal conditions k = 0.024 SI-units.
    // A value of k = 0.01 seems to work OK for the Sod shock tube.
    for ( int ix = ixmin; ix <= ixmax - 1; ++ix ) {
        double dT = Cell[ix+1].gas->T[0] - Cell[ix].gas->T[0];
        double dx = Cell[ix+1].xmid - Cell[ix].xmid;
        // Try quadratic, instead of linear...
        Cell[ix].qstar = -k * dT / dx * fabs(dT / dx);
    }
    // No axial heat transfer to other gas slugs.
    Cell[ixmin-1].qstar = 0.0;
    Cell[ixmax].qstar = 0.0;
    return SUCCESS;
} // end function L_axial_heat_flux()


double min_increment(double y0, double y1, double y2)
// Computes the minimum increment to bring y0 into line
// with its neighbour values y1 and y2.
{
    /* Target values */
    double t1 = y1;               /* no extrapolation */
    double t2 = y1 - (y2 - y1);   /* linear extrapolation */
    /* Increments */
    double d1 = t1 - y0;
    double d2 = t2 - y0;
    /* Select the one with minimum magnitude. */
    if ( fabs(d1) < fabs(d2) ) {
	return d1;
    } else {
	return d2;
    }
}


int GasSlug::adjust_end_cells()
// Adjust the gas properties of the end-of-slug cells
// to kill off glitches which seem to appear at the 
// contact surfaces between gas slugs.
{
    LCell *c, *cn1, *cn2;
    double relax_factor;
    Gas_model *gmodel = get_gas_model_ptr();
    /* 
     * Select a value of the relaxation such that we don't do the 
     * damage all in one go.  
     * Combined with occasional use rather than application every step,
     * this should allow real waves to pass through.
     * However, it doesn't work really well for the start of
     * strong expansions where the pressure and temperature in
     * the first cell *should* drop away really fast.
     */
    relax_factor = 0.1;
    /* 
     * First, do left-end
     */
    c   = &(Cell[ixmin]);      /* end cell    */
    cn1 = &(Cell[ixmin + 1]);  /* neighbour 1 */
    cn2 = &(Cell[ixmin + 2]);  /* neighbour 2 */
    c->gas->T[0] += relax_factor * min_increment(c->gas->T[0], cn1->gas->T[0], cn2->gas->T[0]);
    c->gas->p += relax_factor * min_increment(c->gas->p, cn1->gas->p, cn2->gas->p);
    gmodel->eval_thermo_state_pT(*(c->gas));
    gmodel->eval_transport_coefficients(*(c->gas), gmodel);
    /* c->u += relax_factor * (cn->u - c->u); */
    c->encode_conserved();
    /* 
     * Second, do right-end
     */
    c   = &(Cell[ixmax]);      /* end cell    */
    cn1 = &(Cell[ixmax - 1]);  /* neighbour 1 */
    cn2 = &(Cell[ixmax - 2]);  /* neighbour 2 */
    c->gas->T[0] += relax_factor * min_increment(c->gas->T[0], cn1->gas->T[0], cn2->gas->T[0]);
    c->gas->p += relax_factor * min_increment(c->gas->p, cn1->gas->p, cn2->gas->p);
    gmodel->eval_thermo_state_pT(*(c->gas));
    gmodel->eval_transport_coefficients(*(c->gas), gmodel);
    /* c->u += relax_factor * (cn->u - c->u); */
    c->encode_conserved();

    return SUCCESS;
} // end adjust_end_cells()


//------------------------------------------------------------------------

#define MIN_MOD_LIMIT(c,d)   ( ((c)*(d) <= 0.0) ? 0.0 : SMALLEST(c,d) )


double limit_to_within(double v, double a, double b)
{
    double minimum_ab = MINIMUM(a,b);
    double maximum_ab = MAXIMUM(a,b);
    if ( v < minimum_ab ) v = minimum_ab;
    if ( v > maximum_ab ) v = maximum_ab;
    return v;
}

//---------------------------------------------------------------------------

std::vector<LFlowState> QL, QR;

int set_up_workspace_for_apply_rivp()
{
    printf( "set_up_workspace_forapply_rivp(): NDIM=%d.\n", NDIM);
    Gas_model *gmodel = get_gas_model_ptr();
    while ( QL.size() < NDIM ) {
    	    QL.push_back(LFlowState(gmodel));
    }
    while ( QR.size() < NDIM ) {
	QR.push_back(LFlowState(gmodel));
    }
    return SUCCESS;
} // end set_up_workspace_for_apply_rivp()

int dispose_workspace_for_apply_rivp()
{
    QL.clear();
    QR.clear();
    return SUCCESS;
}

//---------------------------------------------------------------------------

int GasSlug::apply_rivp()
{
    // Apply the Riemann solver to obtain the pressure and
    // velocity at each interface.
    int ix;
    static double del[NDIM], dplus[NDIM], dminus[NDIM];
    static double rhoL, rhoR, eL, eR;
    static double pstar[NDIM], ustar[NDIM];
    static double onedx[NDIM];
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();

    // On first encounter, the gas components inside the QL, QR
    // flow state structures will need to be created.

    // Interpolate the cell average values to obtain
    // LEFT and RIGHT states at each interface.
    //
    // The approach taken is to ignore grid distortions and perform
    // one dimensional projection/interpolation in the ix index-direction.
    // Left cell has index [ix], Right cell has index [ix+1]

    // Always start with first-order interpolation.
    for (ix = ixmin - 1; ix <= ixmax; ++ix) {
	QL[ix].gas->copy_values_from(*(Cell[ix].gas));
	QL[ix].u = Cell[ix].u;
	QR[ix].gas->copy_values_from(*(Cell[ix+1].gas));
	QR[ix].u = Cell[ix + 1].u;
    }

    if ( Xorder == 2 ) {
        // Higher-order interoplation of some quantities.
        // Assume 2 ghost points are available.
        for ( ix = ixmin - 1; ix <= ixmax + 2; ++ix ) {
            onedx[ix] = 1.0 / (Cell[ix].xmid - Cell[ix-1].xmid);
        }
        // Density.
        for ( ix = ixmin - 1; ix <= ixmax + 1; ++ix ) {
            dminus[ix] = (Cell[ix].gas->rho - Cell[ix-1].gas->rho) * onedx[ix];
            dplus[ix] = (Cell[ix+1].gas->rho - Cell[ix].gas->rho) * onedx[ix+1];
            del[ix] = MIN_MOD_LIMIT(dminus[ix], dplus[ix]);
        }
        for ( ix = ixmin-1; ix <= ixmax; ++ix ) {
            rhoL = Cell[ix].gas->rho + del[ix]*(Cell[ix].x - Cell[ix].xmid);
            rhoR = Cell[ix+1].gas->rho - del[ix+1]*(Cell[ix+1].xmid - Cell[ix].x);
            QL[ix].gas->rho = limit_to_within(rhoL, Cell[ix].gas->rho, Cell[ix+1].gas->rho);
            QR[ix].gas->rho = limit_to_within(rhoR, Cell[ix].gas->rho, Cell[ix+1].gas->rho);
        }
	// Individual species mass-fractions.
	for ( size_t isp = 0; isp < nsp; ++isp ) {
	    for ( ix = ixmin - 1; ix <= ixmax + 1; ++ix ) {
		dminus[ix] = (Cell[ix].gas->massf[isp] - Cell[ix-1].gas->massf[isp]) * onedx[ix];
		dplus[ix] = (Cell[ix+1].gas->massf[isp] - Cell[ix].gas->massf[isp]) * onedx[ix + 1];
		del[ix] = MIN_MOD_LIMIT(dminus[ix], dplus[ix]);
	    }
	    for ( ix = ixmin - 1; ix <= ixmax; ++ix ) {
		eL = Cell[ix].gas->massf[isp] + del[ix] * (Cell[ix].x - Cell[ix].xmid);
		eR = Cell[ix+1].gas->massf[isp] - del[ix+1] * (Cell[ix+1].xmid - Cell[ix].x);
		QL[ix].gas->massf[isp] = limit_to_within(eL, Cell[ix].gas->massf[isp],
						       Cell[ix+1].gas->massf[isp]);
		QR[ix].gas->massf[isp] = limit_to_within(eR, Cell[ix].gas->massf[isp],
						       Cell[ix+1].gas->massf[isp]);
	    }
	} // end for isp
        // Axial Velocity
        for ( ix = ixmin-1; ix <= ixmax+1; ++ix ) {
            dminus[ix] = (Cell[ix].u - Cell[ix-1].u) * onedx[ix];
            dplus[ix] = (Cell[ix+1].u - Cell[ix].u) * onedx[ix + 1];
            del[ix] = MIN_MOD_LIMIT(dminus[ix], dplus[ix]);
        }
        for ( ix = ixmin - 1; ix <= ixmax; ++ix ) {
            QL[ix].u = Cell[ix].u + del[ix] * (Cell[ix].x - Cell[ix].xmid);
            QR[ix].u = Cell[ix+1].u - del[ix+1] * (Cell[ix+1].xmid - Cell[ix].x);
        }
        // Internal Energy modes.
	for ( size_t imode = 0; imode < nmodes; ++imode ) {
	    for ( ix = ixmin - 1; ix <= ixmax + 1; ++ix ) {
		dminus[ix] = (Cell[ix].gas->e[imode] - Cell[ix-1].gas->e[imode]) * onedx[ix];
		dplus[ix] = (Cell[ix+1].gas->e[imode] - Cell[ix].gas->e[imode]) * onedx[ix + 1];
		del[ix] = MIN_MOD_LIMIT(dminus[ix], dplus[ix]);
	    }
	    for ( ix = ixmin - 1; ix <= ixmax; ++ix ) {
		eL = Cell[ix].gas->e[imode] + del[ix] * (Cell[ix].x - Cell[ix].xmid);
		eR = Cell[ix+1].gas->e[imode] - del[ix+1] * (Cell[ix+1].xmid - Cell[ix].x);
		QL[ix].gas->e[imode] = limit_to_within(eL, Cell[ix].gas->e[imode],
						       Cell[ix+1].gas->e[imode]);
		QR[ix].gas->e[imode] = limit_to_within(eR, Cell[ix].gas->e[imode],
						       Cell[ix+1].gas->e[imode]);
	    }
	} // end for imode
        // Pressure, Local Speed of Sound and Temperature.
        for ( ix = ixmin-1; ix <= ixmax; ++ix ) {
	    gmodel->eval_thermo_state_rhoe(*(QL[ix].gas));
	    gmodel->eval_thermo_state_rhoe(*(QR[ix].gas));
        }
    } // End of Higher-order interpolation.

    // *****************************
    // * Apply the Riemann solver. *
    // *****************************
    // Apply specified interface velocities if required.
    if (set_left_end_ustar == 1) {
        ustar[ixmin-1] = left_ustar;
    }
    if (set_right_end_ustar == 1) {
        ustar[ixmax] = right_ustar;
    }
    L_rivp(QL, QR, ustar, pstar, ixmin - 1, ixmax,
           set_left_end_ustar, set_right_end_ustar);
    // Save the interface pressures and velocities.
    for (ix = ixmin - 1; ix <= ixmax; ++ix) {
        Cell[ix].pface = pstar[ix];
        Cell[ix].uface = ustar[ix];
    }
    // Save the interface pressures at the ends for interaction with the pistons.
    left_pstar = pstar[ixmin-1];
    right_pstar = pstar[ixmax];

    return SUCCESS;
} // end apply_rivp
//---------------------------------------------------------------------------

int GasSlug::time_derivatives(int time_level)
// Compute the time derivatives for the conserved quantities
// for each active cell.  These are the spatial (RHS) terms 
// in the semi-discrete governing equations.
//
// Input...
// time_level : specifies where the computed derivatives are to be stored.
{
    int ix;
    LCell *C, *Cm1;
    // Interface motions.
    for (ix = ixmin - 1; ix <= ixmax; ++ix) {
        C = &(Cell[ix]);
        C->DxDt[time_level] = C->uface;
    }
    // Cell average properties.
    for (ix = ixmin; ix <= ixmax; ++ix) {
        C = &(Cell[ix]);
        Cm1 = &(Cell[ix - 1]);
        // Mass.
        C->DmDt[time_level] = C->Q_m;
        // Momentum.
        C->DmomDt[time_level] =
            Cm1->pface * Cm1->area - C->pface * C->area +
            C->gas->p * (C->area - Cm1->area) + C->Q_mom;
	// Energy.
        C->DEDt[time_level] =
            Cm1->pface * Cm1->area * Cm1->uface
            - C->pface * C->area * C->uface
            + Cm1->qstar * Cm1->area - C->qstar * C->area + C->Q_E;
	// Length scale for the application of Mirel's theory.
        C->DLDt[time_level] = fabs(C->u);
    } // end for (ix = ...
    return SUCCESS;
} // end function L_time_derivatives

//-----------------------------------------------------------------

int GasSlug::record_state()
// Record the slug state before attempting a time step.
{
    LCell* C;
    // Interfaces.
    for ( int ix = ixmin - 1; ix <= ixmax; ++ix ) {
        C = &(Cell[ix]);
        C->x_old = C->x;
    }
    // Conserved quantities within the cell.
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        C = &(Cell[ix]);
        C->mass_old = C->mass;
        C->moment_old = C->moment;
        C->Energy_old = C->Energy;
        C->L_bar_old = C->L_bar;
    }
    return SUCCESS;
} // end record_state


int GasSlug::restore_state()
// Restore the cell state to that before the attempted time step.
{
    LCell* C;
    // Interfaces.
    for ( int ix = ixmin - 1; ix <= ixmax; ++ix ) {
        C = &(Cell[ix]);
        C->x = C->x_old;
    }
    // Conserved quantities within the cell.
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        C = &(Cell[ix]);
        C->mass = C->mass_old;
        C->moment = C->moment_old;
        C->Energy = C->Energy_old;
        C->L_bar = C->L_bar_old;
    }
    return SUCCESS;
} // end restore_state


int GasSlug::predictor_step()
// Use the time derivatives to advance the conserved quantities forward
// by time step dt.
{
    LCell* C;
    // Interfaces.
    for ( int ix = ixmin - 1; ix <= ixmax; ++ix ) {
        C = &(Cell[ix]);
        C->x = C->x_old + dt * C->DxDt[0];
    }
    // Conserved quantities within the cell.
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        C = &(Cell[ix]);
        double mass = C->mass_old;
        double dm = dt * C->DmDt[0];
        if (mass + dm > MINIMUM_MASS)
            C->mass = mass + dm;
	else
	    C->mass = mass;
        C->moment = C->moment_old + dt * C->DmomDt[0];
        C->Energy = C->Energy_old + dt * C->DEDt[0];
        C->L_bar = C->L_bar_old + dt * C->DLDt[0];
    }
    return SUCCESS;
} // end predictor_step


int GasSlug::corrector_step()
// Use the time derivatives to advance the conserved quantities
// forward by time step dt.
{
    LCell* C;
    // Interfaces.
    for ( int ix = ixmin - 1; ix <= ixmax; ++ix ) {
        C = &(Cell[ix]);
        C->x = C->x_old + dt * 0.5 * (C->DxDt[1] + C->DxDt[0]);
    }
    // Conserved quantities within the cell.
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        C = &(Cell[ix]);
        double mass = Cell[ix].mass_old;
        double dm = dt * 0.5 * (C->DmDt[1] + C->DmDt[0]);
        if ( mass + dm > MINIMUM_MASS )
            C->mass = mass + dm;
	else
	    C->mass = mass;
        C->moment = C->moment_old + dt * 0.5 * (C->DmomDt[1] + C->DmomDt[0]);
        C->Energy = C->Energy_old + dt * 0.5 * (C->DEDt[1] + C->DEDt[0]);
        // length scale
        C->L_bar = C->L_bar_old + dt * 0.5 * (C->DLDt[1] + C->DLDt[0]);
    }
    return SUCCESS;
} // end corrector_step


int GasSlug::check_cells(int js)
// Check cells within a given slug for invalid data.
//
// Input...
// js     : index for the particular slug
//
// Returns the number of bad cells.
{
    int number_bad_cells = 0;
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        LCell* C = &(Cell[ix]);
        LCell* Cm1 = &(Cell[ix - 1]);

        int this_cell_bad = 0;
        double dx = C->x - Cm1->x;
        if (dx <= 0.0) this_cell_bad = 1;
        if (C->mass <= 0.0) this_cell_bad = 1;
        if (C->gas->T[0] <= 1.0) this_cell_bad = 1;
        if (C->gas->p <= 1.0e-6) this_cell_bad = 1;
        if (this_cell_bad == 1) {
            ++number_bad_cells;
            printf("Bad Cell: js=%d, cell=%d\n", js, ix - ixmin);
            printf("        : dx=%e, xL=%e, xR=%e\n", dx, Cm1->x, C->x);
            printf("        : mass=%e, moment=%e, E=%e\n",
                   C->mass, C->moment, C->Energy);
            printf("        : mass_old=%e, moment_old=%e, E_old=%e\n",
                   C->mass_old, C->moment_old, C->Energy_old);
            printf("        : uLface=%e, u=%e, uRface=%e\n",
                   Cm1->uface, C->u, C->uface);
            printf("        : pLface=%e, p=%e, pRface=%e\n",
                   Cm1->pface, C->gas->p, C->pface);
            printf("        : ALface=%e, volume=%e, ARface=%e\n",
                   Cm1->area, C->volume, C->area);
        } // end if
    } // end for
    return number_bad_cells;
} // end function L_check_bad_cells()


int GasSlug::check_cfl()
// Compute the local time step limit for each cell.
// The overall time step is limited by the worst-case cell.
{
    // Some Definitions...
    // ----------------
    // dt       : global time step for the block
    // cfl_target : desired CFL number
    // cfl_min  : approximate minimum CFL number in the block
    // cfl_max  : approximate maximum CFL number in the block
    // dt_allow : allowable time step (i.e. the maximum dt that satisfies the CFL target)
    double signal_time, dt_local, cfl_local;
    double uL_mag, uR_mag, u_mag, a_local, signal_speed, length;
    cfl_min = 1.0e6; /* something outrageously large */
    cfl_max = 0.0;
    dt_allow = 1.0e6;
    for ( int ix = ixmin; ix <= ixmax; ++ix ) {
        // With the Lagrangian formulation,
	// signal speed is essentially the speed of sound.
        a_local = Cell[ix].gas->a;
        signal_speed = a_local;
#       if 0
        /* Select the largest signal speed from 
         * the velocities, also. */
        u_mag = fabs(Cell[ix].u);
        uL_mag = fabs(Cell[ix-1].uface);
        uR_mag = fabs(Cell[ix].uface);
        if (signal_speed < uL_mag)
            signal_speed = uL_mag;
        if (signal_speed < uR_mag)
            signal_speed = uR_mag;
        if (signal_speed < u_mag)
            signal_speed = u_mag;
#       else
	UNUSED_VARIABLE(u_mag);
	UNUSED_VARIABLE(uL_mag);
	UNUSED_VARIABLE(uR_mag);
#       endif
        // Check the INVISCID time step limit.
        length = Cell[ix].x - Cell[ix-1].x;
        signal_time = length / signal_speed;
        // Current (Local) CFL number
        cfl_local = fabs(dt / signal_time);
	// Recommend a time.
        dt_local = cfl_target * signal_time;
	// Search for the worst case.
        if (cfl_local < cfl_min) cfl_min = cfl_local;
        if (cfl_local > cfl_max) cfl_max = cfl_local;
        if (dt_local < dt_allow) dt_allow = dt_local;
	// Some debug for problem situations ...
        if (cfl_local > 10.0) {
            printf("\n-----------------\n");
            printf("CFL = %e , cell[%d]\n", cfl_local, ix);
            printf("rho = %e, u = %e, a = %e, length = %e\n",
                   Cell[ix].gas->rho, Cell[ix].u, 
		   Cell[ix].gas->a, length);
        } // end if
    } // end for ix loop
    if (cfl_max > 0.9) {
        printf("WARNING: large CFL number, cfl_max = %e\n", cfl_max);
    }
    return SUCCESS;
} // end check_cfl


/// \brief Interpolate gas properties at a particular x-location.
///
/// If the x-location falls within the gas slug, the values 
/// are interpolated from the adjacent cells and the function returns 1.
/// If the location does not fall within the bounds of the 
/// gas slug, the values are set to zero and the function returns 0.
/// Look at the code to see what values correspond to what properties.
int GasSlug::interpolate_cell_data(double xloc, LCell& icell)
{
    double alpha, xmid_m1, xmid_p1, xmid, dx_plus, dx_minus;
    LCell* c0;
    LCell* c1;

    int found = 0;
    // We assume that the x-locations of the cells increase with index.
    if (xloc > Cell[ixmin-1].x && xloc <= Cell[ixmax].x) {
        // ...and then find the cell containing the x-location.
        for ( int ix = ixmin; ix <= ixmax; ++ix ) {
            if ( xloc <= Cell[ix].x ) {
		// We have found the cell 
		// Linearly interpolate the quantities.
                xmid = 0.5 * (Cell[ix-1].x + Cell[ix].x);
                if (ix > ixmin) {
                    xmid_m1 = 0.5 * (Cell[ix-2].x + Cell[ix-1].x);
                } else {
                    xmid_m1 = Cell[ixmin-1].x - (xmid - Cell[ixmin-1].x);
                }
                if (ix < ixmax) {
                    xmid_p1 = 0.5 * (Cell[ix].x + Cell[ix+1].x);
                } else {
                    xmid_p1 = Cell[ixmax].x + (Cell[ixmax].x - xmid);
                }
                dx_minus = xmid - xmid_m1;
                dx_plus = xmid_p1 - xmid;

                alpha = xloc - xmid;
                if (alpha >= 0.0) {
                    alpha /= dx_plus; c0 = &Cell[ix]; c1 = &Cell[ix+1];
                } else {
                    alpha /= -1.0 * dx_minus; c0 = &Cell[ix]; c1 = &Cell[ix-1];
                }
		icell.xmid = xloc;
		icell.volume = (1.0-alpha) * c0->volume + alpha * c1->volume;
		icell.gas->average_values_from(*(c0->gas), (1.0-alpha), *(c1->gas), alpha, 0);
		icell.L_bar = (1.0-alpha) * c0->L_bar + alpha * c1->L_bar;
		icell.u = (1.0-alpha) * c0->u + alpha * c1->u;
		if ( viscous_effects ) {
		    icell.shear_stress = (1.0-alpha) * c0->shear_stress + alpha * c1->shear_stress;
		    icell.heat_flux = (1.0-alpha) * c0->heat_flux + alpha * c1->heat_flux;
		} else {
		    icell.shear_stress = 0.0;
		    icell.heat_flux = 0.0;
		}
		icell.entropy = (1.0-alpha) * c0->entropy + alpha * c1->entropy;
                found = 1;
                break;
            }   /* end if ... */
        }   /* end for (ix... */
    }   /* end if ... */

    return found;
} // end function L_interpolate_cell_data


/// \brief Compute the average pressure in a region near the end of the gas slug.
///
/// \param A : pointer to the gas slug data structure.
/// \param which_end : integer to indicate which end of the gas slug we want
/// \param dx : distance over which we will sample the cells
/// \returns the average of the pressures within the cells.
double GasSlug::end_pressure(int which_end, double dx)
{
    int ix, n;
    LCell *c;
    double x, x0, p_avg;
    x0 = 0.0;
    p_avg = 0.0;
    /* gcc could not determine if the following code initialized x0 and p_avg */
    if ( which_end == LEFT ) {
        n = 0;
        for (ix = ixmin; ix <= ixmax; ++ix) {
            c = &(Cell[ix]);
            x = c->xmid;
            if ( n == 0 ) {
                x0 = x;
                p_avg = c->gas->p;
                n = 1;
            } else {
                if ( fabs(x - x0) > dx ) break;
                p_avg += c->gas->p;
                ++n;
            }
        }   /* end for */
    } else {
        /* Assume right-hand end. */
        n = 0;
        for (ix = ixmax; ix >= ixmin; --ix) {
            c = &(Cell[ix]);
            x = c->xmid;
            if ( n == 0 ) {
                x0 = x;
                p_avg = c->gas->p;
                n = 1;
            } else {
                if ( fabs(x - x0) > dx ) break;
                p_avg += c->gas->p;
                ++n;
            }
        }   /* end for */
    }   /* end if */
    p_avg /= n;
    return p_avg;
} // end end_pressure()

/// \brief Compute the average flow properties in a region near the end of the gas slug.
///
/// \param A : pointer to the gas slug data structure.
/// \param which_end : integer to indicate which end of the gas slug we want
/// \param dx : distance over which we will sample the cells (initially)
/// \param total_mass : total mass included in dx region
/// \param Q : store the average flow properties in here
/// \returns : success or failure
int GasSlug::end_properties(int which_end, double dx, double* total_mass, LFlowState& Q)
{
    int ix, n, i;
    LCell *c;
    double x, x0;
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    x0 = 0.0;
    Q.gas->p=0.0;
    Q.gas->T[0]=0.0;
    Q.u=0.0;
    /* gcc could not determine if the following code initialized x0 and p_avg */
    if ( which_end == LEFT ) {
        n = 0;
        for (ix = ixmin; ix <= ixmax; ++ix) {
            c = &(Cell[ix]);
            x = c->xmid;
            if ( n == 0 ) {
                x0 = x;
                Q.gas->p = c->gas->p;
		Q.gas->T[0] = c->gas->T[0];
		Q.u = c->u;
		*total_mass=c->mass;
		for (i=0; i<nsp; i++) {
		    Q.gas->massf[i]=c->gas->massf[i];
		}
                n = 1;
            } else {
                if ( fabs(x - x0) > dx ) break;
                Q.gas->p += c->gas->p;
		Q.gas->T[0] += c->gas->T[0];
		Q.u += c->u;
		*total_mass +=c->mass;
                ++n;
            }
        }   /* end for */
    } else {
        /* Assume right-hand end. */
        n = 0;
        for (ix = ixmax; ix >= ixmin; --ix) {
            c = &(Cell[ix]);
            x = c->xmid;
            if ( n == 0 ) {
                x0 = x;
                Q.gas->p = c->gas->p;
		Q.gas->T[0] = c->gas->T[0];
		Q.u = c->u;
		*total_mass=c->mass;
		for (i=0; i<nsp; i++) {
		    Q.gas->massf[i]=c->gas->massf[i];
		}
                n = 1;
            } else {
                if ( fabs(x - x0) > dx ) break;
                Q.gas->p += c->gas->p;
		Q.gas->T[0] += c->gas->T[0];
		Q.u += c->u;
		*total_mass +=c->mass;
                ++n;
            }
        }   /* end for */
    }   /* end if */
    Q.gas->p /= n;
    Q.gas->T[0] /= n;
    Q.u /= n;
    return SUCCESS;
} // end end_properties()

/** \file init.cxx
 * \ingroup eilmer3
 * \brief Initialisation routines for Multiple-Block Navier-Stokes code.
 *
 * \version 02-Mar-08 Elmer3 port from mbcns2
 * \version 30-Jun-08 Conversion to use of C++ strings and Rowan's ConfigParser.
 */

//-----------------------------------------------------------------

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc_defs.hh"
#include "bc.hh"
#include "init.hh"
#include "diffusion.hh"
#include "visc.hh"

using namespace std;

std::string tostring(int i)
{
    ostringstream ost;
    ost << i;
    return ost.str();
}

/*-----------------------------------------------------------------*/

/// \brief Read simulation config parameters from the INI file.
///
/// These are grouped into global configuration parameters and
/// those which describe the blocks for holding the flow data.
///
/// \param filename : name of the INI parameter file
/// \param master: flag to indicate that this process is master
///
int read_config_parameters(const string filename, int master)
{
    global_data &G = *get_global_data_ptr();
    int jb;

    // Default values for some configuration variables.
    G.dimensions = 2;
    set_radiation_flag(0);

    // --FIX ME--
    // Note to Andrew D. -- adjust the following lines as appropriate.
    // Filter application parameters.
    G.do_filter = 0;
    G.filter_tstart = 0.0;
    G.filter_tend = 0.0;
    G.filter_dt = 0.0;
    G.filter_next_time = 0.0;
    G.filter_alpha = 0.0;
    G.filter_npass = 0;

    // variables for Andrew's time averaging routine.
    G.nav = 0;
    G.tav_0 = 0.0;
    G.tav_f = 0.0;
    G.dtav = 0.0;
    G.tav = 0.0;
    if (G.do_tavg == 1) {
	G.tav = G.tav_0; /* time at which averaging is to commence */
	G.nav = 0;
    }

    // variables for profile recording.
    G.do_record = 0;
    G.block_record = 0;
    G.x_record = 0;
    G.n_record = 0;
    G.step_record = 0;
    G.step = 0;
    G.t_record = 0.0;
    G.t_start = 0.0;
    if ( G.do_record == 1 ) {
	if ( master ) printf( "\ndo_record = %d\n\n", G.do_record );
	G.n_record = 0;
       	G.t_start = 0;
    }

    // variables for time dependent profile input.
    G.do_vary_in = 0;
    G.t_old = 0.0;
    G.ta_old = 0.0;
    G.n_vary = 0;
    if ( G.do_vary_in == 1 ) {
	G.n_vary = 0;
	G.t_old = -100.0; /* something guaranteed to be smaller than the first time */
	G.ta_old = -101.0;
    }
    // error checking
    if ( G.do_vary_in == 1 && G.do_record == 1 ) {
	if ( master ) printf( "\n Cannot simultaneously write and read from profile file\n" );
	exit(BAD_INPUT_ERROR);
    }

    // these variables are all for special cases.
    G.dn_ruptured = 0;
    G.d2_ruptured = 0;
    G.secondary_diaphragm_ruptured = 0;
    /* defaults to no diaphragm (in block -5) */
    G.diaphragm_block = -1;
    G.diaphragm_rupture_time = 0.0;
    G.diaphragm_rupture_diameter = 0.0;
    G.drummond_progressive = 0;

    G.fixed_time_step = false; // Set to false as a default
    G.cfl_count = 10;
    G.print_count = 20;
    G.control_count = 10;

    // At the start of a fresh simulation,
    // we need to set a few items that will be updated later.
    G.sim_time = 0.0;   // Global simulation time.
    G.cfl_max = 1.0;    // Dummy value
    G.cfl_min = 1.0;    // Dummy value
    G.cfl_tiny = 1.0;   // Smallest CFL so far, dummy value
    G.time_tiny = 1.0e6;

    // Most configuration comes from the previously-generated INI file.
    ConfigParser dict = ConfigParser(filename);
    int i_value;
    string s_value, s_value2;
    double d_value;

    dict.parse_string("global_data", "title", G.title, "unknown");
    dict.parse_int("global_data", "dimensions", G.dimensions, 2);
    dict.parse_int("global_data", "control_count", G.control_count, 10);
    if ( get_verbose_flag() ) {
	cout << "title = " << G.title << endl;
	cout << "dimensions = " << G.dimensions << endl;
	cout << "control_count = " << G.control_count << endl;
    }

    dict.parse_string("global_data", "udf_file", G.udf_file, "");
    dict.parse_int("global_data", "udf_source_vector_flag", G.udf_source_vector_flag, 0);

    dict.parse_string("global_data", "gas_model_file", s_value, "gas-model.lua");
    Gas_model *gmodel = set_gas_model_ptr(create_gas_model(s_value));
    if ( get_verbose_flag() ) {
	cout << "gas_model_file = " << s_value << endl;
	cout << "nsp = " << gmodel->get_number_of_species() << endl;
	cout << "nmodes = " << gmodel->get_number_of_modes() << endl;
    }

    dict.parse_int("global_data", "reacting_flag", i_value, 0);
    set_reacting_flag( i_value );
    dict.parse_double("global_data", "reaction_time_start", G.reaction_time_start, 0.0);
    dict.parse_string("global_data", "reaction_update", s_value, "dummy_scheme");
    if( get_reacting_flag() ) set_reaction_update( s_value );
    if ( get_verbose_flag() ) {
	cout << "reacting_flag = " << get_reacting_flag() << endl;
	cout << "reaction_time_start = " << G.reaction_time_start << endl;
	cout << "reaction_update = " << s_value << endl;
    }

    dict.parse_int("global_data", "energy_exchange_flag", i_value, 0);
    set_energy_exchange_flag( i_value );
    dict.parse_string("global_data", "energy_exchange_update", 
		      s_value, "dummy_scheme");
    if( get_energy_exchange_flag() ) set_energy_exchange_update(s_value);
    if( get_verbose_flag() ) {
	cout << "energy_exchange_flag = " << get_energy_exchange_flag() << endl;
	cout << "energy_exchange_update = " << s_value << endl;
    }

    dict.parse_int("global_data", "mhd_flag", i_value, 0);
    set_mhd_flag( i_value );

    dict.parse_int("global_data", "radiation_flag", i_value, 0);
    set_radiation_flag( i_value );
    dict.parse_string("global_data", "radiation_input_file", s_value, "no_file");
    dict.parse_int("global_data", "radiation_update_frequency", i_value, 1);
    if( get_radiation_flag() ) {
    	set_radiation_transport_model( s_value );
    	set_radiation_update_frequency( i_value );
    }
    if ( get_verbose_flag() ) {
	cout << "radiation_flag = " << get_radiation_flag() << endl;
	cout << "radiation_input_file = " << s_value << endl;
	cout << "radiation_update_frequency = " << i_value << endl;
    }

    dict.parse_int("global_data", "viscous_flag", i_value, 0);
    set_viscous_flag( i_value );
    dict.parse_double("global_data", "viscous_delay", G.viscous_time_delay, 0.0);
    dict.parse_double("global_data", "viscous_factor_increment", d_value, 0.01);
    set_viscous_factor_increment( d_value );
    dict.parse_double("global_data", "max_mu_t_factor", G.max_mu_t_factor, 300.0);
    dict.parse_double("global_data", "transient_mu_t_factor", G.transient_mu_t_factor, 1.0);
    dict.parse_int("global_data", "axisymmetric_flag", i_value, 0);
    set_axisymmetric_flag( i_value );
    dict.parse_int("global_data", "turbulence_flag", i_value, 0);
    set_turbulence_flag( i_value );
    // By default, turn off all turbulence models.
    set_k_omega_flag(0);
    set_baldwin_lomax_flag(0);
    dict.parse_double("global_data", "turbulence_prandtl_number", d_value, 0.89);
    set_turbulence_prandtl_number(d_value);
    dict.parse_double("global_data", "turbulence_schmidt_number", d_value, 0.75);
    set_turbulence_schmidt_number(d_value);
    if ( get_turbulence_flag() ) {
	// decide which model
	dict.parse_string("global_data", "turbulence_model", s_value, "none");
	if ( s_value == "baldwin_lomax" ) {
	    set_baldwin_lomax_flag(1);
	} else if ( s_value == "k_omega" ) {
	    set_k_omega_flag(1);
	} else if ( s_value == "spalart_allmaras" ) {
	    cout << "Spalart-Allmaras turbulence model not available." << endl;
	    exit( NOT_IMPLEMENTED_ERROR );
	} else {
	    // assume no turbulence model
	}
    }
    dict.parse_int("global_data", "diffusion_flag", i_value, 0);
    set_diffusion_flag( i_value );
    if ( get_verbose_flag() ) {
	cout << "viscous_flag = " << get_viscous_flag() << endl;
	cout << "viscous_delay = " << G.viscous_time_delay << endl;
	cout << "viscous_factor_increment = " << get_viscous_factor_increment() << endl;
	cout << "max_mu_t_factor = " << G.max_mu_t_factor << endl;
	cout << "transient_mu_t_factor = " << G.transient_mu_t_factor << endl;
	cout << "axisymmetric_flag = " << get_axisymmetric_flag() << endl;
	cout << "turbulence_flag = " << get_turbulence_flag() << endl;
	cout << "k_omega_flag = " << get_k_omega_flag() << endl;
	cout << "baldwin_lomax_flag = " << get_baldwin_lomax_flag() << endl;
	cout << "turbulence_prandtl_number = " << get_turbulence_prandtl_number() << endl;
	cout << "turbulence_schmidt_number = " << get_turbulence_schmidt_number() << endl;
	cout << "diffusion_flag = " << get_diffusion_flag() << endl;
    }

    if( get_diffusion_flag() && ( gmodel->get_number_of_species() > 1 ) ) { 
 	dict.parse_string("global_data", "diffusion_model", s_value, "Stefan-Maxwell");
	set_diffusion_model(s_value);
	if( get_verbose_flag() ) {
	    cout << "diffusion_model = " << s_value << endl;
	}
    }


    dict.parse_int("global_data", "max_invalid_cells", G.max_invalid_cells, 10);
    dict.parse_int("global_data", "flux_calc", i_value, 0);
    set_flux_calculator( i_value );
    dict.parse_double("global_data", "compression_tolerance", d_value, -0.30);
    set_compression_tolerance(d_value);
    dict.parse_double("global_data", "shear_tolerance", d_value, 0.20);
    set_shear_tolerance(d_value);
    dict.parse_string("global_data", "interpolation_type", s_value, "rhoe");
    set_thermo_interpolator( s_value );
    dict.parse_int("global_data", "apply_limiter_flag", i_value, 1);
    set_apply_limiter_flag(i_value);
    dict.parse_int("global_data", "extrema_clipping_flag", i_value, 1);
    set_extrema_clipping_flag(i_value);
    if ( get_verbose_flag() ) {
	cout << "max_invalid_cells = " << G.max_invalid_cells << endl;
	cout << "flux_calc = " << get_flux_calculator() << endl;
	cout << "compression_tolerance = " << get_compression_tolerance() << endl;
	cout << "shear_tolerance = " << get_shear_tolerance() << endl;
	cout << "interpolation_type = " << s_value << endl;
	cout << "apply_limiter_flag = " << get_apply_limiter_flag() << endl;
	cout << "extrema_clipping_flag = " << get_extrema_clipping_flag() << endl;
    }

    dict.parse_int("global_data", "sequence_blocks", i_value, 0);
    G.sequence_blocks = (i_value == 1);
    if ( get_verbose_flag() ) {
	cout << "sequence_blocks = " << G.sequence_blocks << endl;
    }

    // Read a number of gas-states.
    dict.parse_int("global_data", "nflow", G.n_gas_state, 0);
    if ( get_verbose_flag() ) {
	cout << "nflow = " << G.n_gas_state << endl;
    }
    for ( int ig = 0; ig < G.n_gas_state; ++ig ) {
	G.gas_state.push_back(read_flow_condition_from_ini_dict(dict, ig, master));
	if ( get_verbose_flag() ) {
	    cout << "flow condition[" << ig << "]: " << *(G.gas_state[ig]) << endl;
	}
    }

    // Read the parameters for a number of blocks.
    dict.parse_int("global_data", "nblock", G.nblock, 0);
    if ( get_verbose_flag() ) {
	printf( "nblock = %d\n", G.nblock);
    }
    // We keep a record of all of the configuration data for all blocks
    // but, eventually, we may allocate the flow-field data for only a 
    // subset of these blocks. 
    G.bd.resize(G.nblock);

    // Number of pistons
    // FIX-ME code needs to be reworked...
    dict.parse_int("global_data", "npiston", G.npiston, 0);
#   ifdef _MPI
    if ( G.npiston > 0 ) {
	G.npiston = 0;
	cout << "Pistons cannot be used with MPI code, yet." << endl;
    }
#   endif
    G.pistons.resize(G.npiston);
    if ( get_verbose_flag() ) {
	printf( "npiston = %d\n", G.npiston);
    }
    for ( int jp = 0; jp < G.npiston; ++jp ) {
	string piston_label;
	bool piston_cvf, piston_pvf;
	double piston_D, piston_L, piston_m, piston_x0, piston_v0, piston_f;
	string section = "piston/" + tostring(jp);
	dict.parse_string(section, "label", piston_label, "");
	dict.parse_double(section, "D", piston_D, 0.0);
	dict.parse_double(section, "L", piston_L, 0.0);
	dict.parse_double(section, "m", piston_m, 0.0);
	dict.parse_double(section, "x0", piston_x0, 0.0);
	dict.parse_double(section, "v0", piston_v0, 0.0);
	dict.parse_double(section, "f", piston_f, 0.0);
	dict.parse_boolean(section, "const_v_flag", piston_cvf, false);
	dict.parse_boolean(section, "postv_v_flag", piston_pvf, false);
	if ( get_verbose_flag() ) {
	    cout << "piston/" << jp << ": label= " << piston_label << endl;
	    cout << "    L=" << piston_L << ", m=" << piston_m
		 << ", D=" << piston_D << ", x0=" << piston_x0
		 << ", v0=" << piston_v0
		 << ", f=" << piston_f << endl;
	}
	
	G.pistons[jp] = new Piston();
	
	// This variable is used as a default to set the
	// vanishing distance.  We hope nobody ever simulates
	// something of this dimension.
	static const double VERY_LARGE_X = 1.0e6; // m

	// one way of making defaults
	vector<double> bore_resistance_f, bore_resistance_x;
	bore_resistance_x.push_back(0.0);
	bore_resistance_f.push_back(piston_f);
	double rifling_twist = 0.0;
	double rog = 0.0;
	double vanish_at_x = VERY_LARGE_X;
	cout << "Attempting to set piston values from config file..." << endl;
	int status = G.pistons[jp]->set_values(jp, piston_D, piston_L, piston_m, 
					       piston_x0, piston_v0, 
					       rifling_twist, rog, vanish_at_x, 
					       bore_resistance_f, bore_resistance_x);
	G.pistons[jp]->set_const_v_flag(piston_cvf);
	G.pistons[jp]->set_postv_v_flag(piston_pvf);
	if (status != SUCCESS) {
	    cout << "Failure setting piston parameters." << endl;
	    return FAILURE;
	} else {
	    cout << "Success setting piston parameters." << endl;
	}
    }

    dict.parse_int("global_data", "nheatzone", G.n_heat_zone, 0);
    dict.parse_double("global_data", "heat_time_start", G.heat_time_start, 0.0);
    dict.parse_double("global_data", "heat_time_stop", G.heat_time_stop, 0.0);
    dict.parse_double("global_data", "heat_factor_increment", d_value, 0.01);
    set_heat_factor_increment( d_value );
    if ( get_verbose_flag() ) {
	printf("nheatzone = %d\n", G.n_heat_zone);
	printf("heat_time_start = %e\n", G.heat_time_start);
	printf("heat_time_stop = %e\n", G.heat_time_stop);
	printf("heat_factor_increment = %e\n", get_heat_factor_increment() );
    }
    G.heat_zone.resize(G.n_heat_zone);
    for ( int indx = 0; indx < G.n_heat_zone; ++indx ) {
	struct CHeatZone* hzp = &(G.heat_zone[indx]);
	string section = "heat_zone/" + tostring(indx);
	dict.parse_double(section, "qdot", hzp->qdot, 0.0);
	dict.parse_double(section, "x0", hzp->x0, 0.0);
	dict.parse_double(section, "y0", hzp->y0, 0.0);
	dict.parse_double(section, "z0", hzp->z0, 0.0);
	dict.parse_double(section, "x1", hzp->x1, 0.0);
	dict.parse_double(section, "y1", hzp->y1, 0.0);
	dict.parse_double(section, "z1", hzp->z1, 0.0);
	if ( get_verbose_flag() ) {
	    cout << "heat_zone/" << indx << " qdot= " << hzp->qdot << endl;
	    cout << "    point0= " << hzp->x0 << " " << hzp->y0 << " " << hzp->z0 << endl;
	    cout << "    point1= " << hzp->x1 << " " << hzp->y1 << " " << hzp->z1 << endl;
	}
    }

    dict.parse_int("global_data", "nreactionzone", G.n_reaction_zone, 0);
    if ( get_verbose_flag() ) {
	printf("nreactionzone = %d\n", G.n_reaction_zone);
    }
    G.reaction_zone.resize(G.n_reaction_zone);
    for ( int indx = 0; indx < G.n_reaction_zone; ++indx ) {
	struct CReactionZone* rzp = &(G.reaction_zone[indx]);
	string section = "reaction_zone/" + tostring(indx);
	dict.parse_double(section, "x0", rzp->x0, 0.0);
	dict.parse_double(section, "y0", rzp->y0, 0.0);
	dict.parse_double(section, "z0", rzp->z0, 0.0);
	dict.parse_double(section, "x1", rzp->x1, 0.0);
	dict.parse_double(section, "y1", rzp->y1, 0.0);
	dict.parse_double(section, "z1", rzp->z1, 0.0);
	if ( get_verbose_flag() ) {
	    cout << "reaction_zone/" << indx << endl;
	    cout << "    point0= " << rzp->x0 << " " << rzp->y0 << " " << rzp->z0 << endl;
	    cout << "    point1= " << rzp->x1 << " " << rzp->y1 << " " << rzp->z1 << endl;
	}
    }

    dict.parse_int("global_data", "nturbulencezone", G.n_turbulent_zone, 0);
    if ( get_verbose_flag() ) {
	printf("nturbulencezone = %d\n", G.n_turbulent_zone);
    }
    G.turbulent_zone.resize(G.n_turbulent_zone);
    for ( int indx = 0; indx < G.n_turbulent_zone; ++indx ) {
	struct CTurbulentZone* tzp = &(G.turbulent_zone[indx]);
	string section = "turbulence_zone/" + tostring(indx);
	dict.parse_double(section, "x0", tzp->x0, 0.0);
	dict.parse_double(section, "y0", tzp->y0, 0.0);
	dict.parse_double(section, "z0", tzp->z0, 0.0);
	dict.parse_double(section, "x1", tzp->x1, 0.0);
	dict.parse_double(section, "y1", tzp->y1, 0.0);
	dict.parse_double(section, "z1", tzp->z1, 0.0);
	if ( get_verbose_flag() ) {
	    cout << "turbulence_zone/" << indx << endl;
	    cout << "    point0= " << tzp->x0 << " " << tzp->y0 << " " << tzp->z0 << endl;
	    cout << "    point1= " << tzp->x1 << " " << tzp->y1 << " " << tzp->z1 << endl;
	}
    }

    // Now, for the individual block configuration.
    for ( jb = 0; jb < G.nblock; ++jb ) {
        set_block_parameters( jb, dict, master );
    }
    check_connectivity();

    return SUCCESS;
} // end read_config_parameters()


int read_control_parameters( const string filename, int master, int first_time )
{
    int i_value;
    global_data &G = *get_global_data_ptr();
    // Parse the previously-generated INI file.
    ConfigParser dict = ConfigParser(filename);

    dict.parse_int("control_data", "x_order", i_value, 0);
    set_Xorder_flag( i_value );
    dict.parse_int("control_data", "t_order", i_value, 0);
    set_Torder_flag( i_value );
    dict.parse_double("control_data", "dt", G.dt_init, 1.0e-6);
    if ( first_time ) G.dt_global = G.dt_init;
    dict.parse_boolean("control_data", "fixed_time_step", G.fixed_time_step, 0);
    dict.parse_double("control_data", "dt_reduction_factor", G.dt_reduction_factor, 0.2);
    dict.parse_double("control_data", "cfl", G.cfl_target, 0.5);
    dict.parse_int("control_data", "stringent_cfl", i_value, 0);
    set_stringent_cfl_flag( i_value );
    dict.parse_int("control_data", "print_count", G.print_count, 20);
    dict.parse_int("control_data", "cfl_count", G.cfl_count, 10);
    dict.parse_double("control_data", "dt_plot", G.dt_plot, 1.0e-3);
    dict.parse_double("control_data", "dt_history", G.dt_his, 1.0e-3);
    dict.parse_double("control_data", "dt_fstc", G.dt_fstc, 1.0e-3);
    dict.parse_double("control_data", "max_time", G.max_time, 1.0e-3);
    dict.parse_int("control_data", "max_step", G.max_step, 10);
    dict.parse_int("control_data", "halt_now", G.halt_now, 0);
    dict.parse_int("control_data", "implicit_flag", i_value, 0);
    set_implicit_flag( i_value );
    dict.parse_int("control_data", "radiation_update_frequency", i_value, 1);
    set_radiation_update_frequency(i_value);
    if ( first_time && get_verbose_flag() ) {
	cout << "Time-step control parameters:" << endl;
	cout << "    x_order = " << get_Xorder_flag() << endl;
	cout << "    t_order = " << get_Torder_flag() << endl;
	cout << "    dt = " << G.dt_init << endl;
	cout << "    fixed_time_step = " << G.fixed_time_step << endl;
	cout << "    dt_reduction_factor = " << G.dt_reduction_factor << endl;
	cout << "    cfl = " << G.cfl_target << endl;
	cout << "    stringent_cfl = " << get_stringent_cfl_flag() << endl;
	cout << "    print_count = " << G.print_count << endl;
	cout << "    cfl_count = " << G.cfl_count << endl;
	cout << "    dt_plot = " << G.dt_plot << endl;
	cout << "    dt_history = " << G.dt_his << endl;
        cout << "    dt_fstc = " << G.dt_fstc << endl;
	cout << "    max_time = " << G.max_time << endl;
	cout << "    max_step = " << G.max_step << endl;
	cout << "    halt_now = " << G.halt_now << endl;
	cout << "    halt_now = " << G.halt_now << endl;
	cout << "    radiation_update_frequency = " << get_radiation_update_frequency();
	if (get_implicit_flag() == 0) {
	    cout << " (Explicit viscous advancements)" << endl;
	}
	else if (get_implicit_flag() == 1) {
	    cout << " (Point implicit viscous advancements)" << endl;
	}
	else if (get_implicit_flag() == 2) {
	    cout << " (Fully implicit viscous advancements)" << endl;
	}
	else {
	    cout << "\nInvalid implicit flag value: " << get_implicit_flag() << endl;
	    exit( BAD_INPUT_ERROR );
	}
    }
    return SUCCESS;
} // end read_control_parameters()


/// \brief Read simulation config parameters from the INI file.
///
/// At this point, we know the number of blocks in the calculation.
/// Depending on whether we are running all blocks in the one process
/// or we are running a subset of blocks in this process, talking to
/// the other processes via MPI, we need to decide what blocks belong
/// to the current process.
///
/// \param filename : name of the INI file containing the mapping
/// \param master: flag to indicate that this process is master
///
int assign_blocks_to_mpi_rank(const string filename, int master)
{
    global_data &G = *get_global_data_ptr();
    if ( get_verbose_flag() && master ) printf("Assign blocks to processes:\n");
    if ( G.mpi_parallel ) {
	if ( filename.size() > 0 ) {
	    if ( get_verbose_flag() && master ) {
		printf("    MPI parallel, mpimap filename=%s\n", filename.c_str());
		printf("    Assigning specific blocks to specific MPI processes.\n");
	    }
	    G.mpi_rank_for_block.resize(G.nblock);
	    // The mapping comes from the previously-generated INI file.
	    ConfigParser dict = ConfigParser(filename);
	    int nrank = 0;
	    int nblock;
	    int nblock_total = 0;
	    std::vector<int> block_ids, dummy_block_ids;
	    dict.parse_int("global", "nrank", nrank, 0);
	    if ( G.num_mpi_proc != nrank ) {
		if ( master ) {
		    printf("    Error in specifying mpirun -np\n");
		    printf("    It needs to match number of nrank; present values are:\n");
		    printf("    num_mpi_proc= %d nrank= %d\n", G.num_mpi_proc, nrank);
		}
		return FAILURE;
	    }
	    for ( int rank=0; rank < nrank; ++rank ) {
		string section = "rank/" + tostring(rank);
		dict.parse_int(section, "nblock", nblock, 0);
		block_ids.resize(0);
		dummy_block_ids.resize(nblock);
		for ( int i = 0; i < nblock; ++i ) dummy_block_ids[i] = -1;
		dict.parse_vector_of_ints(section, "blocks", block_ids, dummy_block_ids);
		if ( nblock != block_ids.size() ) {
		    if ( master ) {
			printf("    Did not pick up correct number of block_ids:\n");
			printf("        rank=%d, nblock=%d, block_ids.size()=%d\n",
			       rank, nblock, block_ids.size());
		    }
		    return FAILURE;
		}
		for ( int i = 0; i < nblock; ++i ) {
		    int this_block_id = block_ids[i];
		    if ( this_block_id < 0 ) {
			if ( master ) printf("    Error, invalid block id: %d\n", this_block_id);
			return FAILURE;
		    }
		    if ( G.my_mpi_rank == rank ) G.my_blocks.push_back(&(G.bd[this_block_id]));
		    G.mpi_rank_for_block[this_block_id] = rank;
		    nblock_total += 1;
		} // end for i
	    } // end for rank
	    if ( get_verbose_flag() ) {
		printf("    my_rank=%d, block_ids=", G.my_mpi_rank);
		for ( int i=0; i < G.my_blocks.size(); ++i ) {
		    printf(" %d", G.my_blocks[i]->id);
		}
		printf("\n");
	    }
	    if ( master ) {
		if ( nblock_total != G.nblock ) {
		    printf("    Error, total number of blocks incorrect: total=%d G.nblock=%d\n",
			   nblock_total, G.nblock);
		    return FAILURE;
		}
	    }
	} else {
	    if ( get_verbose_flag() && master ) {
		printf("    MPI parallel, No MPI map file specified.\n");
		printf("    Identify each block with the corresponding MPI rank.\n");
	    }
	    if ( G.num_mpi_proc != G.nblock ) {
		if ( master ) {
		    printf("    Error in specifying mpirun -np\n");
		    printf("    It needs to match number of blocks; present values are:\n");
		    printf("    num_mpi_proc= %d nblock= %d\n", G.num_mpi_proc, G.nblock);
		}
		return FAILURE;
	    }
	    G.my_blocks.push_back(&(G.bd[G.my_mpi_rank]));
	    for ( int jb=0; jb < G.nblock; ++jb ) {
		G.mpi_rank_for_block.push_back(jb);
	    }
	}
    } else {
	if ( get_verbose_flag() ) {
	    printf("    Since we are not doing MPI, all blocks in same process.\n");
	}
	for ( int jb=0; jb < G.nblock; ++jb ) {
	    G.my_blocks.push_back(&(G.bd[jb]));
	    G.mpi_rank_for_block.push_back(G.my_mpi_rank);
	} 
    } // endif
    return SUCCESS;
} // end assign_blocks_to_mpi_rank()


/** \brief Use config_parser functions to read the flow condition data. 
 */
CFlowCondition *read_flow_condition_from_ini_dict(ConfigParser &dict, int indx, int master)
{
    double p, u, v, w, Bx, By, Bz, mu_t, k_t, tke, omega, sigma_T, sigma_c;
    string value_string, flow_label;
    std::vector<double> massf, T, vnf;
    CFlowCondition *cfcp;
    Gas_model *gmodel = get_gas_model_ptr();

    string section = "flow/" + tostring(indx);
    dict.parse_string(section, "label", flow_label, "");
    dict.parse_double(section, "p", p, 100.0e3);
    dict.parse_double(section, "u", u, 0.0);
    dict.parse_double(section, "v", v, 0.0);
    dict.parse_double(section, "w", w, 0.0);
    dict.parse_double(section, "Bx", Bx, 0.0);
    dict.parse_double(section, "By", By, 0.0);
    dict.parse_double(section, "Bz", Bz, 0.0);
    dict.parse_double(section, "mu_t", mu_t, 0.0);
    dict.parse_double(section, "k_t", k_t, 0.0);
    dict.parse_double(section, "tke", tke, 0.0);
    dict.parse_double(section, "omega", omega, 1.0);
    dict.parse_double(section, "sigma_T", sigma_T, 0.0);
    dict.parse_double(section, "sigma_c", sigma_c, 1.0);
    int nsp = gmodel->get_number_of_species();
    vnf.resize(nsp);
    for ( int isp = 0; isp < nsp; ++isp ) vnf[isp] = 0.0;
    dict.parse_vector_of_doubles(section, "massf", massf, vnf);
    int nmodes = gmodel->get_number_of_modes();
    vnf.resize(nmodes);
    for ( int imode = 0; imode < nmodes; ++imode ) vnf[imode] = 300.0; 
    dict.parse_vector_of_doubles(section, "T", T, vnf);
    int S = 0;  // shock indicator
    cfcp = new CFlowCondition( gmodel, p, u, v, w, T, massf, flow_label, 
			       tke, omega, mu_t, k_t, S, Bx, By, Bz);
    vnf.resize(0);
    massf.resize(0);
    T.resize(0);
    flow_label.clear();
    return cfcp;
} // end  read_flow_condition_from_ini_dict()


/** \brief Set the parameters for the current block.
 *
 * Copy some of the global parameter values into the block and
 * then read a few parameters for the current block from parameter file.
 * Finally, and set a few more useful parameters from these.
 *
 * Also, Label the block with an integer (usually block number)
 * so that log-file messages can be labelled and special code
 * can be implemented for particular blocks.
 *
 * \param id   : (integer) block identity
 * \param dict : the dictionary associated with the INI-file
 * \param master : int flag to indicate whether this block is
 *                 associated with the master (0) process.
 *
 * \version 07-Jul-97 : Added iturb flag (in the place of nghost)
 * \version 22-Sep-97 : Updated the input format.
 * \version 10-May-2006 : Read from INI dictionary.
 * \version 20-Aug-2006 : Read in wall catalytic b.c.
 * \version 08-Sep-2006 : absorbed impose_global_parameters()
 */
int set_block_parameters(int id, ConfigParser &dict, int master)
{
    global_data &G = *get_global_data_ptr();
    Block &bdp = *get_block_data_ptr(id);
    int indx, iface, other_block, other_face, neighbour_orientation;
    int wc_bc, x_order, sponge_flag, xforce_flag, n_profile;
    string value_string, block_label, filename, wcbc_fname;
    int inflow_condition_id, bc_type_code, is_wall, use_udf_flux, assume_ideal;
    double Twall, Pout, epsilon;
    std::vector<double> f_wall_va;
    std::vector<double> f_wall, vnf;
    std::vector<double> mdot;
    double Twall_i, Twall_f, t_i, t_f;
    vnf.resize(get_gas_model_ptr()->get_number_of_species(), 0.0);

    bdp.id = id;
  
    // Read parameters from INI-file dictionary.
    indx = id;
    string section = "block/" + tostring(indx);

    dict.parse_string(section, "label", block_label, "");
    // Assume all blocks are active. 
    // The active flag will be used to skip over inactive
    // or unused blocks in later sections of code.
    dict.parse_int(section, "active", bdp.active, 1);
    if ( get_verbose_flag() ) {
	cout << section << ":label = " << block_label << endl;
	cout << "    active = " << bdp.active << endl;
    }

    // Number of active cells in block.
    dict.parse_int(section, "nni", bdp.nni, 2);
    dict.parse_int(section, "nnj", bdp.nnj, 2);
    // Allow for ghost cells
    bdp.nidim = bdp.nni + 2 * NGHOST;
    bdp.njdim = bdp.nnj + 2 * NGHOST;
    // Set up min and max indices for convenience in later work.
    // Active cells should then be addressible as
    // get_cell(i,j), imin <= i <= imax, jmin <= j <= jmax.
    bdp.imin = NGHOST;
    bdp.imax = bdp.imin + bdp.nni - 1;
    bdp.jmin = NGHOST;
    bdp.jmax = bdp.jmin + bdp.nnj - 1;
    if ( G.dimensions == 3 ) {
	dict.parse_int(section, "nnk", bdp.nnk, 2);
	bdp.nkdim = bdp.nnk + 2 * NGHOST;
	bdp.kmin = NGHOST;
	bdp.kmax = bdp.kmin + bdp.nnk - 1;
    } else {
	// For purely 2D flow geometry, we keep only one layer of cells.
	bdp.nnk = 1;
	bdp.nkdim = 1;
	bdp.kmin = 0;
	bdp.kmax = 0;
    }
    if ( get_verbose_flag() ) {
	printf( "    nni = %d, nnj = %d, nnk = %d\n", bdp.nni, bdp.nnj, bdp.nnk );
	printf( "    nidim = %d, njdim = %d, nkdim = %d\n", bdp.nidim, bdp.njdim, bdp.nkdim );
    }

    // Rotating frame of reference.
    dict.parse_double(section, "omegaz", bdp.omegaz, 0.0);

    // Boundary condition flags, 
    for ( iface = NORTH; iface <= ((G.dimensions == 3)? BOTTOM : WEST); ++iface ) {
	section = "block/" + tostring(indx) + "/face/" + get_face_name(iface);
	dict.parse_int(section, "bc", bc_type_code, SLIP_WALL);
	dict.parse_int(section, "inflow_condition", inflow_condition_id, 0);
	dict.parse_string(section, "filename", filename, "");
	dict.parse_int(section, "n_profile", n_profile, 1);
	dict.parse_int(section, "x_order", x_order, 0);
	dict.parse_int(section, "sponge_flag", sponge_flag, 0);
	dict.parse_int(section, "xforce_flag", xforce_flag, 0);
	dict.parse_double(section, "Twall", Twall, 300.0);
	dict.parse_double(section, "Pout", Pout, 100.0e3);
	dict.parse_int(section, "is_wall", is_wall, 0);
	dict.parse_int(section, "use_udf_flux", use_udf_flux, 0);
	dict.parse_int(section, "assume_ideal", assume_ideal, 0);
	dict.parse_vector_of_doubles(section, "mdot", mdot, vnf);
	dict.parse_double(section, "epsilon", epsilon, 0.9);
	dict.parse_double(section, "Twall_i", Twall_i, 300.0);
	dict.parse_double(section, "Twall_f", Twall_f, 300.0);
	dict.parse_double(section, "t_i", t_i, 0.0);
	dict.parse_double(section, "t_f", t_f, 0.0);
	dict.parse_int(section, "other_block", other_block, -1);
	dict.parse_string(section, "other_face", value_string, "");
	other_face = get_face_index(value_string);
	dict.parse_int(section, "neighbour_orientation", neighbour_orientation, 0);

	// Wall catalytic boundary condition flags, 
	// f_wall for SUPER_CATALYTIC 
	// wcbc_input_file for CATALYTIC and PARTIALLY_CATALYTIC
	//
	// Rowan or Dan **** FIX-ME *****
	dict.parse_int(section, "wcbc", wc_bc, NON_CATALYTIC);
	dict.parse_string(section, "wcbc_input_file", wcbc_fname, "");
	dict.parse_vector_of_doubles(section, "f_wall", f_wall, vnf);

	bdp.bcp[iface] = create_BC( bdp, iface, bc_type_code, inflow_condition_id, 
				    filename, n_profile, Twall, Pout,
				    x_order, is_wall, use_udf_flux,
				    other_block, other_face, neighbour_orientation,
				    sponge_flag, xforce_flag, mdot, epsilon,
				    wc_bc, wcbc_fname, f_wall,
				    Twall_i, Twall_f, t_i, t_f, assume_ideal);
	if ( get_verbose_flag() ) {
	    cout << "    " << get_face_name(iface) << " face:" << endl;
	    bdp.bcp[iface]->print_info("        ");
	    cout << "        Pout= " << Pout << endl;
	    cout << "        Twall= " << Twall << endl;
	}
    } // end for iface

    // History Cells.
    section = "block/" + tostring(indx);
    dict.parse_int(section, "nhcell", bdp.hncell, 0);
    bdp.hncell = MINIMUM(bdp.hncell, MAX_HNCELL); // limit to array size
    for ( int ih = 0; ih < bdp.hncell; ++ih ) {
	section = "block/" + tostring(indx);
	string key = "history-cell-" + tostring(ih);
	dict.parse_string(section, key, value_string, "0 0 0");
	if ( G.dimensions == 3 ) {
	    sscanf( value_string.c_str(), "%d %d %d", &(bdp.hicell[ih]), &(bdp.hjcell[ih]),
		    &(bdp.hkcell[ih]) );
	} else {
	    sscanf( value_string.c_str(), "%d %d", &(bdp.hicell[ih]), &(bdp.hjcell[ih]) );
	    bdp.hkcell[ih] = 0;
	}
	if ( get_verbose_flag() ) {
	    printf( "    History cell[%d] located at indices [%d][%d][%d]\n",
		    ih, bdp.hicell[ih], bdp.hjcell[ih], bdp.hkcell[ih] );
	}
    }

    return SUCCESS;
} // end set_block_parameters()

//------------------------------------------------------------------

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
#include <stdexcept>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/string_util.hh"
#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "block.hh"
#include "kernel.hh"
#include "cell.hh"
#include "bc.hh"
#include "init.hh"
#include "diffusion.hh"
#include "visc.hh"
#include "flux_calc.hh"
#include "one_d_interp.hh"

using namespace std;

std::map<std::string,update_scheme_t> available_schemes;
int init_available_schemes_map()
{
    // Yes, this has quite a few entries and they're here because
    // I've already made a couple of errors in the input scripts.
    // And, yes, it would be tidier with an initialization list
    // but the Intel C++ compiler doesn't implement such a nicety.
    typedef std::pair<std::string,update_scheme_t> name_scheme_t;
    available_schemes.insert(name_scheme_t("euler",EULER_UPDATE));
    available_schemes.insert(name_scheme_t("Euler",EULER_UPDATE));
    available_schemes.insert(name_scheme_t("pc",PC_UPDATE)); 
    available_schemes.insert(name_scheme_t("PC",PC_UPDATE));
    available_schemes.insert(name_scheme_t("predictor_corrector",PC_UPDATE));
    available_schemes.insert(name_scheme_t("predictor-corrector",PC_UPDATE));
    available_schemes.insert(name_scheme_t("Predictor_corrector",PC_UPDATE));
    available_schemes.insert(name_scheme_t("Predictor-corrector",PC_UPDATE));
    available_schemes.insert(name_scheme_t("midpoint",MIDPOINT_UPDATE)); 
    available_schemes.insert(name_scheme_t("mid-point",MIDPOINT_UPDATE));
    available_schemes.insert(name_scheme_t("mid_point",MIDPOINT_UPDATE)); 
    available_schemes.insert(name_scheme_t("Midpoint",MIDPOINT_UPDATE));
    available_schemes.insert(name_scheme_t("Mid-point",MIDPOINT_UPDATE));
    available_schemes.insert(name_scheme_t("Mid_point",MIDPOINT_UPDATE));
    available_schemes.insert(name_scheme_t("classic-rk3",CLASSIC_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("classic_rk3",CLASSIC_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("Classic-RK3",CLASSIC_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("Classic_RK3",CLASSIC_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("tvd-rk3",TVD_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("tvd_rk3",TVD_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("TVD-RK3",TVD_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("TVD_RK3",TVD_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("denman-rk3",DENMAN_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("denman_rk3",DENMAN_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("Denman-RK3",DENMAN_RK3_UPDATE));
    available_schemes.insert(name_scheme_t("Denman_RK3",DENMAN_RK3_UPDATE));
    return SUCCESS;
}

std::map<std::string,flux_calc_t> available_calculators;
int init_available_calculators_map()
{
    typedef std::pair<std::string,flux_calc_t> name_flux_t;
    available_calculators.insert(name_flux_t("0",FLUX_RIEMANN));
    available_calculators.insert(name_flux_t("riemann",FLUX_RIEMANN));
    available_calculators.insert(name_flux_t("Riemann",FLUX_RIEMANN));
    available_calculators.insert(name_flux_t("1",FLUX_AUSM));
    available_calculators.insert(name_flux_t("ausm",FLUX_AUSM));
    available_calculators.insert(name_flux_t("AUSM",FLUX_AUSM));
    available_calculators.insert(name_flux_t("2",FLUX_EFM));
    available_calculators.insert(name_flux_t("efm",FLUX_EFM));
    available_calculators.insert(name_flux_t("EFM",FLUX_EFM));
    available_calculators.insert(name_flux_t("3",FLUX_AUSMDV));
    available_calculators.insert(name_flux_t("ausmdv",FLUX_AUSMDV));
    available_calculators.insert(name_flux_t("AUSMDV",FLUX_AUSMDV));
    available_calculators.insert(name_flux_t("4",FLUX_ADAPTIVE));
    available_calculators.insert(name_flux_t("adaptive",FLUX_ADAPTIVE));
    available_calculators.insert(name_flux_t("Adaptive",FLUX_ADAPTIVE));
    available_calculators.insert(name_flux_t("5",FLUX_AUSM_PLUS_UP));
    available_calculators.insert(name_flux_t("ausm_plus_up",FLUX_AUSM_PLUS_UP));
    available_calculators.insert(name_flux_t("AUSM_plus_up",FLUX_AUSM_PLUS_UP));
    available_calculators.insert(name_flux_t("6",FLUX_HLLE));
    available_calculators.insert(name_flux_t("hlle",FLUX_HLLE));
    available_calculators.insert(name_flux_t("HLLE",FLUX_HLLE));
    return SUCCESS;
}

std::map<std::string,thermo_interp_t> available_interpolators;
int init_available_interpolators_map()
{
    typedef std::pair<std::string,thermo_interp_t> name_interp_t;
    available_interpolators.insert(name_interp_t("pt",INTERP_PT));
    available_interpolators.insert(name_interp_t("pT",INTERP_PT));
    available_interpolators.insert(name_interp_t("PT",INTERP_PT));
    available_interpolators.insert(name_interp_t("rhoe",INTERP_RHOE));
    available_interpolators.insert(name_interp_t("rhoE",INTERP_RHOE));
    available_interpolators.insert(name_interp_t("RHOE",INTERP_RHOE));
    available_interpolators.insert(name_interp_t("rhop",INTERP_RHOP));
    available_interpolators.insert(name_interp_t("rhoP",INTERP_RHOP));
    available_interpolators.insert(name_interp_t("RHOP",INTERP_RHOP));
    available_interpolators.insert(name_interp_t("rhot",INTERP_RHOT));
    available_interpolators.insert(name_interp_t("rhoT",INTERP_RHOT));
    available_interpolators.insert(name_interp_t("RHOT",INTERP_RHOT));
    return SUCCESS;
}

std::map<std::string,turbulence_model_t> available_turbulence_models;
int init_available_turbulence_models_map()
{
    typedef std::pair<std::string,turbulence_model_t> name_turb_model_t;
    available_turbulence_models.insert(name_turb_model_t("none",TM_NONE));
    available_turbulence_models.insert(name_turb_model_t("None",TM_NONE));
    available_turbulence_models.insert(name_turb_model_t("baldwin_lomax",TM_BALDWIN_LOMAX));
    available_turbulence_models.insert(name_turb_model_t("baldwin-lomax",TM_BALDWIN_LOMAX));
    available_turbulence_models.insert(name_turb_model_t("Baldwin_Lomax",TM_BALDWIN_LOMAX));
    available_turbulence_models.insert(name_turb_model_t("Baldwin-Lomax",TM_BALDWIN_LOMAX));
    available_turbulence_models.insert(name_turb_model_t("k_omega",TM_K_OMEGA));
    available_turbulence_models.insert(name_turb_model_t("k-omega",TM_K_OMEGA));
    available_turbulence_models.insert(name_turb_model_t("K_Omega",TM_K_OMEGA));
    available_turbulence_models.insert(name_turb_model_t("K-Omega",TM_K_OMEGA));
    available_turbulence_models.insert(name_turb_model_t("spalart_allmaras",TM_SPALART_ALLMARAS));
    return SUCCESS;
}

std::map<std::string,bc_t> available_bcs;
int init_available_bcs_map()
{
    typedef std::pair<std::string,bc_t> name_bc_t;
    // We keep the integer values for backward compatibility.
    available_bcs.insert(name_bc_t("adjacent",ADJACENT));
    available_bcs.insert(name_bc_t("ADJACENT",ADJACENT));
    available_bcs.insert(name_bc_t("0",ADJACENT));
    available_bcs.insert(name_bc_t("sup_in",SUP_IN));
    available_bcs.insert(name_bc_t("SUP_IN",SUP_IN));
    available_bcs.insert(name_bc_t("1",SUP_IN));
    available_bcs.insert(name_bc_t("extrapolate_out",EXTRAPOLATE_OUT));
    available_bcs.insert(name_bc_t("EXTRAPOLATE_OUT",EXTRAPOLATE_OUT));
    available_bcs.insert(name_bc_t("2",EXTRAPOLATE_OUT));
    available_bcs.insert(name_bc_t("slip_wall",SLIP_WALL));
    available_bcs.insert(name_bc_t("SLIP_WALL",SLIP_WALL));
    available_bcs.insert(name_bc_t("3",SLIP_WALL));
    available_bcs.insert(name_bc_t("adiabatic",ADIABATIC));
    available_bcs.insert(name_bc_t("ADIABATIC",ADIABATIC));
    available_bcs.insert(name_bc_t("4",ADIABATIC));
    available_bcs.insert(name_bc_t("fixed_T",FIXED_T));
    available_bcs.insert(name_bc_t("FIXED_T",FIXED_T));
    available_bcs.insert(name_bc_t("5",FIXED_T));
    available_bcs.insert(name_bc_t("subsonic_in",SUBSONIC_IN));
    available_bcs.insert(name_bc_t("SUBSONIC_IN",SUBSONIC_IN));
    available_bcs.insert(name_bc_t("6",SUBSONIC_IN));
    available_bcs.insert(name_bc_t("subsonic_out",SUBSONIC_OUT));
    available_bcs.insert(name_bc_t("SUBSONIC_OUT",SUBSONIC_OUT));
    available_bcs.insert(name_bc_t("7",SUBSONIC_OUT));
    available_bcs.insert(name_bc_t("transient_uni",TRANSIENT_UNI));
    available_bcs.insert(name_bc_t("TRANSIENT_UNI",TRANSIENT_UNI));
    available_bcs.insert(name_bc_t("8",TRANSIENT_UNI));
    available_bcs.insert(name_bc_t("transient_prof",TRANSIENT_PROF));
    available_bcs.insert(name_bc_t("TRANSIENT_PROF",TRANSIENT_PROF));
    available_bcs.insert(name_bc_t("9",TRANSIENT_PROF));
    available_bcs.insert(name_bc_t("static_prof",STATIC_PROF));
    available_bcs.insert(name_bc_t("STATIC_PROF",STATIC_PROF));
    available_bcs.insert(name_bc_t("10",STATIC_PROF));
    available_bcs.insert(name_bc_t("fixed_p_out",FIXED_P_OUT));
    available_bcs.insert(name_bc_t("FIXED_P_OUT",FIXED_P_OUT));
    available_bcs.insert(name_bc_t("11",FIXED_P_OUT));
    available_bcs.insert(name_bc_t("transient_T_wall",TRANSIENT_T_WALL));
    available_bcs.insert(name_bc_t("TRANSIENT_T_WALL",TRANSIENT_T_WALL));
    available_bcs.insert(name_bc_t("13",TRANSIENT_T_WALL));
    available_bcs.insert(name_bc_t("surface_energy_balance",SEB));
    available_bcs.insert(name_bc_t("SURFACE_ENERGY_BALANCE",SEB));
    available_bcs.insert(name_bc_t("SEB",SEB));
    available_bcs.insert(name_bc_t("15",SEB));
    available_bcs.insert(name_bc_t("user_defined",USER_DEFINED));
    available_bcs.insert(name_bc_t("USER_DEFINED",USER_DEFINED));
    available_bcs.insert(name_bc_t("16",USER_DEFINED));
    available_bcs.insert(name_bc_t("adjacent_plus_udf",ADJACENT_PLUS_UDF));
    available_bcs.insert(name_bc_t("ADJACENT_PLUS_UDF",ADJACENT_PLUS_UDF));
    available_bcs.insert(name_bc_t("17",ADJACENT_PLUS_UDF));
    available_bcs.insert(name_bc_t("ablating",ABLATING));
    available_bcs.insert(name_bc_t("ABLATING",ABLATING));
    available_bcs.insert(name_bc_t("18",ABLATING));
    available_bcs.insert(name_bc_t("sliding_T",SLIDING_T));
    available_bcs.insert(name_bc_t("SLIDING_T",SLIDING_T));
    available_bcs.insert(name_bc_t("19",SLIDING_T));
    available_bcs.insert(name_bc_t("fstc",FSTC));
    available_bcs.insert(name_bc_t("FSTC",FSTC));
    available_bcs.insert(name_bc_t("20",FSTC));
    available_bcs.insert(name_bc_t("shock_fitting_in",SHOCK_FITTING_IN));
    available_bcs.insert(name_bc_t("SHOCK_FITTING_IN",SHOCK_FITTING_IN));
    available_bcs.insert(name_bc_t("21",SHOCK_FITTING_IN));
    available_bcs.insert(name_bc_t("non_catalytic",NON_CATALYTIC));
    available_bcs.insert(name_bc_t("NON_CATALYTIC",NON_CATALYTIC));
    available_bcs.insert(name_bc_t("22",NON_CATALYTIC));
    available_bcs.insert(name_bc_t("equil_catalytic",EQUIL_CATALYTIC));
    available_bcs.insert(name_bc_t("EQUIL_CATALYTIC",EQUIL_CATALYTIC));
    available_bcs.insert(name_bc_t("23",EQUIL_CATALYTIC));
    available_bcs.insert(name_bc_t("super_catalytic",SUPER_CATALYTIC));
    available_bcs.insert(name_bc_t("SUPER_CATALYTIC",SUPER_CATALYTIC));
    available_bcs.insert(name_bc_t("24",SUPER_CATALYTIC));
    available_bcs.insert(name_bc_t("partially_catalytic",PARTIALLY_CATALYTIC));
    available_bcs.insert(name_bc_t("PARTIALLY_CATALYTIC",PARTIALLY_CATALYTIC));
    available_bcs.insert(name_bc_t("25",PARTIALLY_CATALYTIC));
    return SUCCESS;
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
int read_config_parameters(const string filename, bool master)
{
    global_data &G = *get_global_data_ptr();
    size_t jb;
    init_available_schemes_map();
    init_available_calculators_map();
    init_available_turbulence_models_map();
    init_available_bcs_map();

    // Default values for some configuration variables.
    G.dimensions = 2;
    G.Xorder = 2;
    set_radiation_flag(0);

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

    G.turbulence_model = TM_NONE;
    G.turbulence_prandtl = 0.89;
    G.turbulence_schmidt = 0.75;

    // Most configuration comes from the previously-generated INI file.
    ConfigParser dict = ConfigParser(filename);
    int i_value;
    string s_value, s_value2;
    double d_value;

    dict.parse_string("global_data", "title", G.title, "unknown");
    dict.parse_size_t("global_data", "dimensions", G.dimensions, 2);
    dict.parse_size_t("global_data", "control_count", G.control_count, 10);
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

    dict.parse_int("global_data", "BGK_flag", i_value, 0);
    set_BGK_flag( i_value );

    if (get_BGK_flag() > 0) {

	dict.parse_int("global_data", "velocity_buckets", i_value, 0);
	set_velocity_buckets( i_value );
    
	if (get_velocity_buckets() > 0) {
	    std::vector<Vector3> *vct = get_vcoords_ptr();
	    std::vector<double> tmp;
	    dict.parse_vector_of_doubles("global_data", "vcoords_x", tmp, tmp);
	    for (size_t tid = 0; tid < get_velocity_buckets(); ++tid) {
		(*vct)[tid].x = tmp[tid];
	    }
	    tmp.resize(0);
	    dict.parse_vector_of_doubles("global_data", "vcoords_y", tmp, tmp);
	    for (size_t tid = 0; tid < get_velocity_buckets(); ++tid) {
		(*vct)[tid].y = tmp[tid];
	    }
	    tmp.resize(0);
	    dict.parse_vector_of_doubles("global_data", "vcoords_z", tmp, tmp);
	    for (size_t tid = 0; tid < get_velocity_buckets(); ++tid) {
		(*vct)[tid].z = tmp[tid];
	    }
	    std::vector<double> *vwt = get_vweights_ptr();
	    dict.parse_vector_of_doubles("global_data", "vweights", tmp, tmp);
	    for (size_t tid = 0; tid < get_velocity_buckets(); ++tid) {
		(*vwt)[tid] = tmp[tid];
	    }
	}
	else {
	    cout << "Failure setting BGK velocities." << endl;
	    return FAILURE;
	}
	    
    }

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

    dict.parse_int("global_data", "axisymmetric_flag", i_value, 0);
    set_axisymmetric_flag( i_value );
    if ( get_verbose_flag() ) {
	cout << "axisymmetric_flag = " << get_axisymmetric_flag() << endl;
    }

    dict.parse_int("global_data", "viscous_flag", i_value, 0);
    set_viscous_flag( i_value );
    dict.parse_double("global_data", "viscous_delay", G.viscous_time_delay, 0.0);
    dict.parse_double("global_data", "viscous_factor_increment", d_value, 0.01);
    set_viscous_factor_increment( d_value );
    dict.parse_int("global_data", "viscous_upwinding_flag", i_value, 0);
    set_viscous_upwinding_flag( i_value );
    // FIX-ME 2013-04-23 should probably merge diffusion_model and diffusion_flag
    // as we have done for turbulence_model, below.
    dict.parse_int("global_data", "diffusion_flag", i_value, 0);
    set_diffusion_flag( i_value );
    if ( get_verbose_flag() ) {
	cout << "viscous_flag = " << get_viscous_flag() << endl;
	cout << "viscous_delay = " << G.viscous_time_delay << endl;
	cout << "viscous_factor_increment = " << get_viscous_factor_increment() << endl;
	cout << "viscous_upwinding_flag = " << get_viscous_upwinding_flag() << endl;
	cout << "diffusion_flag = " << get_diffusion_flag() << endl;
    }
    if( get_diffusion_flag() && ( gmodel->get_number_of_species() > 1 ) ) { 
 	dict.parse_string("global_data", "diffusion_model", s_value, "Stefan-Maxwell");
	set_diffusion_model(s_value);
	if( get_verbose_flag() ) {
	    cout << "diffusion_model = " << s_value << endl;
	}
    }

    dict.parse_int("global_data", "shock_fitting_flag", i_value, 0);
    set_shock_fitting_flag( i_value );
    dict.parse_int("global_data", "shock_fitting_decay_flag", i_value, 0);
    set_shock_fitting_decay_flag( i_value );
    dict.parse_double("global_data", "shock_fitting_speed_factor", G.shock_fitting_speed_factor, 1.0);
    dict.parse_int("global_data", "moving_grid_flag", i_value, 0);
    set_moving_grid_flag( i_value );
    dict.parse_int("global_data", "write_vertex_velocities_flag", i_value, 0);
    set_write_vertex_velocities_flag( i_value );
    dict.parse_int("global_data", "adaptive_reconstruction_flag", i_value, 0);
    set_adaptive_reconstruction_flag( i_value );
    if ( get_verbose_flag() ) {
	cout << "shock_fitting_flag = " << get_shock_fitting_flag() << endl;
	cout << "shock_fitting_decay_flag = " << get_shock_fitting_decay_flag() << endl;
	cout << "shock_fitting_speed_factor = " << G.shock_fitting_speed_factor << endl;
	cout << "moving_grid_flag = " << get_moving_grid_flag() << endl;
	cout << "write_vertex_velocities_flag = " << get_write_vertex_velocities_flag() << endl;
	cout << "adaptive_reconstruction_flag = " << get_adaptive_reconstruction_flag() << endl;
    }

    // 2013-apr-23 New specification scheme for turbulence models.
    dict.parse_string("global_data", "turbulence_model", s_value, "none");
    G.turbulence_model = available_turbulence_models[s_value];
    dict.parse_double("global_data", "turbulence_prandtl_number", d_value, 0.89);
    G.turbulence_prandtl = d_value;
    dict.parse_double("global_data", "turbulence_schmidt_number", d_value, 0.75);
    G.turbulence_schmidt = d_value;
    dict.parse_double("global_data", "max_mu_t_factor", G.max_mu_t_factor, 300.0);
    dict.parse_double("global_data", "transient_mu_t_factor", G.transient_mu_t_factor, 1.0);
    if ( G.turbulence_model == TM_SPALART_ALLMARAS )
	throw std::runtime_error("Spalart-Allmaras turbulence model not available.");
    if ( get_verbose_flag() ) {
	cout << "turbulence_model = " << get_name_of_turbulence_model(G.turbulence_model) << endl;
	cout << "turbulence_prandtl_number = " << G.turbulence_prandtl << endl;
	cout << "turbulence_schmidt_number = " << G.turbulence_schmidt << endl;
	cout << "max_mu_t_factor = " << G.max_mu_t_factor << endl;
	cout << "transient_mu_t_factor = " << G.transient_mu_t_factor << endl;
    }

    dict.parse_size_t("global_data", "max_invalid_cells", G.max_invalid_cells, 10);
    dict.parse_string("global_data", "flux_calc", s_value, "adaptive");
    set_flux_calculator(available_calculators[s_value]);
    dict.parse_double("global_data", "compression_tolerance", d_value, -0.30);
    set_compression_tolerance(d_value);
    dict.parse_double("global_data", "shear_tolerance", d_value, 0.20);
    set_shear_tolerance(d_value);
    dict.parse_double("global_data", "M_inf", d_value, 0.01);
    set_M_inf(d_value);
    dict.parse_string("global_data", "interpolation_type", s_value, "rhoe");
    set_thermo_interpolator(available_interpolators[s_value]);
    dict.parse_int("global_data", "apply_limiter_flag", i_value, 1);
    set_apply_limiter_flag(i_value == 1);
    dict.parse_int("global_data", "extrema_clipping_flag", i_value, 1);
    set_extrema_clipping_flag(i_value == 1);
    if ( get_verbose_flag() ) {
	cout << "max_invalid_cells = " << G.max_invalid_cells << endl;
	cout << "flux_calc = " << get_flux_calculator_name(get_flux_calculator()) << endl;
	cout << "compression_tolerance = " << get_compression_tolerance() << endl;
	cout << "shear_tolerance = " << get_shear_tolerance() << endl;
	cout << "M_inf = " << get_M_inf() << endl;
	cout << "interpolation_type = " << get_thermo_interpolator_name(get_thermo_interpolator()) << endl;
	cout << "apply_limiter_flag = " << get_apply_limiter_flag() << endl;
	cout << "extrema_clipping_flag = " << get_extrema_clipping_flag() << endl;
    }

    dict.parse_int("global_data", "filter_flag", i_value, 0);
    set_filter_flag( i_value );
    dict.parse_double("global_data", "filter_tstart", G.filter_tstart, 0.0);
    //set_filter_tstart(d_value);
    dict.parse_double("global_data", "filter_tend", G.filter_tend, 0.0);
    //set_filter_tend(d_value);
    dict.parse_double("global_data", "filter_dt", G.filter_dt, 0.0);
    //set_filter_dt(d_value);
    dict.parse_double("global_data", "filter_mu", G.filter_mu, 0.0);
    //set_filter_mu(d_value);
    dict.parse_size_t("global_data", "filter_npass", G.filter_npass, 0);

    dict.parse_int("global_data", "sequence_blocks", i_value, 0);
    G.sequence_blocks = (i_value == 1);
    if ( get_verbose_flag() ) {
	cout << "sequence_blocks = " << G.sequence_blocks << endl;
    }

    // Read a number of gas-states.
    dict.parse_size_t("global_data", "nflow", G.n_gas_state, 0);
    if ( get_verbose_flag() ) {
	cout << "nflow = " << G.n_gas_state << endl;
    }
    for ( size_t ig = 0; ig < G.n_gas_state; ++ig ) {
	G.gas_state.push_back(read_flow_condition_from_ini_dict(dict, ig, master));
	if ( get_verbose_flag() ) {
	    cout << "flow condition[" << ig << "]: " << *(G.gas_state[ig]) << endl;
	}
    }

    // Read the parameters for a number of blocks.
    dict.parse_size_t("global_data", "nblock", G.nblock, 0);
    if ( get_verbose_flag() ) {
	printf( "nblock = %d\n", static_cast<int>(G.nblock));
    }
    // We keep a record of all of the configuration data for all blocks
    // but, eventually, we may allocate the flow-field data for only a 
    // subset of these blocks. 
    G.bd.resize(G.nblock);

    // Number of pistons
    // FIX-ME code needs to be reworked...
    dict.parse_size_t("global_data", "npiston", G.npiston, 0);
#   ifdef _MPI
    if ( G.npiston > 0 ) {
	G.npiston = 0;
	cout << "Pistons cannot be used with MPI code, yet." << endl;
    }
#   endif
    G.pistons.resize(G.npiston);
    if ( get_verbose_flag() ) {
	printf( "npiston = %d\n", static_cast<int>(G.npiston));
    }
    for ( size_t jp = 0; jp < G.npiston; ++jp ) {
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

    dict.parse_size_t("global_data", "nheatzone", G.n_heat_zone, 0);
    dict.parse_double("global_data", "heat_time_start", G.heat_time_start, 0.0);
    dict.parse_double("global_data", "heat_time_stop", G.heat_time_stop, 0.0);
    dict.parse_double("global_data", "heat_factor_increment", d_value, 0.01);
    set_heat_factor_increment( d_value );
    if ( get_verbose_flag() ) {
	printf("nheatzone = %d\n", static_cast<int>(G.n_heat_zone));
	printf("heat_time_start = %e\n", G.heat_time_start);
	printf("heat_time_stop = %e\n", G.heat_time_stop);
	printf("heat_factor_increment = %e\n", get_heat_factor_increment() );
    }
    G.heat_zone.resize(G.n_heat_zone);
    for ( size_t indx = 0; indx < G.n_heat_zone; ++indx ) {
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

    dict.parse_size_t("global_data", "nreactionzone", G.n_reaction_zone, 0);
    if ( get_verbose_flag() ) {
	printf("nreactionzone = %d\n", static_cast<int>(G.n_reaction_zone));
    }
    G.reaction_zone.resize(G.n_reaction_zone);
    for ( size_t indx = 0; indx < G.n_reaction_zone; ++indx ) {
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

    dict.parse_size_t("global_data", "nturbulencezone", G.n_turbulent_zone, 0);
    if ( get_verbose_flag() ) {
	printf("nturbulencezone = %d\n", static_cast<int>(G.n_turbulent_zone));
    }
    G.turbulent_zone.resize(G.n_turbulent_zone);
    for ( size_t indx = 0; indx < G.n_turbulent_zone; ++indx ) {
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


int read_control_parameters( const string filename, bool master, bool first_time )
// These are read at the start of every time-step and may be used 
// to alter the simulation behaviour during a run.
{
    int i_value;
    std::string s_value;
    global_data &G = *get_global_data_ptr();
    // Parse the previously-generated INI file.
    ConfigParser dict = ConfigParser(filename);

    dict.parse_int("control_data", "x_order", G.Xorder, 2); // default high-order
    // 2013-03-31 change to use an explicitly-named update scheme.
    dict.parse_string("control_data", "gasdynamic_update_scheme", s_value,
		      "predictor-corrector");
    set_gasdynamic_update_scheme(available_schemes[s_value]);
    // To keep backward compatibility with old simulation files,
    // read Torder if it exists and set the equivalent update scheme.
    dict.parse_int("control_data", "t_order", i_value, 0);
    switch ( i_value ) {
    case 1: set_gasdynamic_update_scheme(available_schemes["euler"]); break;
    case 2: set_gasdynamic_update_scheme(available_schemes["predictor-corrector"]); break;
    case 3: set_gasdynamic_update_scheme(available_schemes["denman-rk3"]); break;
    default: /* do nothing */;
    }
    dict.parse_int("control_data", "separate_update_for_viscous_flag", i_value, 0);
    set_separate_update_for_viscous_flag( i_value );
    dict.parse_double("control_data", "dt", G.dt_init, 1.0e-6);
    if ( first_time ) G.dt_global = G.dt_init;
    dict.parse_boolean("control_data", "fixed_time_step", G.fixed_time_step, 0);
    dict.parse_double("control_data", "dt_reduction_factor",
		      G.dt_reduction_factor, 0.2);
    dict.parse_double("control_data", "cfl", G.cfl_target, 0.5);
    dict.parse_int("control_data", "stringent_cfl", i_value, 0);
    set_stringent_cfl_flag( i_value );
    dict.parse_size_t("control_data", "print_count", G.print_count, 20);
    dict.parse_size_t("control_data", "cfl_count", G.cfl_count, 10);
    dict.parse_double("control_data", "dt_shock", G.dt_shock, 1.0e-3);
    dict.parse_double("control_data", "dt_plot", G.dt_plot, 1.0e-3);
    dict.parse_size_t("control_data", "write_at_step", G.write_at_step, 0);
    dict.parse_double("control_data", "dt_history", G.dt_his, 1.0e-3);
    dict.parse_double("control_data", "dt_fstc", G.dt_fstc, 1.0e-3);
    dict.parse_double("control_data", "max_time", G.max_time, 1.0e-3);
    dict.parse_size_t("control_data", "max_step", G.max_step, 10);
    dict.parse_int("control_data", "halt_now", G.halt_now, 0);
    dict.parse_int("control_data", "implicit_flag", i_value, 0);
    set_implicit_flag( i_value );
    dict.parse_int("control_data", "radiation_update_frequency", i_value, 1);
    set_radiation_update_frequency(i_value);
    if ( first_time && get_verbose_flag() ) {
	cout << "Time-step control parameters:" << endl;
	cout << "    x_order = " << G.Xorder << endl;
	cout << "    gasdynamic_update_scheme = " 
	     << get_name_of_gasdynamic_update_scheme(get_gasdynamic_update_scheme())
	     << endl;
	cout << "separate_update_for_viscous_flag = " 
	     << get_separate_update_for_viscous_flag() << endl;
	cout << "    dt = " << G.dt_init << endl;
	cout << "    fixed_time_step = " << G.fixed_time_step << endl;
	cout << "    dt_reduction_factor = " 
	     << G.dt_reduction_factor << endl;
	cout << "    cfl = " << G.cfl_target << endl;
	cout << "    stringent_cfl = " << get_stringent_cfl_flag() << endl;
	cout << "    print_count = " << G.print_count << endl;
	cout << "    cfl_count = " << G.cfl_count << endl;
	cout << "    dt_plot = " << G.dt_plot << endl;
	cout << "    write_at_step = " << G.write_at_step << endl;
	cout << "    dt_shock = " << G.dt_shock << endl;
	cout << "    dt_history = " << G.dt_his << endl;
	cout << "    dt_fstc = " << G.dt_fstc << endl;
	cout << "    max_time = " << G.max_time << endl;
	cout << "    max_step = " << G.max_step << endl;
	cout << "    halt_now = " << G.halt_now << endl;
	cout << "    halt_now = " << G.halt_now << endl;
	cout << "    radiation_update_frequency = " 
	     << get_radiation_update_frequency();
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
	    cout << "\nInvalid implicit flag value: " 
		 << get_implicit_flag() << endl;
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
int assign_blocks_to_mpi_rank(const string filename, bool master)
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
	    size_t nrank = 0;
	    size_t nblock;
	    size_t nblock_total = 0;
	    std::vector<int> block_ids, dummy_block_ids;
	    dict.parse_size_t("global", "nrank", nrank, 0);
	    if ( G.num_mpi_proc != static_cast<int>(nrank) ) {
		if ( master ) {
		    printf("    Error in specifying mpirun -np\n");
		    printf("    It needs to match number of nrank; present values are:\n");
		    printf("    num_mpi_proc= %d nrank= %d\n", G.num_mpi_proc,
			   static_cast<int>(nrank));
		}
		return FAILURE;
	    }
	    for ( size_t rank=0; rank < nrank; ++rank ) {
		string section = "rank/" + tostring(rank);
		dict.parse_size_t(section, "nblock", nblock, 0);
		block_ids.resize(0);
		dummy_block_ids.resize(nblock);
		for ( size_t i = 0; i < nblock; ++i ) dummy_block_ids[i] = -1;
		dict.parse_vector_of_ints(section, "blocks", block_ids, dummy_block_ids);
		if ( nblock != block_ids.size() ) {
		    if ( master ) {
			printf("    Did not pick up correct number of block_ids:\n");
			printf("        rank=%d, nblock=%d, block_ids.size()=%d\n",
			       static_cast<int>(rank), static_cast<int>(nblock),
			       static_cast<int>(block_ids.size()));
		    }
		    return FAILURE;
		}
		for ( size_t i = 0; i < nblock; ++i ) {
		    int this_block_id = block_ids[i];
		    if ( this_block_id < 0 ) {
			if ( master ) printf("    Error, invalid block id: %d\n", this_block_id);
			return FAILURE;
		    }
		    if ( G.my_mpi_rank == static_cast<int>(rank) )
			G.my_blocks.push_back(&(G.bd[this_block_id]));
		    G.mpi_rank_for_block[this_block_id] = static_cast<int>(rank);
		    nblock_total += 1;
		} // end for i
	    } // end for rank
	    if ( get_verbose_flag() ) {
		printf("    my_rank=%d, block_ids=", static_cast<int>(G.my_mpi_rank));
		for ( size_t i=0; i < G.my_blocks.size(); ++i ) {
		    printf(" %d", static_cast<int>(G.my_blocks[i]->id));
		}
		printf("\n");
	    }
	    if ( master ) {
		if ( nblock_total != G.nblock ) {
		    printf("    Error, total number of blocks incorrect: total=%d G.nblock=%d\n",
			   static_cast<int>(nblock_total), static_cast<int>(G.nblock));
		    return FAILURE;
		}
	    }
	} else {
	    if ( get_verbose_flag() && master ) {
		printf("    MPI parallel, No MPI map file specified.\n");
		printf("    Identify each block with the corresponding MPI rank.\n");
	    }
	    if ( G.num_mpi_proc != static_cast<int>(G.nblock) ) {
		if ( master ) {
		    printf("    Error in specifying mpirun -np\n");
		    printf("    It needs to match number of blocks; present values are:\n");
		    printf("    num_mpi_proc= %d nblock= %d\n", 
			   static_cast<int>(G.num_mpi_proc), static_cast<int>(G.nblock));
		}
		return FAILURE;
	    }
	    G.my_blocks.push_back(&(G.bd[G.my_mpi_rank]));
	    for ( size_t jb=0; jb < G.nblock; ++jb ) {
		G.mpi_rank_for_block.push_back(jb);
	    }
	}
    } else {
	if ( get_verbose_flag() ) {
	    printf("    Since we are not doing MPI, all blocks in same process.\n");
	}
	for ( size_t jb=0; jb < G.nblock; ++jb ) {
	    G.my_blocks.push_back(&(G.bd[jb]));
	    G.mpi_rank_for_block.push_back(G.my_mpi_rank);
	} 
    } // endif
    return SUCCESS;
} // end assign_blocks_to_mpi_rank()


/** \brief Use config_parser functions to read the flow condition data. 
 */
CFlowCondition *read_flow_condition_from_ini_dict(ConfigParser &dict, size_t indx, bool master)
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
    size_t nsp = gmodel->get_number_of_species();
    vnf.resize(nsp);
    for ( size_t isp = 0; isp < nsp; ++isp ) vnf[isp] = 0.0;
    dict.parse_vector_of_doubles(section, "massf", massf, vnf);
    size_t nmodes = gmodel->get_number_of_modes();
    vnf.resize(nmodes);
    for ( size_t imode = 0; imode < nmodes; ++imode ) vnf[imode] = 300.0; 
    dict.parse_vector_of_doubles(section, "T", T, vnf);
    int S = 0;  // shock indicator
    cfcp = new CFlowCondition( gmodel, p, u, v, w, T, massf, flow_label, 
			       tke, omega, mu_t, k_t, S, Bx, By, Bz);
    vnf.clear();
    massf.clear();
    T.clear();
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
int set_block_parameters(size_t id, ConfigParser &dict, bool master)
{
    global_data &G = *get_global_data_ptr();
    Block &bd = *get_block_data_ptr(id);
    int indx, iface, other_block, other_face, neighbour_orientation;
    bc_t wc_bc, bc_type_code;
    int x_order, sponge_flag, xforce_flag;
    size_t n_profile;
    string value_string, block_label, filename, wcbc_fname;
    int inflow_condition_id, is_wall, use_udf_flux, assume_ideal;
    double Twall, Pout, epsilon;
    std::vector<double> f_wall_va;
    std::vector<double> f_wall, vnf;
    std::vector<double> mdot;
    double Twall_i, Twall_f, t_i, t_f;
    vnf.resize(get_gas_model_ptr()->get_number_of_species(), 0.0);

    bd.id = id;
  
    // Read parameters from INI-file dictionary.
    indx = id;
    string section = "block/" + tostring(indx);

    dict.parse_string(section, "label", block_label, "");
    // Assume all blocks are active. 
    // The active flag will be used to skip over inactive
    // or unused blocks in later sections of code.
    dict.parse_int(section, "active", bd.active, 1);
    if ( get_verbose_flag() ) {
	cout << section << ":label = " << block_label << endl;
	cout << "    active = " << bd.active << endl;
    }

    // Number of active cells in block.
    dict.parse_size_t(section, "nni", bd.nni, 2);
    dict.parse_size_t(section, "nnj", bd.nnj, 2);
    // Allow for ghost cells
    bd.nidim = bd.nni + 2 * NGHOST;
    bd.njdim = bd.nnj + 2 * NGHOST;
    // Set up min and max indices for convenience in later work.
    // Active cells should then be addressible as
    // get_cell(i,j), imin <= i <= imax, jmin <= j <= jmax.
    bd.imin = NGHOST;
    bd.imax = bd.imin + bd.nni - 1;
    bd.jmin = NGHOST;
    bd.jmax = bd.jmin + bd.nnj - 1;
    if ( G.dimensions == 3 ) {
	dict.parse_size_t(section, "nnk", bd.nnk, 2);
	bd.nkdim = bd.nnk + 2 * NGHOST;
	bd.kmin = NGHOST;
	bd.kmax = bd.kmin + bd.nnk - 1;
    } else {
	// For purely 2D flow geometry, we keep only one layer of cells.
	bd.nnk = 1;
	bd.nkdim = 1;
	bd.kmin = 0;
	bd.kmax = 0;
    }
    if ( get_verbose_flag() ) {
	printf( "    nni = %d, nnj = %d, nnk = %d\n", 
		static_cast<int>(bd.nni), static_cast<int>(bd.nnj),
		static_cast<int>(bd.nnk) );
	printf( "    nidim = %d, njdim = %d, nkdim = %d\n",
		static_cast<int>(bd.nidim), static_cast<int>(bd.njdim),
		static_cast<int>(bd.nkdim) );
    }

    // Rotating frame of reference.
    dict.parse_double(section, "omegaz", bd.omegaz, 0.0);

    // Boundary condition flags, 
    for ( iface = NORTH; iface <= ((G.dimensions == 3)? BOTTOM : WEST); ++iface ) {
	section = "block/" + tostring(indx) + "/face/" + get_face_name(iface);
	dict.parse_string(section, "bc", value_string, "slip_wall");
	// cout << "setting bc_type value_string=" << value_string << endl;
	bc_type_code = available_bcs[value_string];
	dict.parse_int(section, "inflow_condition", inflow_condition_id, 0);
	dict.parse_string(section, "filename", filename, "");
	dict.parse_size_t(section, "n_profile", n_profile, 1);
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
	// cout << "face= " << get_face_name(iface) << " bc=" << get_bc_name(bc_type_code) << endl;

	// Wall catalytic boundary condition flags, 
	// f_wall for SUPER_CATALYTIC 
	// wcbc_input_file for CATALYTIC and PARTIALLY_CATALYTIC
	//
	// Rowan or Dan **** FIX-ME *****
	dict.parse_string(section, "wcbc", value_string, "non_catalytic");
	wc_bc = available_bcs[value_string];
	// cout << "setting wc_bc, value_string=" << value_string << endl;
	dict.parse_string(section, "wcbc_input_file", wcbc_fname, "");
	dict.parse_vector_of_doubles(section, "f_wall", f_wall, vnf);
	// cout << "wc_bc= " << get_bc_name(wc_bc) << endl;

	bd.bcp[iface] = create_BC( &bd, iface, bc_type_code, inflow_condition_id, 
				    filename, n_profile, Twall, Pout,
				    x_order, is_wall, use_udf_flux,
				    other_block, other_face, neighbour_orientation,
				    sponge_flag, xforce_flag, mdot, epsilon,
				    wc_bc, wcbc_fname, f_wall,
				    Twall_i, Twall_f, t_i, t_f, assume_ideal);
	if ( get_verbose_flag() ) {
	    cout << "    " << get_face_name(iface) << " face:" << endl;
	    bd.bcp[iface]->print_info("        ");
	    cout << "        Pout= " << Pout << endl;
	    cout << "        Twall= " << Twall << endl;
	}
    } // end for iface

    // History Cells.
    section = "block/" + tostring(indx);
    dict.parse_size_t(section, "nhcell", bd.hncell, 0);
    for ( size_t ih = 0; ih < bd.hncell; ++ih ) {
	section = "block/" + tostring(indx);
	string key = "history-cell-" + tostring(ih);
	dict.parse_string(section, key, value_string, "0 0 0");
	if ( G.dimensions == 3 ) {
	    unsigned int hicell, hjcell, hkcell;
	    sscanf( value_string.c_str(), "%u %u %u", &hicell, &hjcell, &hkcell );
	    bd.hicell.push_back(hicell);
	    bd.hjcell.push_back(hjcell);
	    bd.hkcell.push_back(hkcell);
	} else {
	    unsigned int hicell, hjcell;
	    sscanf( value_string.c_str(), "%u %u", &hicell, &hjcell );
	    bd.hicell.push_back(hicell);
	    bd.hjcell.push_back(hjcell);
	    bd.hkcell.push_back(0);
	}
	if ( get_verbose_flag() ) {
	    printf( "    History cell[%d] located at indices [%d][%d][%d]\n",
		    static_cast<int>(ih), static_cast<int>(bd.hicell[ih]),
		    static_cast<int>(bd.hjcell[ih]), static_cast<int>(bd.hkcell[ih]) );
	}
    }

    return SUCCESS;
} // end set_block_parameters()

//------------------------------------------------------------------

// l_kernel.cxx

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/kinetics/reaction-update.hh"
#include "l_kernel.hh"
#include "l1d.hh"

SimulationData::SimulationData(std::string config_file_name, int echo_input)
{
    ConfigParser dict = ConfigParser(config_file_name);
    string reaction_scheme_file, gas_model_file;
    if (echo_input == 1) cout << endl << "Reading global_data..." << endl;

    dict.parse_int("global_data", "case_id", test_case, 0);
    L_set_case_id( test_case );
    dict.parse_string("global_data", "gas_model_file", gas_model_file, "gas-model.lua");
    Gas_model *gmodel = set_gas_model_ptr(create_gas_model(gas_model_file));
    dict.parse_string("global_data", "reaction_scheme_file", reaction_scheme_file, "None");
    dict.parse_int("global_data", "reacting_flag", fr_chem, 0);
    if( fr_chem ) set_reaction_update( reaction_scheme_file );
    if (echo_input == 1) {
	cout << "    test_case_id = " << test_case << endl;
	cout << "    gas_model_file = " << gas_model_file << endl;
	cout << "    nsp = " << gmodel->get_number_of_species() << endl;
	cout << "    nmodes = " << gmodel->get_number_of_modes() << endl;
	cout << "    reacting_flag = " << fr_chem << endl;
	cout << "    reaction_scheme_file = " << reaction_scheme_file << endl;
    }
    dict.parse_int("global_data", "nslug", nslug, 0);
    dict.parse_int("global_data", "npiston", npiston, 0);
    dict.parse_int("global_data", "ndiaphragm", ndiaphragm, 0);
    if (echo_input == 1) {
	cout << "    nslug = " << nslug << endl;
	cout << "    npiston = " << npiston << endl;
	cout << "    ndiaphragm = " << ndiaphragm << endl;
    }
    dict.parse_double("global_data", "max_time", max_time, 1.0e-3);
    dict.parse_int("global_data", "max_step", max_step, 100);
    dict.parse_double("global_data", "dt_init", dt_init, 1.0e-9);
    dict.parse_double("global_data", "cfl", CFL, 0.5);
    dict.parse_int("global_data", "x_order", Xorder, 2);
    dict.parse_int("global_data", "t_order", Torder, 2);
    dict.parse_double("global_data", "thermal_damping", k, 0.0);
    if (echo_input == 1) {
	cout << "    max_time = " << max_time << endl;
	cout << "    max_step = " << max_step << endl;
	cout << "    dt_init = " << dt_init << endl;
	cout << "    cfl = " << CFL << endl;
	cout << "    x_order = " << Xorder << endl;
	cout << "    t_order = " << Torder << endl;
	cout << "    thermal_damping = " << k << endl;
    }
    dict.parse_int("global_data", "n_dt_plot", n_dt_plot, 0);
    std::vector<double> vdbl, vdbl_default;
    vdbl_default.resize(n_dt_plot);
    for ( size_t i = 0; i < vdbl_default.size(); ++i ) vdbl_default[i] = 0.0;
    dict.parse_vector_of_doubles("global_data", "t_change", t_change, vdbl_default);
    for ( size_t i = 0; i < vdbl_default.size(); ++i ) vdbl_default[i] = 1.0e-3;
    dict.parse_vector_of_doubles("global_data", "dt_plot", dt_plot, vdbl_default);
    dict.parse_vector_of_doubles("global_data", "dt_his", dt_his, vdbl_default);
    if (echo_input == 1) {
	cout << "    n_dt_plot = " << n_dt_plot << endl;
	cout << "        t_change  dt_plot  dt_his" << endl;
	for ( int i = 0; i < n_dt_plot; ++i ) {
	    cout << "        " << t_change[i] 
		 << " " << dt_plot[i] << " " << dt_his[i] << endl;
	}
    }
    dict.parse_int("global_data", "hloc_n", hnloc, 0);
    vdbl_default.resize(hnloc);
    for ( size_t i = 0; i < vdbl_default.size(); ++i ) vdbl_default[i] = 0.0;
    dict.parse_vector_of_doubles("global_data", "hloc_x", hxloc, vdbl_default);
    if (echo_input == 1) {
	cout << "    hloc_n = " << hnloc << endl;
	cout << "    hloc_x =";
	for ( int i = 0; i < hnloc; ++i ) cout << " " << hxloc[i]; 
	cout << endl;
    }
} // end SimulationData constructor


SimulationData::~SimulationData()
{}


// The managed gas model lives here.
Gas_model *gmodel;

Gas_model *set_gas_model_ptr(Gas_model *gmptr)
{
    return gmodel = gmptr;
}

Gas_model *get_gas_model_ptr()
{
    return gmodel;
}

// The managed reaction update model lives here.
Reaction_update *rupdate;

int set_reaction_update(std::string file_name)
{
    rupdate = create_Reaction_update(file_name, *(get_gas_model_ptr()));
    if ( rupdate != 0 )
	return SUCCESS;
    else
	return FAILURE;
}

Reaction_update *get_reaction_update_ptr()
{
    return rupdate;
}

// The managed energy exchange update model lives here.
Energy_exchange_update *eeupdate;

int set_energy_exchange_update(std::string file_name)
{
    eeupdate = create_Energy_exchange_update(file_name, *(get_gas_model_ptr()));
    if ( eeupdate != 0 )
	return SUCCESS;
    else
	return FAILURE;
}

Energy_exchange_update *get_energy_exchange_update_ptr()
{
    return eeupdate;
}


int L_case_id = 0;

void L_set_case_id( int id ) {
    L_case_id = id;
}

int L_get_case_id( void ) {
    return L_case_id;
}

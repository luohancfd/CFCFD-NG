#include "e3conn.hh"

#include "../../eilmer3/source/kernel.hh"
#include "../../../lib/util/source/useful.h"
#include "initialise.hh"
#include "output.hh"
#include "solid_bc.hh"

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

Wall_model* initialise_wall_model(std::string fname, cht_coupling_t cht_coupling, int start_tindx)
{
    //Assign memory
    Wall_model* my_wall = new Wall_model();
    my_wall->cht_coupling = cht_coupling;
    my_wall->tindx = start_tindx;
    
    //Initialise block 
    my_wall->my_block = initialise_block(fname);
    if (start_tindx > 0){
	initialise_tindx_restart_from_eilmer(*my_wall, start_tindx);
    }
    return my_wall;
}

void get_near_wall_solid_cell_pos(Wall_model &wm, vector<Vector3> &solid_cell_pos)
{
    // Number of nodes in each direction as defined in grid file (1 to max)
    int nni; 
    nni = wm.my_block->nni; 
    
    //Only south wall connected => use j = o for i = 0 to nni-2 (index)
    //Fix me: Need ability to connect any boundary
    int j = 0;
    for (int i = 0; i<(nni-1); i++){
	solid_cell_pos.push_back(Vector3(wm.my_block->block_cells[i][j].pos[0], wm.my_block->block_cells[i][j].pos[1]));
    };
}

int write_solution(Wall_model &wm, double time_elapsed, int print_number)
{
    //Output the solutions
    write_temperatures_to_file(*(wm.my_block), time_elapsed, print_number); 
    return SUCCESS;
}

int gather_near_wall_solid_data(Wall_model &wm, vector<double> &T_near_wall, vector<double> &k_near_wall)
{
    int nni;
    nni = wm.my_block->nni;
    
    //Check size of array    
    if ( T_near_wall.size() != static_cast<size_t>(nni-1) ) {
	std::cout << "Warning T_near_wall vector is not correct size in e3conn";
	exit(1);
    }; 
    
    //Get south boundary near-wall cell data
    int j=0; 
    for (size_t i = 0; i < T_near_wall.size(); ++i ) {
	T_near_wall[i] = wm.my_block->block_cells[i][j].T;
	k_near_wall[i] = wm.my_block->block_cells[i][j].k;
    };
    
    for (size_t i = 0; i < T_near_wall.size(); i++) {
	   
      //      std::cout << "i = " << i << ", T_1st_row_wall = " << T_near_wall[i] << std::endl;
    }
    return SUCCESS;
}


int update_state(Wall_model &wm, double dt, const vector<double> q_wall, const vector<double> T_wall)
{
    //Assuming this performs one time step iteration does not keep track of time.

    BC_E3CONNECTION step_q_connection(q_wall); //create e3 connection bc. Only need for this scope/step    
    BC_E3CONNECTION step_T_connection(T_wall); //create e3 connection bc. Only need for this scope/step
    // Apply bc_e3connection to south wall

    if (wm.cht_coupling == TFS_QWS || wm.cht_coupling == QFS_QWS) {
	wm.my_block->e3connection_type_flag[0] = 0;
	wm.my_block->block_bc[0] = &step_q_connection; //assign address of south wall
	// std::cout << "Connection_type_flag = " << wm.my_block->e3connection_type_flag[0] << std::endl;
    }
    else if (wm.cht_coupling == TFS_TWS || wm.cht_coupling == QFS_TWS) {
	// std::cout << "-------------Testing the BC_E3CONNECTION_TEMPERATURE--------------" << std::endl;
	
	//Temperature flag checked in the boundary condition. But expecting Temperatures passed.
	wm.my_block->e3connection_type_flag[0] = 1;
	wm.my_block->block_bc[0] = &step_T_connection; //assign address of south wall
    }
    else {
	std::cout << "Something went wrong in e3conn.cxx adding BC_E3CONNECTION" << std::endl;
   	std::cout << "Connection_type_flag when wrong = " << wm.my_block->e3connection_type_flag[0] << std::endl;
	exit(1);
    }
    //Print the incoming flux and T_interface to verify it gets passed correctly
    // for (size_t i = 0; i < q_wall.size(); i++) {
    // 	std::cout << std::setprecision(12) << "i = " << i << ", q_flux = " << q_wall[i] << ", T_interface = " << std::setprecision(12) << T_wall[i] << std::endl;
    // }
    // std::cout << "Cell[0][0].T before " << setprecision(12) << wm.my_block->block_cells[0][0].T << std::endl;
    //Calculate fluxes and derivatives fluxes
    wm.my_block->update_boundary_conditions();
    wm.my_block->update_boundary_secondary_interface_temperatures();
    wm.my_block->update_internal_secondary_interface_temperatures();
    wm.my_block->update_vertex_derivatives();
    wm.my_block->update_interface_fluxes();
    wm.my_block->update_cell_energy_derivative();
    
    //Increment solution in time
    wm.my_block->time_update(dt, wm.update_scheme);
    // std::cout << "Cell[0][0].T After " << setprecision(12) << wm.my_block->block_cells[0][0].T << std::endl;
    // std::cout << "dt from eilmer " << setprecision(12)<< dt << std::endl;
    return SUCCESS;
}

void set_wallcon_time_update_scheme(Wall_model &wm, int update_scheme)
{
    wm.update_scheme = update_scheme;
}

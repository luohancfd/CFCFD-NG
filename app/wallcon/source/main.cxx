#include <iostream>

#include "solid_cell.hh" //cells,interfaces,vertecies
#include "solid_bc.hh"
#include "solid_block.hh" //block

#include "output.hh"
#include "initialise.hh"
#include "e3conn.hh"
#include <string>

using namespace std;

int main()
{
    cout << "Start\n" << endl;
    
    // Wall_model * wall_0 = initialise_wall_model("solid.config");
    Wall_model * wall_0 = new Wall_model();

    //This is for the standalne version only
    initialise_stand_alone(*wall_0, "solid.config");
    
    // // Run for n minutes => t0=0, tf= n*60
    double t = wall_0->t;
    double dt = wall_0->dt;
    int when_to_write = wall_0->when_to_write;
    int n = wall_0->n;
    int tindx = wall_0->tindx;

    create_directory("solid/temperature");
    
    cout << "\nNumber of iterations to be performed: " << n << "\n" << endl;
    cout << "Writing every " << when_to_write << " steps. \n" << endl;
    
    
    cout << "Starting iterations..." << endl;
    
    int counter = 0;
    int counter2 = tindx;


    if (tindx == 0) {
	write_solution(*wall_0, t, tindx);	
    };

    string source_file;

    //Iterations
    for (int i = 0; i < (n+1); i++){
	
	if (counter== when_to_write){
	    counter2 += 1;
	    write_solution(*wall_0, t, counter2);    
	    counter = 0;
	}

	//This is done for time varying source terms
	//BCs are update in respective boundary condition files
	//set block iteration counter
	wall_0->my_block->time_varying_iteration = i; //need this to communicate with BCs
	if (wall_0->my_block->time_varying_terms_flag == 1){
	    source_file = wall_0->my_block->time_varying_terms_dir + "/source/" + to_string( i );
	    read_source_terms_from_file(*(wall_0->my_block), source_file);
	}
	// cout<<"Cell[0][0] source " << wall_0->my_block->block_cells[0][0].source << endl;


	//     //Spatial fluxes
        wall_0->my_block->update_boundary_conditions();
        wall_0->my_block->update_boundary_secondary_interface_temperatures();
        wall_0->my_block->update_internal_secondary_interface_temperatures();
        wall_0->my_block->update_vertex_derivatives();
        wall_0->my_block->update_interface_fluxes();
        wall_0->my_block->update_cell_energy_derivative();
        //Time step
        wall_0->my_block->time_update(dt, wall_0->update_scheme);
        t += dt;
	counter += 1;

	

    }
    
    cout << "Finished iterations\n" << endl;
    
    
    cout << "\nFinish" << endl;

    // cout << "cell[0][0].iface[0].print_tangent() ";
    // wall_0->my_block->block_cells[0][0].iface[0]->print_tangent();
    
    // cout << "cell[0][0].iface[1].print_tangent() ";
    // wall_0->my_block->block_cells[0][0].iface[1]->print_tangent();

    // cout << "cell[0][0].iface[2].print_tangent() ";
    // wall_0->my_block->block_cells[0][0].iface[2]->print_tangent();
    
    // cout << "cell[0][0].iface[3].print_tangent() ";
    // wall_0->my_block->block_cells[0][0].iface[3]->print_tangent();


    return 0;
}


#include "initialise.hh"

#include "solid_cell.hh" //cells,interfaces,vertecies
#include "solid_bc.hh"
#include "useful_funcs.hh"
#include "e3conn.hh"

#include <iostream>
#include <stdlib.h>
#include <fstream> //class to read from files
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstdlib>

int set_block_boundary(SolidBlock &blk, std::string line_in) {
    //Parse line_in and create memory in block for valid boundary condition types.
    //Aids in allowing multiblock in future.
    std::vector<std::string> line = strip_string(line_in);
    
    //Valid boundary conditions
    if ( line[0] ==  "BC_FIXED_T" ) {
	if ( line[1] != "T" ) { 
	    std::cout << "Something went wrong adding BC_FIXED_T.";
	    std::cout << "Check spelling and format: BC_FIXED_T T #\n Bailing out" << std::endl;
	    exit(1);
	}
	//Create assign memory for bc.
	blk.block_bc.push_back(new BC_FIXED_T(atof(line[2].c_str())));
    }
    else if (line[0] == "BC_ADIABATIC") {
	// Need to implement som error catching here line[1] will not exsist
	// if ( line[1] != "" ) {
	//     std::cout << "Something went wrong adding BC_ADIABATIC. Check spelling. #\n";
	//     std::cout << "Bailing out\n";
	//     exit(1);
	// }
	blk.block_bc.push_back(new BC_ADIABATIC);
    }
    else if (line[0] == "BC_CONVECTION") {
	if ( !(line[1] == "H" && line[3] == "T") ) {
	    std::cout << "Something went wrong adding BC_CONVECTION.";
	    std::cout << "Check spelling and format: BC_CONVECTION H # T # \n";
	    std::cout << "Bailing out\n";
	    exit(1);
	}
	blk.block_bc.push_back(new BC_CONVECTION( atof(line[2].c_str()), atof(line[4].c_str())));
    }
    else if (line[0] == "BC_RADIATION") {
	if ( !(line[1] == "EPSILON" && line[3] == "T") ) {
	    std::cout << "Something went wrong adding BC_RADIATION.";
	    std::cout << "Check spelling and format: BC_RADIATION EPSILON # T # \n";
	    std::cout << "Bailing out\n";
	    exit(1);
	}
	blk.block_bc.push_back(new BC_RADIATION( atof(line[2].c_str()), atof(line[4].c_str())));
    }
    else if (line[0] == "BC_CONV_RAD") {
	if ( !(line[1] == "EPSILON" && line[3] =="H" && line[5] == "T") ) {
	    std::cout << "Something went wrong adding BC_CONV_RAD.";
	    std::cout << "Check spelling and format: BC_CONV_RAD EPSILON # H # T # \n";
	    std::cout << "Bailing out\n";
	    exit(1);
	}
	blk.block_bc.push_back(new BC_CONV_RAD( atof(line[2].c_str()), atof(line[4].c_str()), atof(line[6].c_str()) ));
    }
    else if ( line[0] == "BC_FIXED_Q" ) {
	blk.block_bc.push_back( new BC_FIXED_Q( atof(line[2].c_str())));
	if ( line[1] != "Q" ) {
	    std::cout << "Something went wrong adding BC_FIXED_Q.";
	    std::cout << "Check spelling and format: BC_FIXED_Q Q #\n";
	    std::cout << "Bailing out\n";
	    exit(1);
	}
    }
    else if ( line[0] == "BC_E3CONNECTION" ) {
	blk.block_bc.push_back( new BC_E3CONNECTION() );
	if ( line[1] != "south" ) {
	    std::cout << "Something went wrong adding BC_E3CONNECTION.";
	    std::cout << "Check spelling and format: BC_E3CONNECTION south\n";
	    std::cout << "Bailing out\n";
	    exit(1);
	}
    }
    //This has been superseeded by eilmer contoling the connection type
    // else if ( line[0] == "BC_E3CONNECTION_TEMPERATURE" ) {
    // 	blk.block_bc.push_back( new BC_E3CONNECTION() );
    // 	blk.e3connection_type_flag[0] = 1;
    // 	if ( line[1] != "south" ) {
    // 	    std::cout << "Something went wrong adding BC_E3CONNECTION_TEMPERATURE.";
    // 	    std::cout << "Check spelling and format: BC_E3CONNECTION_TEMPERATURE south\n";
    // 	    std::cout << "Bailing out\n";
    // 	    exit(1);
    // 	}
    // }
    else if ( line[0] == "BC_TIMEVARYING_Q" ) {
	if ( line[1] != "path" ) {
	    std::cout << "Something went wrong adding BC_TIMEVARYING_Q.";
	    std::cout << "Check spelling and format: BC_TIMEVARYING_Q path \"/path_to_file\"\n";
	    std::cout << "Bailing out\n";
	    exit(1);
	}
	//Read file with fluxes \n delimiter format, no extenstion on file
	std::ifstream user_file(line[2].c_str());
	if (!(user_file.is_open())){
	    std::cout << "In initialise.cxx Failed to open user file defining q at path "<< line[2] << std::endl;
	    exit(1);
	}
	//Create vector to pass to boundary condition. This is not time varying.
	std::vector<double> user_q;
	std::string line;
	while ( !(user_file.eof()) ){
	    getline(user_file,line);
	    user_q.push_back( atof(line.c_str()));
	}
	//Add boundary condition with flux vector
	blk.block_bc.push_back( new BC_TIMEVARYING_Q( user_q )) ;
    }
    else if ( line[0] == "BC_TIMEVARYING" ) {
	if ( line[1] != "path" ) {
	    std::cout << "Something went wrong adding BC_TIMEVARYING.";
	    std::cout << "Check spelling and format: BC_TIMEVARYING path \"/path_to_file\"\n";
	    std::cout << "Bailing out\n";
	    exit(1);
	}
	//Read file with fluxes \n delimiter format, no extenstion on file
	std::ifstream user_file(line[2].c_str());
	if (!(user_file.is_open())){
	    std::cout << "In initialise.cxx Failed to open user for bc_timevarying file at path "<< line[2] << std::endl;
	    exit(1);
	}
	//Create vector to pass to boundary condition. This is not time varying.
	std::vector<double> user_q;
	std::string line;
	while ( !(user_file.eof()) ){
	    getline(user_file,line);
	    user_q.push_back( atof(line.c_str()));
	}
	//Add boundary condition with flux vector
	blk.block_bc.push_back( new BC_TIMEVARYING( user_q )) ;
    }


    else if ( line[0] == "BC_USERDEF_Q" ) {
	if ( line[1] != "path" ) {
	    std::cout << "Something went wrong adding BC_USERDEF_Q.";
	    std::cout << "Check spelling and format: BC_USERDEF_Q path \"/path_to_file\"\n";
	    std::cout << "Bailing out\n";
	    exit(1);
	}
	//Read file with fluxes \n delimiter format, no extenstion on file
	std::ifstream user_file(line[2].c_str());
	if (!(user_file.is_open())){
	    std::cout << "In initialise.cxx Failed to open user file defining q at path "<< line[2] << std::endl;
	    exit(1);
	}
	//Create vector to pass to boundary condition. This is not time varying.
	std::vector<double> user_q;
	std::string line;
	while ( !(user_file.eof()) ){
	    getline(user_file,line);
	    user_q.push_back( atof(line.c_str()));
	}
	//Add boundary condition with flux vector
	blk.block_bc.push_back( new BC_USERDEF_Q( user_q )) ;
    }
    else if ( line[0] == "BC_USERDEF_T" ) {
	if ( line[1] != "path" ) {
	    std::cout << "Something went wrong adding BC_USERDEF_T.";
	    std::cout << "Check spelling and format: BC_USERDEF_T path \"/path_to_file\"\n";
	    std::cout << "Bailing out\n";
	    exit(1);
	}
	//Read file with temperatures \n delimiter format, no extenstion on file
	//path is line[2]
	std::ifstream user_file(line[2].c_str());
	if (!(user_file.is_open())){
	    std::cout << "Failed to open user file defining T" << std::endl;
	    exit(1);
	};
	//Create vector
	std::vector<double> user_T;
	std::string line;
	while ( !(user_file.eof()) ){
	    getline(user_file,line);
	    user_T.push_back( atof(line.c_str()));
	}
	//Add boundary contion
	blk.block_bc.push_back( new BC_USERDEF_T( user_T )) ;
    }

    //Any other case    
    else {
	std::cout << "Wallcon boundary condition undefined: " << line[0] << std::endl;
	std::cout << "Bailing out\n" << std::endl;
	exit(1);
    }

    return 0;
}


int read_source_terms_from_file(SolidBlock &blk, std::string fname){

    // intended file format...
    //   i j x y source
    //   0 0 0.5 0.5 1000.0
    //   ...
    int i,j;
    double epsx,epsy;
    std::vector<std::string> line; //vector holding stripped line
    std::string string_line;

    //Open file
    // std::cout << "Reading source terms" << std::endl;
    std::ifstream myfile (fname.c_str());

    if(!(myfile.is_open())){
	std::cout << "Unable to open source term file at path: "<< fname << std::endl;
	exit(1);
    }

    getline(myfile, string_line); //Throw away first line

    getline(myfile, string_line);
    while ( !( myfile.eof() ) ) {
	line = strip_string(string_line);

	//Check source term location is correct
	i = atoi( line[0].c_str() );
	j = atoi( line[1].c_str() );
	epsx = blk.block_cells[i][j].pos[0] - atof(line[2].c_str());
	epsy = blk.block_cells[i][j].pos[1] - atof(line[3].c_str());
	if ( abs( epsx ) > 0.0000001 || abs( epsy ) > 0.0000001 ){ //this should be accurate enough
	    	    std::cout << "Source position not alinged with cell center at i=" << line[0] << " j=" << line[1] << std::endl;
	    	    std::cout << "Source pos = (" << line[2] << ", " << line[3] << ")" << std::endl;
	    	    std::cout << "Cell pos = (" << blk.block_cells[i][j].pos[0] << ", " << blk.block_cells[i][j].pos[1] << ")" << std::endl;
	    	    std::cout << "Bailing out\n";
	    exit(1);
	}

	//Fill cell source
	blk.block_cells[i][j].source = atof( line[4].c_str() );
	getline( myfile, string_line ); //Get next line here make sure no seg faul if eof
    }
    myfile.close();
    return 0;
}

int read_initial_state_from_file(SolidBlock &blk, std::string fname){
    //Check to if user whishes to prefill initial condition.
    //Fills initial state with energy and temperature.
    //Options to fill other variables are below but note that rho, cp , k will not work as formulation currently uses block properties not cell properties.

    int ii, jj, kk;
    std::string string_line, path;
    std::vector<std::string> line;
	    
    //Open initialisation file
    std::ifstream init_file (fname.c_str());
    if(!(init_file.is_open())){
	std::cout << "Unable to open initialise_from_file file at path: " << path << std::endl;
	exit(1);
    }
		

    //Parse file
    getline(init_file, string_line); //throw first line
    getline(init_file, string_line); //Second line is number of cells
    line = strip_string(string_line);
    
    //Check size of domain (initialise should already have been run
    ii = atoi ( line[0].c_str() );
    jj = atoi ( line[1].c_str() );
    kk = atoi ( line[2].c_str() );
    
    //Check dimensions
    if ( kk > 1 ){
	std::cout << "Initialisation using prefill currently only supports 2D but in file i, j, k = " << ii <<", " << jj << ", " << kk << std::endl;
	exit(1);
    }
    else if ( !( ii == blk.nni - 1 )){
	std::cout << "Initialisation using prefill mismatch in cells. Number i cells in file " << ii << " but in grid "<<  blk.nni - 1 << std::endl;
	exit(1);
		}
    else if ( !( jj == blk.nnj - 1 )){
	std::cout << "Initialisation using prefill mismatch in cells. Number j cells in file " << jj << " but in grid "<<  blk.nnj - 1 << std::endl;
	exit(1);
    }
    
    //Set properties.
    for ( int j = 0; j < jj; j++ ) {
	for ( int i = 0; i < ii; i++ ){
	    //Expected format
	    //pos.x pos.y pos.z  e T

	    getline(init_file, string_line);
	    line = strip_string(string_line);
	    //Only need e and T
	    blk.block_cells[i][j].e = atof( line[3].c_str() );
	    blk.block_cells[i][j].T = atof( line[4].c_str() );
	}
    }
    init_file.close();
    return 0;
}

int read_number_of_nodes_from_file(SolidBlock &blk, std::string filename){
    //Parse in node numbers to block.
    std::string string_line;
    std::vector<std::string> line;

    std::ifstream myfile (filename.c_str());
    if (!(myfile.is_open())) {
		    std::cout << " Unable to open grid file" <<std::endl;
	exit(1);
    }

    getline (myfile,string_line); //read fist line
    line = strip_string(string_line);
    blk.nni = atoi(line[0].c_str()); //assign number of nodes to block
    blk.nnj = atoi(line[1].c_str());
    blk.nnk = atoi(line[2].c_str());
    myfile.close();
    return 0;
}

int set_block_vertices_from_file(SolidBlock &blk, std::string filename){
    //Read in primary verticies from eilmer grid 
    std::string string_line;
    std::vector<std::string> line;
    
    std::ifstream myfile (filename.c_str());//ifstream myfile ("../simple_rect.grid.b0000.t0000");
    if (!(myfile.is_open())) {
		    std::cout << "Unable to open grid file" << std::endl;
	exit(1);
    }

    getline (myfile,string_line); //throw away first line
    
    for (int j= 0; j<blk.nnj; j++) {
	for (int i=0; i<blk.nni; i++){
	    getline (myfile,string_line);
	    line = strip_string(string_line);
	    blk.block_vertices[i][j].pos[0] = atof(line[0].c_str());
	    blk.block_vertices[i][j].pos[1] = atof(line[1].c_str());
	}
    }
    myfile.close();
    return 0;
}

SolidBlock* initialise_block(std::string fname){
    //This function parses the solid.config file intialises block.

    SolidBlock*  block = new SolidBlock();

    std::string string_line;
    std::vector<std::string> line;

    std::string source_terms_path, init_state_path;
    int source_flag = 0;
    int init_state_flag = 0;
    int anisotropic_set = 0;
    int material_set = 0;
    int current_face = 0;

    //Try to open .onfig file
    std::ifstream myfile (fname.c_str());
    if(!(myfile.is_open())){
	std::cout << "Unable to open solid.config file" << std::endl;
	exit(1);
    };

    //Parse through until hit a defblock
    getline(myfile, string_line);
    while ( !(myfile.eof()) ){
	line = strip_string(string_line);

	//Catch an error
	if (line[0]=="material"){
	    std::cout << "You seem to have identified a material outside a block.\n" \
		      << "This was from a previous version of wallcon.\n" \
		      << "Please now identify properties in the defblock section.\n"\
		      << "Bailing out!" <<std::endl;
	    exit(1);
	}

	if (line[0] == "defblock"){
	    //But for now make sure only one block
	    if ( atoi(line[1].c_str() ) > 0 ){
		std::cout << "At the moment wallcon can only support 1 block." << std::endl;
		std::cout << "Use \"defblock 0\"." << std::endl;
		exit(1);
	    }
	    //Set the block id for when implementing multiblocks
	    block->id = atoi( line[1].c_str() );
	    

	    //While still in that same block. Initialise...
	    while ( !(line[0] == "endblock") ){
		getline(myfile, string_line);
		line = strip_string(string_line);
		
		//Catch an error
		if (line[0]=="material"){
		    std::cout << "You seem to have identified a material inside a block.\n" \
			      << "This was from a previous version of wallcon.\n" \
			      << "Please now identify properties in the defblock section.\n"\
			      << "Bailing out!" <<std::endl;
		    exit(1);
		}

				
		//Check anisotropic_flag set
		if (line[0] == "anisotropic_flag"){
		    //Make sure 1 or 0
		    if ( !(line[1] == "1") && !(line[1]== "0" ) ){
			std::cout << "Set \"anisotropic_flag 0\" or \"anisotropic_flag 1\"."<<std::endl;
			exit(1);
		    }
		    block->anisotropic_flag = atoi( line[1].c_str() );
		}

		//Check time_varying_terms_flag
		if (line[0] == "time_varying_terms_flag"){
		    //Make sure 1 or 0
		    if ( !(line[1] == "1") && !(line[1]== "0" ) ){
			std::cout << "wallcon: Set \"time_varying_terms_flag 0\" or \"time_varying_terms_flag 1\"."<<std::endl;
			exit(1);
		    }
		    block->time_varying_terms_flag = atoi( line[1].c_str() );
		}
		//Get time_varying_terms_dir
		if (line[0] == "time_varying_terms_dir"){
		    //Make sure time_varying_terms_flag
		    if ( block->time_varying_terms_flag == 0 ){
			std::cout << "wallcon: time_varying_terms_dir given but time_varying_terms_flag is 0. \n Set this first." << std::endl;
			exit(1);
		    }
		    block->time_varying_terms_dir = line[1].c_str();
		}

		
		
		//Get materials
		if (line[0] == "k"){
		    block->k = atof( line[1].c_str() );
		    material_set += 1;
		}
		if (line[0] == "rho"){
		    block->rho = atof( line[1].c_str() );
		    material_set += 1;
		}
		if (line[0] == "cp"){
		    block->cp = atof( line[1].c_str() );
		    material_set += 1;
		}
		if (line[0] == "k11"){
		    block->k11 = atof( line[1].c_str() );
		    anisotropic_set += 1;
		}
		if (line[0] == "k12"){
		    block->k12 = atof( line[1].c_str() );
		    anisotropic_set += 1;
		}
		if (line[0] == "k22"){
		    block->k22 = atof( line[1].c_str() );
		    anisotropic_set += 1;
		}
		
		//Initial temp
		if ( line[0] == "T_init" ){block->T_init = atof( line[1].c_str() );}
		
		//Check for source terms for each block
		if ( line[0] == "source" ){
		    source_terms_path = line[2];
		    source_flag = 1;   
		}
		
		//Check for initial state block
		if ( line[0] == "initialise_from_file" ){
		    init_state_path = line[2];
		    init_state_flag = 1;
		}
		
		
		//Boundary conditions in order and consecutivley
		if ( line[0] == "defblock_face" && line[1] == "south"){
		    if ( !(current_face == 0) ){
			std::cout << "Define boundaries in order and consecutively, south, east, north, west. \n Bailing out" << std::endl; 
			exit(1);
		    }
		    getline(myfile, string_line);
		    set_block_boundary(*block, string_line);
		    current_face += 1;
		}
		if ( line[0] == "defblock_face" && line[1] == "east"){
		    if ( !(current_face == 1) ) {
			std::cout << "Define boundaries in order and consecutively, south, east, north, west. \n Bailing out" << std::endl;
			exit(1); 
		    }
		    getline(myfile, string_line);
		    set_block_boundary(*block, string_line);
		    current_face += 1;
		}
		if ( line[0] == "defblock_face" && line[1] == "north"){
		    if ( !(current_face == 2) ){
			std::cout << "Define boundaries in order and consecutively, south, east, north, west. \n Bailing out" << std::endl; 
			exit(1);
		    }
		    getline(myfile, string_line);
		    set_block_boundary(*block, string_line);
		    current_face += 1;
		}
		if ( line[0] == "defblock_face" && line[1] == "west"){
		    if ( !(current_face == 3) ) {
			std::cout << "Define boundaries in order, south, east, north, west. \n Bailing out" << std::endl; 
			exit(1);
		    }
		    getline(myfile, string_line);
		    set_block_boundary(*block, string_line);
		    current_face += 1;
		}
	    }
	}
	getline(myfile, string_line);
    }

    //Check materials if anisotropic
    if (block->anisotropic_flag == 0 && material_set < 3){
	std::cout << "Something went wrong. Number of material defined properties < 3." << std::endl;
	exit(1);
    }
    //Check anisotropic condition
    if (block->anisotropic_flag == 1 && (anisotropic_set < 3 || material_set < 2) ){
	std::cout << "\"anisotropic_flag 1\" but k11,k12,k22 not defined correctly." << std::endl;
	exit(1);
    }
    myfile.close();

    //Read in grid
    std::string grid_file;
    std::stringstream block_id_for_file;
    block_id_for_file << block->id;

    //Allocate memory
    read_number_of_nodes_from_file(*block, "solid/grid/b" + block_id_for_file.str() ); 
    block->allocate_memory();
    set_block_vertices_from_file(*block, "solid/grid/b" + block_id_for_file.str() );
    
    //Assign all the properties to the block object block
    block->assign_ifaces_to_block();
    block->assign_cells_to_block();
    block->assign_ifaces_to_cells();
    block->assign_secondary_ifaces();
    block->assign_internal_secondary_geometry_to_vertices();
    block->assign_properties_to_cells();
    block->initialise_cells(block->T_init); //done before boundary conditions
    
    //After we have the block object. Read in coressponding source terms
    if (source_flag == 1){
	read_source_terms_from_file(*block, source_terms_path.c_str() );
	source_flag = 0;
    }

    //After we have the block object. Read in coressponding initial state
    if (init_state_flag == 1){
	read_initial_state_from_file(*block, init_state_path.c_str() );
	init_state_flag = 0;
    }
    return block;
}



int initialise_tindx_solution(Wall_model & wall){
    // Check to see if user wants to run from last solution and load in last solution if so.

    //To Do: Should really cross-reference to eilmer .times file to be sure as eilmer and wallcon write out seperatley

    //read .config file get tindx or default to time = 0.0
    //idea is to set tindx_flag to desired restart or last = -1
    int tindx_flag = 0; //Default initialisation
    std::string string_line;
    std::vector<std::string> line;
    std::string fname = "solid.config";
    std::ifstream myfile (fname.c_str());
    
    if(!(myfile.is_open())){
	std::cout << "Unable to open solid.config file" << std::endl;
	std::cout << "In initialise_tindx_solution" << std::endl;
	exit(1);
    }

    getline(myfile, string_line);
    line = strip_string(string_line);
    
    //Parse file
    if ( line[0] == "defglobal" ){
	while ( !( line[0] == "endglobal") ){
	    getline(myfile, string_line);
	    line = strip_string(string_line);

	    if (line[0] == "tindx") {
		wall.tindx = atoi ( line[1].c_str()); // set restarting tindx
		tindx_flag = 1; //tell initialiser to restart
	    }
	}
    }
    myfile.close();
    

    if (tindx_flag == 1) {
	std::string path_times = "solid/solid.times";
	std::ifstream solid_times(path_times.c_str());
	int tindx_check = 0;

	if( !(solid_times.is_open()) ){
	    std::cout << "wallcon: Could not open solid.times" << std::endl;
	    exit(1);
	}
	    
	//Check tindx solution has been written
	getline(solid_times, string_line);
	while (!( solid_times.eof() )) {
	    line = strip_string(string_line);
	    if (wall.tindx == atoi(line[0].c_str())){
		tindx_check = 1;
	    }
	    getline(solid_times, string_line);	    
	}
	solid_times.close();
		
	//Now check it can be restarted
	if ( !(tindx_check == 1) ){
	    std::cout << "wallcon: The restart tindx " << wall.tindx << " was not in solid.times" << std::endl;
	    std::cout << "\nBailing out!" << std::endl;
	    exit(1);
	}
	
	//If it gets to here. It should exist. Go to the time directory and initialise.
	//To Do: Would have to do each block here for multi-block
	
	stringstream pad_tindx;
	pad_tindx << setfill('0') << setw(4) << wall.tindx; // pad it with zeros
	string path_to_last = "solid/temperature/t" + pad_tindx.str() + "/b0000.t" + pad_tindx.str();
	std::ifstream last_solution( path_to_last.c_str() );
	
	if(!(last_solution.is_open())){
	    std::cout << "wallcon: Unable to open tindx " << wall.tindx << " solution file at path: "<< path_to_last << std::endl;
	    exit(1);
	}
	
	//Parse and set the temperatures time and energies etc
	getline(last_solution, string_line);
	line = strip_string(string_line);
	wall.t = atof ( line[0].c_str() );
	getline(last_solution, string_line);
	getline(last_solution, string_line); //throw away first three lines
	    
	//"pos.x" "pos.y" "pos.z" "volume" "rho" "cp" "k" "source"  "e" "T"
	for (int j=0; j < wall.my_block->nnj-1; j++){
	    for (int i=0; i < wall.my_block->nni-1; i++){
		getline( last_solution, string_line );
		line =strip_string(string_line);
		
		//Only set previous energy and temp so can change the source if need be
		wall.my_block->block_cells[i][j].e = atof( line[8].c_str() );
		wall.my_block->block_cells[i][j].T = atof( line[9].c_str() );
	    }
	}
	write_solution(wall, 1000, 9999);
	last_solution.close();
    }
    return 0;
}


int initialise_tindx_restart_from_eilmer(Wall_model & wall, int start_tindx){

    std::string string_line;
    std::vector<std::string> line;

    //Set wall tindx
    wall.tindx = start_tindx;

    //First check if user is trying to tindx in solid
    std::string fname = "solid.config";
    std::ifstream myfile (fname.c_str());
    if(!(myfile.is_open())){
	std::cout << "Unable to open solid.config file" << std::endl;
	std::cout << "In initialise_tindx_solution" << std::endl;
	exit(1);
    }
    getline(myfile, string_line);
    line = strip_string(string_line);
    //Parse file
    if ( line[0] == "defglobal" ){
	while ( !( line[0] == "endglobal") ){
	    getline(myfile, string_line);
	    line = strip_string(string_line);

	    if (line[0] == "tindx") {
		std::cout << "You are trying to restart using tindx in wallcon but also specified it in Eilmer.\nPlease only specify tindx through Eilmer." << std::endl;
		std::cout << "Bailing out!" << std::endl;
		exit(1);
	    }
	}
    }
    myfile.close();


    //Check start_t is in solid.times 
    std::string path_times = "solid/solid.times";
    std::ifstream solid_times(path_times.c_str());
    int tindx_check = 0;

    if( !(solid_times.is_open()) ){
	std::cout << "wallcon: Could not open solid.times" << std::endl;
	exit(1);
    }
	    
    //Check tindx solution has been written
    getline(solid_times, string_line);
    while (!( solid_times.eof() )) {
	line = strip_string(string_line);
	if (wall.tindx == atoi(line[0].c_str())){
	    tindx_check = 1;
	}
	getline(solid_times, string_line);	    
    }
    solid_times.close();
		
    //Now check it can be restarted
    if ( !(tindx_check == 1) ){
	std::cout << "wallcon: The restart tindx " << wall.tindx << " was not in solid.times" << std::endl;
	std::cout << "\nBailing out!" << std::endl;
	exit(1);
    }
	
    //If it gets to here. It should exist. Go to the time directory and initialise.
    //Set all the formating
    stringstream pad_tindx;
    pad_tindx << setfill('0') << setw(4) << wall.tindx; // pad it with zeros
    string path_to_tindx = "solid/temperature/t" + pad_tindx.str() + "/b0000.t" + pad_tindx.str();
    std::ifstream tindx_solution( path_to_tindx.c_str() );
	
    if(!(tindx_solution.is_open())){
	std::cout << "wallcon: Unable to open tindx " << wall.tindx << " solution file at path: "<< path_to_tindx << std::endl;
	exit(1);
    }
	
    //Parse and set the temperatures time and energies etc
    getline(tindx_solution, string_line);
    line = strip_string(string_line);
    wall.t = atof ( line[0].c_str() );
    getline(tindx_solution, string_line);
    getline(tindx_solution, string_line); //throw away first three lines
	    
    //"pos.x" "pos.y" "pos.z" "volume" "rho" "cp" "k" "source"  "e" "T"
    for (int j=0; j < wall.my_block->nnj-1; j++){
	for (int i=0; i < wall.my_block->nni-1; i++){
	    getline( tindx_solution, string_line );
	    line =strip_string(string_line);
		
	    //Only set previous energy and temp so can change the source if need be
	    wall.my_block->block_cells[i][j].e = atof( line[8].c_str() );
	    wall.my_block->block_cells[i][j].T = atof( line[9].c_str() );
	}
    }
    tindx_solution.close();

    return 0;
}



int initialise_stand_alone(Wall_model &wall, std::string fname){    
    //Initialise the global time and iteration variables if eilmer is not controling
    
    //Read .config file
    std::ifstream myfile (fname.c_str());
    std::string string_line;
    std::vector<std::string> line;
    int count=0;

    if(!(myfile.is_open())){
	std::cout << "Unable to open solid.config file" << std::endl;
	std::cout << "path is " << fname << std::endl;
	exit(1);
    };
 
    wall.t = 0.0; //Starting time 0.0 by default
    wall.tindx = 0; //Starting iterations 0 by default
    wall.update_scheme = 0; //Default to Euler update scheme change this later 

    getline(myfile, string_line);
    while (!(myfile.eof())){
	line = strip_string(string_line);
    	
	//Initialise iterations info
	if ( line[0] == "defglobal" ){
	    while ( !( line[0] =="endglobal") ){
		getline(myfile, string_line);
		line= strip_string(string_line);
		if (line[0] == "dt") {wall.dt = atof( line[1].c_str() );count +=1; continue;}
		if (line[0] == "when_to_write") {wall.when_to_write = atoi( line[1].c_str() );count += 1; continue;}
		if (line[0] == "n_iterations") {wall.n = atoi( line[1].c_str() );count+=1; continue;}
		if (line[0] == "update_scheme") {
		    if (line[1] == "euler"){
			wall.update_scheme = 0;
			continue;
		    }
		    else if (line[1] == "predictor_corrector" || line[1] == "predictor-corrector" || line[1] == "pc"){
			wall.update_scheme = 1;
			continue;
		    }
		    else{
			std::cout << "Wallcon: unrecognised update scheme" << std::endl;
			exit(1);
		    }
		}
	    }
	}
	getline(myfile, string_line);
    }
    
    if (!(count == 3)) {
	std::cout << "Please initialise global time variables first: dt, when_to_write, n_interations." << std::endl;
	exit(1);
    }

    myfile.close();
    //Init block
    wall.my_block = initialise_block(fname);
    //INcase restart is specified
    initialise_tindx_solution(wall);
    return 0;
}

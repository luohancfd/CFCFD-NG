#ifndef INITIALISE_HH
#define INITIALISE_HH

#include "solid_block.hh"
#include "e3conn.hh" //need the wall model structure from here
#include <string>

SolidBlock * initialise_block(std::string fname);

int set_block_boundary(SolidBlock &blk, std::string line_in);
int read_source_terms_from_file(SolidBlock &blk, std::string fname);
int read_initial_state_from_file(SolidBlock &blk, std::string fname);
int read_number_of_nodes_from_file(SolidBlock &blk, std::string filename);
int set_block_vertices_from_file(SolidBlock &blk, std::string filname);
SolidBlock* initialise_block(std::string fname);
int initialise_tindx_solution(Wall_model &wall);
int initialise_stand_alone(Wall_model &wall, std::string fname);

#endif

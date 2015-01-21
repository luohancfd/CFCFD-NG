#ifndef EILMER3_CONN_HH
#define EILMER3_CONN_HH

#include "../../../lib/geometry2/source/geom.hh"
#include "../../eilmer3/source/cht_coupling_enum.hh"
#include "solid_block.hh"
#include <vector>

//Structure to contain SolidBlock
struct Wall_model {
    SolidBlock* my_block;
    
    //Global variables to control iterations for both stand alone and connection
    double t;
    double dt;
    int n;
    int when_to_write;
    int tindx;
    int eilmer_tindx; //Incase tindx from eilmer is passed
    cht_coupling_t cht_coupling;
    int update_scheme;
};

Wall_model* initialise_wall_model(std::string fname, cht_coupling_t cht_coupling, int start_tindx);

void get_near_wall_solid_cell_pos(Wall_model &wm, std::vector<Vector3> &solid_cell_pos);
int write_solution(Wall_model &wm, double time_elapsed, int print_number);
int gather_near_wall_solid_data(Wall_model &wm, std::vector<double> &T_near_wall, std::vector<double> &k_near_wall);
int update_state(Wall_model &wm, double dt, const std::vector<double> q_wall, const std::vector<double> T_wall);
int initialise_tindx_restart_from_eilmer(Wall_model & wall, int start_tindx);
void set_wallcon_time_update_scheme(Wall_model &wm, int update_scheme);


#endif


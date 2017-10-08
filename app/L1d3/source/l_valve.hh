// l_diaph.hh
// Diaphragm house-keeping.

#ifndef L_VALVE_HH
#define L_VALVE_HH

#include <stdio.h>

// Valves are used to control the interaction of the gas-slugs.
class ValveData {
public:
    double sim_time;        /* current simulation time */

    // Status information.
    int is_open;         /* 1=open, 0=closed, 2=opening     */
    //double P_open;        /*  pressure that valve opens in Pa    */
    double open_period;     /* time it takes for valev to open */
    double open_time;       /* time afer which valve starts to open */

    //int already_blended;    /* 1=indicates that the slugs have been blended 
	//		       after valve opens   */
    //double blend_dx;        /* distance over which slug data is blended */
    //double blend_delay;     /* time after valve opening that blending is applied */

    // Neighbour information.
    int left_slug_id, right_slug_id;
    int left_slug_end_id, right_slug_end_id;
    double left_slug_dx, right_slug_dx;
    // Neighbouring gas slug identifiers.
    // xxxx-slug_id    : number of the adjoining gas slug
    // xxxx-slug_end_id: LEFT or RIGHT end adjoins
    // xxxx-slug_dx    : the sampling distance for gas pressure

    ValveData(int indx, std::string config_file_name, int echo_input=0);
    ValveData(const ValveData& vd);
    ~ValveData();
    int read_state(FILE* infile);
    int write_state(FILE* outfile);
}; // end class ValveData


#endif

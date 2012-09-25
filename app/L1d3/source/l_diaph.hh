// l_diaph.hh
// Diaphragm house-keeping.

#ifndef L_DIAPH_HH
#define L_DIAPH_HH

#include <stdio.h>

// Diaphragms are used to control the interaction of the gas-slugs.
class DiaphragmData {
public:
    double sim_time;        /* current simulation time */

    // Status information.
    int is_burst;           /* 1=burst, 0=intact       */
    double P_burst;         /* burst pressure in Pa    */
    double hold_period;     /* time for which diaphragm will hold */
                            /* after the arrival of the rupture   */
                            /* trigger.                           */
    double trigger_time;    /* time at which rupture is triggered */

    int already_blended;    /* 1=indicates that the slugs have been blended 
			       after diaphragn burst   */
    double blend_dx;        /* distance over which slug data is blended */
    double blend_delay;     /* time after rupture that blending is applied */

    // Neighbour information.
    int left_slug_id, right_slug_id;
    int left_slug_end_id, right_slug_end_id;
    double left_slug_dx, right_slug_dx;
    // Neighbouring gas slug identifiers.
    // xxxx-slug_id    : number of the adjoining gas slug
    // xxxx-slug_end_id: LEFT or RIGHT end adjoins
    // xxxx-slug_dx    : the sampling distance for gas pressure

    DiaphragmData(int indx, std::string config_file_name, int echo_input=0);
    ~DiaphragmData();
    int read_state(FILE* infile);
    int write_state(FILE* outfile);
}; // end class DiaphragmData


#endif

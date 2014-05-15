// l_piston.hh

#ifndef L_PISTON_HH
#define L_PISTON_HH

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "l1d.hh"
#include <iostream>
#include <fstream>

// We may have zero or more pistons pushing against the gas.
class PistonData {
public:
    // Piston properties. 
    double mass;        /* piston mass in kg                */
    double mass_old;    /* old piston mass in kg (for decay)*/
    double diam;        /* effective diameter in m          */
    double area;        /* area of front/back face in m**2  */
    double length;      /* piston length, m                 */
    double front_seal_f;/* friction coefficient             */ 
    double front_seal_area; /* m**2 for gas-facing area     */
    double back_seal_f; /* friction coefficient             */ 
    double back_seal_area;  /* m**2 for gas-facing area     */
    int with_brakes;    /* 1=brakes available, 0=no brakes  */
    double p_restrain;  /* pressure up to which piston is   */
                        /* restrained, Pa                   */
    double x_buffer;    /* x-position of the buffer, m      */
    double V_buffer;    /* speed at impact, m/s             */
    double f_decay;     /* Mass decay time-constant, 1/s    */
    double mass_limit;  /* Mass decay limit, kg             */
    int apply_decay;    /* mass decay flag                  */

    // Piston state 
    int brakes_on;      /* Flag for the piston brakes       */
    int brakes_on_old;
    int is_restrain;    /* Flag for initial restraint       */
                        /* 1 = restrained, 0 = free         */	
                        /* 2 = read-in trajectory from "piston_trajectory.dat"  */						
    int is_restrain_old;
                        /* 1 = restrained, 0 = free         */
    int hit_buffer_count; /* recording buffer strikes       */
    int hit_buffer_count_old;
    int on_buffer;      /* flag for buffer constraint       */
    int on_buffer_old;
                        /* 0 = not constrained by buffer    */
                        /* 1 = is constrained by  buffer    */
    double x0;          /* initial position of centroid     */
    double V0;          /* initial velocity, m/s            */

    double sim_time;    /* present simulation time, s       */
    double x;           /* present position, m              */
    double V;           /* present velocity, m/s            */
    double x_old, V_old; /* values at start of time-step    */

    double Pf;          /* pressure on front face, N/m**2   */
    double Pb;          /* pressure on back face, N/m**2    */
    double Friction;    /* friction force, N                */
    std::vector<double> Ptime;  /* piston timing, m            */
    std::vector<double> Pvel;   /* piston velocity, m/s        */
    std::vector<double> Pacc;   /* piston acceleration, m/s^2  */	
		
    // Time derivatives. 
    double dt;          /* current time step, s             */
    double DxDt[NL];    /* velocity, m/s                    */
    double DVDt[NL];    /* acceleration, m/s**2             */

    // Neighbour information.
    int left_slug_id, right_slug_id;
    int left_slug_end_id, right_slug_end_id;
    // Neighbouring gas slug identifiers.
    // xxxx-slug_id    : number of the adjoining gas slug
    // xxxx-slug_end_id: LEFT or RIGHT end adjoins

    PistonData(int indx, double dt_init, std::string config_file_name, int echo_input=0);
    PistonData(const PistonData& pd);
    ~PistonData();
    int read_state(FILE* infile);
    int write_state(FILE* outfile);
    int record_state(void);
    int restore_state(void);
    int time_derivatives(int time_level, double sim_time);
    int predictor_step(void);
    int corrector_step(void);
}; // end class PistonData

#endif

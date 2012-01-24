/**
 * \file rsp.hh
 *
 * \brief Header for RSP functions and structures
 * 
 * \author DF Potter
 * \version 19-Nov-2006, c++ version for cfcfd2 repository.
 */

#ifndef RSP_HH
#define RSP_HH

/* Internally used structures and functions only declared here */

struct pg_int_state {
double pstar;
double pstar_d;
double astar;
double v_rwt;
double v_cs;
double v_sw;
};

double shock_wave( struct L_flow_state * QI,
                   double ustar,
                   struct pg_int_state * QIstar_pg,
	           int left );

int rarefaction( struct L_flow_state * QI,
                 double ustar,
		 int left,
                 struct pg_int_state * QIstar_pg );

int perfect_gas ( struct L_flow_state * QL,
                  struct L_flow_state * QR,
                  struct pg_int_state * S_pg );

double shock_energy ( double Us,
                      double v_i,
                      double P_i,
                      double rho_i,
                      double P_star,
                      struct L_flow_state * Qi );

double shock_jump ( double P_star,
                    struct L_flow_state * QI,
                    int shock_right,
                    struct L_flow_state * QIstar );

int rare_jump ( double P_star,
                struct L_flow_state * QI,
                int rare_left,
	        struct L_flow_state * QIstar);

double f_zero ( double P_star,
                struct L_flow_state * QL,
                struct L_flow_state * QR );

int final_states ( double P_star,
                   struct L_flow_state * QL,
                   struct L_flow_state * QR,
                   struct L_flow_state * QLstar,
                   struct L_flow_state * QRstar,
                   double * Us_L,
                   double * Us_R );

int step_through_rw( struct L_flow_state * QI,
                     struct L_flow_state * QF,
                     struct ersp_solution * base_sol,
                     int * count,
                     int left,
                     double x_mid,
                     double time);

int sign_check  ( double a,
                  double b );

int fs2bs ( double * x_mid,
            double * time,
            double * u,
            struct L_flow_state * Q,
            struct ersp_solution * base_sol );

#endif

/* end file erf_header.hh */

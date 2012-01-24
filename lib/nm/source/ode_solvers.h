/** \file ode_solvers.h
 *
 * \date 15-Apr-04
 * \author Rowan J Gollan
 * \version 1.0
 *
 */

#ifndef ODE_HEADER_ALREADY_INCLUDED

#include "newt_raph.h"
#include "vecmat.h"
#include "../../util/source/useful.h"

/* ODE debug */
#define DBG_ODE  0


/* Control the verbosity of the warnings */
#define OFW   0
/* OFW =   0 : no warning when ode solver fails
 *         1 : warning on final failure
 */

/* Control the timestep selection after failure */
#define STRICT_TIMESTEP 1

/* STRICT_TIMESTEP =   0 : not strict (limited by MIN_FACTOR)
 *                     1 : strict, as predicted by selection routine
 */

#define SIMPLE_EULER  1
#define SAIM          2
#define QSS           3
#define IMP_EULER     4
#define VODE          5 
#define RK1           6
#define RK2           7
#define RKF           8

#define ALLOW_TIMESTEP_TO_CHANGE 0
#define KEEP_TIMESTEP_CONSTANT 1
#define USE_SUGGESTED_DT 2

#define ZERO_EPS 1.0e-50
#define MIN_CONC 1.0e-30
#define ZERO_TIME 

#define MAX_FACTOR 100.0
#define MIN_FACTOR 1.0e-4
#define DT_INCREASE_FACTOR 1.1
#define DT_REDUCTION_FACTOR 0.333

/*---- Related to the alpha-QSS method ---*/
#define QSS_MAX_CORR 4
#define QSS_C 1.1
#define QSS_EPS2 1.0e-5
#define QSS_EPS1 (QSS_C * QSS_EPS2) 

#define SIMPLE 1
#define MOTT 2

#define STEP_SELECTION SIMPLE

#define SUCCESS 0
#define SYS_FAIL    -1
#define NUM_FAIL    -2

/*---- Related to the rkf2 method ----*/
/* c1t = 0.5*(1 - 1.0/sqrt(3.0)) */
#define c1t 0.211324865
/* c2t = 0.5*(1 + 1.0/sqrt(3.0)) */
#define c2t 0.788675134
/* c1y = 0.5*(0.5 - 1.0/sqrt(3.0)) */
#define c1y -0.038675134
/* c2y = 0.5*(0.5 + 1.0/sqrt(3.0)) */
#define c2y 0.538675134

struct odeOpt {
    int dt_control;
    int first_step;
    int select_timestep;
    double total_t_interval;
    double tolerance;
};

struct odeMem {

    int ndim;
    struct odeOpt *ode_opt;

    /* Function pointers */
    int (*ode_fun) (int ndim, double *y, double *ydot, void *fdata, double *Q, double *L, int option);
    int (*ode_Jfun) (int ndim, double *y, double **dfdy, void *Jdata);
    double (*ode_ts_select) (int ndim, double *y, double *y_dot, void *fdata, double *Q, double *L);
    int (*ode_step) (struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg);
    int (*ode_system_test) (int ndim, double *yout, void *fdata);
 

    /* Data structure pointers */
    /* Additional data that may be required by user-defined functions */
    void *fdata;
    void *Jdata;
    void *tsdata;

    /* Working arrays */
    double *y_dot, *y0, *y_p, *y_c, *y_save;
    double *q0, *q_p, *q_tilde;
    double *p0, *p_p, *p_bar;
    double *a0, *a_bar;
    double *k0, *k1, *k2;
    double *k3, *k4, *k5, *k6, *R;
    double *k1_old, *k2_old;
    double *tmp_a, *tmp_b, *tmp_c, *tmp_d;
    double **dfdy, **hdfdy, **I, **lhs; 

    /* A pointer needed to use the Newton-Raphson solver */
    struct nrMem *nrP;
    double dt; /* A !current! record of the time-step */
               /* Do NOT use this unless you have set it */
               /* yourself previously */

    /* Error logging */
    FILE *fp;
};

struct odeResults {
    int no_steps;
    double dt_final;
};

#ifdef TEST_ODE
struct testOde {
    double **A;
};
#endif

/* ode_solvers.c */

struct odeMem *ode_init(int ndim, int method, FILE *errFp,
			int (*ode_fun)(int, double *, double *, void *, double *, double *, int),
			int (*ode_Jfun)(int, double *, double **, void *),
			double (*ode_ts_select)(int, double *, double *, void *, double *, double *),
			int (*ode_system_test) (int, double *, void *),
			void *fdata, void *Jdata, void *tsdata, struct odeOpt *ode_opt);

int ode_solve(struct odeMem *ode_mem, double *yin, double *yout, double t0, double tf, double dt_guess, struct odeResults *results);
double ode_set_timestep(struct odeMem *ode_mem, double *yin, double t0, double tf, double dt_guess);
int ode_time_stepping(struct odeMem *ode_mem, double *yin, double *yout, double t0, double tf, double dt, struct odeResults *ode_results);
double dt_based_on_success(int dt_control, double dt_used, double dt_sugg, double tdiff);
double dt_based_on_failure(double dt_used, double dt_sugg, double tdiff);
int ode_set_new_data(struct odeMem *ode_mem, void *fdata, void *Jdata, void *tsdata);
int euler_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg);
int saim_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg);
int qss_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg);
inline int real_p(int ndim, double *y, double *p);
inline int calculate_alpha(int ndim, double dt, double *p, double *a);
inline int calculate_p_bar(int ndim, double *P0, double *P_p, double *P_bar);
inline int calculate_q_tilde(int ndim, double *Q0, double *Q_p, double *A_bar,
			     double *Q_tilde);
inline int test_converged(int ndim, double *Y_c, double *Y_p);
inline double qss_ts_select(int ndim, double dt_old, double *Y_c, double *Y_p, int converged);
int imp_euler_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg);
int nrf(int ndim, double *y, double *G, void *fdata);
int nrJac(int ndim, double *y, double **dGdy, void *Jdata);
int vode_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg);
int rk1_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg);
int rk2_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg);
int rkf_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg);


#ifdef TEST_ODE
int f1(int ndim, double *y, double *ydot, void *ode_data, double *C, double *D, int option);
int Jf1(int ndim, double *y, double **dfdy, void *Jode_data);
int f2(int ndim, double *y, double *ydot, void *ode_data, double *C, double *D, int option);
int Jf2(int ndim, double *y, double **dfdy, void *Jode_data);
#endif

#define ODE_HEADER_ALREADY_INCLUDED
#endif

/** \file ode_solvers.c
 * \ingroup nm
 * \brief An interface for some ODE solvers.
 *
 * This module provides an interface to some numerical ODE solvers.
 * The interface allows all ODE solvers to be called in a similar fashion.
 *
 * \date 12-Aug-04
 * \author Rowan J Gollan
 * \version 1.0
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ode_solvers.h"

#define y0   (ode_mem->y0)
#define y_dot (ode_mem->y_dot)
#define y_p  (ode_mem->y_p)
#define y_c  (ode_mem->y_c)
#define y_save (ode_mem->y_save)
#define q0   (ode_mem->q0)
#define q_p  (ode_mem->q_p)
#define q_tilde  (ode_mem->q_tilde)
#define p0   (ode_mem->p0)
#define p_p  (ode_mem->p_p)
#define p_bar (ode_mem->p_bar)
#define a0   (ode_mem->a0)
#define a_bar (ode_mem->a_bar)
#define k0   (ode_mem->k0)
#define k1   (ode_mem->k1)
#define k2   (ode_mem->k2)
#define k3   (ode_mem->k3)
#define k4   (ode_mem->k4)
#define k5   (ode_mem->k5)
#define k6   (ode_mem->k6)
#define k1_old   (ode_mem->k1_old)
#define k2_old   (ode_mem->k2_old)
#define tmp_a (ode_mem->tmp_a)
#define tmp_b (ode_mem->tmp_b)
#define tmp_c (ode_mem->tmp_c)
#define tmp_d (ode_mem->tmp_d)
#define R (ode_mem->R)
#define hdfdy (ode_mem->hdfdy)
#define I    (ode_mem->I)
#define fRhs (ode_mem->ode_fun)
#define lhs  (ode_mem->lhs)
#define Jac  (ode_mem->ode_Jfun)
#define ts_select (ode_mem->ode_ts_select)
#define step (ode_mem->ode_step)
#define system_test (ode_mem->ode_system_test)
#define fData (ode_mem->fdata)
#define JData (ode_mem->Jdata)
#define tsData (ode_mem->tsdata)
#define fp    (ode_mem->fp)


struct odeMem* ode_init(int ndim, int method, FILE *errFp,
			int (*ode_fun) (int, double*, double*, void*, double*, double*, int),
			int (*ode_Jfun) (int, double*, double**, void*),
			double (*ode_ts_select) (int, double *, double *, void *, double *, double *),
			int (*ode_system_test) (int, double *, void *),
			void *fdata, void *Jdata, void *tsdata, struct odeOpt *ode_opt)
{
#   define TOL 1.0e-8
    struct odeMem *ode_mem;
    ode_mem = (struct odeMem *) malloc(sizeof(struct odeMem));

    if (errFp == NULL) {
	fp = stdout;
    } else {
	fp = errFp;
    }

    ode_mem->ndim = ndim;
    ode_mem->ode_opt = ode_opt;

    /* Allocate memory for temporary work arrays
     * (within the ode_mem structure).
     */ 

    y_dot = mem_alloc_vector(ndim);
    y0 = mem_alloc_vector(ndim);
    y_c = mem_alloc_vector(ndim);
    y_p = mem_alloc_vector(ndim);
    y_save = mem_alloc_vector(ndim);
    q0 = mem_alloc_vector(ndim);
    q_p = mem_alloc_vector(ndim);
    q_tilde = mem_alloc_vector(ndim);
    p0 = mem_alloc_vector(ndim);
    p_p = mem_alloc_vector(ndim);
    p_bar = mem_alloc_vector(ndim);
    a0 = mem_alloc_vector(ndim);
    a_bar = mem_alloc_vector(ndim);
    k0 = mem_alloc_vector(ndim);
    k1 = mem_alloc_vector(ndim);
    k2 = mem_alloc_vector(ndim);
    k3 = mem_alloc_vector(ndim);
    k4 = mem_alloc_vector(ndim);
    k5 = mem_alloc_vector(ndim);
    k6 = mem_alloc_vector(ndim);
    k1_old = mem_alloc_vector(ndim);
    k2_old = mem_alloc_vector(ndim);
    R = mem_alloc_vector(ndim);
    tmp_a = mem_alloc_vector(ndim);
    tmp_b = mem_alloc_vector(ndim);
    tmp_c = mem_alloc_vector(ndim);
    tmp_d = mem_alloc_vector(ndim);
    ode_mem->dfdy = mem_alloc_sqmat(ndim);
    hdfdy = mem_alloc_sqmat(ndim);
    I = mem_alloc_sqmat(ndim);
    lhs = mem_alloc_sqmat(ndim);

    /* Set the function pointers */
    fRhs = ode_fun;
    Jac = ode_Jfun;
    ts_select = ode_ts_select;
    system_test = ode_system_test;

    if ( method == SIMPLE_EULER ) { 
 	step = euler_step; 
    } else if ( method == SAIM ) {
	step = saim_step;
    } else if ( method == QSS ) {
	step = qss_step;
    } else if ( method == IMP_EULER ) {
	step = imp_euler_step;
	/* We need to do some extra things to use the 
	 * Newton-Raphson solver.
	 */
	ode_mem->nrP = newt_raph_init(ndim, nrf, nrJac, ode_mem, ode_mem, TOL);
	eye(I, ndim); /* Set I to the identity matrix for all time. */

    } else if ( method == VODE ) {
	step = vode_step;
    } else if ( method == RK1 ) {
	step = rk1_step;
	eye(I, ndim);
    } else if ( method == RK2 ) {
	step = rk2_step;
    } else if ( method == RKF ) {
	step = rkf_step;
    } else {
	fprintf(fp, "Integration method: %d is not known \n", method);
	return (NULL);
    }

    /* Set the pointers for any data structures required
     * by the user-supplied functions.
     */

    fData = fdata;
    JData = Jdata;
    tsData = tsdata;

#   if DBG_ODE >= 1
    printf("DBG_ODE-1: ode_init() has proceeded without problems.\n");
#   endif

#   undef TOL
    return (struct odeMem *) ode_mem;
}

int ode_solve(struct odeMem *ode_mem, double *yin, double *yout, double t0,
	      double tf, double dt_guess, struct odeResults *results)
{

#   define MAX_ATTEMPT 4

    int  ndim, flag, attempt;
    double dt, dt_temp, tdiff;

    ndim = ode_mem->ndim;
    ode_mem->ode_opt->first_step = 1;

    /* 1. Check the time bounds */

    tdiff = tf - t0;
    if ( tdiff <= 0.0 ) {
	fprintf(fp, "The finishing time is earlier than the starting time\n");
	fprintf(fp, "for the ODE integration. t0 = %e, tf = %e \n", t0, tf);
	return (SYS_FAIL);
    }

    ode_mem->ode_opt->total_t_interval = tdiff;

    /* 2. Set what to do about timestep */
    
    dt = ode_set_timestep(ode_mem, yin, t0, tf, dt_guess);

#   if DBG_ODE >= 1
    printf("DBG_ODE-1: t0= %g  tf= %g  tdiff= %g \n", t0, tf, tdiff);
    printf("DBG_ODE-1: dt= %g (selected by ode_set_timestep() )\n", dt);
#   endif

    /* 3. Try to solve the system */
    /* Loop through the ODE time-stepping,
     * starting with a given dt.
     */

    /* Keep a copy because yin will change */
    copy_vector(yin, y_save, ndim);

#   if DBG_ODE >= 2
    printf("DBG_ODE-2: Check that the y_save vector is correct \n");
    print_vector(ndim, yin, "yin");
    print_vector(ndim, y_save, "y_save");
#   endif

    for (attempt = 0; attempt <= MAX_ATTEMPT; attempt++) {
#       if DBG_ODE >= 2
        printf("DBG_ODE-2: Before ode_time_stepping(), attempt= %d, dt= %g\n", attempt, dt);
	print_vector(ndim, yin, "yin");
        print_vector(ndim, yout, "yout");
#       endif

	flag = ode_time_stepping(ode_mem, yin, yout, t0, tf, dt, results);

#       if DBG_ODE >= 2
        printf("DBG_ODE-2: After ode_time_stepping(), attempt= %d, flag= %d\n", attempt, flag);
	print_vector(ndim, yin, "yin");
        print_vector(ndim, yout, "yout");
#       endif

	if (flag == SUCCESS) {
#           if DBG_ODE >= 2
	    printf("DBG_ODE-2: On attempt= %d, the ode_time_stepping was successful.\n", attempt);
#           endif
	    break;
	}
	else {
	    /* We might have failed due to the system test failure OR
	     * the ODE time-stepping algorithm failed.
	     */

	    /* If we are going to try again, we need the original yin */
	    copy_vector(y_save, yin, ndim);

#           if DBG_ODE >= 2
            printf("DBG_ODE-2: Attempt= %d failed.  Check that yin is reset. \n", attempt);
	    print_vector(ndim, yin, "yin");
	    print_vector(ndim, y_save, "y_save");
#           endif

	    /* Do something about the timestep before we try again. */
	    if (attempt == 0) {
		/* Maybe we can use the suggestion of the user-defined
		 * timestep selection algorithm
		 */
#               if DBG_ODE >= 2
		printf("DBG_ODE-2: On attempt= %d, we'll try to use some smarts to select dt.\n", attempt);
#               endif

		fRhs(ndim, yin, y_dot, fData, q0, p0, 1);
		if (ts_select != NULL) {
		    dt_temp = ts_select(ndim, yin, y_dot, fData, q0, p0);
#                   if DBG_ODE >= 2
		    printf("DBG_ODE-2: dt_temp= %g (from ts_select) \n", dt_temp);
#                   endif
		    /* Now sometimes the selection can be quite harsh */
#                   if (STRICT_TIMESTEP == 1)
		    dt = dt_temp;
#                       if DBG_ODE >= 2
		    printf("DBG_ODE-2: STRICT_TIMESTEP enabled, dt = %g\n", dt);
#                       endif
#                   elif
		    if (dt_temp  < (MIN_FACTOR * dt) ) {
			/* Instead we will try to limit the reduction to MIN_FACTOR */
			dt = dt * MIN_FACTOR;
#                       if DBG_ODE >= 2
			printf("DBG_ODE-2: dt is too small, limited by MIN_FACTOR to= %g \n", dt);
#                       endif
		    } else {
			dt = dt_temp;
#                       if DBG_ODE >= 2
			printf("DBG_ODE-2: dt seems reasonale= %g \n", dt);
#                       endif
		    }
#                   endif 
		}
		else {
		    dt *= DT_REDUCTION_FACTOR;
#                   if DBG_ODE >= 2
		    printf("DBG_ODE-2: No specific ts_select routine, therefore= %g (from ts_select) \n", dt);
#                   endif
		}
		    
	    } else {
		/* We just reduce by a given factor */
		dt *= DT_REDUCTION_FACTOR;
#               if DBG_ODE >= 2
		printf("DBG_ODE-2: On attempt= %d, we just reduce by the reduction factor, dt= %g.\n", attempt, dt);
#               endif

	    }

	}
	
    }

    if (attempt > MAX_ATTEMPT) {
#       if (OFW >= 1)
	printf("The ODE integrator tried %d times to reduce the timestep but still failed.\n", MAX_ATTEMPT);
#       endif
	return (flag);
    }
	
    /* 6. If we made it this far, we return a success! */

#   undef MAX_ATTEMPT

    return(SUCCESS);
}

double ode_set_timestep(struct odeMem *ode_mem, double *yin, double t0,
			double tf, double dt_guess)
{

#   define DT_DUMMY 1.0e-6

    int ndim;
    double dt, tdiff;


    ndim = ode_mem->ndim;
    tdiff = tf - t0;

    /* Check that a dt_guess has been given
     * By convention, any value greater than zero
     * will be assumed to be a supplied guess.
     * Users would usually set dt_guess to -1.0
     * to force timestep selection.
     */
    if (dt_guess > 0.0) {
	dt = dt_guess;
    } else {
	/* We need to do some sort of timestep selection.*/
	if (ts_select != NULL) {
	    fRhs(ndim, yin, y_dot, fData, q0, p0, 1);
	    dt = ts_select(ndim, yin, y_dot, fData, q0, p0);
	} else {
            /* We need to find a better solution here but at the moment
	     * it is too much implementation effort to try and code a
	     * generic all encompassing method to select a time step
	     * for an ODE (that is a black art in itself).
	     */
	    dt = DT_DUMMY;
	}
    }

    if (tdiff < dt) {
	/* The timestep is greater than the desired integration 
	 * interval, therefore set dt to this interval.
	 */
	dt = tdiff;
    }

#   undef DT_DUMMY

    return dt;

}

int ode_time_stepping(struct odeMem *ode_mem, double *yin, double *yout,
		       double t0, double tf, double dt, struct odeResults *ode_results)
{

#   define MAX_ITER 4
    int i, j, ndim, flag, steps, test_result;
    double time, dt_store, dt_sugg, frac, t_old, delta_y, tdiff;

    ndim = ode_mem->ndim;
    time = t0;
    t_old = t0;
    copy_vector(yin, yout, ndim);
    steps = 0;
    tdiff = tf - t0;
    dt_store = tdiff;

    while (time < tf) { /* Not ideal to test floating points - but logic below should help out... */
	copy_vector(yout, yin, ndim);
	t_old = time;

	for (i = 0; i <= MAX_ITER; ++i) {

#           if DBG_ODE >= 3
	    printf("DBG_ODE-3: Before taking an individual step, i= %d, dt= %g\n", i, dt);
	    print_vector(ndim, yin, "yin");
	    print_vector(ndim, yout, "yout");
#           endif

	    flag = step(ode_mem, dt, yin, yout, &dt_sugg);

#           if DBG_ODE >= 3
	    printf("DBG_ODE-3: After taking an individual step, flag= %d, i= %d, dt_sugg= %g\n", flag, i, dt_sugg);
	    print_vector(ndim, yin, "yin");
	    print_vector(ndim, yout, "yout");
	    for (j = 0; j < ndim; j++) {
		printf("yout[%d] - yin[%d] = %e \n", j, j, yout[j] - yin[j]);
	    }
#           endif

	    if (flag == SUCCESS) {
		time += dt;
		if (dt == tdiff) {
		    /* We can leave the routine and NOT update dt_store */
		    break;
		}

#               if DBG_ODE >= 3
		printf("DBG_ODE-3: Success. time= %g  dt= %e \n", time, dt);
#               endif

		dt_store = dt;

		tdiff = tf - time;
		dt = dt_based_on_success(ode_mem->ode_opt->dt_control, dt, dt_sugg, tdiff);
#               if DBG_ODE >= 3
		printf("DBG_ODE-3: New dt (based on success) dt= %e \n", dt);
#               endif

		break;
	    } else { /* The step failed */
		if (dt != tdiff) {
		    dt_store = dt;
		}
		dt = dt_based_on_failure(dt, dt_sugg, tdiff);
	    }
	}
	steps++;

	if (i > MAX_ITER) { /* Then we failed on every attempt */
	    return (NUM_FAIL);
	}

    } /* end of while for time-stepping */

    /* We need to do our linear interpolation before we leave */

    if ( time > tf ) {
	/* Would expect it to be larger.
	 * Therefore, linearly interpolate to get result.
	 */

	frac = ( tf - t_old ) / ( time - t_old );
	for (i = 0; i < ndim; i++) {
	    delta_y = yout[i] - yin[i];
	    yout[i] = yin[i] + frac * delta_y;
	}

    }

    /* If we've made it this far, we've had 'numerical' success, in as much
     * as the stepping scheme has converged.
     * BUT now we need to test that the ODE system passes its test based on
     * the constraints applied by a user-defined function.
     */
    if (system_test != NULL) {
	test_result = system_test(ndim, yout, fData);

    } else {
	test_result = SUCCESS;
    }

    /* 5. Record the results */

    ode_results->no_steps = steps;
    ode_results->dt_final = dt_store;

#   undef MAX_ITER

    return (test_result);

}

double dt_based_on_success(int dt_control, double dt_used, double dt_sugg, double tdiff)
{

    double dt;

    if (dt_sugg > 0.0) {
	/* We have been given a dt_sugg to use by the stepping scheme.
	 * But it reamins to be seen whether or not we should use it.
	 */
	if (dt_sugg < dt_used) {
	    dt = dt_used;
	} else if ( dt_sugg > (MAX_FACTOR * dt_used) ) {
	    /* This prevents jumping by too many orders of magnitude */
	    dt = dt_used * DT_INCREASE_FACTOR;
	} else {
	    /* For all other cases, we can use the suggestion */
	    dt = dt_sugg;
	}
    } else {
	if (dt_control == ALLOW_TIMESTEP_TO_CHANGE) {
	    dt = dt_used * DT_INCREASE_FACTOR;
	} else {
	    /* For all other cases, leave the timestep unchanged. */
	    dt = dt_used;
	}
    }

    /* No point stepping larger than tdiff */
    if (dt > tdiff) {
	dt = tdiff;
    }

    return dt;
}

double dt_based_on_failure(double dt_used, double dt_sugg, double tdiff)
{
    double dt;

    if (dt_sugg > 0.0) {
	/* We have been given a dt_sugg to use by the stepping scheme.
	 * But it remains to be seen if it is realistic.
	 */
	if (dt_sugg > dt_used) {
	    /* There's no point in trying a larger timestep */
	    dt = dt_used * DT_REDUCTION_FACTOR;
	} else if (dt_sugg < (MIN_FACTOR * dt_used) ) {
	    /* This is a warning if the timestep is being altered
	     * too dramatically.
	     */
	    printf("Previous timestep = %e, new suggested timestep = %e \n", dt_used, dt_sugg);
	    /* Instead we will try to limit the reduction to MIN_FACTOR */
	    dt = dt_used * MIN_FACTOR;
	} else {
	    dt = dt_sugg;
	}
    } else { /* There is no suggestion */
	dt = dt_used * DT_REDUCTION_FACTOR;
    }

    /* finally, there's no need to step greater than tdiff */
    if (dt > tdiff)
	dt = tdiff;

    return dt;
}

int ode_set_new_data(struct odeMem *ode_mem, void *fdata,
		      void *Jdata, void *tsdata)
{    
    fData = fdata;
    JData = Jdata;
    tsData = tsdata;

    return 0;
}
    

/*---------------------------------------------------------------------*/
/*                       Various ODE solvers                           */
/*---------------------------------------------------------------------*/


int euler_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg)
{

    int i, ndim;
    double h;    
    ndim = ode_mem->ndim;

    /* 1. Evaluate the function derivative */

    if (ode_mem->ode_opt->first_step == 1) {
	ode_mem->ode_opt->first_step = 0;
	fRhs(ndim, yin, y_dot, fData, q0, p0, 1);
    } else {
	fRhs(ndim, yin, y_dot, fData, q0, p0, 0);
    }
    
    h = dt;

    /* There's no point stepping larger than the required time. */

    for (i = 0; i < ndim; i++) {
	yout[i] = yin[i] + h * y_dot[i];
#       if DBG_ODE >= 4
	printf("DBG_ODE-4: euler_step() dt = %g: yout[%d] = %g  yin[%d] = %g  y_dot[%d] = %g \n",
		dt, i, yout[i], i, yin[i], i, y_dot[i]);
#       endif
    }

    *dt_sugg = 0.0;

    return(SUCCESS);
}

int saim_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg)
{
    UNUSED_VARIABLE(ode_mem);
    UNUSED_VARIABLE(dt);
    UNUSED_VARIABLE(yin);
    UNUSED_VARIABLE(yout);
    UNUSED_VARIABLE(dt_sugg);

    return(SUCCESS);
}


int qss_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg)
{
    int corr, i, ndim, converged;
    double h;

    ndim = ode_mem->ndim;
    copy_vector(yin, y0, ndim);
    converged = NUM_FAIL;
    /* First derivative evaluation! */
    /* Now it is possible to estimate a timestep (if necessary)
     * and perform the predictor step. */
    if (ode_mem->ode_opt->first_step == 1) {
	ode_mem->ode_opt->first_step = 0;
	fRhs(ndim, yin, y_dot, fData, q0, p0, 1);
    } else {
	fRhs(ndim, yin, y_dot, fData, q0, p0, 0);
    }

    h = dt;

    /* --- Predictor Step --- */

    real_p(ndim, yin, p0);
    calculate_alpha(ndim, h, p0, a0);

    for (i = 0; i < ndim; ++i) {
	y_p[i] = y0[i] + ( (h * y_dot[i]) / (1.0 + a0[i] * h * p0[i]));
    }

    
    /* once off to get things started */
    copy_vector(y_p, y_c, ndim);

    /* --- Corrector step(s) --- */
    for (corr = 0; corr <= QSS_MAX_CORR; ++corr) {
	copy_vector(y_c, y_p, ndim); /* Our old corrector is the new predictor */
	fRhs(ndim, y_p, y_dot, fData, q_p, p_p, 0);
	real_p(ndim, y_p, p_p);
	calculate_p_bar(ndim, p0, p_p, p_bar);
	calculate_alpha(ndim, h, p_bar, a_bar);
	calculate_q_tilde(ndim, q0, q_p, a_bar, q_tilde);

	/* Actual corrector */
	for (i = 0; i < ndim; ++i) {
	    y_c[i] = y0[i] + (( h * ( q_tilde[i] - p_bar[i] * y0[i] ) ) / ( 1.0 + a_bar[i] * h * p_bar[i] ));
            if (y_c[i] < MIN_CONC)
		y_c[i] = 0.0;
	}

	converged = test_converged(ndim, y_c, y_p);

	if (converged == SUCCESS) {
	    copy_vector(y_c, yout, ndim);
	    /* printf("Iterations converged BEFORE max_count\n"); */
	    *dt_sugg = qss_ts_select(ndim, h, y_c, y_p, SUCCESS);
	    return (SUCCESS);
	}
    }
    
    /* If we got this far, we've exceeded maximum iterations without
     * converging, therefore, we've failed.
     */
    copy_vector(y_c, yout, ndim);
    *dt_sugg = qss_ts_select(ndim, h, y_c, y_p, NUM_FAIL);
    
    return (NUM_FAIL);
}

inline int real_p(int ndim, double *y, double *p)
{
    int i;
    
    /* We use the ZERO_EPS in case one of the mole fractions
     * is 0.
     */

    for (i = 0; i < ndim; ++i) {
	p[i] = p[i] / (y[i] + ZERO_EPS);
    }
    return 0;
}

inline int calculate_alpha(int ndim, double dt, double *p, double *a)
{
    int i;
    double r, temp_a, temp_b;

    for(i = 0; i < ndim; ++i) {
	r = 1.0 / ( (p[i] * dt) + ZERO_EPS);
	temp_a = 180.0 * r * r * r + 60.0 * r * r + 11.0 * r + 1;
	temp_b = 360.0 * r * r * r + 60.0 * r * r + 12.0 * r + 1;
	a[i] = temp_a / temp_b;
    }

    return 0;
}
	
inline int calculate_p_bar(int ndim, double *P0, double *P_p, double *P_bar)
{
    int i;
 
    for(i = 0; i < ndim; ++i) {
	P_bar[i] = 0.5 * ( P0[i] + P_p[i] );
    }

    return 0;
}

inline int calculate_q_tilde(int ndim, double *Q0, double *Q_p, double *A_bar,
			     double *Q_tilde)
{
    int i;
    
    for (i = 0; i < ndim; ++i) {
	Q_tilde[i] = A_bar[i] * Q_p[i] + (1.0 - A_bar[i]) * Q0[i];
    }
    
    return 0;
}

inline int test_converged(int ndim, double *Y_c, double *Y_p)
{
    int i;
    int flag = 0;
    double test;
    
    test = 0.0;

    for (i = 0; i < ndim; ++i) {
	if (Y_c[i] < MIN_CONC) 
	    continue;
	test = fabs(Y_c[i] - Y_p[i]);
	if (test > (QSS_EPS1 * Y_c[i])) {
	    ++flag;
	}
    }

    if (flag == 0)
	return (SUCCESS);
    else
	return (NUM_FAIL);
}


inline double qss_ts_select(int ndim, double dt_old, double *Y_c, double *Y_p, int converged)
{

#   if (STEP_SELECTION == MOTT)
    int i;
    double dt_new, sigma, test;

    sigma = 0.0;

    for (i = 0; i < ndim; ++i) {
	if (Y_c[i] < MIN_CONC)
	    continue;
	test = ( fabs(Y_c[i] - Y_p[i]) ) / (QSS_EPS2 * Y_c[i]);
	if (test > sigma)
	    sigma = test;
    }

    /* In the following, we proivde some sanity checks on the
     * timestep selection.  Mott's recipe doesn't seem to 
     * work straight from the book.
     */


    /* Check that sigma is not still zero */
    /*    printf("In here, converged = %d\n", converged); */
    if (sigma <= 0.0) {
	if (converged == SUCCESS) {
	    dt_new = dt_old * DT_INCREASE_FACTOR;
	} else {
	    dt_new = dt_old * DT_REDUCTION_FACTOR;
	}
    } else {
	dt_new = dt_old * ( (1.0 / sqrt(sigma)) + 0.005);
        /* Check we've moved in the right direction */
        if (converged == SUCCESS) {
	    if ( dt_new < dt_old ) {
		dt_new = dt_old * DT_INCREASE_FACTOR;
	    } else {
		if ( dt_new > dt_old ) {
		    dt_new = dt_old * DT_REDUCTION_FACTOR;
		}
	    }
	}
    }

    return dt_new;

#   else

    double dt_new;

    UNUSED_VARIABLE(ndim);
    UNUSED_VARIABLE(Y_c);
    UNUSED_VARIABLE(Y_p);

    if (converged == SUCCESS)
	dt_new = dt_old * DT_INCREASE_FACTOR;
    else
	dt_new = dt_old * DT_REDUCTION_FACTOR;

    return dt_new;

#   endif

}


int imp_euler_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg)
{

    /* Make sure ode_mem has current values for these things
     * because they are accessed by the function which
     * the Newton-Raphson solver tries to solve.
     */

    ode_mem->dt = dt;
    copy_vector(yin, y0, ode_mem->ndim);

    newt_raph(ode_mem->nrP, yin, yout);

    *dt_sugg = 0.0;

    return(SUCCESS);
}

int nrf(int ndim, double *y, double *G, void *fdata)
{

    int i;
    double h;
    struct odeMem *ode_mem;

    ode_mem = (struct odeMem *) fdata;
    h = ode_mem->dt;

    /* Evalute all ydots */
    if (ode_mem->ode_opt->first_step == 1) {
	ode_mem->ode_opt->first_step = 0;
	fRhs(ndim, y, y_dot, fData, q0, p0, 1);
    } else {
	fRhs(ndim, y, y_dot, fData, q0, p0, 0);
    }


    for (i = 0; i < ndim; i++) {
	G[i] = y[i] - y0[i] - h * y_dot[i];
    }

    return 0;

}


int nrJac(int ndim, double *y, double **dGdy, void *Jdata)
{

    struct odeMem *ode_mem;
    double h;
    double **dfdy;

    ode_mem = (struct odeMem *) Jdata;
    h = ode_mem->dt;
    dfdy = ode_mem->dfdy;

    Jac(ndim, y, dfdy, JData);
    msys_scalmult(h, dfdy, hdfdy, ndim);
    msys_subtract(dGdy, I, hdfdy, ndim);

    return 0;

}
    

int vode_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg)
{
    UNUSED_VARIABLE(ode_mem);
    UNUSED_VARIABLE(dt);
    UNUSED_VARIABLE(yin);
    UNUSED_VARIABLE(yout);
    UNUSED_VARIABLE(dt_sugg);
    return(SUCCESS);
}

int rk1_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg)
{
    int ndim;
    double **dfdy;

    UNUSED_VARIABLE(dt_sugg);

    ndim = ode_mem->ndim;
    dfdy = ode_mem->dfdy;

    if (ode_mem->ode_opt->first_step == 1) {
	ode_mem->ode_opt->first_step = 0;
	fRhs(ndim, yin, k0, fData, q0, p0, 1);
    } else {
	fRhs(ndim, yin, k0, fData, q0, p0, 0);
    }
    
    Jac(ndim, yin, dfdy, JData);

    msys_scalmult(dt/2.0, dfdy, hdfdy, ndim);
    msys_subtract(lhs, I, hdfdy, ndim);
   
    msys_solve(lhs, k1, k0, ndim);

    scale_vector(ndim, dt, k1, tmp_a);

    add_vectors(ndim, yin, tmp_a, yout);

    *dt_sugg = 0.0;
    return (SUCCESS);

}


int rk2_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg)
{
    int ndim, converged, count, MAX_COUNT;
    double TOL;
    double d1, d2;

    MAX_COUNT = 100;
    TOL = 1.0e-7;
    ndim = ode_mem->ndim;

    if (ode_mem->ode_opt->first_step == 1) {
	ode_mem->ode_opt->first_step = 0;
	fRhs(ndim, yin, y_dot, fData, q0, p0, 1);
    } else {
	fRhs(ndim, yin, y_dot, fData, q0, p0, 0);
    }

    copy_vector(k0, k1, ndim);
    copy_vector(k0, k2, ndim);
    copy_vector(k0, k1_old, ndim);
    copy_vector(k0, k2_old, ndim);
    /* Apply fixed-point iteration to get k1, k2 */
    converged = 0;
    count = 0;

    while (converged != 1 && count < MAX_COUNT) {
	scale_vector(ndim, 0.25, k1, tmp_a);
	scale_vector(ndim, c1y,  k2, tmp_b);
	add_vectors(ndim, tmp_a, tmp_b, tmp_c);
	scale_vector(ndim, dt, tmp_c, tmp_c);
	add_vectors(ndim, yin, tmp_c, tmp_d);
	fRhs(ndim, tmp_d, k1, fData, q0, p0, 0);
	
	scale_vector(ndim, c2y, k1, tmp_a);
	scale_vector(ndim, 0.25,  k2, tmp_b);
	add_vectors(ndim, tmp_a, tmp_b, tmp_c);
	scale_vector(ndim, dt, tmp_c, tmp_c);
	add_vectors(ndim, yin, tmp_c, tmp_d);
	fRhs(ndim, tmp_d, k2, fData, q0, p0, 0);

	subtract_vectors(ndim, k1, k1_old, tmp_a);
	add_scalar_to_vector(ndim, 1.0, k1, tmp_b);
	element_divide(ndim, tmp_a, tmp_b, tmp_c);
	abs_vector(ndim, tmp_c, tmp_d);
	d1 = max_vector(ndim, tmp_d);

	subtract_vectors(ndim, k2, k2_old, tmp_a);
	add_scalar_to_vector(ndim, 1.0, k2, tmp_b);
	element_divide(ndim, tmp_a, tmp_b, tmp_c);
	abs_vector(ndim, tmp_c, tmp_d);
	d2 = max_vector(ndim, tmp_d);

	if (d1 < TOL && d2 < TOL)
	    converged = 1;

	count++;
	copy_vector(k1, k1_old, ndim);
	copy_vector(k2, k2_old, ndim);
    }

    if (count >= MAX_COUNT)
	return (NUM_FAIL);
    
    add_vectors(ndim, k1, k2, tmp_a);
    scale_vector(ndim, 0.5*dt, tmp_a, tmp_b);
    add_vectors(ndim, yin, tmp_b, yout);

    *dt_sugg = 0.0;

    return (SUCCESS);
}


int rkf_step(struct odeMem *ode_mem, double dt, double *yin, double *yout, double *dt_sugg)
{
    int i, ndim;
    double max_R, tol, h, delta;

    ndim = ode_mem->ndim;
    tol = ode_mem->ode_opt->tolerance;

    if (tol <= 0.0)
	tol = 1.0e-8;

    /* 1st function evaulation */

    if (ode_mem->ode_opt->first_step == 1) {
	ode_mem->ode_opt->first_step = 0;
	fRhs(ndim, yin, y_dot, fData, q0, p0, 1);
    } else {
	fRhs(ndim, yin, y_dot, fData, q0, p0, 0);
    }
    
    h = dt;

    for (i = 0; i < ndim; i++) {
	k1[i] = h * y_dot[i];
	tmp_a[i] = yin[i] + (1.0/4.0) * k1[i];
    }

    /* 2nd function evaluation */
    fRhs(ndim, tmp_a, y_dot, fData, q0, p0, 0);

    for (i = 0; i < ndim; i++) {
	k2[i] = h * y_dot[i];
	tmp_a[i] = yin[i] + (3.0/32.0) * k1[i] + (9.0/32.0)*k2[i];
    }

    /* 3rd function evaluation */
    fRhs(ndim, tmp_a, y_dot, fData, q0, p0, 0);
    
    for (i = 0; i < ndim; i++) {
	k3[i] = h * y_dot[i];
	tmp_a[i] = yin[i] + (1932.0/2197.0) * k1[i]
	    - (7200.0/2197.0) * k2[i] + (7296.0/2197.0) * k3[i];
    }

    /* 4th function evaulation */
    fRhs(ndim, tmp_a, y_dot, fData, q0, p0, 0);

    for (i = 0; i < ndim; i++) {
	k4[i] = h * y_dot[i];
	tmp_a[i] = yin[i] + (439.0/216.0) * k1[i]
	    - 8.0 * k2[i] + (3680.0/513.0) * k3[i]
	    - (845.0/4104.0) * k4[i];
    }

    /* 5th function evaluation */
    fRhs(ndim, tmp_a, y_dot, fData, q0, p0, 0);

    for (i = 0; i < ndim; i++) {
	k5[i] = h * y_dot[i];
	tmp_a[i] = yin[i] - (8.0/27.0) * k1[i]
	    + 2.0 * k2[i] - (3544.0/2565.0) * k3[i]
	    + (1859.0/4104.0) * k4[i] - (11.0/40.0) * k5[i];
    }

    /* 6th function evaluation */
    fRhs(ndim, tmp_a, y_dot, fData, q0, p0, 0);

    for (i = 0; i < ndim; i++) {
	k6[i] = h * y_dot[i];
	R[i] = fabs( (1.0/360.0) * k1[i] - (128.0/4275.0) * k3[i]
		     - (2197.0/75240.0) * k4[i] + (1.0/50.0) * k5[i]
		     + (2.0/55.0) * k6[i] ) / h;
    }

    max_R = R[0];
    for (i = 1; i < ndim; i++) {
	if (R[i] > max_R)
	    max_R = R[i];
    }

    delta = 0.84 * pow(tol / max_R, 0.25);

    /* New timestep selection */
    if (delta <= 0.1)
	*dt_sugg = 0.1 * h;
    else if (delta >= 4.0) 
	*dt_sugg = 4.0 * h;
    else
	*dt_sugg = delta * h;

    if (max_R <= tol) { /* Our timestep was OK! */
	for (i = 0; i < ndim; i++) {
	    yout[i] = yin[i] + (25.0/216.0) * k1[i] 
		+ (1408.0/2565.0) * k3[i] + (2197.0/4104.0) * k4[i]
		- (1.0/5.0) * k5[i];
	}
	return (SUCCESS);
    }

    /* If our timestep wasn't OK, we return FAILURE.
     * The outer program loop will decide if we've had
     * too many attempts at timestep reduction.
     */

    return (NUM_FAIL);
}


/*--------------------------------------------------------------------*/
/*                   Testing routines                                 */
/*--------------------------------------------------------------------*/

#ifdef TEST_ODE

int f1(int ndim, double *y, double *ydot, void *ode_data, double *C, double *D, int option)
{

    struct testOde *td;
    td = (struct testOde *) ode_data;

    ydot[0] = td->A[0][0] * y[0] + td->A[0][1] * y[1];
    ydot[1] = td->A[1][0] * y[0] + td->A[1][1] * y[1];

    return 0;

}

int Jf1(int ndim, double *y, double **dfdy, void *Jode_data) 
{ 
    struct testOde *td;
    td = (struct testOde *) Jode_data;

    dfdy[0][0] = td->A[0][0];  dfdy[0][1] = td->A[0][1];
    dfdy[1][0] = td->A[1][0];  dfdy[1][1] = td->A[1][1];

    return 0;
}

int f2(int ndim, double *y, double *ydot, void *ode_data, double *C, double *D, int option)
{

    ydot[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    ydot[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * pow(y[1], 2);
    ydot[2] = 3.0e7 * pow(y[1], 2);

    return 0;

}

int Jf2(int ndim, double *y, double **dfdy, void *Jode_data)
{

    dfdy[0][0] = -0.04;  dfdy[0][1] = 1.0e4 * y[2];  dfdy[0][2] = 1.0e4 * y[1];
    dfdy[1][0] = 0.04;   dfdy[1][1] = -1.0e4 * y[2] - 2.0 * 3.0e7 * y[1]; dfdy[1][2] = -1.0e4 * y[1];
    dfdy[2][0] = 0.0;    dfdy[2][1] = 2.0 * 3.0e7 * y[1]; dfdy[2][2] = 0.0;

    return 0;

}


int main()
{
    int ndim;
    double *yin, *yout, *anl;
    double lambda1, lambda2;
    double t0, tf, dt_guess;
    double t1, t2, t3;

    struct odeResults test_results;
    struct odeMem *odeP;
    struct odeOpt ode_optP;
    struct testOde test_data;


    printf("======================================================\n");
    printf(" Test for: ode_solvers.c                              \n");
    printf("======================================================\n\n");

    printf("--- Test case 1: MECH3750 Assignment Question ---\n");
    printf("\n");

    ndim = 2;
    lambda1 = -0.5;
    lambda2 = -5.0;

    yin = mem_alloc_vector(ndim);
    yout = mem_alloc_vector(ndim);
    test_data.A = mem_alloc_sqmat(ndim);

    test_data.A[0][0] = lambda1;  test_data.A[0][1] = 0.0;
    test_data.A[1][0] = 1.0;      test_data.A[1][1] = lambda2;
    dt_guess = 0.01;
    t0 = 0.0;
    tf = 1.0;
    yin[0] = 1.0; yin[1] = 0.0;

    printf("Solve the 2 x 2 system of ordinary differential equations:\n");
    printf("y' = A y where\n");
    print_sqmat(ndim, test_data.A, "A");
    
    ode_optP.dt_control = KEEP_TIMESTEP_CONSTANT;
    odeP = ode_init(ndim, SIMPLE_EULER, NULL, f1, NULL, NULL, NULL, &(test_data), NULL, NULL, &(ode_optP));
    ode_solve(odeP, yin, yout, t0, tf, dt_guess, &test_results);
    printf("Using the Euler method with timestep = %g : \n", dt_guess);
    print_vector(ndim, yout, "y(x = 1.0)");
    
    odeP = ode_init(ndim, IMP_EULER, NULL, f1, Jf1, NULL, NULL, &(test_data), &(test_data), NULL, &(ode_optP));
    yin[0] = 1.0; yin[1] = 0.0; /* Restart the initial condition */
    ode_solve(odeP, yin, yout, t0, tf, dt_guess, &test_results);
    printf("\nUsing the implicit Euler method with timestep %g :\n", dt_guess);
    print_vector(ndim, yout, "y(x = 1.0)");

    odeP = ode_init(ndim, RK2, NULL, f1, Jf1, NULL, NULL, &(test_data), &(test_data), NULL, &(ode_optP));
    yin[0] = 1.0; yin[1] = 0.0; /* Restart the initial condition */
    ode_solve(odeP, yin, yout, t0, tf, dt_guess, &test_results);
    printf("\nUsing the semi-implicit two-stage Runge-Kutta method with timestep %g :\n", dt_guess);
    print_vector(ndim, yout, "y(x = 1.0)");

    odeP = ode_init(ndim, RK1, NULL, f1, Jf1, NULL, NULL, &(test_data), &(test_data), NULL, &(ode_optP));
    yin[0] = 1.0; yin[1] = 0.0; /* Restart the initial condition */
    ode_solve(odeP, yin, yout, t0, tf, dt_guess, &test_results);
    printf("\nUsing the semi-implicit one-stage Runge-Kutta method with timestep %g :\n", dt_guess);
    print_vector(ndim, yout, "y(x = 1.0)");

    odeP = ode_init(ndim, RKF, NULL, f1, Jf1, NULL, NULL, &(test_data), &(test_data), NULL, &(ode_optP));
    yin[0] = 1.0; yin[1] = 0.0; /* Restart the initial condition */
    odeP->ode_opt->tolerance = 1.0e-8;
    ode_solve(odeP, yin, yout, t0, tf, dt_guess, &test_results);
    printf("\nUsing the Runge-Kutta-Fehlberg method with error control (tol= %g) and timestep  %g :\n",
	   odeP->ode_opt->tolerance, dt_guess);
    print_vector(ndim, yout, "y(x = 1.0)");


    printf("The analytical solution is..\n");
    anl = mem_alloc_vector(ndim);
    anl[0] = 0.6065306;  anl[1] = 0.1332873;
    print_vector(ndim, anl, "analytical");
    printf("DONE\n\n");
    free_sqmat(test_data.A, ndim);
    free(yin);
    free(yout);
    free(anl);


    printf("--- Test case 2: from CVODE example ---\n");
    printf("\n");

    ndim = 3;
    ode_optP.dt_control = KEEP_TIMESTEP_CONSTANT;
    yin = mem_alloc_vector(ndim);
    yout = mem_alloc_vector(ndim);
    t0 = 0.0;
    t1 = 4.0e-1;
    t2 = 4.0;
    t3 = 4.0e1;
    dt_guess = 0.0002;

    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;

    printf("Solve the following chemical kinetics problem..\n");
    printf("dy1/dt = -0.04 * y1 + 1.0e4*y2*y3\n");
    printf("dy2/dt = 0.04 * y1 - 1.e4*y2*y3 - 3.0e7*(y2)^2\n");
    printf("dy3/dt = 3.0e7*(y2)^2\n\n");
    printf("on the interval t = 0.0 to t = 4.0e10 with\n\n");
    print_vector(ndim, yin, "y0");

    odeP = ode_init(ndim, SIMPLE_EULER, NULL, f2, NULL, NULL, NULL, NULL, NULL, NULL, &(ode_optP));
    printf("Using the Euler method with timestep = %g : \n", dt_guess);
    
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t1, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t1, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t2, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t2, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t3, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t3, yout[0], yout[1], yout[2]);


    odeP = ode_init(ndim, IMP_EULER, NULL, f2, Jf2, NULL, NULL, NULL, NULL, NULL, &(ode_optP));
    printf("\nUsing the implicit Euler method with timestep = %g : \n", dt_guess);
    
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t1, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t1, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t2, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t2, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t3, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t3, yout[0], yout[1], yout[2]);

    odeP = ode_init(ndim, RK2, NULL, f2, Jf2, NULL, NULL, NULL, NULL, NULL, &(ode_optP));
    printf("\nUsing the semi-implicit two-stage Runge-Kutta method with timestep = %g : \n", dt_guess);
    
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t1, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t1, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t2, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t2, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t3, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t3, yout[0], yout[1], yout[2]);

    odeP = ode_init(ndim, RK1, NULL, f2, Jf2, NULL, NULL, NULL, NULL, NULL, &(ode_optP));
    printf("\nUsing the semi-implicit one-stage Runge-Kutta method with timestep = %g : \n", dt_guess);
    
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t1, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t1, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t2, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t2, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t3, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t3, yout[0], yout[1], yout[2]);

    odeP = ode_init(ndim, RKF, NULL, f2, Jf2, NULL, NULL, NULL, NULL, NULL, &(ode_optP));
    printf("\nUsing the Runge-Kutta-Fehlberg method with timestep = %g : \n", dt_guess);
    odeP->ode_opt->tolerance = 1.0e-8;
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t1, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t1, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t2, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t2, yout[0], yout[1], yout[2]);
    yin[0] = 1.0; yin[1] = 0.0; yin[2] = 0.0;
    ode_solve(odeP, yin, yout, t0, t3, dt_guess, &test_results);
    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t3, yout[0], yout[1], yout[2]);


    printf("\nThe solution given in the CVODE report is:\n\n");
    
    printf("At t = 4.0000e-01 \t y =  9.851641e-01 \t 3.386242e-05 \t 1.480205e-02 \n");
    printf("At t = 4.0000e+00 \t y =  9.055097e-01 \t 2.240338e-05 \t 9.446793e-02 \n");
    printf("At t = 4.0000e+01 \t y =  7.158016e-01 \t 9.185045e-06 \t 2.841893e-01 \n");

    printf("DONE\n\n");

    printf("======================================================\n");
    printf(" End of test for: ode_solvers.c                       \n");
    printf("======================================================\n");

    return 0;

}

#endif







 

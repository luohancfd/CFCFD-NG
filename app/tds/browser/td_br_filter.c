/** \file td_br_filter.c
 * \ingroup tds
 * \brief Fast(er) filter functions for td_browser.tcl.
 * \author Peter Jacobs 
 * \author Josh Corbett 2004 : Heat transfer filter
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <tcl.h>
#include <blt.h>


/** \brief Filter a vector of data using a moving average.
 *
 * The filtered data is put back into the original BLT vector.
 *
 * \param interp : pointer to the Tcl interpreter
 * \param vname  : name of the BLT vector containing the data
 * \param nhalf  : half-width (in sample points) of the filter
 *
 * \returns -1 if there was a problem, 0 if all went well.
 */
int cfilter( Tcl_Interp *interp, char *vname, int nhalf ) {
    Blt_Vector *vPtr;
    double *tv, sum;
    int    N, count, i, j;

    printf( "cfilter: vname  %s, nhalf  %d\n", vname, nhalf );

    if ( Blt_VectorExists( interp, vname ) ) {
        if ( Blt_GetVector( interp, vname, &vPtr ) != TCL_OK ) {
	    printf( "cfilter: cannot get vector %s\n", vname );
            return -1;
        }
        N = vPtr->numValues;
        tv = malloc (N * sizeof(double) );
        if ( tv == NULL ) {
	    printf( "cfilter: could not allocate work space.\n" );
            return -1;
        }

        /* Start-up phase: compute the first average. */
        sum = 0.0;
        count = 0;
        for ( i = 0; i <= nhalf; ++i ) {
	    sum += vPtr->valueArr[i];
            ++count;
        }
        tv[0] = sum / count;

        /* Do the remaining points one at a time, keeping a running total. */
        for ( j = 1; j < N; ++j ) {
	    if ( j + nhalf < N ) {
	        sum += vPtr->valueArr[j+nhalf];
                ++count;
	    }
            if ( j >= nhalf ) {
	        sum -= vPtr->valueArr[j-nhalf];
                --count;
            }
            tv[j] = sum / count;
        } /* end for */

        /* Copy back into the original vector such that notifications occur. */
        if ( Blt_ResetVector( vPtr, tv, N, N, TCL_DYNAMIC ) != TCL_OK ) {
	    printf( "cfilter: resetting vector failed.\n" );
            return -1;
        }
    } else {
        return -1;
    }
    return 0;
} /* end cfilter() */


/** \brief Compute heat transfer(W/m*) from raw emf (V).
 *
 * The heat-tnasfer data is put back into the original BLT vector.
 *
 * \param interp : pointer to the Tcl interpreter
 * \param vname  : name of the BLT vector containing the surface temperature history
 * \param dt     : time interval between samples
 * \param rhock  : thermal product (units: J/(m^2 degree-K s^0.5)
 * \param E      : ambient voltage of external amplifier (V)
 * \param a      : transducer sensitivity (1/K)
 *
 * \returns -1 if there was a problem, 0 if all went well.
 */
int c_heat_transfer_emf( Tcl_Interp *interp, char *vname, double dt,
		     double rhock, double E, double a) {
    Blt_Vector *vPtr;
    double *tv, *g, *gc, *qdot;
    int    N, j, k;
    double sqrtpi, accumulator;
    sqrtpi = sqrt(M_PI);

    printf( "c_heat_transfer_emf: vname  %s, dt %e, rhock %e, E %e, a %e\n", 
	    vname, dt, rhock, E, a );

    if ( Blt_VectorExists( interp, vname ) ) {
        if ( Blt_GetVector( interp, vname, &vPtr ) != TCL_OK ) {
	    printf( "c_heat_transfer_emf: cannot get vector %s\n", vname );
            return -1;
        }
        N = vPtr->numValues; /*Number of points in vector*/
        tv = malloc (N * sizeof(double) );
        if ( tv == NULL ) {	/*Holds imported emf files*/
	    printf( "c_heat_transfer_emf: could not allocate space for tv vector.\n" );
            return -1;
        }
        g = malloc (N * sizeof(double) );
        if ( g == NULL ) {	/*Holds implulse response*/
	    printf( "c_heat_transfer_emf: could not allocate space for g vector.\n" );
            return -1;
        }
        gc = malloc (N * sizeof(double) );
        if ( gc == NULL ) {	/*Holds manipulated impulse response*/
	    printf( "c_heat_transfer_emf: could not allocate space for gc vector.\n" );
            return -1;
        }
        qdot = malloc (N * sizeof(double) );
        if ( qdot == NULL ) {	/*Holds calculated heat-transfer*/
	    printf( "c_heat_transfer_emf: could not allocate space for qdot vector.\n" );
            return -1;
        }

        /* Convert emf to temperature. */
        for ( j = 0; j < N; ++j ) {
	    tv[j] = vPtr->valueArr[j] / (a * E);
        }
	
	/* Calculate Impulse Response, gir */  
	for ( j = 0; j < N; j++ ) {
	    g[j] = 2.0/(rhock*sqrtpi) * (sqrt((j+1)*dt) - sqrt(j*dt));
	}
      
	/*
	 * Run Deconfilter 
	 */  
	/* Normalise first component of temperature. */
	for ( j = 0; j < N ; j++ ) {
	    tv[j] /= g[0];
	}    
	/* Set up and normalise g */
	for ( j = 0; j < N-1 ; j++ ) {
	    gc[j] = -1.0 * g[j+1] / g[0];
	}
      /* Filter */    
	for ( j = 0; j < N ; j++ ) {
	    accumulator = tv[j];
	    for ( k = 0; k < N-1; k++ ) { 
		if (j > k) accumulator += gc[k] * qdot[j-k-1];
	    }
	    qdot[j] = accumulator;
	}

        /* Copy back into the original vector such that notifications occur. */
        if ( Blt_ResetVector( vPtr, qdot, N, N, TCL_DYNAMIC ) != TCL_OK ) {
	    printf( "c_heat_transfer_emf: resetting vector failed.\n" );
            return -1;
        }
	free(g);
	free(gc);
	free(tv);
	/* qdot freed by the BLT module when it is done (?check?) */
    } else {
        return -1;
    }
    return 0;
} /* end c_heat_transfer_emf() */


/** \brief Compute heat transfer(W/m*) from surface temperature (K).
 *
 * The heat-tnasfer data is put back into the original BLT vector.
 *
 * \param interp : pointer to the Tcl interpreter
 * \param vname  : name of the BLT vector containing the surface temperature history
 * \param dt     : time interval between samples
 * \param rhock  : thermal product (units: J/(m^2 degree-K s^0.5)
 * \param E      : ambient voltage of external amplifier (V)
 * \param a      : transducer sensitivity (1/K)
 *
 * \returns -1 if there was a problem, 0 if all went well.
 */
int c_heat_transfer_temp( Tcl_Interp *interp, char *vname, double dt,
		     double rhock, double E, double a) {
    Blt_Vector *vPtr;
    double *tv, *g, *gc, *qdot;
    int    N, j, k;
    double sqrtpi, accumulator;
    sqrtpi = sqrt(M_PI);

    printf( "c_heat_transfer_temp: vname  %s, dt %e, rhock %e, E %e, a %e\n", 
	    vname, dt, rhock, E, a );

    if ( Blt_VectorExists( interp, vname ) ) {
        if ( Blt_GetVector( interp, vname, &vPtr ) != TCL_OK ) {
	    printf( "c_heat_transfer_temp: cannot get vector %s\n", vname );
            return -1;
        }
        N = vPtr->numValues; /*Number of points in vector*/
        tv = malloc (N * sizeof(double) );
        if ( tv == NULL ) {	/* Holds imported temp */
	    printf( "c_heat_transfer_temp: could not allocate space for tv vector.\n" );
            return -1;
        }
        g = malloc (N * sizeof(double) );
        if ( g == NULL ) { /* Holds impulse response */
	    printf( "c_heat_transfer_temp: could not allocate space for g vector.\n" );
            return -1;
        }
        gc = malloc (N * sizeof(double) );
        if ( gc == NULL ) { /* Holds manipulated impulse response */
	    printf( "c_heat_transfer_temp: could not allocate space for gc vector.\n" );
            return -1;
        }
        qdot = malloc (N * sizeof(double) );
        if ( qdot == NULL ) { /* Holds calculated heat-transfer*/
	    printf( "c_heat_transfer_temp: could not allocate space for qdot vector.\n" );
            return -1;
        }

        /* Import temperature. */
        for ( j = 0; j < N; ++j ) {
	    tv[j] = vPtr->valueArr[j];
        }
	
	/* Calculate Impulse Response, gir */  
	for ( j = 0; j < N; j++ ) {
	    g[j] = 2.0/(rhock*sqrtpi) * (sqrt((j+1)*dt) - sqrt(j*dt));
	}
      
	/*
	 * Run Deconfilter 
	 */  
	/* Normalise first component of temperature. */
	for ( j = 0; j < N ; j++ ) {
	    tv[j] /= g[0];
	}    
	/* Set up and normalise g */
	for ( j = 0; j < N-1 ; j++ ) {
	    gc[j] = -1.0 * g[j+1] / g[0];
	}    
      /* Run filter */
	for ( j = 0; j < N ; j++ ) {
	    accumulator = tv[j];
	    for ( k = 0; k < N-1; k++ ) { 
		if (j > k) accumulator += gc[k] * qdot[j-k-1];
	    }
	    qdot[j] = accumulator;
	}

        /* Copy back into the original vector such that notifications occur. */
        if ( Blt_ResetVector( vPtr, qdot, N, N, TCL_DYNAMIC ) != TCL_OK ) {
	    printf( "c_heat_transfer_temp: resetting vector failed.\n" );
            return -1;
        }
	free(g);
	free(gc);
	free(tv);
	/* qdot freed by the BLT module when it is done (?check?) */
    } else {
        return -1;
    }
    return 0;
} /* end c_heat_transfer_temp() */


/** \brief Compute surface temperature (K) from raw emf (V).
 *
 * The temperature data is put back into the original BLT vector.
 *
 * \param interp : pointer to the Tcl interpreter
 * \param vname  : name of the BLT vector containing the surface temperature history
 * \param dt     : time interval between samples
 * \param rhock  : thermal product (units: J/(m^2 degree-K s^0.5)
 * \param E      : ambient voltage of external amplifier (V)
 * \param a      : transducer sensitivity (1/K)
 *
 * \returns -1 if there was a problem, 0 if all went well.
 */
int c_emf_temp( Tcl_Interp *interp, char *vname, double dt,
		     double rhock, double E, double a) {
    Blt_Vector *vPtr;
    double *tv;
    int    N, j;
    
    
    
    printf( "c_emf_temp: vname  %s, dt %e, rhock %e, E %e, a %e\n", 
	    vname, dt, rhock, E, a );

    if ( Blt_VectorExists( interp, vname ) ) {
        if ( Blt_GetVector( interp, vname, &vPtr ) != TCL_OK ) {
	    printf( "c_heat_transfer_emf: cannot get vector %s\n", vname );
            return -1;
        }
        N = vPtr->numValues; /*Number of points in vector*/
        tv = malloc (N * sizeof(double) );
        if ( tv == NULL ) { /* Holds temperature data */
	    printf( "c_emf_temp: could not allocate space for tv vector.\n" );
            return -1;
        }
        
        /* Convert emf to temperature. */
        for ( j = 0; j < N; ++j ) {
	    tv[j] = vPtr->valueArr[j] / (a * E);
        }
		
        /* Copy back into the original vector such that notifications occur. */
        if ( Blt_ResetVector( vPtr, tv, N, N, TCL_DYNAMIC ) != TCL_OK ) {
	    printf( "c_emf_temp: resetting vector failed.\n" );
            return -1;
        }
	
    } else {
        return -1;
    }
    return 0;
} /* end c_emf_temp() */


/** \brief Compute surface temperature from heat transfer history.
 *
 * The surface-temperature data is put back into the original BLT vector.
 *
 * \param interp : pointer to the Tcl interpreter
 * \param vname  : name of the BLT vector containing the heat-transfer history
 * \param dt     : time interval between samples
 * \param rhock  : thermal product (units: J/(m^2 degree-K s^0.5)
 * \param E      : ambient voltage of external amplifier (V)
 * \param a      : transducer sensitivity (1/K)
 *
 * \returns -1 if there was a problem, 0 if all went well.
 */
int c_surface_temperature( Tcl_Interp *interp, char *vname, double dt,
			   double rhock, double E, double a ) {
    Blt_Vector *vPtr;
    double *tv, *g, *qdot;
    int    N, j, k;
    double sqrtpi, accumulator;
    sqrtpi = sqrt(M_PI);

    printf( "c_surface_temperature: vname  %s, dt %e, rhock %e, E %e, a %e\n", 
	    vname, dt, rhock, E, a );

    if ( Blt_VectorExists( interp, vname ) ) {
        if ( Blt_GetVector( interp, vname, &vPtr ) != TCL_OK ) {
	    printf( "c_surface_temperature: cannot get vector %s\n", vname );
            return -1;
        }
        N = vPtr->numValues; /*Number of points in vector*/
        tv = malloc (N * sizeof(double) );
        if ( tv == NULL ) {	/* Holds calculated temperature data */
	    printf( "c_surface_temperature: could not allocate space for tv vector.\n" );
            return -1;
        }
        g = malloc (N * sizeof(double) );
        if ( g == NULL ) {	/* Holds impulse response */
	    printf( "c_surface_temperature: could not allocate space for g vector.\n" );
            return -1;
        }
        qdot = malloc (N * sizeof(double) );
        if ( qdot == NULL ) { /* Holds input heat-transfer */
	    printf( "c_surface_temperature: could not allocate space for qdot vector.\n" );
            return -1;
        }

	/* Get qdot from BLT vector. */
        for ( j = 0; j < N; ++j ) {
	    qdot[j] = vPtr->valueArr[j];
        }
	
	/* Calculate Impulse Response, gir */  
	for ( j = 0; j < N; j++ ) {
	    g[j] = 2.0/(rhock*sqrtpi) * (sqrt((j+1)*dt) - sqrt(j*dt));
	}    
	      
	/*
	 * Run confilter 
	 */  
	for ( j = 0; j < N ; j++ ) {
	    accumulator = qdot[0] * g[j];
	    for ( k = 0; k < N-1; k++ ) { 
		if (j > k) accumulator += qdot[k+1] * g[j-k-1];
	    }
	    tv[j] = accumulator;
	}

        /* Copy back into the original vector such that notifications occur. */
        if ( Blt_ResetVector( vPtr, tv, N, N, TCL_DYNAMIC ) != TCL_OK ) {
	    printf( "c_surface_temperature: resetting vector failed.\n" );
            return -1;
        }
	free(g);
	free(qdot);
	/* qdot freed by the BLT module when it is done (?check?) */
    } else {
        return -1;
    }
    return 0;
} /* end c_surface_temperature() */

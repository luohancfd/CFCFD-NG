/** \file moc_gasdynamic.c
 * \ingroup imoc
 * \brief Basic Gas-Dynamic relations.
 */

#include <stdio.h>
#include <math.h>
#include "moc_gasdynamic.h"

/*------------------------------------------------------------------*/

/* @function */
double T0_T(double M, double g) {
   /**
    Purpose: Isentropic flow relations. <BR>
    Input  : <BR>
    M : Mach number <BR>
    g : ratio of specific heats <BR>
    Output : <BR>
    Returns ratio of total T over static T. <BR>
    (Available from the Tcl interpreter.) 
    */
   double t0_t;
   t0_t = 1.0 + 0.5 * (g - 1.0) * M * M;
   return t0_t;
} /* end function T0_T */

/*------------------------------------------------------------------*/

/* @function */
double P0_P(double M, double g) {
   /**
    Purpose: Isentropic flow relations. <BR>
    Input  : <BR>
    M : Mach number <BR>
    g : ratio of specific heats <BR>
    Output : <BR>
    Returns ratio of total P over static P. <BR>
    (Available from the Tcl interpreter.) 
    */
   double p0_p;
   p0_p = pow( T0_T(M,g), g/(g-1.0) );
   return p0_p;
} /* end function P0_P */

/*------------------------------------------------------------------*/

/* @function */
double NuFromM(double M, double g) {
   /**
    Purpose: Calculate the Prandtl-Meyer function given Mach number. <BR>
             See equation 4.21b in Liepmann and Roshko 1957. <BR>
    Input  : <BR>
    M : Mach number <BR>
    g : ratio of specific heats <BR>
    Output : <BR>
    Returns Prandtl-Meyer function (Nu) or 0.0 on failure. <BR>
    Note that Nu = 0.0 when M = 1.0 <BR>
    (Available from the Tcl interpreter.) 
    */
   double  temp1, temp2, Nu;

   temp1 = sqrt((g + 1.0) / (g - 1.0));
   temp2 = M * M - 1.0;
   if (temp2 < 0.0) {
      printf( "NuFromM received a subsonic Mach number: %g\n", M );
      Nu = 0.0;
   } else {
      temp2 = sqrt(temp2);
      Nu = temp1 * atan(temp2 / temp1) - atan(temp2);
   } /* end if */

   return Nu;
} /* end function NuFromM */

/*------------------------------------------------------------------*/

/* @function */
double MFromNu( double Nu, double g ) {
   /**
    Purpose: Compute Mach number given Prandtl Meyer function. <BR>
    Input  : <BR>
    Nu : Parndtl-Meyer function (radians) <BR>
    g  : ratio of specific heats <BR>
    Output : <BR>
    Returns Mach number (>= 1.0) or a value of 0.0 on failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   double Nu_0, Nu_1, Nu_2, M_0, M_1, M_2;
   double f_0, f_1, f_2, slope;
   int    count;

   if ( Nu < 0.0 ) {
      return 0.0;
   } /* end if */

   /*
    * Generate an initial guess and, it it is good, return it.
    */
   M_0  = MFromNu_approximate( Nu, g );
   Nu_0 = NuFromM( M_0, g);
   f_0  = Nu - Nu_0;
   if ( fabs(f_0) < 1.0e-6 ) {
      return M_0;
   } /* end if */

   /*
    * Make some improvements using the secant method.
    */
   M_1  = M_0 * 1.001;
   Nu_1 = NuFromM( M_1, g);
   f_1  = Nu - Nu_1;

   count = 0;
   do {
      ++count;
      slope = (f_1 - f_0) / (M_1 - M_0);
      M_2   = M_1 - f_1 / slope;
      Nu_2  = NuFromM( M_2, g );
      f_2   = Nu - Nu_2;

      M_0 = M_1; Nu_0 = Nu_1; f_0 = f_1;
      M_1 = M_2; Nu_1 = Nu_2; f_1 = f_2;
   } while ( fabs(f_2) > 1.0e-6  && count < 30 );

   if ( fabs(f_2) > 1.0e-6 ) {
      printf( "MFromNu: Warning, iteration did not converge.\n");
   } /* end if */

   return M_2;
} /* end function MFromNu */

/*------------------------------------------------------------------*/

/* @function */
double MFromNu_approximate( double Nu, double g ) {
   /**
    Purpose: Compute Mach number given Prandtl Meyer function. <BR>
    Use the polynomial approximation from <BR>
    S. M. Fraser (1975) <BR>
    Calculation of Mach number from given turning angle in 
    supersonic flow.   <BR>
    The Aeronautical Journal, February 1975. <BR>
    For large Mach numbers, use an asymptotic expansion. <BR>
    Input  : <BR>
    Nu : Parndtl-Meyer function (radians) <BR>
    g  : ratio of specific heats <BR>
    Output : <BR>
    Returns Mach number (>= 1.0) or a value of 0.0 on failure.
    */
   double Nu_d, M, Nu_max, bigG;
   double myPi = 3.1415927;

   Nu_d = Nu * 180.0 / myPi;
   if (Nu_d < 0.0) {
      M = 0.0;
   } else if ( Nu_d < 5.0 ) {
      M = 1.0 + 7.932e-2 * pow(Nu_d, 2.0/3.0) * (
         1.0 + Nu_d * (3.681e-2 + Nu_d * (-5.99e-3 + Nu_d * 5.719e-4)) );
   } else if ( Nu_d < 65.0 ) {
      M = 1.071 + Nu_d * (3.968e-2 + Nu_d * (-4.615e-4 + 
         Nu_d * (1.513e-5 + Nu_d * (-1.840e-7 + Nu_d * 1.186e-9))));
   } else {
      /* Use an asymptotic expansion for large M. */
      bigG = sqrt((g + 1.0) / (g - 1.0));
      Nu_max = myPi * 0.5 * (bigG - 1.0);
      M = (1.0 - bigG * bigG) / (Nu - Nu_max);
   } /* end if */
   
   return M;
} /* end function MFromNu_approximate */

/*------------------------------------------------------------------*/

/* @function */
double MachAngle( double M ) {
   /**
    Purpose: Compute Mach Angle from Mach number <BR>
    Input  : <BR>
    M : Mach number (M >= 1.0) <BR>
    Output : <BR>
    Returns Mach Angle or a value of 0.0 on failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   double mu;

   if ( M <= 1.0 ) {
      printf( "MachAngle: subsonic Mach number: %g\n", M );
      mu = 0.0;
   } else {
      mu = asin( 1.0 / M );
   } /* end if */

   return mu;
} /* end function MachAngle */

/*------------------------------------------------------------------*/

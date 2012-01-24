/** \file moc_wall.c
 * \ingroup imoc
 * \brief Wall functions for Method-of-Characteristics program.
 *
 * This set of functions forms part of the computational and
 * data storage kernel for the MOC program.
 * Interaction with the user is handled by the set of
 * Tcl/Tk functions.
 *
 * \author PA Jacobs
 *
 * \version 04-Jan-2000 : Initial hack
 * \version 26-Jan-2000 : Add more functions to access wall data.
 *
 */

/*-----------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "moc_bezier.h"
#include "moc_wall.h"

/*-----------------------------------------------------*/

/*
 * We shall store the wall data in an array of structures.
 */
#define NWALL     2
#define MAX_POINT 10

struct point_3D WallPoint[NWALL][MAX_POINT];
struct bezier_3_poly_data WallBezier[NWALL];
int    npoint[NWALL];
int    can_use_now[NWALL];

/*-----------------------------------------------------*/

/* @function */
int InitWall( void ) {
   /**
    Purpose: Initialize the wall data structures. <BR>
    Input:   None <BR>
    Output:  Returns 0 indicating success. <BR>
    (Available from the Tcl interpreter.) 
    */
   int result_flag, iw;

   for ( iw = 0; iw < NWALL; ++iw ) {
      npoint[iw] = 0;
      can_use_now[iw] = 0;
   } /* end for iw */

   result_flag = 0;
   for ( iw = 0; iw < NWALL; ++iw ) {
      result_flag = alloc_bezier_3_poly( &(WallBezier[iw]), MAX_POINT-1 );
      if (result_flag != 0) {
         printf( "Bezier allocation failed for iw = %d\n", iw );
         break;
      } /* end if */
   } /* end for iw */

   return result_flag;
} /* end function InitWall */

/*-----------------------------------------------------*/

/* @function */
int WallIsPresent( int iw ) {
   /**
    Returns 1 if the wall is present and can be used. <BR>
    (Available from the Tcl interpreter.) 
    */
   if ( iw >= 0 && iw < NWALL ) {
      return can_use_now[iw];
   } else {
      return 0;
   } /* end if */
} /* end function */

/*-----------------------------------------------------*/

/* @function */
int WallGetNumberOfPoints( int iw ) {
   /**
    Returns number of points defining the wall. <BR>
    (Available from the Tcl interpreter.) 
    */
   if ( iw >= 0 && iw < NWALL ) {
      return npoint[iw];
   } else {
      return 0;
   } /* end if */
} /* end function */

/*-----------------------------------------------------*/

/* @function */
double WallGetPointX( int iw, int ip ) {
   /**
    Returns the X-coordinate of the ip-th point defining the wall. <BR>
    0 <= ip < npoint[iw]. <BR>
    (Available from the Tcl interpreter.) 
    */
   double x;

   if ( iw >= 0 && iw < NWALL ) {
      if (ip < npoint[iw]) {
         x = WallPoint[iw][ip].x;
      } else {
         x = 0.0;
      } /* end if */
   } else {
      x = 0.0;
   } /* end if */

   return x;
} /* end function */

/*-----------------------------------------------------*/

/* @function */
double WallGetPointY( int iw, int ip ) {
   /**
    Returns the Y-coordinate of the ip-th point defining the wall. <BR>
    0 <= ip < npoint[iw]. <BR>
    (Available from the Tcl interpreter.) 
    */
   double y;

   if ( iw >= 0 && iw < NWALL ) {
      if (ip < npoint[iw]) {
         y = WallPoint[iw][ip].y;
      } else {
         y = 0.0;
      } /* end if */
   } else {
      y = 0.0;
   } /* end if */

   return y;
} /* end function */

/*-----------------------------------------------------*/

/* @function */
int WallDeletePoints( int iw ) {
   /** 
    Purpose: Set all of the data back to scratch. <BR>
    Input: iw is the index of the Wall to be reset. <BR>
    Output: Returns 0 indicating success or -1 on failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   int result_flag;

   if ( iw >= 0 && iw < NWALL ) {
      npoint[iw] = 0;
      can_use_now[iw] = 0;
      result_flag = 0;
   } else {
      printf( "Invalid index specified for wall: %d\n", iw );
      result_flag = -1;
   } /* end if */

   return result_flag;
} /* end function WallDeletePoints */

/*-----------------------------------------------------*/

/* @function */
int WallAddPoint( int iw, double x, double y ) {
   /**
    Purpose: Add an interpolation point to the specified wall. <BR>
    Input: <BR>
    iw   : specifies the wall <BR>
    x, y : coordinates of the point. <BR>
    Output: <BR>
    Returns the index of the point if successfully added, 
    or -1 on failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   int ip, result_flag, bez_result;

   if ( iw < 0 || iw >= NWALL ) {
      printf( "Invalid index specified for wall: %d\n", iw );
      result_flag = -1;
      return result_flag;
   } /* end if */

   if ( npoint[iw] < MAX_POINT ) {
      ip = npoint[iw];
      WallPoint[iw][ip].x = x;
      WallPoint[iw][ip].y = y;
      WallPoint[iw][ip].z = 0.0;
      ++(npoint[iw]);
      bez_result = bezier_3_spline( &(WallBezier[iw]), 
                   npoint[iw]-1, &(WallPoint[iw][0]) );
      if ( bez_result == 0 ) {
         can_use_now[iw] = 1;
      } else if ( bez_result == -2 ) {
         /* Assume that the cause is too few points (< 2) */
         can_use_now[iw] = 0;
      } else {
         printf( "Problem making spline: result flag = %d\n", bez_result );
         can_use_now[iw] = 0;
         return -1;
      } /* end if */
   } else {
      printf( "Point not added; too many already\n" );
      return -1;
   } /* end if */

   return npoint[iw];
} /* end function WallAddPoint */

/*-----------------------------------------------------*/

/* @function */
double WallPosX( int iw, double t ) {
   /**
    Purpose: Evaluate the X-value of the spline for parameter value t. <BR>
    Input  : <BR>
    iw   : specified wall <BR>
    t    : parametric distance 0 <= t <= 1 <BR>
    Output : <BR>
    Returns the X-value.  
    If anything goes wrong, a value of zero is returned. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct point_3D pos;
   double x;
   int    bez_result;

   if ( iw < 0 || iw >= NWALL ) {
      printf( "Invalid index specified for wall: %d\n", iw );
      return 0.0;
   } /* end if */

   if ( can_use_now[iw] == 1 ) {
      bez_result = eval_bezier_3_poly( &(WallBezier[iw]), t, &pos );
      if ( bez_result == 0 ) {
         x = pos.x;
      } else { 
         printf( "Bezier evaluation failed.\n");
         x = 0.0;
      } /* end if */
   } else {
      x = 0.0;
   } /* end if */

   return x;
} /* end function WallPosX */

/* @function */
double WallPosY( int iw, double t ) {
   /**
    Purpose: Evaluate the Y-value of the spline for parameter value t. <BR>
    Input  : <BR>
    iw   : specified wall <BR>
    t    : parametric distance 0 <= t <= 1 <BR>
    Output : <BR>
    Returns the Y-value.  
    If anything goes wrong, a value of zero is returned. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct point_3D pos;
   double y;
   int    bez_result;

   if ( iw < 0 || iw >= NWALL ) {
      printf( "Invalid index specified for wall: %d\n", iw );
      return 0.0;
   } /* end if */

   if ( can_use_now[iw] == 1 ) {
      bez_result = eval_bezier_3_poly( &(WallBezier[iw]), t, &pos );
      if ( bez_result == 0 ) {
         y = pos.y;
      } else { 
         printf( "Bezier evaluation failed.\n");
         y = 0.0;
      } /* end if */
   } else {
      y = 0.0;
   } /* end if */

   return y;
} /* end function WallPosY */

/*---------------------------------------------------------*/

/* @function */
double WallSlope( int iw, double t ) {
   /**
    Purpose: Evaluate the slope (in the x,y-plane) 
             of the spline for parameter value t. <BR>
    Input  : <BR>
    iw   : specified wall <BR>
    t    : parametric distance 0 <= t <= 1 <BR>
    Output :  <BR>
    Returns the slope.  
    If anything goes wrong, a value of zero is returned. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct point_3D deriv;
   double slope;
   int    bez_result;

   if ( iw < 0 || iw >= NWALL ) {
      printf( "Invalid index specified for wall: %d\n", iw );
      return 0.0;
   } /* end if */

   if ( can_use_now[iw] == 1 ) {
      bez_result = deriv_bezier_3_poly( &(WallBezier[iw]), t, &deriv );
      if ( bez_result == 0 ) {
         slope = deriv.y / deriv.x;
      } else {
         printf( "Bezier derivative failed.\n" );
         slope = 0.0;
      } /* end if */
   } else {
      slope = 0.0;
   } /* end if */

   return slope;
} /* end function WallSlope */

/*--------------------------------------------------------------*/

/* @function */
double WallFindT( int iw, 
   double xp, double yp, double cos_th, double sin_th ) {
   /**
    Purpose: Given a line (specified by a point and a slope)
             find the intersection point along the spline.  
             It is assumed that the wall is not too far from 
             being a straight line, itself. <BR>
    Input  : <BR>
    iw    : specified wall <BR>
    xp,yp : a point on the straight line. <BR>
    sin_th
    cos_th: Slope of the line specified as sin(theta)
            and cos(theta). <BR>
    Output :  <BR>
    Returns the (parametric) t-coordinatei of the
    intersection, 0.0 <= t <= 1.0.
    If anything goes wrong a value of -1.0 is returned. <BR>
    (Available from the Tcl interpreter.) 
    */
   double t, tL, tR, tL_new, tR_new;
   double xL, yL, xR, yR, xb, yb;
   double dx, dy, dxpL, dypL;
   double denominator, beta, alpha;
   double xi, yi, distance_error, distance_tolerance;
   int    count;

   if ( iw < 0 || iw >= NWALL ) {
      printf( "Invalid Wall index %d.\n", iw );
      return -1.0;
   } /* end if */

   if ( can_use_now[iw] != 1 ) {
      printf( "Cannot use Wall Bezier %d.\n", iw);
      return -1.0;
   } /* end if */

   tL_new = 0.0; xL = WallPosX(iw, tL_new); yL = WallPosY(iw, tL_new);
   tR_new = 1.0; xR = WallPosX(iw, tR_new); yR = WallPosY(iw, tR_new);
   distance_tolerance = 1.0e-6 * sqrt( (xR-xL)*(xR-xL) + (yR-yL)*(yR-yL) );
   alpha  = 0.5;
   count  = 0;
   do {
      ++count; 

      tL = tL_new; xL = WallPosX(iw, tL); yL = WallPosY(iw, tL);
      tR = tR_new; xR = WallPosX(iw, tR); yR = WallPosY(iw, tR);
      dx   = xR - xL; dy   = yR - yL;
      dxpL = xp - xL; dypL = yp - yL;

      denominator = dx * sin_th - dy * cos_th;
      if ( fabs(denominator) <= 1.0e-12 ) {
         printf( "WallFindT: Lines are parallel.\n");
         return -1.0;
      } /* end if */

      beta = (dxpL * sin_th - dypL * cos_th) / denominator;
      xi = xL + beta * dx; yi = yL + beta * dy;
      if ( beta < 0.0 || beta > 1.0) {
         printf( "WallFindT: Step %d\n", count );
         printf( "Intersection is outside current range.\n");
         printf( "Here, beta=%g, intersection=(%g,%g)\n",
            beta, xi, yi );
         printf( "Current range: tL=%g, tR=%g, L(%g,%g)->R(%g,%g)\n",
            tL, tR, xL, yL, xR, yR );
         printf( "p=(%g,%g) + cos_th=%g, sin_th=%g\n", 
            xp, yp, cos_th, sin_th );
      } /* end if */

      t = tL + beta * (tR - tL);
      xb = WallPosX(iw, t); yb = WallPosY(iw, t);

      distance_error = sqrt( (xb-xi)*(xb-xi) + (yb-yi)*(yb-yi) );

      if ( beta < 0.0 ) {
          /* Seems to be left of current range; move tL only
           and reduce the rate at which future ranges are contracted. */
         tL_new = tL - (1.0 + fabs(beta)) * (tR - tL);
         if ( tL_new < 0.0 ) tL_new = 0.0;
         tR_new = tR;
         alpha *= 0.5;
      } else if ( beta > 1.0 ) {
          /* Seems to be right of current range; move tR only
           and reduce the rate at which future ranges are contracted. */
         tL_new = tL;
         tR_new = tR + (1.0 + beta) * (tR - tL);
         if ( tR_new > 1.0 ) tR_new = 1.0;
         alpha *= 0.5;
      } else {
         /* In between the current range; contract the range. */
         tL_new = tL + alpha * beta * (tR - tL);
         tR_new = tR - alpha * (1.0 - beta) * (tR - tL);
      } /* end if */

   } while ( distance_error > distance_tolerance && 
             count < 30 && 
             fabs(tR - tL) > 1.0e-6 );

   if ( distance_error > distance_tolerance ) {
      printf( "WallFindT: Did not get close to wall point: %g\n",
         distance_error );
      printf( "           returning with t=%g\n", t);
   } /* end if */

   return t;
} /* end function WallFindT */

/*-----------------------------------------------------*/

/* @function */
int SaveWall( int iw, char *FileName ) {
   /**
    Purpose: Write the interpolation points for the specified wall. <BR>
    Input  :  <BR>
    iw       : index specifying the wall <BR>
    FileName : pointer to the filename string. <BR>
    Output :  <BR>
    Returns 0 if successful, -1 on failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   FILE *fp;
   int ip, result_flag;

   if ( iw < 0 || iw >= NWALL ) {
      printf( "Invalid index specified for wall: %d\n", iw );
      result_flag = -1;
      return result_flag;
   } /* end if */

   fp = fopen( FileName, "w" );
   if ( fp == NULL ) {
      printf( "Couldn't open file %s\n", FileName );
      return -1; 
   } /* end if */

   fprintf( fp, "%d  : number of points \n", npoint[iw] );
   for (ip = 0; ip < npoint[iw]; ++ip) {
      fprintf( fp, "%e %e  : x y\n", WallPoint[iw][ip].x, WallPoint[iw][ip].y );
   } /* end for */
 
   fclose( fp );
   return 0;
} /* end function SaveWall */

/* @function */
int LoadWall( int iw, char *FileName ) {
   /**
    Purpose: Read the interpolation points for the specified wall
             and set up the spline ready for evaluation. <BR>
    Input  :  <BR>
    iw       : index specifying the wall <BR>
    FileName : pointer to the filename string. <BR>
    Output :  <BR>
    Returns the number of points if successful, -1 on failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   FILE *fp;
   int ip, arg_count, failed_flag, bez_result, result;
   char textstring[132];
   double x, y;

   if ( iw < 0 || iw >= NWALL ) {
      printf( "Invalid index specified for wall: %d\n", iw );
      return -1;
   } /* end if */

   fp = fopen( FileName, "r" );
   if ( fp == NULL ) {
      printf( "Couldn't open file %s\n", FileName );
      return -1; 
   } /* end if */

   fgets( textstring, sizeof(textstring), fp );
   arg_count = sscanf( textstring, "%d", &(npoint[iw]) );
   if ( arg_count != 1 ) {
      printf( "Failed to read number of wall points.\n");
      return -1;
   } /* end if */

   failed_flag = 0;
   for (ip = 0; ip < npoint[iw]; ++ip) {
      fgets( textstring, sizeof(textstring), fp );
      arg_count = sscanf( textstring, "%lf %lf", &x, &y );
      if ( arg_count != 2 ) {
         printf( "Failed reading wall point %d.\n", ip );
         failed_flag = 1;
         break;
      } /* end if */
      WallPoint[iw][ip].x = x;
      WallPoint[iw][ip].y = y;
   } /* end for */
 
   fclose( fp );

   if (failed_flag == 0) {
      bez_result = bezier_3_spline( &(WallBezier[iw]), 
                   npoint[iw]-1, &(WallPoint[iw][0]) );
      if ( bez_result == 0 ) {
         can_use_now[iw] = 1;
      } else if ( bez_result == -2 ) {
         /* Assume that the cause is too few points (< 2) */
         can_use_now[iw] = 0;
         failed_flag = -1;
      } else {
         printf( "Problem making spline: result flag = %d\n", bez_result );
         can_use_now[iw] = 0;
         failed_flag = -1;
      } /* end if */
   } /* end if */

   if (failed_flag == 0) {
      result = npoint[iw];
   } else {
      result = -1;
   } /* end if */

   return result;
} /* end function LoadWall */

/*-----------------------------------------------------*/


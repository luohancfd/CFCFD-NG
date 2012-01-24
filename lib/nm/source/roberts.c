/** \file roberts.c
 * \ingroup nm
 * \brief Grid point distribution and stretching in 1D.
 *
 * \author PA Jacobs
 *
 */


/*---------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../util/source/compiler.h"
#include "roberts.h"

/*---------------------------------------------------------------*/


int distribute_points (double x1, double x2, int n, double *x,
		       double beta_end1, double beta_end2)


/*
 * Purpose...
 * -------
 * Given two end points, distribute a set of points 
 * between these end points.
 *
 * Input...
 * -----
 * x1      : end point 1
 * x2      : end point 2
 * n       : number of intervals (the number of points, including
 *           end points is n+1)
 * end1    =0 : points are not clustered to end 1
 *         =1 : points ARE clustered to end 1
 * end2    =0 : points are not clustered to end 2
 *         =1 : points ARE clustered to end 2
 * beta    : grid stretching parameter
 *         1 < beta < +inf : points are clustered 
 *         The closer to 1, the more the clustering.
 *         beta < 1        : no clustering
 *
 * Output...
 * ------
 * x[] : array of distributed points, 0...n
 *
 */
 

{  /* Begin distribute_points() */

int    i, end1 = 0, end2 = 0;
double del_eta, etabar, beta = 0.0;

if (beta_end1 > 1.0 && beta_end1 == beta_end2) {
    /* Cluster the same to each end. */ 
    distribute_points_1(x1,x2,n,x,1,1,beta_end1);
} else if (beta_end1 > 1.0 && beta_end2 > 1.0) { 
    /* Cluster with differing amounts to each end. */
    distribute_points_2( x1,x2,n,x,beta_end1,beta_end2);
} else if (beta_end1 <= 1.0 && beta_end2 <= 1.0) {
    /* Just distribute uniformly. */
    del_eta = 1.0 / n;
    for (i = 0; i <= n; ++i)  {
        etabar = del_eta * i;
        x[i] = (1.0 - etabar) * x1 + etabar * x2;
    }  
} else {
    /* Cluster to only one end. */
    if (beta_end1 <= 1.0) {end1 = 0; end2 = 1; beta = beta_end2;}
    if (beta_end2 <= 1.0) {end1 = 1; end2 = 0; beta = beta_end1;}
    distribute_points_1 (x1,x2,n,x,end1,end2,beta);
}

return 0;
}  /* End of distribute_points() */



/*---------------------------------------------------------------*/

int distribute_points_1 (double x1, double x2, int n, double x[],
			 int  end1,  int end2, double beta)

/*
 * Purpose...
 * -------
 * Given two parameter values, distribute a set of points
 * between them.
 *
 * Input...
 * -----
 * x1  : parameter value  1
 * x2  : parameter value  2
 * n   : number of intervals (the number of points, including
 *       end points is n+1)
 * end1    =0 : points are not clustered to end 1
 *         =1 : points ARE clustered to end 1
 * end2    =0 : points are not clustered to end 2
 *         =1 : points ARE clustered to end 2
 * beta    : grid stretching parameter
 *         1 < beta < +inf : points are clustered 
 *         The closer to 1, the more the clustering.
 *         beta < 1        : no clustering
 *
 * Output...
 * ------
 * x[] : array of distributed values, 0...n
 *
 */
 

{  /* Begin distribute_points_1() */

int    i, reverse, cluster;
double alpha, del_eta, eta, etabar;

/*
 * decide on stretching parameters for Robert's transform.
 */
alpha = 0.5;
reverse = 0;
cluster = 1;

if (end1 == 0 && end2 == 0) cluster = 0;
if (beta < 1.0) cluster = 0;

if (end1 == 1 && end2 == 1) alpha = 0.5;

if (end1 == 1 && end2 == 0) 
   {
   reverse = 1;
   alpha   = 0.0;
   }

if (end1 == 0 && end2 == 1)
   {
   reverse = 0;
   alpha   = 0.0;
   }

/*
 * Compute the grid points.
 */

del_eta = 1.0 / n;

for (i = 0; i <= n; ++i)
   {
   /* Compute the intermediate parameter. */
   eta = del_eta * i;

   if (reverse == 1) eta = 1.0 - eta;

   /* Cluster the points */
   if (cluster == 1)
      etabar = roberts (eta, alpha, beta);
   else 
      etabar = eta;

   if (reverse == 1) etabar = 1.0 - etabar;

   /* Compute the parameter value. */
   x[i] = (1.0 - etabar) * x1 + etabar * x2;

   }  /* end of for(i=... */

return 0;
}  /* End of distribute_points_1() */



/*---------------------------------------------------------------*/


int distribute_points_2 (double x1, double x2, int n, double x[],
			 double beta_end1, double beta_end2)

/*
 * Purpose...
 * -------
 * Given two parameter values, distribute a set of points
 * between them.
 *
 * Input...
 * -----
 * x1  : parameter value  1
 * x2  : parameter value  2
 * n   : number of intervals (the number of points, including
 *       end points is n+1)
 *      Grid Stretching parameters beta < 1 : no clustering
 * beta_end1 : 1 < beta < +inf : points are clustered 
 *              The closer to 1, the more the clustering.
 * beta_end2 : 1 < beta < +inf : points are clustered 
 *             The closer to 1, the more the clustering.
 *
 * Output...
 * ------
 * x[] : array of distributed values, 0...n
 *
 */
 

{  /* Begin distribute_points_2() */

int    i;
double del_eta, eta, etabar;
double theta,pi,frac,sign;


/* Need to check if Beta's > 1.0 */
if (beta_end1 <= 1.0 || beta_end2 <= 1.0) {
    printf("Clustering Beta value less than 1 \n");
    exit(-1);
}

/*
 * Compute the grid points.
 */

del_eta = 1.0 / n;

for (i = 0; i <= n; ++i)
   {
   /* Compute the intermediate parameter. */
   eta = del_eta * i;

   pi = 3.141592654;
   theta = (eta - 0.5) * pi;
   if (theta != 0.0) sign = theta / fabs(theta);
   else              sign = 1.0;
   frac = (pow(sin(fabs(theta)),(1.0/1.0)) * sign + 1.0)*0.5;
   etabar = (1.0 - frac) * roberts_rev (eta,beta_end1) + 
               frac      * roberts (eta,0.0,beta_end2);

   /* Compute the parameter value. */
   x[i] = (1.0 - etabar) * x1 + etabar * x2;

   }  /* end of for(i=... */

return 0;
}  /* End of distribute_points_2() */

/*---------------------------------------------------------------*/


int distribute_points_3 (double x1, double x2,int n, double x[],
			 double beta, double yc)

/*
 * Purpose...
 * -------
 * Given two parameter values, distribute a set of points
 * between them where the points are clustered around yc.
 * beta is the clustering paramater which should be greater
 * than 0. The closer beta is to 0, the less the clustering.
 * Shouldn't go much higher than say 5.
 *
 * Input...
 * -----
 * n   : number of intervals (the number of points, including
 *       end points is n+1)
 *      Grid Stretching parameters beta < 1 : no clustering
 *
 * Output...
 * ------
 * x[] : array of distributed values, 0...n
 *
 */
 

{  /* Begin distribute_points_2() */

int    i;
double del_eta, eta, etabar;
double B;


/* Need to check if Beta in range */
if (beta <= 0.0 ) {
    printf("Clustering Beta for distribute 3 out of range \n");
    exit(-1);
}

/*
 * Compute the grid points.
 */

B = (1.0/(2.0*beta)) * log((1.0 + (exp(beta) - 1.0)*yc)/
                          (1.0 + (exp(-beta) - 1.0)*yc));

del_eta = 1.0 / n;

for (i = 0; i <= n; ++i)
   {
    
   /* Compute the intermediate parameter. */
   eta = del_eta * i;

   etabar = yc * (1.0 + (sinh(beta * (eta - B))/sinh(beta * B)));

   /* Compute the parameter value. */
   x[i] = (1.0 - etabar) * x1 + etabar * x2;

   }  /* end of for(i=... */

return 0;
}  /* End of distribute_points_3() */




/*---------------------------------------------------------------*/

double roberts (double eta, double alpha, double beta)


/*
 * Input ...
 * -----
 * eta   : unstretched coordinate, 0 <= eta <= 1.0
 * beta  : stretching factor (more stretching as beta --> 1.0)
 * alpha : location of stretching
 *         alpha = 0.5 : clustering of nodes at both extremes of eta
 *         alpha = 0.0 : nodes will be clustered near eta=1.0
 *
 * Output ...
 * ------
 * roberts() : stretched coordinate, 0 <= roberts <= 1.0
 */

{
double lambda, etabar;

lambda = (beta + 1.0) / (beta - 1.0);
lambda = pow ( lambda, ((eta - alpha)/(1.0 - alpha)) );
etabar = (beta + 2.0 * alpha) * lambda - beta + 2.0 * alpha;
etabar = etabar / ((2.0 * alpha + 1.0) * (1.0 + lambda));

return etabar;
}

/*---------------------------------------------------------------*/

double roberts_rev (double eta, double beta)


/*
 * Input ...
 * -----
 * eta   : unstretched coordinate, 0 <= eta <= 1.0
 * beta  : stretching factor (more stretching as beta --> 1.0)
 *
 * Output ...
 * ------
 * roberts() : stretched coordinate, 0 <= roberts <= 1.0
 */

{
double lambda, etabar;

lambda = (beta + 1.0) / (beta - 1.0);
lambda = pow ( lambda, (1.0 - eta ) );
etabar = (2.0 * beta) / ( 1.0 + lambda );
etabar = (1.0 - beta + etabar);

return etabar;
}


/*=================================================================*/

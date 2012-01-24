/** \file moc_unitproc.c
 * \ingroup imoc
 * \brief Basic Unit Processes.
 *
 * Basic unit processes implemented in C.
 * Each function computes the data for a new node based on information 
 * from other nodes.
 * The functions are computationally intensive and access the internals of
 * the node data structure directly.  
 * Thus, we believe that the implementation is faster and tidier than a pure
 * Tcl implementation -- at least we hope so.
 */

#include <stdio.h>
#include <math.h>
#include "moc_kernel.h"
#include "moc_gasdynamic.h"
#include "moc_wall.h"
#include "moc_unitproc.h"

/*------------------------------------------------------------------*/

/**
 Most of these functions involve some iteration.
 The parameters for convergence check are: <BR>
 max_iteration = 15 <BR>
 position_tolerance = 1.0e-5
 */
int    max_iteration = 15;
double position_tolerance = 1.0e-5;

/*------------------------------------------------------------------*/

/* @function */
int InteriorNode( int node1, int node2, int node4 ) {
   /**
    Purpose: Calculate an interior point from two initial points. <BR>
    Input  :  <BR>
    node1 : index of initial point along C- characteristic <BR>
    node2 : index of initial point along C+ characteristic <BR>
    node4 : index of solution point (may have a value of -1) <BR>
    If -1 is specified as the index for node4, a new node will be
    created for the solution point. <BR>
    Output :  <BR>
    Returns the index of the solution point or a value of -1
    if there has been a failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct NodeData *n1;
   struct NodeData *n2;
   struct NodeData *n4;
   int    node4_created;
   double x1, y1, pm1, mu1, th1, M1;
   double x2, y2, pm2, mu2, th2, M2;
   double x4, y4, pm4, mu4, th4, M4;
   double xStream, yStream, dot_product;
   double lengthCplus, lengthCminus;
   int    directionCplus, directionCminus;
   double integralCplus, integralCminus;
   double axiterm1, axiterm2, axiterm4;
   double cosCplus, sinCplus, cosCminus, sinCminus;
   double lambdaCminus, numerator, denominator;
   double x4_old, y4_old, dx, dy, change_in_position;
   int    iteration_count;

   /*
    * First, get pointers to the node data.
    */
   if ( node1 == node2 ) {
      printf( "InteriorNode: error, same index given for node1 and node2.\n" );
      return -1;
   } /* end if */

   n1 = GetNodePtr( node1 );
   if ( n1 == NULL ) {
      printf( "InteriorNode: error: node1 with id=%d doesn't exist.\n", node1 );
      return -1;
   } /* end if */

   n2 = GetNodePtr( node2 );
   if ( n2 == NULL ) {
      printf( "InteriorNode: error: node2 with id=%d doesn't exist.\n", node2 );
      return -1;
   } /* end if */

   node4_created = 0;
   if ( node4 == -1 ) {
      node4 = CreateNode( -1 );
      node4_created = 1;
   } /* end if */

   n4 = GetNodePtr( node4 );
   if ( n4 == NULL ) {
      node4 = CreateNode( node4 );
      if ( node4 == -1 ) {
         printf( "InteriorNode: error: couldn't create node %d.\n", node4 );
         return -1;
      } else {
         n4    = GetNodePtr( node4 );
      } /* end if */
   } /* end if */

   /*
    * Make some local copies of the node data.
    */
   x1 = n1->X; y1 = n1->Y; pm1 = n1->Nu; th1 = n1->Theta; M1 = n1->Mach;
   x2 = n2->X; y2 = n2->Y; pm2 = n2->Nu; th2 = n2->Theta; M2 = n2->Mach;
   mu1 = MachAngle( M1 );
   mu2 = MachAngle( M2 );
   /* Use point 2 to get the streamline direction cosines. */
   xStream = cos(th2);
   yStream = sin(th2);

   /*
    * Guess at some of the solution point properties.
    * The position will be way off but it is used as part of
    * the convergence check a little further on.
    */
   x4  = 0.5 * (x1  + x2);
   y4  = 0.5 * (y1  + y2);
   th4 = 0.5 * (th1 + th2);
   mu4 = 0.5 * (mu1 + mu2);

   /*
    * Compute the solution point position and flow properties.
    */
   iteration_count = 0;
   do {
      ++iteration_count;
      x4_old = x4; y4_old = y4;

      /*
       * Locate solution point by assuming straight-line segments.
       */
      sinCminus = 0.5 * ( sin(th1 - mu1) + sin(th4 - mu4) );
      cosCminus = 0.5 * ( cos(th1 - mu1) + cos(th4 - mu4) );
      sinCplus  = 0.5 * ( sin(th2 + mu2) + sin(th4 + mu4) );
      cosCplus  = 0.5 * ( cos(th2 + mu2) + cos(th4 + mu4) );

      numerator = (x2 - x1) * sinCplus - (y2 - y1) * cosCplus;
      denominator = cosCminus * sinCplus - sinCminus * cosCplus;
      if ( fabs(denominator) <= 1.0e-12 ) {
         printf( "InteriorNode: characteristics are parallel.\n" );
         if ( node4_created ) DeleteNode( node4 );
         return -1;
      } /* end if */
      lambdaCminus = numerator / denominator;
      x4 = x1 + lambdaCminus * cosCminus;
      y4 = y1 + lambdaCminus * sinCminus;
      dx = x4 - x4_old; dy = y4 - y4_old;
      change_in_position = sqrt( dx * dx + dy * dy );

      /* Lengths of the characteristic segments */
      dx = x4 - x1; dy = y4 - y1;
      lengthCminus = sqrt( dx * dx + dy * dy );
      dot_product = dx * xStream + dy * yStream;
      if ( dot_product < 0.0 ) {
         directionCminus = -1;
      } else {
         directionCminus = +1;
      } /* end if */

      dx = x4 - x2; dy = y4 - y2;
      lengthCplus  = sqrt( dx * dx + dy * dy );
      dot_product = dx * xStream + dy * yStream;
      if ( dot_product < 0.0 ) {
         directionCplus = -1;
      } else {
         directionCplus = +1;
      } /* end if */

      /*
       * Update flow properties at solution point
       * First, assume 2D planar geometry then add
       * axisymmetric contributions if flag is set.
       */
      pm4 = 0.5 * (pm1 + pm2) + 0.5 * (th1 - th2);
      th4 = 0.5 * (pm1 - pm2) + 0.5 * (th1 + th2);
      if ( GetAxiFlag() == 1 ) {
         if ( y1 < 1.0e-6 && y2 < 1.0e-6 ) {
            printf( "InteriorNode: both initial nodes are too close to axis.\n" );
            if (node4_created ) DeleteNode( node4 );
            return -1;
         } /* end if */

         /* Axisymmetric components. */
         if ( y1 < 1.0e-6 ) {
            axiterm1 = sin(mu2) * sin(th2) / y2;
         } else {
            axiterm1 = sin(mu1) * sin(th1) / y1;
         } /* end if */

         if ( y2 < 1.0e-6 ) {
            axiterm2 = sin(mu1) * sin(th1) / y1;
         } else {
            axiterm2 = sin(mu2) * sin(th2) / y2;
         } /* end if */

         axiterm4 = sin(mu4) * sin(th4) / y4;
         integralCminus = 0.5 * directionCminus * lengthCminus * (axiterm1 + axiterm4);
         integralCplus  = 0.5 * directionCplus * lengthCplus  * (axiterm2 + axiterm4);

         /* Include axisymmetric components. */
         pm4 += 0.5 * (integralCminus + integralCplus);
         th4 += 0.5 * (integralCminus - integralCplus);
      } /* end if */
 
   } while ( iteration_count < max_iteration && 
           change_in_position > position_tolerance );

   /*
    * Save the solution-point properties and connect the 
    * node into the characteristic mesh.
    */
   M4 = MFromNu( pm4, GetGamma() );
   n4->X = x4; n4->Y = y4; n4->Nu = pm4; n4->Theta = th4; n4->Mach = M4;
   if ( x4 > x2 ) {
      n4->CPlusUp    = node2;
      n2->CPlusDown  = node4;
   } else {
      n4->CPlusDown  = node2;
      n2->CPlusUp    = node4;
   } /* end if */
   if ( x4 > x1 ) {
      n4->CMinusUp   = node1;
      n1->CMinusDown = node4;
   } else {
      n4->CMinusDown = node1;
      n1->CMinusUp   = node4;
   } /* end if */

   /*
    * Assuming a successful calculation, 
    * return the index of the solution node.
    */
   return node4;
} /* end function InteriorNode */

/*------------------------------------------------------------------*/

/* @function */
int InsertNode( int node1, int node2, int node4, double alpha ) {
   /**
    Purpose: Insert a node (node4) in between two initial 
        nodes (node1 and node2). <BR>
        If node1 and node2 are adjacent nodes along a characteristic line,
        node4 will be connected in between. <BR>
    Input  :  <BR>
    node1 : index of initial point 1 <BR>
    node2 : index of initial point 2 <BR>
    node4 : index of solution point (may have a value of -1) <BR>
    If -1 is specified as the index for node4, a new node will be
    created for the solution point. <BR>
    alpha : fraction that node4 is like node2;
        n4.value = alpha n2.value + (1-alpha) n1.value <BR>
    Output :  <BR>
    Returns the index of the solution point or a value of -1
    if there has been a failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct NodeData *n1;
   struct NodeData *n2;
   struct NodeData *n4;
   int    node4_created;

   /*
    * First, get pointers to the node data.
    */
   if ( node1 == node2 ) {
      printf( "InteriorNode: error, same index given for node1 and node2.\n" );
      return -1;
   } /* end if */

   n1 = GetNodePtr( node1 );
   if ( n1 == NULL ) {
      printf( "InteriorNode: error: node1 with id=%d doesn't exist.\n", node1 );
      return -1;
   } /* end if */

   n2 = GetNodePtr( node2 );
   if ( n2 == NULL ) {
      printf( "InteriorNode: error: node2 with id=%d doesn't exist.\n", node2 );
      return -1;
   } /* end if */

   node4_created = 0;
   if ( node4 == -1 ) {
      node4 = CreateNode( -1 );
      node4_created = 1;
   } /* end if */

   n4 = GetNodePtr( node4 );
   if ( n4 == NULL ) {
      node4 = CreateNode( node4 );
      if ( node4 == -1 ) {
         printf( "InteriorNode: error: couldn't create node %d.\n", node4 );
         return -1;
      } else {
         n4    = GetNodePtr( node4 );
      } /* end if */
   } /* end if */

   /*
    * Enforce a 0.0..1.0 range for alpha
    */
   if (alpha > 1.0) alpha = 1.0;
   if (alpha < 0.0) alpha = 0.0;

   /*
    * Linearly Interpolate all node properties.
    */
   n4->X     = (1.0 - alpha) * n1->X     + alpha * n2->X;
   n4->Y     = (1.0 - alpha) * n1->Y     + alpha * n2->Y;
   n4->Nu    = (1.0 - alpha) * n1->Nu    + alpha * n2->Nu;
   n4->Theta = (1.0 - alpha) * n1->Theta + alpha * n2->Theta;
   n4->P0    = (1.0 - alpha) * n1->P0    + alpha * n2->P0;
   n4->T0    = (1.0 - alpha) * n1->T0    + alpha * n2->T0;
   n4->Mach  = MFromNu( n4->Nu, GetGamma() );

   /*
    * Connect into the mesh only if nodes 1 and 2 are adjacent.
    */
   if ( n1->CPlusDown == n2->CPlusUp && n1->CPlusDown != -1 ) {
      n4->CPlusUp    = node1;
      n1->CPlusDown  = node4;
      n2->CPlusUp    = node4;
      n4->CPlusDown  = node2;
   } else if ( n1->CPlusUp == n2->CPlusDown && n1->CPlusUp != -1 ) {
      n4->CPlusUp    = node2;
      n2->CPlusDown  = node4;
      n1->CPlusUp    = node4;
      n4->CPlusDown  = node1;
   } else if ( n1->CMinusDown == n2->CMinusUp && n1->CMinusDown != -1 ) {
      n4->CMinusUp   = node1;
      n1->CMinusDown = node4;
      n2->CMinusUp   = node4;
      n4->CMinusDown = node2;
   } else if ( n1->CMinusUp == n2->CMinusDown && n1->CMinusUp != -1 ) {
      n4->CMinusUp   = node2;
      n2->CMinusDown = node4;
      n1->CMinusUp   = node4;
      n4->CMinusDown = node1;
   } else if ( n1->CZeroDown == n2->CZeroUp && n1->CZeroDown != -1 ) {
      n4->CZeroUp    = node1;
      n1->CZeroDown  = node4;
      n2->CZeroUp    = node4;
      n4->CZeroDown  = node2;
   } else if ( n1->CZeroUp == n2->CZeroDown && n1->CZeroUp != -1 ) {
      n4->CZeroUp    = node2;
      n2->CZeroDown  = node4;
      n1->CZeroUp    = node4;
      n4->CZeroDown  = node1;
   } /* end if */

   /*
    * Assuming a successful calculation, 
    * return the index of the solution node.
    */
   return node4;
} /* end function InsertNode */

/*------------------------------------------------------------------*/

/* @function */
int CMinusWallNode( int iw, int node1, int node4 ) {
   /**
    Purpose: Calculate a wall point from one initial (C-) point. <BR>
    Input  : <BR>
    iw    : Index of selected wall. <BR>
    node1 : index of initial point along C- characteristic <BR>
    node4 : index of solution point (may have a value of -1) <BR>
    If -1 is specified as the index for node4, a new node will be
    created for the solution point. <BR>
    Output : <BR>
    Returns the index of the solution point or a value of -1
    if there has been a failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct NodeData *n1;
   struct NodeData *n4;
   int    node4_created;
   double x1, y1, pm1, mu1, th1, M1;
   double x4, y4, pm4, mu4, th4, M4;
   double lengthCminus, t;
   double xStream, yStream, dot_product;
   int    directionCminus;
   double integralCminus;
   double axiterm1, axiterm4;
   double cosCminus, sinCminus;
   double x4_old, y4_old, dx, dy, change_in_position;
   int    iteration_count;

   /*
    * First, get pointers to the node data.
    */
   if ( WallIsPresent(iw) != 1 ) {
      printf( "CMinusWallNode: Wall %d is not present.\n", iw );
      return -1;
   } /* end if */

   n1 = GetNodePtr( node1 );
   if ( n1 == NULL ) {
      printf( "CMinusWallNode: error: node1 with id=%d doesn't exist.\n", node1 );
      return -1;
   } /* end if */

   node4_created = 0;
   if ( node4 == -1 ) {
      node4 = CreateNode( -1 );
      node4_created = 1;
   } /* end if */

   n4 = GetNodePtr( node4 );
   if ( n4 == NULL ) {
      node4 = CreateNode( node4 );
      if ( node4 == -1 ) {
         printf( "CMinusWallNode: error: couldn't create node %d.\n", node4 );
         return -1;
      } else {
         n4    = GetNodePtr( node4 );
      } /* end if */
   } /* end if */

   /*
    * Make some local copies of the node data.
    */
   x1 = n1->X; y1 = n1->Y; pm1 = n1->Nu; th1 = n1->Theta; M1 = n1->Mach;
   mu1 = MachAngle( M1 );
   /* Use point 1 to get the streamline direction cosines. */
   xStream = cos(th1);
   yStream = sin(th1);

   /*
    * Guess at some of the solution point properties.
    * The position will be way off but it is used as part of
    * the convergence check a little further on.
    */
   x4 = x1; y4 = y1; th4 = th1; mu4 = mu1;
 
   /*
    * Compute the solution point position and flow properties.
    */
   iteration_count = 0;
   do {
      ++iteration_count;
      x4_old = x4; y4_old = y4;

      /*
       * Locate solution point by assuming straight-line segments.
       */
      sinCminus = 0.5 * ( sin(th1 - mu1) + sin(th4 - mu4) );
      cosCminus = 0.5 * ( cos(th1 - mu1) + cos(th4 - mu4) );

      t = WallFindT( iw, x1, y1, cosCminus, sinCminus );
      x4 = WallPosX( iw, t );
      y4 = WallPosY( iw, t );
      if ( GetDebugLevel() >= 1 ) {
          printf( "CMinusWallNode: Find t: iw=%d x1=%g y1=%g cos=%g sin=%g\n",
             iw, x1, y1, cosCminus, sinCminus );
          printf( "CMinusWallNode: Iter=%d, t=%g, x4=%g, y4=%g\n", 
             iteration_count, t, x4, y4 );
      } /* end if */
      dx = x4 - x4_old; dy = y4 - y4_old;
      change_in_position = sqrt( dx * dx + dy * dy );

      /* Lengths of the characteristic segments */
      dx = x4 - x1; dy = y4 - y1;
      lengthCminus = sqrt( dx * dx + dy * dy );
      dot_product = dx * xStream + dy * yStream;
      if ( dot_product < 0.0 ) {
         directionCminus = -1;
      } else {
         directionCminus = +1;
      } /* end if */

      /*
       * Update flow properties at solution point
       * First, assume 2D planar geometry then add
       * axisymmetric contributions if flag is set.
       */
      th4 = atan( WallSlope( iw, t ) );
      pm4 = pm1 - (th4 - th1);
      if ( GetAxiFlag() == 1 ) {
         if ( y1 < 1.0e-6 ) {
            printf( "CMinusWallNode: point 1 is too close to axis: y=%g", y1 );
            if ( node4_created ) DeleteNode( node4 );
            return -1;
         } /* end if */

         /* Axisymmetric components. */
         axiterm1 = sin(mu1) * sin(th1) / y1;
         if ( y4 < 1.0e-6 ) {
            axiterm4 = axiterm1;
         } else {
            axiterm4 = sin(mu4) * sin(th4) / y4;
         } /* end if */
         integralCminus = 0.5 * directionCminus * lengthCminus * (axiterm1 + axiterm4);

         /* Include axisymmetric components. */
         pm4 += integralCminus;
      } /* end if */
 
   } while ( iteration_count < max_iteration && 
           change_in_position > position_tolerance );

   /*
    * Save the solution-point properties and connect the 
    * node into the characteristic mesh.
    */
   M4 = MFromNu( pm4, GetGamma() );
   n4->X = x4; n4->Y = y4; n4->Nu = pm4; n4->Theta = th4; n4->Mach = M4;
   if ( x4 > x1 ) {
      n4->CMinusUp   = node1;
      n1->CMinusDown = node4;
   } else {
      n4->CMinusDown = node1;
      n1->CMinusUp   = node4;
   } /* end if */

   /*
    * Assuming a successful calculation, 
    * return the index of the solution node.
    */
   return node4;
} /* end function CMinusWallNode */

/*------------------------------------------------------------------*/

/* @function */
int CPlusWallNode( int iw, int node2, int node4 ) {
   /**
    Purpose: Calculate a wall point from one upstream (C+) point. <BR>
    Input  : <BR>
    iw    : index of the wall <BR>
    node2 : index of initial point along C+ characteristic <BR>
    node4 : index of solution point (may have a value of -1) <BR>
    If -1 is specified as the index for node4, a new node will be
    created for the solution point. <BR>
    Output : <BR>
    Returns the index of the solution point or a value of -1
    if there has been a failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct NodeData *n2;
   struct NodeData *n4;
   int    node4_created;
   double x2, y2, pm2, mu2, th2, M2;
   double x4, y4, pm4, mu4, th4, M4;
   double lengthCplus, t;
   double xStream, yStream, dot_product;
   int    directionCplus;
   double integralCplus;
   double axiterm2, axiterm4;
   double cosCplus, sinCplus;
   double x4_old, y4_old, dx, dy, change_in_position;
   int    iteration_count;

   /*
    * First, get pointers to the node data.
    */
   if ( WallIsPresent(iw) != 1 ) {
      printf( "CPlusWallNode: Wall %d is not present.\n", iw );
      return -1;
   } /* end if */

   n2 = GetNodePtr( node2 );
   if ( n2 == NULL ) {
      printf( "CPlusWallNode: error: node2 with id=%d doesn't exist.\n", node2 );
      return -1;
   } /* end if */

   node4_created = 0;
   if ( node4 == -1 ) {
      node4 = CreateNode( -1 );
      node4_created = 1;
   } /* end if */

   n4 = GetNodePtr( node4 );
   if ( n4 == NULL ) {
      node4 = CreateNode( node4 );
      if ( node4 == -1 ) {
         printf( "CPlusWallNode: error: couldn't create node %d.\n", node4 );
         return -1;
      } else {
         n4    = GetNodePtr( node4 );
      } /* end if */
   } /* end if */

   /*
    * Make some local copies of the node data.
    */
   x2 = n2->X; y2 = n2->Y; pm2 = n2->Nu; th2 = n2->Theta; M2 = n2->Mach;
   mu2 = MachAngle( M2 );
   /* Use point 2 to get the streamline direction cosines. */
   xStream = cos(th2);
   yStream = sin(th2);

   /*
    * Guess at some of the solution point properties.
    * The position will be way off but it is used as part of
    * the convergence check a little further on.
    */
   x4  = x2; y4  = y2; th4 = th2; mu4 = mu2;

   /*
    * Compute the solution point position and flow properties.
    */
   iteration_count = 0;
   do {
      ++iteration_count;
      x4_old = x4; y4_old = y4;

      /*
       * Locate solution point by assuming straight-line segments.
       */
      sinCplus  = 0.5 * ( sin(th2 + mu2) + sin(th4 + mu4) );
      cosCplus  = 0.5 * ( cos(th2 + mu2) + cos(th4 + mu4) );

      t = WallFindT( iw, x2, y2, cosCplus, sinCplus );
      x4 = WallPosX( iw, t );
      y4 = WallPosY( iw, t );
      dx = x4 - x4_old; dy = y4 - y4_old;
      change_in_position = sqrt( dx * dx + dy * dy );

      /* Lengths of the characteristic segments */
      dx = x4 - x2; dy = y4 - y2;
      lengthCplus  = sqrt( dx * dx + dy * dy );
      dot_product = dx * xStream + dy * yStream;
      if ( dot_product < 0.0 ) {
         directionCplus = -1;
      } else {
         directionCplus = +1;
      } /* end if */

      /*
       * Update flow properties at solution point
       * First, assume 2D planar geometry then add
       * axisymmetric contributions if flag is set.
       */
      th4 = atan( WallSlope( iw, t ) );
      pm4 = pm2 + (th4 - th2);
      if ( GetAxiFlag() == 1 ) {
         if ( y4 < 1.0e-6 ) {
            printf( "CPlusWallNode: point 4 is too close to axis: y=%g", y4 );
            if ( node4_created ) DeleteNode( node4 );
            return -1;
         } /* end if */

         /* Axisymmetric components. */
         axiterm4 = sin(mu4) * sin(th4) / y4;
         if ( y2 < 1.0e-6 ) {
            axiterm2 = axiterm4;
         } else {
            axiterm2 = sin(mu2) * sin(th2) / y2;
         } /* end if */
         integralCplus  = 0.5 * directionCplus * lengthCplus * (axiterm2 + axiterm4);

         /* Include axisymmetric components. */
         pm4 += integralCplus;
      } /* end if */
 
   } while ( iteration_count < max_iteration && 
           change_in_position > position_tolerance );

   /*
    * Save the solution-point properties and connect the 
    * node into the characteristic mesh.
    */
   M4 = MFromNu( pm4, GetGamma() );
   n4->X = x4; n4->Y = y4; n4->Nu = pm4; n4->Theta = th4; n4->Mach = M4;
   if ( x4 > x2 ) {
      n4->CPlusUp    = node2;
      n2->CPlusDown  = node4;
   } else {
      n4->CPlusDown  = node2;
      n2->CPlusUp    = node4;
   } /* end if */

   /*
    * Assuming a successful calculation, 
    * return the index of the solution node.
    */
   return node4;
} /* end function CPlusWallNode */

/*------------------------------------------------------------------*/

/* @function */
int CPlusFreeBndyNode( int node0, int node2, int node4 ) {
   /**
    Purpose: Calculate a free-boundary point from one point (node0) 
             already on the boundary and one point (node2)
             on a C+ characteristic. <BR>
    Input  : <BR>
    node0 : index of initial point along C0 streamline <BR>
    node2 : index of initial point along C+ characteristic <BR>
    node4 : index of solution point (may have a value of -1) <BR>
    If -1 is specified as the index for node4, a new node will be
    created for the solution point. <BR>
    Output : <BR>
    Returns the index of the solution point or a value of -1
    if there has been a failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct NodeData *n0;
   struct NodeData *n2;
   struct NodeData *n4;
   int    node4_created;
   double x0, y0, pm0, mu0, th0, M0;
   double x2, y2, pm2, mu2, th2, M2;
   double x4, y4, pm4, mu4, th4, M4;
   double lengthCplus, lengthCzero;
   double xStream, yStream, dot_product;
   int    directionCplus;
   double integralCplus;
   double axiterm2, axiterm4;
   double cosCplus, sinCplus, cosCzero, sinCzero;
   double lambdaCzero, numerator, denominator;
   double x4_old, y4_old, dx, dy, change_in_position;
   int    iteration_count;

   /*
    * First, get pointers to the node data.
    */
   if ( node0 == node2 ) {
      printf( "CPlusFreeBndyNode: error, same index given for node0 and node2.\n" );
      return -1;
   } /* end if */

   n0 = GetNodePtr( node0 );
   if ( n0 == NULL ) {
      printf( "CPlusFreeBndyNode: error: node0 with id=%d doesn't exist.\n", node0 );
      return -1;
   } /* end if */

   n2 = GetNodePtr( node2 );
   if ( n2 == NULL ) {
      printf( "CPlusFreeBndyNode: error: node2 with id=%d doesn't exist.\n", node2 );
      return -1;
   } /* end if */

   node4_created = 0;
   if ( node4 == -1 ) {
      node4 = CreateNode( -1 );
      node4_created = 1;
   } /* end if */

   n4 = GetNodePtr( node4 );
   if ( n4 == NULL ) {
      node4 = CreateNode( node4 );
      if ( node4 == -1 ) {
         printf( "CPlusFreeBndyNode: error: couldn't create node %d.\n", node4 );
         return -1;
      } else {
         n4    = GetNodePtr( node4 );
      } /* end if */
   } /* end if */

   /*
    * Make some local copies of the node data.
    */
   x0 = n0->X; y0 = n0->Y; pm0 = n0->Nu; th0 = n0->Theta; M0 = n0->Mach;
   x2 = n2->X; y2 = n2->Y; pm2 = n2->Nu; th2 = n2->Theta; M2 = n2->Mach;
   mu0 = MachAngle( M0 );
   mu2 = MachAngle( M2 );
   /* Use point 0 to get the streamline direction cosines. */
   xStream = cos(th0);
   yStream = sin(th0);

   /*
    * Guess at some of the solution point properties.
    * The position will be way off but it is used as part of
    * the convergence check a little further on.
    */
   x4  = 0.5 * (x0 + x2);
   y4  = 0.5 * (y0 + y2);
   th4 = th0;
   mu4 = mu0;

   /*
    * Compute the solution point position and flow properties.
    */
   iteration_count = 0;
   do {
      ++iteration_count;
      x4_old = x4; y4_old = y4;

      /*
       * Locate solution point by assuming straight-line segments.
       */
      sinCzero = 0.5 * ( sin(th0) + sin(th4) );
      cosCzero = 0.5 * ( cos(th0) + cos(th4) );
      sinCplus = 0.5 * ( sin(th2 + mu2) + sin(th4 + mu4) );
      cosCplus = 0.5 * ( cos(th2 + mu2) + cos(th4 + mu4) );

      numerator = (x2 - x0) * sinCplus - (y2 - y0) * cosCplus;
      denominator = cosCzero  * sinCplus - sinCzero  * cosCplus;
      if ( fabs(denominator) <= 1.0e-12 ) {
         printf( "CPlusFreeBndyNode: characteristics are parallel.\n" );
         if ( node4_created ) DeleteNode( node4 );
         return -1;
      } /* end if */
      lambdaCzero  = numerator / denominator;
      x4 = x0 + lambdaCzero * cosCzero;
      y4 = y0 + lambdaCzero * sinCzero;
      dx = x4 - x4_old; dy = y4 - y4_old;
      change_in_position = sqrt( dx * dx + dy * dy );

      /* Lengths of the characteristic segments */
      dx = x4 - x0; dy = y4 - y0;
      lengthCzero = sqrt( dx * dx + dy * dy );

      dx = x4 - x2; dy = y4 - y2;
      lengthCplus = sqrt( dx * dx + dy * dy );
      dot_product = dx * xStream + dy * yStream;
      if ( dot_product < 0.0 ) {
         directionCplus = -1;
      } else {
         directionCplus = +1;
      } /* end if */

      /*
       * Update flow properties at solution point
       * First, assume 2D planar geometry then add
       * axisymmetric contributions if flag is set.
       */
      pm4 = pm0;
      th4 = th2 + (pm4 - pm2);
      if ( GetAxiFlag() == 1 ) {
         if ( y0 < 1.0e-6 && y2 < 1.0e-6 ) {
            printf( "CPlusFreeBndyNode: both initial nodes are too close to axis.\n" );
            if ( node4_created ) DeleteNode( node4 );
            return -1;
         } /* end if */

         if ( y2 < 1.0e-6 ) {
            axiterm2 = sin(mu0) * sin(th0) / y0;
         } else {
            axiterm2 = sin(mu2) * sin(th2) / y2;
         } /* end if */

         axiterm4 = sin(mu4) * sin(th4) / y4;
         integralCplus  = 0.5 * directionCplus * lengthCplus  * (axiterm2 + axiterm4);

         /* Include axisymmetric components. */
         th4 -= integralCplus;
      } /* end if */
 
   } while ( iteration_count < max_iteration && 
           change_in_position > position_tolerance );

   /*
    * Save the solution-point properties and connect the 
    * node into the characteristic mesh.
    */
   M4 = MFromNu( pm4, GetGamma() );
   n4->X = x4; n4->Y = y4; n4->Nu = pm4; n4->Theta = th4; n4->Mach = M4;
   if ( x4 > x2 ) {
      n4->CPlusUp    = node2;
      n2->CPlusDown  = node4;
   } else {
      n4->CPlusDown  = node2;
      n2->CPlusUp    = node4;
   } /* end if */
   if ( x4 > x0 ) {
      n4->CZeroUp   = node0;
      n0->CZeroDown = node4;
   } else {
      n4->CZeroDown = node0;
      n0->CZeroUp   = node4;
   } /* end if */

   /*
    * Assuming a successful calculation, 
    * return the index of the solution node.
    */
   return node4;
} /* end function CPlusFreeBndyNode */


/*------------------------------------------------------------------*/

/* @function */
int CMinusFreeBndyNode( int node0, int node1, int node4 ) {
   /**
    Purpose: Calculate a free-boundary point from one point (node0) 
             already on the boundary and one point (node1)
             on a C- characteristic. <BR>
    Input  :  <BR>
    node0 : index of initial point along C0 streamline <BR>
    node1 : index of initial point along C- characteristic <BR>
    node4 : index of solution point (may have a value of -1) <BR>
    If -1 is specified as the index for node4, a new node will be
    created for the solution point. <BR>
    Output :  <BR>
    Returns the index of the solution point or a value of -1
    if there has been a failure. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct NodeData *n1;
   struct NodeData *n0;
   struct NodeData *n4;
   int    node4_created;
   double x1, y1, pm1, mu1, th1, M1;
   double x0, y0, pm0, mu0, th0, M0;
   double x4, y4, pm4, mu4, th4, M4;
   double lengthCzero, lengthCminus;
   double xStream, yStream, dot_product;
   int    directionCminus;
   double integralCminus;
   double axiterm1, axiterm4;
   double cosCzero, sinCzero, cosCminus, sinCminus;
   double lambdaCminus, numerator, denominator;
   double x4_old, y4_old, dx, dy, change_in_position;
   int    iteration_count;

   /*
    * First, get pointers to the node data.
    */
   if ( node1 == node0 ) {
      printf( "CMinusFreeBndyNode: error, same index given for node1 and node0.\n" );
      return -1;
   } /* end if */

   n1 = GetNodePtr( node1 );
   if ( n1 == NULL ) {
      printf( "CMinusFreeBndyNode: error: node1 with id=%d doesn't exist.\n", node1 );
      return -1;
   } /* end if */

   n0 = GetNodePtr( node0 );
   if ( n0 == NULL ) {
      printf( "CMinusFreeBndyNode: error: node0 with id=%d doesn't exist.\n", node0 );
      return -1;
   } /* end if */

   node4_created = 0;
   if ( node4 == -1 ) {
      node4 = CreateNode( -1 );
      node4_created = 1;
   } /* end if */

   n4 = GetNodePtr( node4 );
   if ( n4 == NULL ) {
      node4 = CreateNode( node4 );
      if ( node4 == -1 ) {
         printf( "CMinusFreeBndyNode: error: couldn't create node %d.\n", node4 );
         return -1;
      } else {
         n4    = GetNodePtr( node4 );
      } /* end if */
   } /* end if */

   /*
    * Make some local copies of the node data.
    */
   x1 = n1->X; y1 = n1->Y; pm1 = n1->Nu; th1 = n1->Theta; M1 = n1->Mach;
   x0 = n0->X; y0 = n0->Y; pm0 = n0->Nu; th0 = n0->Theta; M0 = n0->Mach;
   mu1 = MachAngle( M1 );
   mu0 = MachAngle( M0 );
   /* Use point 0 to get the streamline direction cosines. */
   xStream = cos(th0);
   yStream = sin(th0);

   /*
    * Guess at some of the solution point properties.
    * The position will be way off but it is used as part of
    * the convergence check a little further on.
    */
   x4  = 0.5 * (x1 + x0);
   y4  = 0.5 * (y1 + y0);
   th4 = th0;
   mu4 = mu0;

   /*
    * Compute the solution point position and flow properties.
    */
   iteration_count = 0;
   do {
      ++iteration_count;
      x4_old = x4; y4_old = y4;

      /*
       * Locate solution point by assuming straight-line segments.
       */
      sinCminus = 0.5 * ( sin(th1 - mu1) + sin(th4 - mu4) );
      cosCminus = 0.5 * ( cos(th1 - mu1) + cos(th4 - mu4) );
      sinCzero  = 0.5 * ( sin(th0) + sin(th4) );
      cosCzero  = 0.5 * ( cos(th0) + cos(th4) );

      numerator = (x0 - x1) * sinCzero - (y0 - y1) * cosCzero;
      denominator = cosCminus * sinCzero - sinCminus * cosCzero;
      if ( fabs(denominator) <= 1.0e-12 ) {
         printf( "CMinusFreeBndyNode: characteristics are parallel.\n" );
         if ( node4_created ) DeleteNode( node4 );
         return -1;
      } /* end if */
      lambdaCminus = numerator / denominator;
      x4 = x1 + lambdaCminus * cosCminus;
      y4 = y1 + lambdaCminus * sinCminus;
      dx = x4 - x4_old; dy = y4 - y4_old;
      change_in_position = sqrt( dx * dx + dy * dy );

      /* Lengths of the characteristic segments */
      dx = x4 - x1; dy = y4 - y1;
      lengthCminus = sqrt( dx * dx + dy * dy );
      dot_product = dx * xStream + dy * yStream;
      if ( dot_product < 0.0 ) {
         directionCminus = -1;
      } else {
         directionCminus = +1;
      } /* end if */

      dx = x4 - x0; dy = y4 - y0;
      lengthCzero  = sqrt( dx * dx + dy * dy );

      /*
       * Update flow properties at solution point
       * First, assume 2D planar geometry then add
       * axisymmetric contributions if flag is set.
       */
      pm4 = pm0;
      th4 = th1 + (pm1 - pm4);
      if ( GetAxiFlag() == 1 ) {
         if ( y1 < 1.0e-6 && y0 < 1.0e-6 ) {
            printf( "CMinusFreeBndyNode: both initial nodes are too close to axis.\n" );
            if ( node4_created ) DeleteNode( node4 );
            return -1;
         } /* end if */

         /* Axisymmetric components. */
         if ( y1 < 1.0e-6 ) {
            axiterm1 = sin(mu0) * sin(th0) / y0;
         } else {
            axiterm1 = sin(mu1) * sin(th1) / y1;
         } /* end if */

         axiterm4 = sin(mu4) * sin(th4) / y4;
         integralCminus = 0.5 * directionCminus * lengthCminus * (axiterm1 + axiterm4);

         /* Include axisymmetric components. */
         th4 += integralCminus;
      } /* end if */
 
   } while ( iteration_count < max_iteration && 
           change_in_position > position_tolerance );

   /*
    * Save the solution-point properties and connect the 
    * node into the characteristic mesh.
    */
   M4 = MFromNu( pm4, GetGamma() );
   n4->X = x4; n4->Y = y4; n4->Nu = pm4; n4->Theta = th4; n4->Mach = M4;
   if ( x4 > x0 ) {
      n4->CZeroUp    = node0;
      n0->CZeroDown  = node4;
   } else {
      n4->CZeroDown  = node0;
      n0->CZeroUp    = node4;
   } /* end if */
   if ( x4 > x1 ) {
      n4->CMinusUp   = node1;
      n1->CMinusDown = node4;
   } else {
      n4->CMinusDown = node1;
      n1->CMinusUp   = node4;
   } /* end if */

   /*
    * Assuming a successful calculation, 
    * return the index of the solution node.
    */
   return node4;
} /* end function CMinusFreeBndyNode */

/*------------------------------------------------------------------*/

/* @function */
int AddStreamNode( int node0, int node1, int node2, int node4,
                   int test_only ) {
   /**
    Purpose: Calculate a new streamline node, extending the streamline
             to the line joining nodeA and nodeB. <BR>
    Input  : <BR>
    node0 : index of initial point on the streamline <BR>
    node1 : index of first initial interpolation point <BR>
    node2 : index of second initial interpolation point <BR>
    node4 : index of solution point (may have a value of -1) <BR>
    If -1 is specified as the index for node4, a new node will be
    created for the solution point. <BR>
    test_only : flag to indicate whether we want to test for intersection only
                or if we actually want to add the node to the streamline <BR>
                test_only == 0 : add the node <BR>
                test_only == 1 : test for intersection only <BR>
    Output : <BR>
    if test_only == 0 :
    returns the index of the solution point or a value of 0
    if there has been a failure. <BR>
    if test_only == 1 : 
    returns 1 if intersection occurred between nodes 1 and 2.
    A value of 0 indicates that intersection did not occur 
    between nodes 1 and 2. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct NodeData *n0;
   struct NodeData *n1;
   struct NodeData *n2;
   struct NodeData *n4;
   int    node4_created;
   double x0, y0, pm0, th0, M0;
   double x1, y1, pm1, th1, M1;
   double x2, y2, pm2, th2, M2;
   double x4, y4, pm4, th4, M4;
   double alpha, dx14, dy14;
   double cosCzero, sinCzero, dx12, dy12;
   double lambda12, numerator, denominator;
   double x4_old, y4_old, dx, dy, change_in_position;
   int    iteration_count, valid_intersection;

   /*
    * First, get pointers to the node data.
    */
   if ( node1 == node2 ) {
      printf( "AddStreamNode: error, same index given for node1 and node2.\n" );
      return -1;
   } /* end if */

   n0 = GetNodePtr( node0 );
   if ( n0 == NULL ) {
      printf( "AddStreamNode: error: node0 with id=%d doesn't exist.\n", node0 );
      return -1;
   } /* end if */

   n1 = GetNodePtr( node1 );
   if ( n1 == NULL ) {
      printf( "AddStreamNode: error: node1 with id=%d doesn't exist.\n", node1 );
      return -1;
   } /* end if */

   n2 = GetNodePtr( node2 );
   if ( n2 == NULL ) {
      printf( "AddStreamNode: error: node2 with id=%d doesn't exist.\n", node2 );
      return -1;
   } /* end if */

   /*
    * Make some local copies of the node data.
    */
   x0 = n0->X; y0 = n0->Y; pm0 = n0->Nu; th0 = n0->Theta; M0 = n0->Mach;
   x1 = n1->X; y1 = n1->Y; pm1 = n1->Nu; th1 = n1->Theta; M1 = n1->Mach;
   x2 = n2->X; y2 = n2->Y; pm2 = n2->Nu; th2 = n2->Theta; M2 = n2->Mach;

   /*
    * Guess at some of the solution point properties.
    * The position will be way off but it is used as part of
    * the convergence check a little further on.
    */
   x4  = 0.5 * (x1  + x2);
   y4  = 0.5 * (y1  + y2);
   th4 = 0.5 * (th1 + th2);

   /*
    * Compute the solution point position and flow properties.
    */
   iteration_count = 0;
   do {
      ++iteration_count;
      x4_old = x4; y4_old = y4;

      /*
       * Locate solution point by assuming straight-line segments.
       */
      dx12     = x2 - x1;
      dy12     = y2 - y1;
      sinCzero = 0.5 * ( sin(th0) + sin(th4) );
      cosCzero = 0.5 * ( cos(th0) + cos(th4) );

      numerator = (x0 - x1) * sinCzero - (y0 - y1) * cosCzero;
      denominator = dx12 * sinCzero - dy12 * cosCzero;
      if ( fabs(denominator) <= 1.0e-12 ) {
         printf( "AddStreamNode: line segments are parallel.\n" );
         return -1;
      } /* end if */
      lambda12 = numerator / denominator;
      x4 = x1 + lambda12 * dx12;
      y4 = y1 + lambda12 * dy12;
      dx = x4 - x4_old; dy = y4 - y4_old;
      change_in_position = sqrt( dx * dx + dy * dy );

      /*
       * Update flow properties at solution point
       * Use linear interpolation between nodes 1 and 2.
       */
      dx14 = x4 - x1; dy14 = y4 - y1;
      if ( fabs(dx12) > 1.0e-6 ) {
         /* Assume that 1->2 is not vertical. */
         alpha = dx14 / dx12;
      } else {
         /* Assume that 1->2 is vertical. */
         alpha = dy14 / dy12;
      } /* end if */
      pm4 =  (1.0 - alpha) * pm1 + alpha * pm2;
      th4 =  (1.0 - alpha) * th1 + alpha * th2;
 
   } while ( iteration_count < max_iteration && 
           change_in_position > position_tolerance );

   if ( test_only == 0 ) {
      /*
       * Save the solution-point properties and connect the 
       * node into the characteristic mesh.
       */
      node4_created = 0;
      if ( node4 == -1 ) {
         node4 = CreateNode( -1 );
         node4_created = 1;
      } /* end if */

      n4 = GetNodePtr( node4 );
      if ( n4 == NULL ) {
         node4 = CreateNode( node4 );
         if ( node4 == -1 ) {
            printf( "AddStreamNode: error: couldn't create node %d.\n", node4 );
            return -1;
         } else {
            n4    = GetNodePtr( node4 );
         } /* end if */
      } /* end if */

      M4 = MFromNu( pm4, GetGamma() );
      n4->X = x4; n4->Y = y4; n4->Nu = pm4; n4->Theta = th4; n4->Mach = M4;
      if ( x4 > x0 ) {
         n4->CZeroUp    = node0;
         n0->CZeroDown  = node4;
      } else {
         n4->CZeroDown  = node0;
         n0->CZeroUp    = node4;
      } /* end if */

      /*
       * Assuming a successful calculation, 
       * return the index of the solution node.
       */
      return node4;
   } else {
      /*
       * Send back a flag to indicate whether or not the intersection
       * lies between nodes 1 and 2.
       */
      valid_intersection = ( alpha >= 0.0 && alpha <= 1.0 );
      return valid_intersection;
   } /* end if */

} /* end function AddStreamNode */

/*------------------------------------------------------------------*/

/* @function */
int StepStreamNode( int node0, int node4, double dL ) {
   /**
    Purpose: Calculate a new streamline node, extending the streamline
             by length dL <BR>
    Input  :  <BR>
    node0 : index of initial point on the streamline <BR>
    node4 : index of solution point (may have a value of -1) <BR>
    dL    : step-size along streamline;
            A positive value will step downstream while a negative
            value will step upstream. <BR>
    If -1 is specified as the index for node4, a new node will be
    created for the solution point. <BR>
    Output :  <BR>
    Returns the index of the solution point or a value of -1
    if there has been a failure.  One possible failure is that
    there are no nodes close enough to include in the interpolation phase. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct NodeData *n0;
   struct NodeData *np;
   struct NodeData *n4;
   int    node4_created;
   double x0, y0, pm0, th0, M0;
   double x4, y4, pm4, th4, M4;
   double cosCzero, sinCzero;
   #define MAX_NEAR 30
   int    j, count, near_node_count, near_node_array[MAX_NEAR];
   double x[MAX_NEAR], y[MAX_NEAR];
   double R, mu, sum_Xi, Xi[MAX_NEAR], r[MAX_NEAR], w[MAX_NEAR];
   double th[MAX_NEAR], pm[MAX_NEAR];
   double x4_old, y4_old, dx, dy, change_in_position;
   int    iteration_count;

   /*
    * First, get pointers to the node data.
    */
   n0 = GetNodePtr( node0 );
   if ( n0 == NULL ) {
      printf( "StepStreamNode: error: node0 with id=%d doesn't exist.\n", node0 );
      return -1;
   } /* end if */

   node4_created = 0;
   if ( node4 == -1 ) {
      node4 = CreateNode( -1 );
      node4_created = 1;
   } /* end if */

   n4 = GetNodePtr( node4 );
   if ( n4 == NULL ) {
      node4 = CreateNode( node4 );
      if ( node4 == -1 ) {
         printf( "StepStreamNode: error: couldn't create node %d.\n", node4 );
         return -1;
      } else {
         n4    = GetNodePtr( node4 );
      } /* end if */
   } /* end if */

   /*
    * Make some local copies of the node data.
    */
   x0 = n0->X; y0 = n0->Y; pm0 = n0->Nu; th0 = n0->Theta; M0 = n0->Mach;

   /*
    * Guess at some of the solution point properties.
    * The position may be a little off but it is used to 
    * do the search for local nodes and as part of
    * the convergence check a little further on.
    * Hopefully, flow angle does not change dramatically
    * along the streamline.
    */
   x4 = x0 + dL * cos(th0); 
   y4 = y0 + dL * sin(th0); 
   th4 = th0; pm4 = pm0;

   R = 0.9 * dL; /* Radius of influence. */
   mu = 2.0; /* smoothing parameter for Shepard interpolation. */
   near_node_count = FindNodesNear( x4, y4, R, near_node_array, MAX_NEAR );
   if ( GetDebugLevel() >= 1 ) {
      printf( "StepStreamNode: %d near-by nodes found.\n", near_node_count);
   } /* end if */
   count = 0;
   for (j = 0; j < near_node_count; ++j) {
      np    = GetNodePtr( near_node_array[j] );
      if ( np != n4 && np != n0 ) {
         x[count]  = np->X;
         y[count]  = np->Y;
         pm[count] = np->Nu;
         th[count] = np->Theta;
         ++count;
      } /* end if */
   } /* end for */
   near_node_count = count;
   if ( GetDebugLevel() >= 1 ) {
      printf( "StepStreamNode: %d near-by nodes after excluding n4, n0.\n", 
         near_node_count );
   } /* end if */
   if (near_node_count == 0) {
      printf( "StepStreamNode: No near-by nodes.\n" );
      if ( node4_created ) DeleteNode( node4 );
      return -1;
   } /* end if */

   /*
    * Compute the solution point position and flow properties.
    */
   iteration_count = 0;
   do {
      ++iteration_count;
      x4_old = x4; y4_old = y4;

      /*
       * Interpolate the flow properties by Shepard interpolation
       * over the near-by nodes.
       */
      sum_Xi = 0.0;
      for (j = 0; j < near_node_count; ++j) {
         dx   = x4 - x[j];
         dy   = y4 - y[j];
         r[j] = sqrt(dx * dx + dy * dy);
         Xi[j] = pow(1.0 - r[j]/R, mu);
         sum_Xi += Xi[j];
      } /* end for */
      for (j = 0; j < near_node_count; ++j) {
         w[j] = Xi[j] / sum_Xi;
      } /* end for */
      pm4 = 0.0;
      th4 = 0.0;
      for (j = 0; j < near_node_count; ++j) {
         pm4 += w[j] * pm[j];
         th4 += w[j] * th[j];
      } /* end for */

      /*
       * Locate solution point by using average slope.
       */
      sinCzero = 0.5 * ( sin(th0) + sin(th4) );
      cosCzero = 0.5 * ( cos(th0) + cos(th4) );

      x4 = x0 + cosCzero * dL;
      y4 = y0 + sinCzero * dL;
      dx = x4 - x4_old; dy = y4 - y4_old;
      change_in_position = sqrt( dx * dx + dy * dy );
 
   } while ( iteration_count < max_iteration && 
           change_in_position > position_tolerance );

   /*
    * Save the solution-point properties and connect the 
    * node into the characteristic mesh.
    */
   M4 = MFromNu( pm4, GetGamma() );
   n4->X = x4; n4->Y = y4; n4->Nu = pm4; n4->Theta = th4; n4->Mach = M4;
   if ( x4 > x0 ) {
      n4->CZeroUp    = node0;
      n0->CZeroDown  = node4;
   } else {
      n4->CZeroDown  = node0;
      n0->CZeroUp    = node4;
   } /* end if */

   /*
    * Assuming a successful calculation, 
    * return the index of the solution node.
    */
   return node4;
} /* end function StepStreamNode */

/*------------------------------------------------------------------*/

/* @function */
int InterpolateNode( double x_point, double y_point, double R, int node4 ) {
   /**
    Purpose: Locate a new node at coordinates (x,y), interpolating
       the node's properties from other near-by nodes. <BR>
    Input  :  <BR>
    x_point, 
    y_point  : coordinates of the new node <BR>
    R     : radius-of-influence for the Shepard interpolation <BR>
    node4 : index of solution point (may have a value of -1) <BR>
    If -1 is specified as the index for node4, a new node will be
    created for the solution point. <BR>
    Output :  <BR>
    Returns the index of the solution point or a value of -1
    if there has been a failure.  One possible failure is that
    there are no nodes close enough to include in the interpolation. <BR>
    (Available from the Tcl interpreter.) 
    */
   struct NodeData *np;
   struct NodeData *n4;
   int    node4_created;
   double x4, y4, pm4, th4, M4;
   double dx, dy;
   #define MAX_NEAR 30
   int    j, count, near_node_count, near_node_array[MAX_NEAR];
   double x[MAX_NEAR], y[MAX_NEAR];
   double mu, sum_Xi, Xi[MAX_NEAR], r[MAX_NEAR], w[MAX_NEAR];
   double th[MAX_NEAR], pm[MAX_NEAR];

   /*
    * First, get pointers to the node data.
    */
   node4_created = 0;
   if ( node4 == -1 ) {
      node4 = CreateNode( -1 );
      node4_created = 1;
   } /* end if */

   n4 = GetNodePtr( node4 );
   if ( n4 == NULL ) {
      node4 = CreateNode( node4 );
      if ( node4 == -1 ) {
         printf( "InterpolateNode: error: couldn't create node %d.\n", node4 );
         return -1;
      } else {
         n4    = GetNodePtr( node4 );
      } /* end if */
   } /* end if */

   mu = 2.0; /* smoothing parameter for Shepard interpolation. */
   near_node_count = FindNodesNear( x_point, y_point, R, near_node_array, MAX_NEAR );
   if ( GetDebugLevel() >= 1 ) {
      printf( "InterpolateNode: %d near-by nodes found.\n", near_node_count);
   } /* end if */
   count = 0;
   for (j = 0; j < near_node_count; ++j) {
      np    = GetNodePtr( near_node_array[j] );
      if ( np != n4 ) {
         x[count]  = np->X;
         y[count]  = np->Y;
         pm[count] = np->Nu;
         th[count] = np->Theta;
         ++count;
      } /* end if */
   } /* end for */
   near_node_count = count;
   if ( GetDebugLevel() >= 1 ) {
      printf( "InterpolateNode: %d near-by nodes after excluding n4.\n", 
         near_node_count );
   } /* end if */
   if (near_node_count == 0) {
      printf( "InterpolateNode: No near-by nodes.\n" );
      if ( node4_created ) DeleteNode( node4 );
      return -1;
   } /* end if */

   /*
    * The location of the solution point is given.
    */
   x4 = x_point; y4 = y_point;

   /*
    * Interpolate the flow properties by Shepard interpolation
    * over the near-by nodes.
    */
   sum_Xi = 0.0;
   for (j = 0; j < near_node_count; ++j) {
      dx   = x4 - x[j];
      dy   = y4 - y[j];
      r[j] = sqrt(dx * dx + dy * dy);
      Xi[j] = pow(1.0 - r[j]/R, mu);
      sum_Xi += Xi[j];
   } /* end for */
   for (j = 0; j < near_node_count; ++j) {
      w[j] = Xi[j] / sum_Xi;
   } /* end for */
   pm4 = 0.0;
   th4 = 0.0;
   for (j = 0; j < near_node_count; ++j) {
      pm4 += w[j] * pm[j];
      th4 += w[j] * th[j];
   } /* end for */

   /*
    * Save the solution-point properties.
    */
   M4 = MFromNu( pm4, GetGamma() );
   n4->X = x4; n4->Y = y4; n4->Nu = pm4; n4->Theta = th4; n4->Mach = M4;

   /*
    * Assuming a successful calculation, 
    * return the index of the solution node.
    */
   return node4;
} /* end function InterpolateNode */

/*------------------------------------------------------------------*/

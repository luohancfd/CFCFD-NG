/** \file ssp3d.c
 * \ingroup plot
 * \brief Three-dimensional plotting functions -- need fixing.
 *
 * A set of three-dimensional routines for the ssp package.
 *
 * The routines here are based on the conventions in
 * David F. Rogers and J. Alan Adams
 * Mathematical Elements for Computer Graphics (2nd Ed.)
 * McGraw Hill Publishing Company 1990
 *
 * \author PA Jacobs
 * \version 1.0, November 1992
 * \version 1.1,   25-Apr-93, Added ssp_ApplyProjectionMatrix()
 *                   Filled out the z-components of the trimetric
 *                   projection matrix.
 *                   Polygons in 3D.
 * \version 1.11,  04-Jun-93, Orthographic projection
 * \version 1.2,   08-Jun-93, 3-D axes with tick marks
 *
 */

/*
 * ------------------
 * Global Definitions
 * ------------------
 */

#include "../../util/source/compiler.h"
#include "../../geometry/source/geom.h"
#include "ssp.h"

#include <stdio.h>
#include <math.h>

#if (STDLIBH)
#  include <stdlib.h>
#  include <stdarg.h>
#endif

/*
 * -----------
 * Global Data
 * -----------
 */

struct point_3D origin3D;

double proj_matrix[3][3];

/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetOrigin3D (double x, double y, double z)

#else

int ssp_SetOrigin3D (x, y, z)
double x, y, z;

#endif

/*
 * Purpose...
 * -------
 * Set the origin in three-space for rotations and projections.
 * This point in three-space will correspond to the current origin
 * on the plotted page or screen.
 *
 * Input...
 * -----
 * x, y, z : location of the new origin in user coordinates
 *
 */

{  /* begin ssp_SetOrigin3D() */

origin3D.x = x;
origin3D.y = y;
origin3D.z = z;

return (0);
}  /* end of ssp_SetOrigin3D() */

/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetTrimetricProjection (double phi, double theta)

#else

int ssp_SetTrimetricProjection (phi, theta)
double phi, theta;

#endif

/*
 * Purpose...
 * -------
 * Set the projection matrix for a trimetric projection from three-space
 * onto the (x,y) plotting page.
 *
 * This projection can be specified in three stages...
 * (1) rotate about the y-axis by phi radians,
 * (2) rotate about the x-axis by theta radians,
 * (3) use parallel projection onto the z-plane (which is
 *     the plotting page)
 * Note that these rotations are performed about the three-space origin.
 *
 * Input...
 * -----
 * phi   : angle of rotation around the y-axis (radians)
 * theta : angle of rotation around the x-axis (radians)
 *
 */

{  /* begin ssp_SetTrimetricProjection() */
double sin_phi, cos_phi, sin_th, cos_th;

sin_phi = sin(phi);
cos_phi = cos(phi);
sin_th  = sin(theta);
cos_th  = cos(theta);

proj_matrix[0][0] = cos_phi;
proj_matrix[0][1] = sin_phi * sin_th;
proj_matrix[0][2] = -sin_th * cos_th;

proj_matrix[1][0] = 0.0;
proj_matrix[1][1] = cos_th;
proj_matrix[1][2] = sin_th;

proj_matrix[2][0] = sin_phi;
proj_matrix[2][1] = -cos_phi * sin_th;
proj_matrix[2][2] = cos_phi * cos_th;

return (0);
}  /* end of ssp_SetTrimetricProjection() */


/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetIsometricProjection (void)

#else

int ssp_SetIsometricProjection ()

#endif

/*
 * Purpose...
 * -------
 * Set the projection matrix for one of the isometric projections.
 * (The one that seemed to be most natural to me, anyway.)
 *
 */

{  /* begin ssp_SetIsometricProjection() */
double phi_y, theta_x;

phi_y = 45.0 / 180.0 * ssp_PI;
theta_x = 35.26 / 180.0 * ssp_PI;

ssp_SetTrimetricProjection (phi_y, theta_x);

return (0);
}  /* end of ssp_SetIsometricProjection() */


/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetOrthographicProjection (int iplane)

#else

int ssp_SetOrthographicProjection (iplane)
int iplane;

#endif

/*
 * Purpose...
 * -------
 * Set the projection matrix for one of the orthographic projections.
 *
 * Input...
 * -----
 * iplane  : integer selection for the type of projection
 *           0 : XY plane projection
 *           1 : ZY plane projection
 *           2 : XZ plane projection
 *
 */

{  /* begin ssp_OrthographicProjection() */
double phi_y, theta_x;

if (iplane == 2)
   { /* XZ plane projection */
   phi_y = 0.0;
   theta_x = -0.5 * ssp_PI;
   }
else if (iplane == 1)
   { /* ZY plane projection */
   phi_y = 0.5 * ssp_PI;
   theta_x = 0.0;
   }
else
   { /* XY plane projection */
   phi_y = 0.0;
   theta_x = 0.0;
   }

ssp_SetTrimetricProjection (phi_y, theta_x);

return (0);
}  /* end of ssp_SetOrthographicProjection() */


/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_Plot3D (double x, double y, double z)

#else

int ssp_Plot3D (x, y, z)
double x, y, z;

#endif

/*
 * Purpose...
 * -------
 * Draw (a projection of) a line from the previous point in three-space
 * to the currently specified point in three-space.
 * The current projection matrix is used to transform to 2D plotting
 * coordinates about the currently specified origin in three-space.
 *
 * Note...
 * ----
 * Before using this routine, the projection matrix and the three-space
 * origin must have valid settings.
 *
 * Input...
 * -----
 * x,y,z  : coordinates for the new position (user coordinates)
 *
 */

{  /* begin ssp_Plot3D() */
double x_2D, y_2D;

/* Transform to two-space */
x_2D =   (x - origin3D.x) * proj_matrix[0][0]
       + (y - origin3D.y) * proj_matrix[1][0]
       + (z - origin3D.z) * proj_matrix[2][0];
y_2D =   (x - origin3D.x) * proj_matrix[0][1]
       + (y - origin3D.y) * proj_matrix[1][1]
       + (z - origin3D.z) * proj_matrix[2][1];

/* Use the standard ssp 2D plotting routine. */
ssp_Plot (x_2D, y_2D);

return (0);
}  /* end of ssp_Plot3D() */

/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_Move3D (double x, double y, double z)

#else

int ssp_Move3D (x, y, z)
double x, y, z;

#endif

/*
 * Purpose...
 * -------
 * Move to (the projection of) the currently specified point in three-space.
 * The current projection matrix is used to transform to 2D plotting
 * coordinates about the currently specified origin in three-space.
 *
 * Note...
 * ----
 * Before using this routine, the projection matrix and the three-space
 * origin must have valid settings.
 *
 * Input...
 * -----
 * x,y,z  : coordinates for the new position (user coordinates)
 *
 */

{  /* begin ssp_Move3D() */
double x_2D, y_2D;

/* Transform to two-space */
x_2D =   (x - origin3D.x) * proj_matrix[0][0]
       + (y - origin3D.y) * proj_matrix[1][0]
       + (z - origin3D.z) * proj_matrix[2][0];
y_2D =   (x - origin3D.x) * proj_matrix[0][1]
       + (y - origin3D.y) * proj_matrix[1][1]
       + (z - origin3D.z) * proj_matrix[2][1];

/* Use the standard ssp 2D plotting routine. */
ssp_Move (x_2D, y_2D);

return (0);
}  /* end of ssp_Move3D() */

/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_ApplyProjectionMatrix (struct point_3D *pt1,
                               struct point_3D *pt2)

#else

int ssp_ApplyProjectionMatrix (pt1, pt2)
struct point_3D *pt1;
struct point_3D *pt2;

#endif

/*
 * Purpose...
 * -------
 * Apply the progection matrix to get the x,y,z coordinates in
 * viewing space (e.g. after ofsetting for the zero and rotating
 * for trimetric projection)
 *
 * Note...
 * ----
 * Before using this routine, the projection matrix and the three-space
 * origin must have valid settings.
 *
 * Input...
 * -----
 * *pt1   : pointer to the coordinates in model space
 * *pt2   : pointer to the coordinates in viewing space
 *
 */

{  /* begin ssp_ApplyProjectionMatrix() */
double x, y, z;

x = pt1->x;
y = pt1->y;
z = pt1->z;

pt2->x =   (x - origin3D.x) * proj_matrix[0][0]
         + (y - origin3D.y) * proj_matrix[1][0]
         + (z - origin3D.z) * proj_matrix[2][0];

pt2->y =   (x - origin3D.x) * proj_matrix[0][1]
         + (y - origin3D.y) * proj_matrix[1][1]
         + (z - origin3D.z) * proj_matrix[2][1];

pt2->z =   (x - origin3D.x) * proj_matrix[0][2]
         + (y - origin3D.y) * proj_matrix[1][2]
         + (z - origin3D.z) * proj_matrix[2][2];

return (0);
}  /* end of ssp_ApplyProjectionMatrix() */

/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_Polygon3D (int n, struct point_3D vtx[],
                   int grey_fill, int colour_fill, int border_flag)

#else

int ssp_Polygon3D (n, vtx, grey_fill, colour_fill, border_flag)
int    n;
struct point_3D vtx[];
int    grey_fill, colour_fill, border_flag;

#endif

/*
 * Purpose...
 * -------
 * Draw a polygon in 3D space and optionally fill it.
 *
 * Input...
 * -----
 * n         : number of vertices (maximum 10)
 * vtx[]     : coordinates of the vertices (in plot units)
 * grey_fill : = 0, don't fill the polygon with grey
 *             = 1, fill the interior of the polygon with grey
 * colour_fill : = 0, don't fill the polygon with colour
 *               = 1, fill the interior of the polygon with colour
 * border_flag : = 0, don't draw the border
 *               = 1, draw the border with the currently selected pen
 *
 */

{  /* begin ssp_Polygon3D() */
int    i;
#define  VTX_DIM  10
struct point_3D vtx_view[VTX_DIM];
double x[VTX_DIM], y[VTX_DIM];

/*
 * Convert the vertex coordinates to viewing space.
 */
for (i = 0; i < n; ++i)
   {
   ssp_ApplyProjectionMatrix ( &(vtx[i]), &(vtx_view[i]) );
   x[i] = vtx_view[i].x;
   y[i] = vtx_view[i].y;
   }

/*
 * Use the 2D function to draw the polygon on the page.
 */
ssp_Polygon (n, x, y, grey_fill, colour_fill, border_flag);

return 0;
}  /* end of ssp_Polygon3D() */

/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_axes3D (double xaxis, double yaxis, double zaxis,
                double xtick, double ytick, double ztick,
                double ticksize)

#else

int ssp_axes3D (xaxis, yaxis, zaxis,
                xtick, ytick, ztick,
                ticksize)
double xaxis, yaxis, zaxis;
double xtick, ytick, ztick;
double ticksize;

#endif

{  /* Begin ssp_axes3D() */
double xx, yy, zz;

/*
 * Draw the main lines for the axes.
 */
ssp_Move3D (0.0, 0.0, 0.0);
ssp_Plot3D (xaxis, 0.0, 0.0);

ssp_Move3D (0.0, 0.0, 0.0);
ssp_Plot3D (0.0, yaxis, 0.0);

ssp_Move3D (0.0, 0.0, 0.0);
ssp_Plot3D (0.0, 0.0, zaxis);

/*
 * Now draw the tick marks.
 */
xx = 0.0;
yy = 0.0;
zz = 0.0;
for (xx = 0.0; xx <= 1.001 * xaxis; xx += xtick)
   {
   ssp_Move3D (xx, yy+ticksize, zz);
   ssp_Plot3D (xx, yy, zz);
   ssp_Plot3D (xx, yy, zz+ticksize);
   }

xx = 0.0;
yy = 0.0;
zz = 0.0;
for (yy = 0.0; yy <= 1.001 * yaxis; yy += ytick)
   {
   ssp_Move3D (xx+ticksize, yy, zz);
   ssp_Plot3D (xx, yy, zz);
   ssp_Plot3D (xx, yy, zz+ticksize);
   }

xx = 0.0;
yy = 0.0;
zz = 0.0;
for (zz = 0.0; zz <= 1.001 * zaxis; zz += ztick)
   {
   ssp_Move3D (xx, yy+ticksize, zz);
   ssp_Plot3D (xx, yy, zz);
   ssp_Plot3D (xx+ticksize, yy, zz);
   }


return 0;
}  /* end of ssp_axes3D() */

/*-------------------------------------------------------------------*/


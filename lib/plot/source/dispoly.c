/** \file dispoly.c
 * \ingroup plot
 * \brief Display 3D polygons -- needs fixing.
 *
 * Read a file containing 3D points and an associated
 * file containing indices for polygon vertices and then
 * plot an orthographic projection of the polygon set.
 *
 * \author PA Jacobs
 *
 * \version 1.0  :    19-Apr-93 : First Cut
 * \version 1.1  :    25-Apr-93 : Depth sorting
 *                    Plot facets using the polygon functions.
 * \version 1.11 :    31-May-93 : rename points to vertices
 * \version 1.2  :    01-Jun-93 : allocate memory for the big arrays
 * \version 1.3  :    02-Jun-93 : fix problem in depth calc. for triangles
 *                    and lines.
 * \version 1.4  :    04-Jun-93 : Added Orthographic projection
 * \version 2.0  :    06-Jun-93 : added contouring
 * \version 2.1  :    08-Jun-93 : better contouring for quadrilaterals
 * \version 2.11 :    09-Jun-93 : activate ssp_ContourQuad3D()
 * \version 2.2  :    17-Aug-93 : fixed bug in the basic contouring routine
 *                    Added a section to plot line segments.
 * \version 2.21 :    20-Oct-93 : Adjust default contouring limits and reset
 *                    the x,y or z ranges to something finite is
 *                    any one (but only one) is near zero.
 *                    Coloured contours.
 * \version 2.3  :    24-Oct-93 : Added colour map down the RHS of the page.
 * \version 2.31 :    10-Nov-93 : Changed size of written text at top of page.
 *                    Changed the limiting length of the x-axis.
 * \version 3.0  :    01-Jan-94 : Added colour-filling to show contour levels.
 * \version 3.1  :    11-Jan-94 : Put the value filename onto the plot.
 *                    Get the line filename early on.
 *                    Add counter for the polygons.
 * \version 3.2  :    12-Jun-94 : Grey-scale shading as well as colour-filling
 *                    Move the plot to the right by 5mm to cope with
 *                    Ghostscript
 * \version 3.3  :    10-Jul-94 : Illumination model for grey-shading.
 * \version 3.31 :    12-Jul-94 : brighter colours and lighter illumination model
 *                    thinner lines for facet boundaries
 * \version 3.32 :    20-Aug-94 : Some variable name changes to suit the IBM RISC
 *                    compiler.
 * \version 3.33 :    31-Aug-94 : More compact colour postscript.
 * \version 3.34 :    12-Mar-95 : Do not plot ploygons that have out-of-range
 *                    vertex indices.
 * \version 3.35 :    25-Apr-95 : Change the name hue to hue_d to keep the BCOS2
 *                    linker happy.  It doesn't like global data (in
 *                    separate modules) having the same name.
 *
 */

/*-------------------------------------------------------*/

#include <math.h>
#include <stdio.h>

#include "compiler.h"
#include "cmath.h"

#if (STDLIBH)
#  include <stdlib.h>
#  include <string.h>
#endif

#if (CONIOH)
#  include <conio.h>
#endif

#include "geom.h"
#include "ssp.h"

#if PROTO
   double reflected_light_intensity ( void );
#else
   double reflected_light_intensity ();
#endif

/*-------------------------------------------------------*/

/*
 * Put most of the storage as global data so that we can
 * use subroutines without explicitly passing data
 */

#define  NCHAR  256
char   txt[NCHAR], VertexFileName[32], PolyFileName[32];
char   PlotFileName[32], ValueFileName[32], LineFileName[32];
FILE   *vertex_file, *poly_file, *value_file, *line_file;

struct point_3D *vtx, pt2, *pvtx;
struct point_3D poly_vtx[4];
int    n_vertex, i_vertex, n_value;
static double *z_vtx, *v_vtx;

struct poly_3D *polygn, *ppolygn;
int    n_poly, i_poly, ip, A, B, C, D;
double *z_depth;
int    *z_rank;

int    count, climit;
char   msg[64];

double xstart, xend, xtick;
double ystart, yend, ytick;
double zstart, zend, ztick;
double vstart, vend, vtick;
double minscale, xscale, yscale, zscale;
double xrange, yrange, zrange, vrange;
double level, delA, delB, delC, delD;
int    missed;

int    screen, pfile;
double xplotsize, yplotsize, xorig, yorig;
double symbolsize;
double xaxislength, yaxislength, zaxislength;
int    change, trueshape;
int    do_contours, fill_facets, vertex_symbols;
int    do_lines, facet_boundaries, do_colours;
int    do_filled_colours, do_filled_grey;
double xx, yy, zz, vv, sc, ang;
int    isymb, iproj, finished, i;
double xx1, yy1, zz1, xx2, yy2, zz2;
double hue_d, bright, sat;

double Intensity;


/*-------------------------------------------------------*/

main ()
{

printf ("\n");
printf ("----------------\n");
printf ("Display Polygons\n");
printf ("----------------\n");
printf ("Rev. 3.35, 25-Apr-95\n");
printf ("\n");


/*
 * Set the options.
 */
printf ("Enter a 1 or a 0 for the following options...\n");
printf ("(a) Vertex Symbols : ");
scanf ("%d", &vertex_symbols );
printf ("(b) Facet Boundaries : ");
scanf ("%d", &facet_boundaries );
printf ("(c) Fill Facets (with illumination model) : ");
scanf ("%d", &fill_facets );
printf ("(d) Contours : ");
scanf ("%d", &do_contours );
printf ("(e) Colours for contours : ");
scanf ("%d", &do_colours );
printf ("(f) Filled-colour mapping (instead of contour lines) : ");
scanf ("%d", &do_filled_colours );
printf ("(f) Filled-grey-scale mapping : ");
scanf ("%d", &do_filled_grey );
printf ("(g) Extra Line Segments : ");
scanf ("%d", &do_lines );

/*
 * Read the Vertex file
 */

printf ("\nNow the file names...\n");
printf ("Enter Vertex File Name : ");
scanf ("%s", VertexFileName);
if ((vertex_file = fopen(VertexFileName, "r")) == NULL)
   {
   printf ("Could not open Vertex file: %s\n", VertexFileName);
   exit (-1);
   }

fgets (txt, NCHAR, vertex_file);
sscanf (txt, "%d", &n_vertex);
printf ("Number of vertices: %d\n", n_vertex);
vtx = NULL;
vtx = (struct point_3D *) malloc (n_vertex * sizeof(struct point_3D));
if (vtx == NULL)
   {
   printf ("Could not allocate memory for %d vertices.\n", n_vertex);
   exit (-1);
   }
z_vtx = NULL;
z_vtx = (double *) malloc (n_vertex * sizeof(double));
if (z_vtx == NULL)
   {
   printf ("Could not allocate memory for %d vertex depths.\n", n_vertex);
   exit (-1);
   }

for (i = 0; i < n_vertex; ++i)
   {
   pvtx = &vtx[i];

   fgets (txt, NCHAR, vertex_file);
   sscanf (txt, "%lf %lf %lf", &xx, &yy, &zz );
#  if 0
   printf ("vertex[%d]: (%f, %f, %f)\n", i, xx, yy, zz );
#  endif

   pvtx->x = xx;
   pvtx->y = yy;
   pvtx->z = zz;
   }

/*
 * Read the Value file only if we are going to contour the facets.
 */

if ( do_contours || do_filled_colours || do_filled_grey )
   {
   printf ("Enter Value File Name  : ");
   scanf ("%s", ValueFileName);
   if ((value_file = fopen(ValueFileName, "r")) == NULL)
      {
      printf ("Could not open Value file: %s\n", ValueFileName);
      exit (-1);
      }

   fgets (txt, NCHAR, value_file);
   sscanf (txt, "%d", &n_value);
   if (n_value != n_vertex)
      {
      printf ("Mismatch %d vertices and %d values.\n", n_vertex, n_value);
      exit (-1);
      }
   v_vtx = NULL;
   v_vtx = (double *) malloc (n_vertex * sizeof(double));
   if (v_vtx == NULL)
      {
      printf ("Could not allocate memory for %d values.\n", n_vertex);
      exit (-1);
      }

   for (i = 0; i < n_vertex; ++i)
      {
      fgets (txt, NCHAR, value_file);
      sscanf (txt, "%lf",  &vv );
#     if 0
      printf ("value[%d]: %f\n", i, vv );
#     endif
      v_vtx[i] = vv;
      }
   }  /* end of reading the value-file */

/*
 * Read the Polygon file
 */

printf ("Enter Poly File Name   : ");
scanf ("%s", PolyFileName);
if ((poly_file = fopen(PolyFileName, "r")) == NULL)
   {
   printf ("Could not open Polygon file: %s\n", PolyFileName);
   exit (-1);
   }

fgets (txt, NCHAR, poly_file);
sscanf (txt, "%d", &n_poly);
printf ("Number of polygons: %d\n", n_poly);
polygn = NULL;
polygn = (struct poly_3D *) malloc (n_poly * sizeof(struct poly_3D));
if (polygn == NULL)
   {
   printf ("Could not allocate memory for %d polygons.\n", n_poly);
   exit (-1);
   }
z_depth = NULL;
z_depth = (double *) malloc (n_poly * sizeof(double));
if (z_depth == NULL)
   {
   printf ("Could not allocate memory for %d polygon depths.\n", n_poly);
   exit (-1);
   }
z_rank = NULL;
z_rank = (int *) malloc (n_poly * sizeof(int));
if (z_rank == NULL)
   {
   printf ("Could not allocate memory for %d polygon rank.\n", n_poly);
   exit (-1);
   }

for (i = 0; i < n_poly; ++i)
   {
   ppolygn = &(polygn[i]);

   fgets (txt, NCHAR, poly_file);
   sscanf (txt, "%d %d %d %d", &A, &B, &C, &D );
#  if 0
   printf ("polygn[%d]: (%d, %d, %d, %d)\n", A, B, C, D );
#  endif

   /* Weed out invalid vertex indices. */
   if ( A >= n_vertex || B >= n_vertex || C >= n_vertex || D >= n_vertex )
      {
      A = -1; B = -1; C = -1; D = -1;
      printf ("Warning: polygon[%d] has invalid vertices.\n", i);
      }

   ppolygn->A = A;
   ppolygn->B = B;
   ppolygn->C = C;
   ppolygn->D = D;
   }

/*
 * Clean-up after the file reading.
 */
fclose (vertex_file);
fclose (poly_file);
if ( do_contours || do_filled_colours ) fclose (value_file);

/*
 * Get the line file name now, open it but don't read it until later.
 */
if ( do_lines )
   {
   printf ("Enter Line-Segment File Name : ");
   scanf ("%s", LineFileName);
   if ((line_file = fopen(LineFileName, "r")) == NULL)
      {
      printf ("Could not open %s.", LineFileName);
      finished = 1;
      }
   }  /* if ( do_lines... */

/*
 ********************************
 * Set the Plotting Parameters. *
 ********************************
 *
 * Default plotting parameters.
 * Plotting units are millimetres (on an A4 page).
 */
xplotsize = 100.0;
yplotsize = 100.0;
xorig = 45.0;
yorig = 50.0;
trueshape = 1;
iproj = 1;

printf ("\n");
printf ("Projection 0=XY, 1=ZY, 2=XZ, -1=Isometric : ");
scanf ("%d", &iproj);

/*
 * Find extreme values and decide on scales.
 */
printf ("Find extremes...\n");
xstart = vtx[0].x;
xend   = vtx[0].x;
ystart = vtx[0].y;
yend   = vtx[0].y;
zstart = vtx[0].z;
zend   = vtx[0].z;
if ( do_contours || do_filled_colours || do_filled_grey )
   {
   vstart = v_vtx[0];
   vend   = v_vtx[0];
   }
for (i_vertex = 1; i_vertex < n_vertex; ++i_vertex)
   {
   if (vtx[i_vertex].x < xstart) xstart = vtx[i_vertex].x;
   if (vtx[i_vertex].x > xend  ) xend   = vtx[i_vertex].x;
   if (vtx[i_vertex].y < ystart) ystart = vtx[i_vertex].y;
   if (vtx[i_vertex].y > yend  ) yend   = vtx[i_vertex].y;
   if (vtx[i_vertex].z < zstart) zstart = vtx[i_vertex].z;
   if (vtx[i_vertex].z > zend  ) zend   = vtx[i_vertex].z;
   if ( do_contours || do_filled_colours || do_filled_grey )
      {
      if (v_vtx[i_vertex] < vstart) vstart = v_vtx[i_vertex];
      if (v_vtx[i_vertex] > vend  ) vend   = v_vtx[i_vertex];
      }
   }

/*
 * Set up nominal ranges.
 */
yrange = yend - ystart;
ytick  = yrange / 5.0;
xrange = xend - xstart;
xtick  = xrange / 5.0;
zrange = zend - zstart;
ztick  = zrange / 5.0;
if ( do_contours || do_filled_colours || do_filled_grey )
   {
   vrange = vend - vstart;
   vtick  = vrange / 10.0;
   /* Bring the ends in by half an interval */
   vstart += 0.5 * vtick;
   vend   -= 0.5 * vtick;
   }

printf ("Nominal ranges...\n");
printf ("xstart = %e, xend = %e, xtick = %e\n", xstart, xend, xtick);
printf ("ystart = %e, yend = %e, ytick = %e\n", ystart, yend, ytick);
printf ("zstart = %e, zend = %e, ztick = %e\n", zstart, zend, ztick);
if ( do_contours || do_filled_colours || do_filled_grey )
   printf ("vstart = %e, vend = %e, vtick = %e\n", vstart, vend, vtick);
if ( fabs(xrange) < 1.0e-10 )
   {
   printf ("Warning: X-range too small, reset to 0.1\n");
   xend = xstart + 0.1;
   xtick = 0.1;
   }
if ( fabs(yrange) < 1.0e-10 )
   {
   printf ("Warning: Y-range too small, reset to 0.1\n");
   yend = ystart + 0.1;
   ytick = 0.1;
   }
if ( fabs(zrange) < 1.0e-10 )
   {
   printf ("Warning: Z-range too small, reset to 0.1\n");
   zend = zstart + 0.1;
   ztick = 0.1;
   }
if ( (do_contours || do_filled_colours || do_filled_grey)
     && fabs(vrange) < 1.0e-10 )
   {
   printf ("Warning: V-range too small, set to 0.1\n");
   vend = vstart + 0.1;
   vtick = 0.1;
   }

/*
 * Alter the nominal ranges if desired.
 */
printf ("Change ranges (1/0)? ");
scanf ("%d", &change);
if (change == 1)
   {
   printf ("xstart xend xtick: ");
   scanf ("%lf %lf %lf", &xstart, &xend, &xtick);
   printf ("ystart yend ytick: ");
   scanf ("%lf %lf %lf", &ystart, &yend, &ytick);
   printf ("zstart zend ztick: ");
   scanf ("%lf %lf %lf", &zstart, &zend, &ztick);
   xrange = xend - xstart;
   yrange = yend - ystart;
   zrange = zend - zstart;
   if ( do_contours || do_filled_colours || do_filled_grey )
      {
      printf ("vstart vend vtick: ");
      scanf ("%lf %lf %lf", &vstart, &vend, &vtick);
      vrange = vend - vstart;
      }
   }

/*
 * Finalize the Plot scales.
 */
printf ("trueshape (1/0): ");
scanf ("%d", &trueshape);
xscale = xplotsize / xrange;
yscale = yplotsize / yrange;
zscale = xplotsize / zrange;
minscale = xscale;
if (yscale < minscale) minscale = yscale;
if (zscale < minscale) minscale = zscale;
if (trueshape == 1)
   {
   xscale = minscale;
   yscale = minscale;
   zscale = minscale;
   }

/*
 * Scale the spatial variables.
 */
printf ("Scale the vertices...\n");
for (i_vertex = 0; i_vertex < n_vertex; ++i_vertex)
   {
   vtx[i_vertex].x = (vtx[i_vertex].x - xstart) * xscale;
   vtx[i_vertex].y = (vtx[i_vertex].y - ystart) * yscale;
   vtx[i_vertex].z = (vtx[i_vertex].z - zstart) * zscale;
   }

/*
 ****************
 * Do the plot. *
 ****************
 */

printf ("\nPlot to the screen, file (1/0 1/0): ");
scanf ("%d %d", &screen, &pfile);
if (pfile == 1)
   {
   printf ("Enter Plot File Name : ");
   scanf ("%s", PlotFileName);
   }
else
   {
   strcpy(PlotFileName, "noplot");
   }


ssp_BeginPlot (screen, pfile, PlotFileName, "ps", "landscape");

/*
 * Set the 2D plotting-page origin.
 */
ssp_SetOrigin (xorig, yorig);

/*
 * Set the 3D origin and type of projection.
 */
ssp_SetOrigin3D (0.0, 0.0, 0.0);
if (iproj == -1)
   ssp_SetIsometricProjection ();
else
   ssp_SetOrthographicProjection (iproj);

/*
 * Sort the facets into z-depth order (minimum z first).
 * Apply the viewing transform first.
 */
ssp_Message ("Depth sort...");
for (i_vertex = 0; i_vertex < n_vertex; ++i_vertex)
   {
   ssp_ApplyProjectionMatrix( &(vtx[i_vertex]), &pt2 );
   z_vtx[i_vertex] = pt2.z;
   }

for (i_poly = 0; i_poly < n_poly; ++i_poly)
   {
   z_rank[i_poly]  = i_poly;

   A = polygn[i_poly].A;
   B = polygn[i_poly].B;
   C = polygn[i_poly].C;
   D = polygn[i_poly].D;

   if ( A < 0 || B < 0 )
      { /* do nothing */ ; }
   else if (C < 0)
      { /* Line only. */
      z_depth[i_poly] = 0.25 * (z_vtx[A] + z_vtx[B]);
      }
   else if (D < 0)
      { /* Triangular facet. */
      z_depth[i_poly] = 0.25 * (z_vtx[A] + z_vtx[B] + z_vtx[C]);
      }
   else
      { /* Quadrilateral facet. */
      z_depth[i_poly] = 0.25 * (z_vtx[A] + z_vtx[B] +
                                z_vtx[C] + z_vtx[D]);
      }
   }  /* for(ipoly...  */

indexx (n_poly, z_depth, z_rank);
ssp_Message ("Depth sort finished.");

/*
 * Write some blurb across the top.
 */
ssp_SetLineThickness (0.2);
ssp_PlotText (-15.0, 30.0+yplotsize, 2.5, 0.0, VertexFileName, 16);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotText (xx+5.0, yy, 2.5, 0.0, PolyFileName, 16);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotText (xx+5.0, yy, 2.5, 0.0, ValueFileName, 16);

ssp_PlotText (-15.0, 30.0+yplotsize-5.0, 2.5, 0.0, "x1 x2 dx", 8);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, xstart, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, xend, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, xtick, 2);

ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotText (xx+2.0, yy, 2.5, 0.0, "y1 y2 dy", 8);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, ystart, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, yend, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, ytick, 2);

ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotText (xx+2.0, yy, 2.5, 0.0, "z1 z2 dz", 8);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, zstart, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, zend, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, ztick, 2);

if ( do_contours || do_filled_colours || do_filled_grey )
   {
   ssp_PlotText (-15.0, 30.0+yplotsize-10.0, 2.5, 0.0, "v1 v2 dv", 8);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, vstart, 2);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, vend, 2);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, vtick, 2);
   /*
    * Put a colour map in here by writing the contour
    * values down the right-hand side of the page in the
    * appropriate colours.
    */
   if ( do_colours || do_filled_colours || do_filled_grey )
      {
      for (i = 15; i >= 0; --i)
         {
         level = vstart + i * (vend - vstart) / 15.0;
         if ( do_filled_grey )
            {
            hue_d = 0.0;
            bright = (level - vstart) / (vend - vstart);
            sat = 0.0;  /* zero saturation should force grey-scale */
	    ssp_SetHSBColour (hue_d, sat, bright);
            }
         else
            {
	    hue_d = ssp_MapToColour (level, vstart, vend);
            bright = 1.0;
            sat = 1.0;
            ssp_SetHSBColour (hue_d, sat, bright);
            }
         xx = 170.0;
         yy = 0.0 + 5.0 * i;
         ssp_PlotENumber (xx, yy, 2.5, 0.0, level, 3);
         }
      ssp_SetNoColour ();
      }
   }


/*
 * Plot Axes.
 */
ssp_SetLineThickness (0.35);
xaxislength = xrange * xscale;
yaxislength = yrange * yscale;
zaxislength = zrange * zscale;
ssp_axes3D (xaxislength, yaxislength, zaxislength,
            xtick*xscale, ytick*yscale, ztick*zscale,
            yaxislength/20.0);

/*
 * Remove 2D clipping for this 3D plotting.
 *
 * ssp_SetClipLimits (0.0, 0.0, xplotsize, yplotsize);
 * ssp_SetClipOn ();
 */

symbolsize = yplotsize / 70.0;   /* plotting symbols for the vertices */
isymb = CIRCLE_SYM;
ssp_SetFillGrey (0.5);               /* grey level for shading */

if (vertex_symbols == 1)
   {
   ssp_SetLineThickness (0.18);
   for (i_vertex = 0; i_vertex < n_vertex; ++i_vertex)
      {
      ssp_Move3D (vtx[i_vertex].x, vtx[i_vertex].y, vtx[i_vertex].z);
      ssp_Where (&xx, &yy, &sc, &ang);
      ssp_PlotSymbol (xx, yy, symbolsize, isymb);
      }
   }  /* if (do_symbols... */

if ( facet_boundaries || do_contours || do_filled_colours || fill_facets
     || do_filled_grey )
   {
   ssp_SetLineThickness (0.25);    /* Facet boundaries and contours */

   if ( fill_facets )              /* Fill in */
      ssp_SetFillGrey (0.5);       /* shaded */
   else
      ssp_SetFillGrey (0.0);       /* black */

   count = 0;
   climit = 100;
   for (i_poly = 0; i_poly < n_poly; ++i_poly)
      {
      ++count;
      if (count >= climit)
         {
         sprintf ( msg, "Polygon %d", i_poly );
         ssp_Message ( msg );
         count = 0;
         }

      ip = z_rank[i_poly];
      A = polygn[ip].A;
      B = polygn[ip].B;
      C = polygn[ip].C;
      D = polygn[ip].D;

      if ( A < 0 || B < 0 )
         { /* do nothing */ ; }

      else if (C < 0)
         { /* Line only. */
         if (facet_boundaries)
            {
            ssp_Move3D (vtx[A].x, vtx[A].y, vtx[A].z);
            ssp_Plot3D (vtx[B].x, vtx[B].y, vtx[B].z);
            }
         /* contouring and filling means nothing for a line... */
         } /* end of Line only */

      else if (D < 0)
         { /* Plot a triangular facet. */
         poly_vtx[0].x=vtx[A].x; poly_vtx[0].y=vtx[A].y; poly_vtx[0].z=vtx[A].z;
         poly_vtx[1].x=vtx[B].x; poly_vtx[1].y=vtx[B].y; poly_vtx[1].z=vtx[B].z;
         poly_vtx[2].x=vtx[C].x; poly_vtx[2].y=vtx[C].y; poly_vtx[2].z=vtx[C].z;
         if (facet_boundaries || fill_facets)
            {
            Intensity = reflected_light_intensity();
            ssp_SetFillGrey ( Intensity );
            ssp_Polygon3D (3, poly_vtx, fill_facets, 0, facet_boundaries);
            }
         if ( do_contours )
            {
            for (level = vstart; level <= vend+0.001*vtick; level += vtick)
               {
               delA = v_vtx[A] - level;
               delB = v_vtx[B] - level;
               delC = v_vtx[C] - level;
               missed = (delA < 0.0 && delB < 0.0 && delC < 0.0) ||
                        (delA > 0.0 && delB > 0.0 && delC > 0.0);
               if (do_colours && !missed)
                  {
                  hue_d = ssp_MapToColour (level, vstart, vend);
                  bright = 1.0;
                  sat = 1.0;
                  ssp_SetHSBColour (hue_d, sat, bright);
                  }
               if (!missed)
                  ssp_ContourTriangle3D ( &vtx[A], &vtx[B], &vtx[C],
                                          v_vtx[A], v_vtx[B], v_vtx[C],
                                          level );
               }
            if (do_colours) ssp_SetNoColour();
            }  /* if (do_contours)... */

         if ( do_filled_colours || do_filled_grey )
            {
            ssp_ColourFillTriangle3D ( &vtx[A], &vtx[B], &vtx[C],
                                       v_vtx[A], v_vtx[B], v_vtx[C],
                                       vstart, vend, do_filled_grey );
            }  /* if ( do_filled_colours or grey ... */

         } /* end of triangular facet */

      else
         { /* Plot a quadrilateral facet. */
         poly_vtx[0].x=vtx[A].x; poly_vtx[0].y=vtx[A].y; poly_vtx[0].z=vtx[A].z;
         poly_vtx[1].x=vtx[B].x; poly_vtx[1].y=vtx[B].y; poly_vtx[1].z=vtx[B].z;
         poly_vtx[2].x=vtx[C].x; poly_vtx[2].y=vtx[C].y; poly_vtx[2].z=vtx[C].z;
         poly_vtx[3].x=vtx[D].x; poly_vtx[3].y=vtx[D].y; poly_vtx[3].z=vtx[D].z;
         if (facet_boundaries || fill_facets)
            {
            Intensity = reflected_light_intensity();
            ssp_SetFillGrey ( Intensity );
            ssp_Polygon3D (4, poly_vtx, fill_facets, 0, facet_boundaries);
            }
         if ( do_contours )
            {
            for (level = vstart; level <= vend; level += vtick)
               {
               delA = v_vtx[A] - level;
               delB = v_vtx[B] - level;
               delC = v_vtx[C] - level;
               delD = v_vtx[D] - level;
               missed = (delA < 0.0 && delB < 0.0 && delC < 0.0 && delD < 0.0) ||
                        (delA > 0.0 && delB > 0.0 && delC > 0.0 && delD > 0.0);
               if (do_colours && !missed)
                  {
                  hue_d = ssp_MapToColour (level, vstart, vend);
                  bright = 1.0;
                  sat = 1.0;
		  ssp_SetHSBColour (hue_d, sat, bright);
                  }
               if (!missed)
                  ssp_ContourQuad3D ( &vtx[A], &vtx[B], &vtx[C], &vtx[D],
                                      v_vtx[A], v_vtx[B], v_vtx[C], v_vtx[D],
                                      level );
               }
            if (do_colours) ssp_SetNoColour();
            }  /* if (do_contours)... */

         if ( do_filled_colours || do_filled_grey )
            {
            ssp_ColourFillQuad3D ( &vtx[A], &vtx[B], &vtx[C], &vtx[D],
                                   v_vtx[A], v_vtx[B], v_vtx[C], v_vtx[D],
                                   vstart, vend, do_filled_grey );
            }  /* if ( do_filled_colours or grey ... */

         } /* end of quadrilateral facet */

      }  /* for (i_poly = ... */

   }  /* if (facet_boundaries || ... */


/*
 * Deal with the line-segments file.
 * This file is used if a separate set of line segments is to
 * be added to the plot after the contouring and/or shading has
 * been completed.
 */

if (do_lines == 1)
   {
   ssp_SetLineThickness (0.35);
   finished = 0;

   while (finished != 1)
      {
      /* Attempt to read the coordinates. */
      if ( fscanf(line_file, "%lf %lf %lf %lf %lf %lf",
                  &xx1, &yy1, &zz1, &xx2, &yy2, &zz2) == EOF ) finished = 1;
      /*
       * printf ("pnt1 = (%f, %f, %f), pnt2 = (%f, %f, %f)\n",
       *          xx1, yy1, zz1, xx2, yy2, zz2);
       */
      if (finished != 1)
         {
         /*
          * Scale the line segment and plot it.
          */
         xx1 = (xx1 - xstart) * xscale;
         yy1 = (yy1 - ystart) * yscale;
         zz1 = (zz1 - zstart) * zscale;

         xx2 = (xx2 - xstart) * xscale;
         yy2 = (yy2 - ystart) * yscale;
         zz2 = (zz2 - zstart) * zscale;

         ssp_Move3D (xx1, yy1, zz1);
         ssp_Plot3D (xx2, yy2, zz2);
         }  /* if (finished != 1... */
      }     /* while (!finished...  */

   fclose (line_file);
   }        /* if (do_lines...      */


/* ssp_SetClipOff (); */
ssp_EndPlot ();

return 0;
}  /* end of dispoly.c */

/*---------------------------------------------------------------*/

#if PROTO

double reflected_light_intensity ( void )

#else

double reflected_light_intensity ()

#endif
/*
 * Lighting model for computing the reflectd light intensity
 * from a facet.
 * The input to the routing (via global data) includes vtx[A]
 * vtx[B] and vtx[C].
 * Yeah, this is a bit daggy but the Topspeed compiler can't cope
 * with too many passed data structures -- we keep getting an error
 * message about exceeding the dynamic pool limit.
 */
{
struct vector_3D light_vector, sideBC, sideBA, unitN;
double Intensity, Int_a, k_a, Int_l, k_d, dd, konst, costheta;

/*
 * Set the parameters for the illumination model.
 */
light_vector.x = 0.0;    /* this must be a unit vector */
light_vector.y = 0.0;
light_vector.z = 1.0;

Int_a = 0.2;
k_a   = 0.5;
Int_l = 1.0;
k_d   = 0.80;   /* originally 0.5 but result was all quite dark */
                /* 0.85 seems to be too light */
konst = 1.0;

sideBA.x = vtx[A].x - vtx[B].x;
sideBA.y = vtx[A].y - vtx[B].y;
sideBA.z = vtx[A].z - vtx[B].z;

sideBC.x = vtx[C].x - vtx[B].x;
sideBC.y = vtx[C].y - vtx[B].y;
sideBC.z = vtx[C].z - vtx[B].z;

vector_prod_3D ( &sideBC, &sideBA, &unitN );
normalize_3D ( &unitN );

costheta = scalar_prod_3D ( &unitN, &light_vector );
if ( costheta < 0.0 ) costheta = -costheta;  /* treat both sides the same */

Intensity = Int_a * k_a + Int_l * k_d * costheta / konst;
if ( Intensity > 1.0 ) Intensity = 1.0;

return (Intensity);
}  /* end of reflected_light_intensity() */

/*---------------------------------------------------------------*/

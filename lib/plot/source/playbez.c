/* playbez.c
 * Play with the control points of Bezier curves.
 * -----------------------------------------------------------------
 *
 * Written by ... P. A. Jacobs
 * ----------     293 Lillian Avenue
 *                Salisbury, Qld 4107.
 *
 * Purpose...
 * -------
 * This interactive (menu driven) program can be used to fit
 * Bezier polylines to arbitrary sets of (x,y) data points.
 * Data for Bezier polylines can also be read and written in
 * a format suitable for the Multi-Block CNS solver.
 *
 * Version...
 * -------
 * 1.0  : 11-Jul-92 : First cut
 * 1.1  : 10-Nov-92
 *        added scaling to the data to improve the curve fitting
 *        performance
 *        The scaling is defined as...
 *        x_normalized = (x_actual - x_offset) * x_factor
 *        y_normalized = (y_actual - y_offset) * y_factor
 * 1.2  : 22-Jan-93
 *        * increased the size of the data arrays
 *        * added a penalty term to the objective function
 *          (for slope discontinuity)
 * 1.21 : 31-Jan-93
 *        * remove the penalty term as it is incompatible with
 *          abutting convex segments
 *        * add print for largest distance error
 *        * add a standard cubic-spline fitting routine which
 *          converts its results to Bezier segments.
 * 1.30 : 25-May-93
 *        * changed location_3D to point_3D
 * 1.31 : 15-May-94
 *        * increase data arrray dimension to 95
 */

/*-----------------------------------------------------------------*/

#include "cmath.h"
#include "compiler.h"

#include <stdio.h>
#include <math.h>
#if (STDLIBH)
#  include <stdlib.h>
#endif

#include "geom.h"
#include "ssp.h"
#include "bezier.h"
#include "playbez.h"

/*-----------------------------------------------------------------*/

/*
 * Global Data
 */
double xstart, xend, xtick;
double ystart, yend, ytick;
double minscale, xscale, yscale, xrange, yrange;
double xplotsize, yplotsize, xorig, yorig;
double xaxislength, yaxislength;
int    change, trueshape;
double xx, yy, sc, ang;
double g_max_dist;

#define  NCHAR  132

#if (TURBO_C || TOPSPEED_C)
/*
 * Limited memory.
 */
#define  NB_MAX      4
#define  ISEG_MAX    5

#define  ND_MAX      2
#define  DATA_MAX    95

#else
/*
 * Lots of memory.
 */
#define  NB_MAX      20
#define  ISEG_MAX    40

#define  ND_MAX      20
#define  DATA_MAX    200

#endif

struct bezier_3_poly_data  bp[NB_MAX];
int    valid_bezier[NB_MAX];

struct ssp_data_vector dv[ND_MAX];
int    valid_data[ND_MAX];

int    target_data_set;          /* index of the target data set   */
int    fitted_bezier;            /* index of bezier to fit to data */
double penalty_function ();      /* measure of fit                 */


/*-----------------------------------------------------------------*/

/*
 * The program proper.
 */

main ()
{
int    screen, pfile;
char   PlotFileName[32];
FILE   *pfp;
int    flag, jbp, jd, ix;

printf ("\n\n");
printf ("----------------------------\n");
printf ("Adjust Bezier Control Points\n");
printf ("----------------------------\n");
printf ("Rev. 1.21, 31-Jan-93\n");


printf ("\n");
printf ("\nAllocate data space for the Bezier curves.\n");
for (jbp = 0; jbp < NB_MAX; ++jbp)
   {
   valid_bezier[jbp] = 0;
   flag = alloc_bezier_3_poly (&bp[jbp], ISEG_MAX);
   if (flag != 0)
      {
      printf ("Bezier allocation failure, %d\n", jbp);
      exit (-1);
      }
   }

printf ("Allocate space for the data point arrays\n");
for (jd = 0; jd < ND_MAX; ++jd)
   {
   ssp_AllocateVector (&dv[jd], DATA_MAX);
   for (ix = 0; ix < DATA_MAX; ++ix)
      {
      dv[jd].x[ix] = 0.0;
      dv[jd].value[ix] = 0.0;
      }
   dv[jd].imin = 0;
   dv[jd].imax = DATA_MAX - 1;
   valid_data[jd] = 0;
   dv[jd].xscale = 1.0;
   dv[jd].vscale = 1.0;
   }

/*
 * Put in default values.
 * Fundamental plotting units are mm on an A4 page.
 */
xplotsize = 240.0;
yplotsize = 150.0;
xorig = 35.0;
yorig = 25.0;

xstart = 0.0;
xend = 1.0;
xtick = (xend - xstart) / 5.0;
ystart = 0.0;
yend = 1.0;
ytick = (yend - ystart) / 5.0;

xrange = xend - xstart;
xscale = xplotsize / xrange;
yrange = yend - ystart;
yscale = yplotsize / yrange;
trueshape = 1;
minscale = xscale;
if (yscale < minscale) minscale = yscale;
if (trueshape == 1)
   {
   xscale = minscale;
   yscale = minscale;
   }


/*
 * Start up the graphics system.
 */
printf ("\nPlot to the screen, file (1/0 1/0): ");
scanf ("%d %d", &screen, &pfile);
if (pfile == 1)
   {
   printf ("Enter Plot File Name : ");
   scanf ("%s", PlotFileName);
   if ((pfp = fopen(PlotFileName, "w")) == NULL)
      {
      printf ("Could not open PLOT file: %s\n", PlotFileName);
      exit (-1);
      }
   }


ssp_BeginPlot (screen, pfile, PlotFileName, "ps", "landscape");

ssp_SetClipLimits (0.0, 0.0, 300.0, 210.0);
ssp_SetClipOn ();

initialize_picture ();

root_menu ();

ssp_EndPlot ();
return 0;
}

/*-----------------------------------------------------------------*/

/*
 * Menu display routines.
 */

#if (PROTO)
   int root_menu (void)
#else
   int root_menu ()
#endif

{
char  line[132];

while (1)  /* infinite loop */
   {
   ssp_ClearText ();
   ssp_Message (
      "ROOT MENU: (F)ile (E)dit (P)lot fit(B)ezier fit(S)pline (Q)uit");
   scanf ("%s", line);

   if (line[0] == 'F' || line[0] == 'f')
      {
      file_menu ();
      }
   else if (line[0] == 'E' || line[0] == 'e')
      {
      edit_menu ();
      }
   else if (line[0] == 'P' || line[0] == 'p')
      {
      plot_menu ();
      }
   else if (line[0] == 'B' || line[0] == 'b')
      {
      fit_bezier ();
      }
   else if (line[0] == 'S' || line[0] == 's')
      {
      fit_spline ();
      }
   else if (line[0] == 'Q' || line[0] == 'q')
      {
      return (0);
      }
   else
      {
      /* No valid response */
      }
   }

return (0);  /* should never reach this point */
}  /* end of root_menu() */


#if (PROTO)
   int file_menu (void)
#else
   int file_menu ()
#endif
{
char line[132];

ssp_ClearText ();
ssp_Message
   ("FILE: (1)save_bezier (2)read_bezier (3)save_data (4)read_data (Q)uit");
scanf ("%s", line);

if (line[0] == '1')
   {
   save_bezier ();
   }
else if (line[0] == '2')
   {
   read_bezier ();
   }
else if (line[0] == '3')
   {
   save_data_set ();
   }
else if (line[0] == '4')
   {
   read_data_set ();
   }
else if (line[0] == 'Q' || line[0] == 'q')
   {
   return (0);
   }
else
   {
   /* No valid response */
   }

/*
 * Drop back into the ROOT menu.
 */

return 0;
}  /* end of file_menu() */


#if (PROTO)
   int edit_menu (void)
#else
   int edit_menu ()
#endif
{
char line[132];

ssp_ClearText ();
ssp_Message ("EDIT: (B)ezier_polyline; (D)ata_set; (Q)uit");
scanf ("%s", line);

if (line[0] == 'B' || line[0] == 'b')
   {
   edit_bezier ();
   }
else if (line[0] == 'D' || line[0] == 'd')
   {
   edit_data_set ();
   }
else if (line[0] == 'Q' || line[0] == 'q')
   {
   return (0);
   }
else
   {
   /* No valid response */
   }

/*
 * Drop back into the ROOT menu.
 */

return 0;
}  /* end of edit_menu() */


#if (PROTO)
   int plot_menu (void)
#else
   int plot_menu ()
#endif
{
char line[132];

ssp_ClearText ();
ssp_Message ("PLOT: (R)edraw_all; (C)hange_view; (Q)uit");
scanf ("%s", line);

if (line[0] == 'R' || line[0] == 'r')
   {
   redraw_all ();
   }
else if (line[0] == 'C' || line[0] == 'c')
   {
   change_view ();
   }
else if (line[0] == 'Q' || line[0] == 'q')
   {
   return (0);
   }
else
   {
   /* No valid response */
   }

/*
 * Drop back into the ROOT menu.
 */

return 0;
}  /* end of plot_menu() */

/*-----------------------------------------------------------------*/

/*
 * Data editing routines.
 */


#if (PROTO)
   int edit_bezier (void)
#else
   int edit_bezier ()
#endif
{
int    jbp, iseg, jcp;
struct point_3D  point;

/*
 * Assume that we have come from the graphics environment and
 * wish to use test only here.
 */
ssp_TextOn ();
printf ("Edit a Bezier polyline.\n\n");

for (jbp = 0; jbp < NB_MAX; ++jbp)
   {
   printf ("Polyline [%d] is ", jbp);
   if (valid_bezier[jbp] == 1) 
      printf ("valid\n");
   else
      printf ("NOT valid\n");
   }

printf ("\nPolyline index (0...%d): ", NB_MAX - 1);
scanf ("%d", &jbp);
if (jbp < 0 || jbp >= NB_MAX) goto BailOut;

printf ("Polyline [%d], number of Bezier segments = %d\n",
	jbp, bp[jbp].n);

for (iseg = 0; iseg < bp[jbp].n; ++iseg)
   {
   printf ("seg[%d] (%f %f) (%f %f)\n", iseg,
      bp[jbp].b3_seg[iseg].B[0].x, bp[jbp].b3_seg[iseg].B[0].y,
      bp[jbp].b3_seg[iseg].B[1].x, bp[jbp].b3_seg[iseg].B[1].y );
   printf ("        (%f %f) (%f %f)\n",
      bp[jbp].b3_seg[iseg].B[2].x, bp[jbp].b3_seg[iseg].B[2].y,
      bp[jbp].b3_seg[iseg].B[3].x, bp[jbp].b3_seg[iseg].B[3].y );
   }


printf
   ("\nNow, edit the control points... (remember iseg < 0 exits)\n");

Ask_for_data:
printf ("\n\nEnter [iseg] [jcp] (x y): ");
scanf ("%d %d %lf %lf", &iseg, &jcp, &(point.x), &(point.y) );
if (iseg >= 0)
   {
   if (iseg >= 0 && iseg < bp[jbp].n
       && jcp >= 0 && jcp <= 3 
       && valid_bezier[jbp] == 1 )
      {
      bp[jbp].b3_seg[iseg].B[jcp].x = point.x;
      bp[jbp].b3_seg[iseg].B[jcp].y = point.y;
      goto Ask_for_data;
      }
   }

ssp_WaitPrompt ();

BailOut:
/*
 * Return to the previous menu and the graphic environment.
 */
ssp_TextOff ();
return (0);
}


#if (PROTO)
   int edit_data_set (void)
#else
   int edit_data_set ()
#endif
{


return (0);
}


/*-----------------------------------------------------------------*/

/*
 * Input-Output routines.
 */


#if (PROTO)
   int save_data_set (void)
#else
   int save_data_set ()
#endif
{

return (0);
}


#if (PROTO)
   int read_data_set (void)
#else
   int read_data_set ()
#endif
{
char   line[NCHAR], filename[40];
FILE   *datafile;
int    ix, jd, nd;
struct point_3D  point;
double x_factor, y_factor, x_offset, y_offset;

/*
 * Assume that we have come from the graphics environment and
 * wish to use test only here.
 */
ssp_TextOn ();
printf ("Read a Data Set.\n\n");

for (jd = 0; jd < ND_MAX; ++jd)
   {
   printf ("Data Set [%d] is ", jd);
   if (valid_data[jd] == 1) 
      printf ("valid\n");
   else
      printf ("NOT valid\n");
   }

printf ("\nData Set index (0...%d): ", ND_MAX - 1);
scanf ("%d", &jd);
if (jd < 0 || jd >= ND_MAX) goto BailOut;

printf ("File Name: ");
scanf ("%s", filename);
datafile = NULL;
datafile = fopen (filename, "r");
if (datafile == NULL) goto BailOut;

printf ("Scale factors and offsets (xf, xo, yf, yo) : ");
scanf ("%lf %lf %lf %lf", &x_factor, &x_offset, &y_factor, &y_offset);

fgets (line, NCHAR, datafile);
sscanf (line, "%d", &nd);
if (nd > DATA_MAX) nd = DATA_MAX;
dv[jd].imin = 0;
dv[jd].imax = nd - 1;
#if (1)
printf ("Data Set [%d], n = %d\n", jd, (dv[jd].imax + 1) );
#endif

for (ix = dv[jd].imin; ix <= dv[jd].imax; ++ix)
   {
   fgets (line, NCHAR, datafile);
   sscanf (line, "%lf %lf", &(point.x), &(point.y) );
#  if (1)
   printf ("point [%d]: (%f %f)\n", ix, point.x, point.y);
#  endif
   dv[jd].x[ix] = x_factor * (point.x - x_offset);
   dv[jd].value[ix] = y_factor * (point.y - y_offset);
   }

fclose (datafile);

valid_data[jd] = 1;

ssp_WaitPrompt ();

BailOut:
/*
 * Return to the previous menu and the graphic environment.
 */
ssp_TextOff ();

return (0);
}


#if (PROTO)
   int save_bezier (void)
#else
   int save_bezier ()
#endif
/*
 * Purpose ...
 * -------
 * Write the data for a specified Bezier polyline out
 * to a file.
 *
 */
{
char   filename[40];
FILE   *bezierfile;
int    jbp, iseg;
double x_factor, y_factor, x_offset, y_offset;

/*
 * Assume that we have come from the graphics environment and
 * wish to use test only here.
 */
ssp_TextOn ();
printf ("Write a Bezier polyline.\n\n");

for (jbp = 0; jbp < NB_MAX; ++jbp)
   {
   printf ("Polyline [%d] is ", jbp);
   if (valid_bezier[jbp] == 1)
      printf ("valid\n");
   else
      printf ("NOT valid\n");
   }

printf ("\nPolyline index (0...%d): ", NB_MAX - 1);
scanf ("%d", &jbp);
if (jbp < 0 || jbp >= NB_MAX) goto BailOut;

printf ("File Name: ");
scanf ("%s", filename);
bezierfile = NULL;
bezierfile = fopen (filename, "w");
if (bezierfile == NULL) goto BailOut;

printf ("Scale factors and offsets (xf, xo, yf, yo) : ");
scanf ("%lf %lf %lf %lf", &x_factor, &x_offset, &y_factor, &y_offset);

fprintf (bezierfile, "%d  number of Bezier segments\n", bp[jbp].n);

for (iseg = 0; iseg < bp[jbp].n; ++iseg)
   {
   fprintf (bezierfile, "%f %f %f %f %f %f %f %f\n",
      bp[jbp].b3_seg[iseg].B[0].x / x_factor + x_offset,
      bp[jbp].b3_seg[iseg].B[0].y / y_factor + y_offset,
      bp[jbp].b3_seg[iseg].B[1].x / x_factor + x_offset,
      bp[jbp].b3_seg[iseg].B[1].y / y_factor + y_offset,
      bp[jbp].b3_seg[iseg].B[2].x / x_factor + x_offset,
      bp[jbp].b3_seg[iseg].B[2].y / y_factor + y_offset,
      bp[jbp].b3_seg[iseg].B[3].x / x_factor + x_offset,
      bp[jbp].b3_seg[iseg].B[3].y / y_factor + y_offset );
   }

fclose (bezierfile);

ssp_WaitPrompt ();

BailOut:
/*
 * Return to the previous menu and the graphic environment.
 */
ssp_TextOff ();
return (0);
}  /* end of save_bezier() */


#if (PROTO)
   int read_bezier (void)
#else
   int read_bezier ()
#endif
/*
 * Purpose ...
 * -------
 * Read the data for a specified Bezier polyline.
 *
 */
{
char   line[NCHAR], filename[40];
FILE   *bezierfile;
int    jbp, iseg;
struct point_3D  loc0, loc1, loc2, loc3;
double x_factor, y_factor, x_offset, y_offset;

/*
 * Assume that we have come from the graphics environment and
 * wish to use test only here.
 */
ssp_TextOn ();
printf ("Read a Bezier polyline.\n\n");

for (jbp = 0; jbp < NB_MAX; ++jbp)
   {
   printf ("Polyline [%d] is ", jbp);
   if (valid_bezier[jbp] == 1) 
      printf ("valid\n");
   else
      printf ("NOT valid\n");
   }

printf ("\nPolyline index (0...%d): ", NB_MAX - 1);
scanf ("%d", &jbp);
if (jbp < 0 || jbp >= NB_MAX) goto BailOut;

printf ("File Name: ");
scanf ("%s", filename);
bezierfile = NULL;
bezierfile = fopen (filename, "r");
if (bezierfile == NULL) goto BailOut;

printf ("Scale factors and offsets (xf, xo, yf, yo) : ");
scanf ("%lf %lf %lf %lf", &x_factor, &x_offset, &y_factor, &y_offset);

fgets (line, NCHAR, bezierfile);
sscanf (line, "%d", &(bp[jbp].n));
#if (1)
printf ("Bezier polyline [%d], n = %d\n", jbp, bp[jbp].n);
#endif

for (iseg = 0; iseg < bp[jbp].n; ++iseg)
   {
   fgets (line, NCHAR, bezierfile);
   sscanf (line, "%lf %lf %lf %lf %lf %lf %lf %lf",
      &(loc0.x), &(loc0.y), &(loc1.x), &(loc1.y),
      &(loc2.x), &(loc2.y), &(loc3.x), &(loc3.y) );
   loc0.z = 0.0; loc1.z = 0.0; loc2.z = 0.0; loc3.z = 0.0;
#  if (1)
   printf ("seg[%d]: (%f %f) (%f %f) (%f %f) (%f %f)\n", iseg,
       loc0.x, loc0.y, loc1.x, loc1.y,
       loc2.x, loc2.y, loc3.x, loc3.y );
#  endif
   loc0.x = x_factor * (loc0.x - x_offset);
   loc0.y = y_factor * (loc0.y - y_offset);
   loc0.z = 0.0;
   loc1.x = x_factor * (loc1.x - x_offset);
   loc1.y = y_factor * (loc1.y - y_offset);
   loc1.z = 0.0;
   loc2.x = x_factor * (loc2.x - x_offset);
   loc2.y = y_factor * (loc2.y - y_offset);
   loc2.z = 0.0;
   loc3.x = x_factor * (loc3.x - x_offset);
   loc3.y = y_factor * (loc3.y - y_offset);
   loc3.z = 0.0;
   segment_bezier_3_poly (&bp[jbp], &loc0, &loc1, &loc2, &loc3, iseg);
   }

fclose (bezierfile);

#if (0)
printf ("Bezier Polyline %d: normalize\n", jbp);
#endif
normalize_bezier_3_poly (&bp[jbp]);

valid_bezier[jbp] = 1;

ssp_WaitPrompt ();

BailOut:
/*
 * Return to the previous menu and the graphic environment.
 */
ssp_TextOff ();
return (0);
}  /* end of read_bezier() */


/*-----------------------------------------------------------------*/

/*
 * Plotting routines.
 */

#if (PROTO)
   initialize_picture (void)
#else
   initialize_picture ()
#endif
/*
 * Purpose...
 * -------
 * Set up the scales, axes and labels for the plot display.
 */
{
ssp_SetOrigin (xorig, yorig);
ssp_PlotText (50.0, yplotsize+15.0, 8.0, 0.0, "Bezier Data", 10);

xaxislength = xrange * xscale;
if (xplotsize > xaxislength) xaxislength = xplotsize;
ssp_PlotXAxis (0.0, 0.0, xaxislength, 0.0, xscale,
   xstart, xtick, "x", -1);

yaxislength = yrange * yscale;
if (yplotsize > yaxislength) yaxislength = yplotsize;
ssp_PlotYAxis (0.0, 0.0, yaxislength, 0.0, yscale, 
   ystart, ytick, "y", -1);

return (0);
}


#if (PROTO)
   redraw_all (void)
#else
   redraw_all ()
#endif
/*
 * Purpose ...
 * -------
 * Redraw all of the valid Bezier polylines and the data points.
 *
 */
{
int    jbp, iseg, np, j, ix, jd;
double t, dt;
struct point_3D point;

/*
 * Assume that we are in graphics mode.  Start with a fresh screen.
 */
ssp_StartNewPlot ();
initialize_picture ();


/*
 * Bezier Polylines...
 */

for (jbp = 0; jbp < NB_MAX; ++jbp)
   {
   if (valid_bezier[jbp] == 1)
      {
      /*
       * We have some valid data, go ahead and plot the polyline.
       * Plot centred symbols for the control points. 
       * Remember, 4 points per segment.
       */
      for (iseg = 0; iseg < bp[jbp].n; ++iseg)
	 {
	 /* Control point 0. */
	 xx = xscale * (bp[jbp].b3_seg[iseg].B[0].x - xstart);
	 yy = yscale * (bp[jbp].b3_seg[iseg].B[0].y - ystart);
	 ssp_PlotSymbol (xx, yy, 5.0, CIRCLE_SYM);
	 /* Control point 1. */
	 xx = xscale * (bp[jbp].b3_seg[iseg].B[1].x - xstart);
	 yy = yscale * (bp[jbp].b3_seg[iseg].B[1].y - ystart);
	 ssp_PlotSymbol (xx, yy, 5.0, CIRCLE_SYM);
	 /* Control point 2. */
	 xx = xscale * (bp[jbp].b3_seg[iseg].B[2].x - xstart);
	 yy = yscale * (bp[jbp].b3_seg[iseg].B[2].y - ystart);
	 ssp_PlotSymbol (xx, yy, 5.0, CIRCLE_SYM);
	 /* Control point 3. */
	 xx = xscale * (bp[jbp].b3_seg[iseg].B[3].x - xstart);
	 yy = yscale * (bp[jbp].b3_seg[iseg].B[3].y - ystart);
	 ssp_PlotSymbol (xx, yy, 5.0, CIRCLE_SYM);
	 }  /* for (iseg = 0...     */

      /*
       * Plot line segments along the Bezier polyline.
       */
      np = 50;
      dt = 1.0 / np;
      t = 0.0;
      eval_bezier_3_poly ( &bp[jbp], t, &point );
      ssp_Move (xscale * (point.x - xstart), 
		yscale * (point.y - ystart) );
      for (j = 1; j <= np; ++j)
	 {
	 t = dt * j;
	 eval_bezier_3_poly ( &bp[jbp], t, &point);
	 ssp_Plot (xscale * (point.x - xstart), 
		   yscale * (point.y - ystart) );
	 }

      }     /* if (valid_bezier ... */
   }        /* for (jbp = 0...      */


/*
 * Data Sets ...
 */

for (jd = 0; jd < ND_MAX; ++jd)
   {
   if (valid_data[jd] == 1)
      {
      /*
       * We have some valid data, go ahead and plot the points.
       */
      for (ix = dv[jd].imin; ix <= dv[jd].imax; ++ix)
	 {
	 xx = xscale * (dv[jd].x[ix] - xstart);
	 yy = yscale * (dv[jd].value[ix] - ystart);
	 ssp_PlotSymbol (xx, yy, 5.0, PLUS_SYM);
	 }

      }     /* if (valid_data ... */
   }        /* for (jd = 0...      */

return (0);
}


#if (PROTO)
   change_view (void)
#else
   change_view ()
#endif
/*
 * Purpose...
 * -------
 * Get a new set of window limits.
 * Reset the plotting scales.
 *
 */
{   /* begin change_view() */

/*
 * Assume that we have come from the graphics environment and
 * wish to use test only here.
 */
ssp_TextOn ();
printf ("Change the view.\n\n");

Ask_for_data:
printf ("Present window : \n");
printf ("xstart = %f, xend = %f, ystart = %f, yend = %f\n",
	xstart, xend, ystart, yend);
printf ("New window : \n");
printf ("xstart  xend  ystart  yend : ");
scanf ("%lf %lf %lf %lf", &xstart, &xend, &ystart, &yend);

if (xstart >= xend || ystart >= yend)
   {
   printf ("Invalid window -- try again.\n");
   goto Ask_for_data;
   }

xtick = (xend - xstart) / 5.0;
ytick = (yend - ystart) / 5.0;

xrange = xend - xstart;
xscale = xplotsize / xrange;
yrange = yend - ystart;
yscale = yplotsize / yrange;
trueshape = 1;
minscale = xscale;
if (yscale < minscale) minscale = yscale;
if (trueshape == 1)
   {
   xscale = minscale;
   yscale = minscale;
   }


ssp_WaitPrompt ();

BailOut:
/*
 * Return to the previous menu and the graphic environment.
 */
ssp_TextOff ();
return (0);
}  /* end of change_view() */


/*-------------------------------------------------------------*/

/*
 * Best-fit routines.
 */

#if (PROTO)
   int fit_spline (void)
#else
   int fit_spline ()
#endif
/*
 * Purpose...
 * -------
 * Fit a cubic spline to a set of data points and convert the
 * result to a Bezier Polyline.
 *
 */
{
int    jbp, jd, iseg;
double xd[DATA_MAX], yd[DATA_MAX], wd[DATA_MAX];
double xs[ISEG_MAX+1], ys[ISEG_MAX+1], s1, s2, sums, distance;
double dx;
int    iflag, last;
double bs[ISEG_MAX+1], cs[ISEG_MAX+1], ds[ISEG_MAX+1];
int    ndata, nseg, ix, i;

/*
 * Assume that we have come from the graphics environment and
 * wish to use test only here.
 */
ssp_TextOn ();
printf ("\nFit a Cubic Spline to a data set.\n\n");

for (jbp = 0; jbp < NB_MAX; ++jbp)
   {
   printf ("Polyline [%d] is ", jbp);
   if (valid_bezier[jbp] == 1)
      printf ("valid\n");
   else
      printf ("NOT valid\n");
   }

printf ("\nPolyline index (0...%d): ", NB_MAX - 1);
scanf ("%d", &jbp);
if (jbp < 0 || jbp >= NB_MAX) goto BailOut;

printf ("\n");
for (jd = 0; jd < ND_MAX; ++jd)
   {
   printf ("Data Set [%d] is ", jd);
   if (valid_data[jd] == 1)
      printf ("valid\n");
   else
      printf ("NOT valid\n");
   }

printf ("\nData Set index (0...%d): ", ND_MAX - 1);
scanf ("%d", &jd);
if (jd < 0 || jd >= ND_MAX) goto BailOut;

if (valid_data[jd] == 1 && valid_bezier[jbp] == 1)
   {
   /*
    * Go ahead and try to fit the line.
    */

   /* number of cubic segments, (nseg+1) knots */
   nseg = bp[jbp].n;

   /* The knots are the end-points of the Bezier segments. */
   for (iseg = 0; iseg < nseg; ++iseg)
      {
      xs[iseg] = bp[jbp].b3_seg[iseg].B[0].x;
      ys[iseg] = bp[jbp].b3_seg[iseg].B[0].y;
      }
   xs[nseg] = bp[jbp].b3_seg[nseg-1].B[3].x;
   ys[nseg] = bp[jbp].b3_seg[nseg-1].B[3].y;

   ndata = dv[jd].imax - dv[jd].imin + 1;
   i = 0;
   for (ix = dv[jd].imin; ix <= dv[jd].imax; ++ix)
      {
      xd[i] = dv[jd].x[ix];
      yd[i] = dv[jd].value[ix];
      wd[i] = 1.0;
      ++i;
      }

   /*
    * This initial guess for the slopes assumes that
    * the x-coordinates are in ascending order.
    */
   s1 = (yd[1] - yd[0]) / (xd[1] - xd[0]);
   s2 = (yd[ndata-1] - yd[ndata-2]) / (xd[ndata-1] - xd[ndata-2]);

   sums = 1.0e-9;

   printf ("Call the spline-fitting routine...\n");
   fitspl (ndata, xd, yd, wd, nseg+1, xs, ys, &s1, &s2, &sums, &iflag);

   printf ("%s\n", cmathmsg(FITSPL_C, iflag) );
   printf ("Weighted sum of residuals = %e\n", sums);

   printf ("End slopes: s1 = %e, s2 = %e\n", s1, s2);
   printf ("Spline knots...\n");
   for (iseg = 0; iseg <= nseg; ++iseg)
      {
      printf ("x,y[%d] = %e %e\n", iseg, xs[iseg], ys[iseg]);
      }

   /*
    * Get the cubic-spline coefficients and compute
    * largest error.
    */
   spline (nseg+1, 1, 1, s1, s2, xs, ys, bs, cs, ds, &iflag);
   g_max_dist = 0.0;
   last = 1;
   for (i = 0; i < ndata; ++i)
      {
      distance =
         fabs( yd[i] - seval(nseg+1, xd[i], xs, ys, bs, cs, ds, &last) );
      if (distance > g_max_dist) g_max_dist = distance;
      }
   printf ("Maximum distance error = %e\n", g_max_dist);

   /* Convert back to Bezier Control Points. */
   for (iseg = 0; iseg < nseg; ++iseg)
      {
      dx = xs[iseg+1] - xs[iseg];

      bp[jbp].b3_seg[iseg].B[0].x = xs[iseg];
      bp[jbp].b3_seg[iseg].B[0].y = ys[iseg];

      bp[jbp].b3_seg[iseg].B[1].x = xs[iseg] + dx / 3.0;
      bp[jbp].b3_seg[iseg].B[1].y = dx * bs[iseg] / 3.0 + ys[iseg];

      bp[jbp].b3_seg[iseg].B[2].x = xs[iseg] + dx * 2.0 / 3.0;
      bp[jbp].b3_seg[iseg].B[2].y =
         dx / 3.0 * (dx * cs[iseg] + 2.0 * bs[iseg]) + ys[iseg];

      bp[jbp].b3_seg[iseg].B[3].x = xs[iseg+1];
      bp[jbp].b3_seg[iseg].B[3].y = ys[iseg+1];
      }

   }  /* end of fitting line */

ssp_WaitPrompt ();

BailOut:
/*
 * Return to the previous menu and the graphic environment.
 */
ssp_TextOff ();
return (0);
}  /* end of fit_spline() */


#if (PROTO)
   int fit_bezier (void)
#else
   int fit_bezier ()
#endif
/*
 * Purpose...
 * -------
 * Fit a Bezier Polyline to a set of data points.
 *
 */
{
int    jbp, jd, iseg, jp;
int    npar, convge, nfe, maxfe, numres, flag;
#define   NPAR_MAX   2 * ISEG_MAX
double par[NPAR_MAX], dpar[NPAR_MAX];
double reqmin, reltol, abstol, fmin;
double dx, dy, length;

/*
 * Assume that we have come from the graphics environment and
 * wish to use test only here.
 */
ssp_TextOn ();
printf ("\nFit a Bezier polyline to a data set.\n\n");

for (jbp = 0; jbp < NB_MAX; ++jbp)
   {
   printf ("Polyline [%d] is ", jbp);
   if (valid_bezier[jbp] == 1)
      printf ("valid\n");
   else
      printf ("NOT valid\n");
   }

printf ("\nPolyline index (0...%d): ", NB_MAX - 1);
scanf ("%d", &jbp);
if (jbp < 0 || jbp >= NB_MAX) goto BailOut;

printf ("\n");
for (jd = 0; jd < ND_MAX; ++jd)
   {
   printf ("Data Set [%d] is ", jd);
   if (valid_data[jd] == 1)
      printf ("valid\n");
   else
      printf ("NOT valid\n");
   }

printf ("\nData Set index (0...%d): ", ND_MAX - 1);
scanf ("%d", &jd);
if (jd < 0 || jd >= ND_MAX) goto BailOut;

if (valid_data[jd] == 1 && valid_bezier[jbp] == 1)
   {
   /*
    * Go ahead and try to fit the line.
    */
   target_data_set = jd;
   fitted_bezier   = jbp;

   npar = 2 * bp[jbp].n;
   for (iseg = 0; iseg < bp[jbp].n; ++iseg)
      {
      jp = 2 * iseg;

      /* Initial segments are straight lines. */
      par[jp] = 0.0;
      par[jp+1] = 0.0;

      /* Perturbations are a fixed fraction of the segment lengths */
      dx = bp[jbp].b3_seg[iseg].B[3].x - bp[jbp].b3_seg[iseg].B[0].x;
      dy = bp[jbp].b3_seg[iseg].B[3].y - bp[jbp].b3_seg[iseg].B[0].y;
      length = sqrt ( dx * dx + dy * dy );
      dpar[jp] = 0.5 * length;
      dpar[jp+1] = 0.5 * length;
      }
   reqmin = 1.0e-6;
   convge = 5;
   maxfe  = 800;
   reltol = 0.0;
   abstol = 0.0;

   printf ("\nCalling function minimizer...\n");
   nelmin (penalty_function, npar, par, &fmin, reqmin, dpar,
	   convge, &nfe, maxfe, &numres, &flag, reltol, abstol);
   printf ("%s\n", cmathmsg(NELMIN_C, flag) );
   printf ("nfe = %d, numres = %d\n", nfe, numres);
   printf ("max error (dist) = %e\n", g_max_dist);

   }  /* end of fitting line */

ssp_WaitPrompt ();

BailOut:
/*
 * Return to the previous menu and the graphic environment.
 */
ssp_TextOff ();
return (0);
}  /* end of fit_bezier() */


#if (PROTO)
   double penalty_function ( int npar, double par[] )
#else
   double penalty_function (npar, par)
   int    npar;
   double par[];
#endif
/*
 * Purpose...
 * -------
 * Evaluate the measure of fit for the selected Bezier polyline
 * to the target data.
 * The parameters are the distances of control points B[1] and B[2]
 * normal to the line joining B[0] and B[3].
 * Control (end) points B[0] and B[3] are assumed given for each
 * segment in the Bezier polyline while B[1] and B[2] are located
 * on normals positioned 1/3 and 2/3 of the total distance from B[0].
 *
 * 22-Jan-93 -- add a penalty term for slope discontinuity
 *
 */
{  /* begin penalty_function() */
int    jbp, iseg, jd, ns, js, ix;
struct point_3D B0, B1, B2, B3, B13, B23;
struct point_3D dxyzdt;
double left_slope, right_slope;
double dxdt_0[ISEG_MAX], dydt_0[ISEG_MAX];
double dxdt_1[ISEG_MAX], dydt_1[ISEG_MAX];
/*
 * NS_MAX is the number of sample points that will be along each
 * Bezier polyline and then compared with the specified data points.
 * NS_MAX needs to be set large enough to get points close to each
 * data point (but not too large).
 */
#define  NS_MAX  201
struct point_3D sample[NS_MAX];
double lambda_13, lambda_23;
double dy, dx, length, cosine, sine;
double t, dt, sum, distance, min_dist, max_dist;

/* Debug...
gotoxy (1,1);
printf ("penalty function...\n");
*/

jbp = fitted_bezier;
jd  = target_data_set;

if (npar != 2 * bp[jbp].n)
   {
   printf ("Wrong number of parameters: npar = %d, iseg = %d\n",
	    npar, bp[jbp].n);
   return (0.0);
   }

/*
 * Set up the Bezier polyline, one segment at a time.
 */
for (iseg = 0; iseg < bp[jbp].n; ++iseg)
   {
   /* Unpack the parameter vector. */
   lambda_13 = par[iseg * 2];
   lambda_23 = par[iseg * 2 + 1];

   /* Debug...
   printf ("iseg[%d], l_13 = %f, l_23 = %f\n",
	   iseg, lambda_13, lambda_23);
   */

   /* The present control points for this segment. */
   B0.x = bp[jbp].b3_seg[iseg].B[0].x;
   B0.y = bp[jbp].b3_seg[iseg].B[0].y;
   /*
   B1.x = bp[jbp].b3_seg[iseg].B[1].x;
   B1.y = bp[jbp].b3_seg[iseg].B[1].y;
   B2.x = bp[jbp].b3_seg[iseg].B[2].x;
   B2.y = bp[jbp].b3_seg[iseg].B[2].y;
   */
   B3.x = bp[jbp].b3_seg[iseg].B[3].x;
   B3.y = bp[jbp].b3_seg[iseg].B[3].y;

   /* The 1/3 and 2/3 points on the B[0] -- B[3] line. */
   dx = B3.x - B0.x;
   dy = B3.y - B0.y;
   length = sqrt( dx * dx + dy * dy );
   sine = dy / length;
   cosine = dx / length;
   B13.x = B0.x + dx / 3.0;
   B13.y = B0.y + dy / 3.0;
   B23.x = B0.x + 2.0 * dx / 3.0;
   B23.y = B0.y + 2.0 * dy / 3.0;

   /* The new control points. */
   B1.x = B13.x - sine * lambda_13;
   B1.y = B13.y + cosine * lambda_13;
   B2.x = B23.x - sine * lambda_23;
   B2.y = B23.y + cosine * lambda_23;

   /* Update the Bezier line. */
   bp[jbp].b3_seg[iseg].B[1].x = B1.x;
   bp[jbp].b3_seg[iseg].B[1].y = B1.y;
   bp[jbp].b3_seg[iseg].B[2].x = B2.x;
   bp[jbp].b3_seg[iseg].B[2].y = B2.y;

   /* Compute the derivatives at each end of the Bezier segment */
   deriv_bezier_3 ( &(bp[jbp].b3_seg[iseg]), 0.0, &dxyzdt);
   dxdt_0[iseg] = dxyzdt.x;
   dydt_0[iseg] = dxyzdt.y;
   deriv_bezier_3 ( &(bp[jbp].b3_seg[iseg]), 1.0, &dxyzdt);
   dxdt_1[iseg] = dxyzdt.x;
   dydt_1[iseg] = dxyzdt.y;

   }  /* for (iseg = 0;... */

/*
 * Normalize the new line and then sample it.
 */
normalize_bezier_3_poly (&bp[jbp]);
ns = NS_MAX;
dt = 1.0 / (ns - 1);
for (js = 0; js < ns; ++js)
   {
   t = js * dt;
   eval_bezier_3_poly ( &bp[jbp], t, &sample[js] );
   }

/*
 * Evaluate the measure of fit by first evaluating the error
 * with respect to the specified data...
 */
sum = 0.0;
max_dist = 0.0;  /* the largest error (distance) */
for (ix = dv[jd].imin; ix <= dv[jd].imax; ++ix)
   {
   /* Find the distance of closest approach for each data point. */
   min_dist = 1.0e12;  /* VERY big */
   for (js = 0; js < ns; ++js)
      {
      dx = dv[jd].x[ix] - sample[js].x;
      dy = dv[jd].value[ix] - sample[js].y;
      distance = sqrt( dx * dx + dy * dy );
      if (distance < min_dist) min_dist = distance;
      }

   /* Add it to the overall measure of fit. */
   sum += min_dist;

   /* Record the largest error (distance) */
   if (min_dist > max_dist) max_dist = min_dist;
   }  /* for (ix = ... */

/*
 * Add an extra penalty for slope mismatch between adjacent
 * Bezier segments...
 */
#if (0)
/* Cut out this section. */
for (iseg = 0; iseg < bp[jbp].n - 1; ++iseg)
   {
   if ( fabs(dxdt_1[iseg]) >= fabs(dydt_1[iseg]) )
      {
      /* Usual slope ... */
      left_slope = dydt_1[iseg] / dxdt_1[iseg];
      right_slope = dydt_0[iseg+1] / dxdt_0[iseg+1];
      }
   else
      {
      /* Inverse slope ... */
      left_slope = dxdt_1[iseg] / dydt_1[iseg];
      right_slope = dxdt_0[iseg+1] / dydt_0[iseg+1];
      }
   /* Add the penalty term multiplied by a weight... */
   sum += 1.0 * fabs(left_slope - right_slope);
   }  /* for (iseg = 0...  */
#endif


/* Debug...
printf ("sum = %f\n", sum);
*/

g_max_dist = max_dist;  /* copy to the global variable */

return (sum);
}  /* end of penalty_function() */



/*----------------------- end of playbez.c ------------------------*/
/* particle.c
 * Display particle positions.
 *------------------------------------------------------------
 *
 * Purpose...
 * -------
 * Read a file containing the (x,y) positions and the types
 * of a set of np particles and then plot them.
 *
 * Version...
 * -------
 * 1.0  : 13-Jan-94 : First attempt
 * 1.1  : 02-Apr-94 : smaller particles as requested by Nick Robinson
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

/*-------------------------------------------------------*/

main ()
{
#define  NCHAR  256
char   txt[NCHAR];
char   FileName[32];
FILE   *pfile;

int    np, ip, *itype;
double *x, *y;

double xplotsize, yplotsize, xorig, yorig;
double xaxislength, yaxislength;
double xstart, xend, xtick, xrange;
double ystart, yend, ytick, yrange;
int    trueshape, change;
double xscale, yscale, minscale;
int    screen, plot_in_file;
double xx, yy, sc, ang;

printf ("\n");
printf ("-----------------\n");
printf ("Display Particles\n");
printf ("-----------------\n");
printf ("Rev. 1.1, 02-Apr-94\n");
printf ("\n");

/*
 * Read the Particle file
 */

printf ("Enter Particle File Name : ");
scanf ("%s", FileName);
if ((pfile = fopen(FileName, "r")) == NULL)
   {
   printf ("Could not open file: %s\n", FileName);
   exit (-1);
   }

fgets (txt, NCHAR, pfile);
sscanf (txt, "%d", &np);
printf ("Number of particles: %d\n", np);
x = NULL;
x = (double *) malloc (np * sizeof(double));
if (x == NULL)
   {
   printf ("Could not allocate memory for %d particles.\n", np);
   exit (-1);
   }
y = NULL;
y = (double *) malloc (np * sizeof(double));
if (y == NULL)
   {
   printf ("Could not allocate memory for %d particles.\n", np);
   exit (-1);
   }
itype = NULL;
itype = (int *) malloc (np * sizeof(int));
if (itype == NULL)
   {
   printf ("Could not allocate memory for %d particles.\n", np);
   exit (-1);
   }

for (ip = 0; ip < np; ++ip)
   {
   fgets (txt, NCHAR, pfile);
   sscanf (txt, "%lf %lf %d", &(x[ip]), &(y[ip]), &(itype[ip]) );
#  if 0
   printf ("particle[%d]: (%f, %f, %d)\n", ip, x[ip], y[ip]);
#  endif
   }

/*
 * Clean-up after the file reading.
 */
fclose (pfile);

/*
 ********************************
 * Set the Plotting Parameters. *
 ********************************
 *
 * Default plotting parameters.
 * Plotting units are millimetres (on an A4 page).
 */
xplotsize = 200.0;
yplotsize = 120.0;
xorig = 40.0;
yorig = 20.0;
trueshape = 1;

/*
 * Find extreme values and decide on scales.
 */
printf ("Find extremes...\n");
xstart = x[0];
xend   = x[0];
ystart = y[0];
yend   = y[0];
for (ip = 1; ip < np; ++ip)
   {
   if (x[ip] < xstart) xstart = x[ip];
   if (x[ip] > xend  ) xend   = x[ip];
   if (y[ip] < ystart) ystart = y[ip];
   if (y[ip] > yend  ) yend   = y[ip];
   }

/*
 * Set up nominal ranges.
 */
yrange = yend - ystart;
ytick  = yrange / 5.0;
xrange = xend - xstart;
xtick  = xrange / 5.0;

printf ("Nominal ranges...\n");
printf ("xstart = %e, xend = %e, xtick = %e\n", xstart, xend, xtick);
printf ("ystart = %e, yend = %e, ytick = %e\n", ystart, yend, ytick);
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
   xrange = xend - xstart;
   yrange = yend - ystart;
   }

/*
 * Finalize the Plot scales.
 */
printf ("trueshape (1/0): ");
scanf ("%d", &trueshape);
xscale = xplotsize / xrange;
yscale = yplotsize / yrange;
minscale = xscale;
if (yscale < minscale) minscale = yscale;
if (trueshape == 1)
   {
   xscale = minscale;
   yscale = minscale;
   }

/*
 * Scale the spatial variables.
 */
printf ("Scale the coordinates...\n");
for (ip = 0; ip < np; ++ip)
   {
   x[ip] = (x[ip] - xstart) * xscale;
   y[ip] = (y[ip] - ystart) * yscale;
   }

/*
 ****************
 * Do the plot. *
 ****************
 */

printf ("\nPlot to the screen, file (1/0 1/0): ");
scanf ("%d %d", &screen, &plot_in_file);
if (plot_in_file == 1)
   {
   printf ("Enter Plot File Name : ");
   scanf ("%s", FileName);
   }
else
   {
   strcpy(FileName, "noplot");
   }


ssp_BeginPlot (screen, plot_in_file, FileName, "ps", "landscape");

/*
 * Set the 2D plotting-page origin.
 */
ssp_SetOrigin (xorig, yorig);

ssp_SetLineThickness (0.35);
ssp_PlotText (-15.0, yplotsize+12.0, 2.5, 0.0, "x1 x2 dx", 8);
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

ssp_SetLineThickness (0.5);
xaxislength = xrange * xscale;
if (xplotsize > xaxislength) xaxislength = xplotsize;
ssp_PlotXAxis (0.0, 0.0, xaxislength, 0.0, xscale,
   xstart, xtick, "x, m", -1);
yaxislength = yrange * yscale;
if (yplotsize > yaxislength) yaxislength = yplotsize;
ssp_PlotYAxis (0.0, 0.0, yaxislength, 0.0, yscale,
   ystart, ytick, "y, m", -1);

ssp_SetClipLimits (0.0, 0.0, xplotsize, yplotsize);
ssp_SetClipOn ();

for (ip = 0; ip < np; ++ip)
   {
   ssp_PlotSymbol ( x[ip], y[ip], 0.8, itype[ip] );
   }

ssp_SetClipOff ();
ssp_EndPlot ();

}  /* end of dispoly.c */


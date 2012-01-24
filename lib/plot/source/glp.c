/** \file glp.c
 * \ingroup plot
 * \brief Generic Line Plotter -- no longer used; GNUPlot is better.
 *
 * Read a GENERIC data file containing data in the form
 * (x, f1, f2, ... fn) and produce line plots.
 * The data format may be deduced from the section of
 * code shown below which reads the file.
 * Note that several blocks of data may be included in the
 * one file.
 *
 * \author PA Jacobs
 * \version 1.0  :    Oct-91 : First Cut
 * \version 1.01 : 04-Jan-95 : DIMVAR increased to 15
 *
 */

/*-------------------------------------------------------*/

#include <math.h>
#include <stdio.h>

#include "compiler.h"

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
#define  BLOCKS  10
int    ivar, nvar, ix, iy, nnx[BLOCKS], nblock, jb;
#define  DIMVAR  15
char   VarName[DIMVAR][32];
double **Var[BLOCKS];

#define  NCHAR  256
char   title[NCHAR];
char   txt[NCHAR], DataFileName[32], PlotFileName[32];
FILE   *dfp;
int    screen, pfile;

struct ssp_data_vector v[BLOCKS];
double xstart, xend, xtick;
double ystart, yend, ytick;
double vstart, vend, vtick;
double minscale, xscale, yscale, xrange, yrange, vrange;
double xplotsize, yplotsize, xorig, yorig;
double symbolsize;
double xaxislength, yaxislength;
int    change, trueshape, do_symbols, do_lines;
double xx, yy, sc, ang;
int    isymb;


printf ("\n");
printf ("---------------------\n");
printf ("Generic Line Plotting\n");
printf ("---------------------\n");
printf ("Rev. 1.01, 04-Jan-95\n");
printf ("\n");

printf ("Enter Data File Name : ");
scanf ("%s", DataFileName);
if ((dfp = fopen(DataFileName, "r")) == NULL)
   {
   printf ("Could not open data file: %s\n", DataFileName);
   exit (-1);
   }

/*
 * ***********************
 * * Read the data file. *
 * ***********************
 */

fgets (title, NCHAR, dfp);
printf ("\nTitle: %s\n", title);

fgets (txt, NCHAR, dfp);
sscanf (txt, "%d", &nvar);
printf ("Number of variables: %d\n", nvar);
if (nvar > DIMVAR)
   {
   printf ("DIMVAR is too small; rebuild code.\n");
   exit (-1);
   }

for (ivar = 0; ivar < nvar; ++ivar)
   {
   fgets (txt, NCHAR, dfp);
   sscanf (txt, "%s", VarName[ivar]);
   printf ("Variable[%d]: %s\n", ivar, VarName[ivar]);
   }

/* NOTE : we start reading file without discarding line ends. */
fscanf (dfp, "%d", &nblock);
printf ("Number of blocks: %d\n", nblock);
if (nblock > BLOCKS)
   {
   printf ("BLOCKS is too small; rebuild code.\n");
   exit (-1);
   }

for (jb = 0; jb < nblock; ++jb)
   {
   fscanf (dfp, "%d", &(nnx[jb]) );
   printf ("nnx = %d\n", nnx[jb] );
   
   printf ("Allocating memory for data block[%d]...\n", jb);
   if ((Var[jb] = (double **) malloc (nvar * sizeof(double *)))
      == NULL)
      {
      printf ("Allocation failure\n");
      exit (-1);
      }
   for (ivar = 0; ivar < nvar; ++ivar)
      {
      if ( (Var[jb][ivar] = 
	   (double *) malloc(nnx[jb] * sizeof(double))) 
           == NULL )
         {
         printf ("Allocation failure: %d\n", ivar);
         exit (-1);
         }
      }

   for (ix = 0; ix < nnx[jb]; ++ix)
      for (ivar = 0; ivar < nvar; ++ivar)
         fscanf (dfp, "%lf", &Var[jb][ivar][ix]);

   }   /* for (jb... */

/*
 *******************************************
 * Decide which variable is to be plotted. *
 *******************************************
 */

Make_a_Plot:
printf ("\n-------------------------------\n");

PickVariable:
for (ivar = 0; ivar < nvar; ++ivar)
   printf ("%d=%s ", ivar, VarName[ivar]);
printf ("\nWhich variable (-1(exit), 0 ... %d): ", nvar-1);
scanf ("%d", &ivar);
if (ivar < 0) exit (0);
if (ivar >= nvar) goto PickVariable;

/*
 * Data space for the contouring routine.
 * We assume that the first two variables are x,y.
 */
printf ("Copy data...\n");
for (jb = 0; jb < nblock; ++jb)
   {
   ssp_AllocateVector (&v[jb], nnx[jb]);
   for (ix = 0; ix < nnx[jb]; ++ix)
      {
      v[jb].x[ix] = Var[jb][0][ix];
      v[jb].value[ix] = Var[jb][ivar][ix];
      }
   v[jb].imin = 0;
   v[jb].imax = nnx[jb] - 1;
   }


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
yorig = 16.0;
trueshape = 1;


/*
 * Find extreme values and decide on scales.
 */
printf ("Find extremes...\n");
for (jb = 0; jb < nblock; ++jb)
   ssp_VectorExtremes ( &(v[jb]) );

/*
 * Set up nominal ranges.
 */
ystart = v[0].vmin;
yend = v[0].vmax;
xstart = v[0].xmin;
xend = v[0].xmax;

if (nblock > 1)
   {
   for (jb = 1; jb < nblock; ++jb)
      {
      if (v[jb].vmin < ystart) ystart = v[jb].vmin;
      if (v[jb].vmax > yend) yend = v[jb].vmax;
      if (v[jb].xmin < xstart) xstart = v[jb].xmin;
      if (v[jb].xmax > xend) xend = v[jb].xmax;
      }
   }

yrange = yend - ystart;
ytick = yrange / 10.0;
xrange = xend - xstart;
xtick = xrange / 10.0;

printf ("Nominal ranges...\n");
printf ("xstart = %e, xend = %e, xtick = %e\n", xstart, xend, xtick);
printf ("ystart = %e, yend = %e, ytick = %e\n", ystart, yend, ytick);

if (fabs(ytick) < 1.0e-10 || fabs(xtick) < 1.0e-10 )
   {
   printf ("Data looks strange!\n");
   exit (1);
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
 * Plot scales.
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
printf ("Scale the spatial data...\n");
for (jb = 0; jb < nblock; ++jb)
   {
   for (ix = 0; ix < nnx[jb]; ++ix)
      {
      v[jb].x[ix] = (v[jb].x[ix] - xstart) * xscale;
      v[jb].value[ix] = (v[jb].value[ix] - ystart) * yscale;
      }
   }

/*
 ****************
 * Do the plot. *
 ****************
 */

/*
 * Symbols or lines.
 */
printf ("Plot Symbols, Lines (1/0 1/0): ");
scanf ("%d %d", &do_symbols, &do_lines);


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


/* ssp_BeginPlot (screen, pfile, PlotFileName, "hpgl", "landscape"); */
ssp_BeginPlot (screen, pfile, PlotFileName, "ps", "landscape");
ssp_SetOrigin (xorig, yorig);

ssp_SetLineThickness (0.35);
ssp_PlotText (2.0, yplotsize+7.0, 3.0, 0.0, DataFileName, 16);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotText (xx+5.0, yy, 3.0, 0.0, VarName[ivar], 16);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotText (xx+5.0, yy, 3.0, 0.0, title, 32);

ssp_PlotText (2.0, yplotsize+2.0, 3.0, 0.0, "x1 x2 dx", 8);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 3.0, 0.0, xstart, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 3.0, 0.0, xend, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 3.0, 0.0, xtick, 2);

ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotText (xx+2.0, yy, 3.0, 0.0, "y1 y2 dy", 8);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 3.0, 0.0, ystart, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 3.0, 0.0, yend, 2);
ssp_Where (&xx, &yy, &sc, &ang);
ssp_PlotENumber (xx+2.0, yy, 3.0, 0.0, ytick, 2);


ssp_SetLineThickness (0.5);
xaxislength = xrange * xscale;
if (xplotsize > xaxislength) xaxislength = xplotsize;
ssp_PlotXAxis (0.0, 0.0, xaxislength, 0.0, xscale, 
   xstart, xtick, VarName[0], -1);
yaxislength = yrange * yscale;
if (yplotsize > yaxislength) yaxislength = yplotsize;
ssp_PlotYAxis (0.0, 0.0, yaxislength, 0.0, yscale, 
   ystart, ytick, VarName[ivar], -1);

ssp_SetClipLimits (0.0, 0.0, xplotsize, yplotsize);
ssp_SetClipOn ();

symbolsize = yplotsize / 70.0;

for (jb = 0; jb < nblock; ++jb)
   {
   if (jb == 0)
      isymb = CIRCLE_SYM;
   else if (jb == 1)
      isymb = CROSS_SYM;
   else
      isymb = DIAMOND_SYM;

   if (do_symbols == 1)
      {
      ssp_SetLineThickness (0.18);
      ssp_PlotSymbol (v[jb].x[0], v[jb].value[0], symbolsize, isymb);
      for (ix = 1; ix < nnx[jb]; ++ix)
         {
         ssp_PlotSymbol (v[jb].x[ix], v[jb].value[ix], 
                        symbolsize, isymb);
         }
      }

   if (do_lines == 1)
      {
      ssp_SetLineThickness (0.35);
      ssp_Move (v[jb].x[0], v[jb].value[0]);
      for (ix = 1; ix < nnx[jb]; ++ix)
         {
         ssp_Plot (v[jb].x[ix], v[jb].value[ix]); 
         }
      }

   }  /* for (jb = 0; ...  */


ssp_SetClipOff ();
ssp_EndPlot (); 
goto Make_a_Plot;

}  /* end of glp.c */


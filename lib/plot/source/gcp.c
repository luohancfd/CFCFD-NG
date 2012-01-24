/** \file gcp.c
 * \ingroup plot
 * \brief Generic Contour Plotter -- no longer used.
 *
 * \author PA Jacobs
 */

/*
 *
 * Purpose...
 * -------
 * Read a GENERIC data file containing data in the form
 * (x, y, f1, f2, ... fn) and produce contour plots.
 * The data format may be deduced from the section of 
 * code shown below which reads the file.
 * Note that several blocks of data may be included in 
 * the one file.
 *
 * As an extra feature, geometry records may also be read
 * after the blocks of data.  There are an arbitrary number
 * of records of the form (x1, y1, x2, y2).
 *
 * Revisions...
 * ---------
 * 1.0 :    Oct-91 : First cut
 * 2.0 : 30-Sep-92 : geometry records added
 *                   optional edge plotting
 * 2.1 : 23-Jan-93 : mirror image in y=0 may also be plotted
 * 2.2 : 20-Mar-93 : increased the number of blocks
 * 2.21: 16-Aug-93 : fixed bug in the basic contouring
 * 2.22: 20-Oct-93 : Adjust contour limits, add colour contours.
 * 2.23: 10-Nov-93 : Changed size of written text to suit glyphs.
 * 2.24: 12-Nov-93 : Reduced the number of tickmarks and mirrored
 *                   the colour contouring.  Also added a colour
 *                   table to the right-hand side.
 * 2.3 : 01-Jan-94 : Added colour-filling option (as an alternative
 *                   to contouring)
 * 2.4 : 11-Jan-94 : Make the user interrogation a little easier.
 * 3.0 : 12-Jan-94 : Added vector plotting as originally coded by
 *                   Keith Weinman.
 * 3.1 : 02-Mar-94 : locate tails of the vectors at the data points
 * 3.2 : 12-Jun-94 : filled grey-scale plots (for interferograms)
 *                   move the plot to the right by 5mm to cope with
 *                   Ghostscript
 * 3.21: 12-Jul-94 : Thinner contour lines and boundary lines.
 * 3.22: 31-Aug-94 : More compact colour postscript.
 * 3.23: 04-Jan-95 : DIMVAR is now 15 to cope with new data from L1d
 *
 */

/*-------------------------------------------------------*/

#include <math.h>
#include <stdio.h>

#include "compiler.h"

#if (STDLIBH)
#  include <stdlib.h>
#endif

#if (CONIOH)
#  include <conio.h>
#endif

#include "geom.h"
#include "ssp.h"

/*-------------------------------------------------------*/

main ()
{
#define  BLOCKS  20
int    ivar, ivar2, nvar, ix, iy, nnx[BLOCKS], nny[BLOCKS];
int    jb, nblock;
#define  DIMVAR  15
char   VarName[DIMVAR][32];
double ***Var[BLOCKS];

#define  NCHAR  256
char   title[NCHAR];
char   txt[NCHAR], DataFileName[32], PlotFileName[32];
FILE   *dfp, *pfp;
int    screen, pfile;

int    icontour, imesh, iedge, nlevel, i;
int    do_vectors, do_colours, do_filled_colours, go_again;
int    do_filled_grey;
int    add_labels, land_scape;
double level, u, v;
double hue, bright, sat;
struct ssp_data_array fn[BLOCKS];
double xstart, xend, xtick;
double ystart, yend, ytick;
double vstart, vend, vtick;
double minscale, xscale, yscale, xrange, yrange, vrange;
double xplotsize, yplotsize, xorig, yorig;
double xaxislength, yaxislength;
int    change, trueshape, mirror_in_yaxis;
double xx, yy, sc, ang;
double y_temp;
double vector_scale, mag, xs, ys, len, headlen, width;
int    count, climit;
int    finished;
double x1, y1, x2, y2;

printf ("\n");
printf ("------------------------\n");
printf ("Generic Contour Plotting\n");
printf ("------------------------\n");
printf ("Rev. 3.3, 12-Dec-96\n");

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
   fscanf (dfp, "%d %d", &nny[jb], &nnx[jb] );
   printf ("nny = %d,  nnx = %d\n", nny[jb], nnx[jb]);

   printf ("Allocating memory for block[%d]...\n", jb);
   if ((Var[jb] = (double ***) malloc (nvar * sizeof(double **)))
      == NULL)
      {
      printf ("Allocation failure\n");
      exit (-1);
      }
   for (ivar = 0; ivar < nvar; ++ivar)
      {
      if ((Var[jb][ivar] =
	   (double **) malloc (nnx[jb] * sizeof(double *))) == NULL)
         {
         printf ("Allocation failure: %d\n", ivar);
         exit (-1);
         }
      for (ix = 0; ix < nnx[jb]; ++ix)
         {
         if ((Var[jb][ivar][ix] =
	      (double *) malloc (nny[jb] * sizeof(double))) == NULL)
	    {
	    printf ("Allocation failure: %d %d\n", ivar, ix);
	    exit (-1);
	    }
         }
      }

   /*
    * Read the field data for this block.
    */
   for (ix = 0; ix < nnx[jb]; ++ix)
      for (iy = 0; iy < nny[jb]; ++iy)
         {
         for (ivar = 0; ivar < nvar; ++ivar)
	    {
	    fscanf (dfp, "%lf", &Var[jb][ivar][ix][iy]);
	    }
         }

   }   /* for (jb = 0; ... */

/*
 *******************************************
 * Decide which variable is to be plotted. *
 *******************************************
 */

Make_a_Plot:
printf ("\n-------------------------------\n");

printf ("Enter a 1 or 0 for the following options...\n");
printf ("Add labels     (1/0) : ");
scanf ("%d", &add_labels);
printf ("Plot landscape (1/0) : ");
scanf ("%d", &land_scape);
printf ("Plot Vectors   (1/0) : ");
scanf ("%d", &do_vectors);
printf ("Filled-Colours (1/0) : ");
scanf ("%d", &do_filled_colours);
printf ("Filled-Grey    (1/0) : ");
scanf ("%d", &do_filled_grey);
printf ("Contours       (1/0) : ");
scanf ("%d", &icontour);
printf ("Colours        (1/0) : ");
scanf ("%d", &do_colours);
printf ("Add mesh       (1/0) : ");
scanf ("%d", &imesh);
printf ("Mesh-Edges     (1/0) : ");
scanf ("%d", &iedge);
printf ("Mirror-Image   (1/0) : ");
scanf ("%d", &mirror_in_yaxis);

printf ("\nYou have the following variables available...\n");
for (ivar = 0; ivar < nvar; ++ivar)
   printf ("%d=%s ", ivar, VarName[ivar]);
printf ("\n");

if ( do_vectors )
   {
   /*
    * Pick variables for the u,v components.
    */
   PickVariable1:
   printf ("Variable for u-component (-1(exit), 0 ... %d): ", nvar-1);
   scanf ("%d", &ivar);
   if (ivar < 0) exit (0);
   if (ivar >= nvar) goto PickVariable1;
   PickVariable2:
   printf ("Variable for v-component (-1(exit), 0 ... %d): ", nvar-1);
   scanf ("%d", &ivar2);
   if (ivar2 < 0) exit (0);
   if (ivar2 >= nvar) goto PickVariable2;
   }
else
   {
   /*
    * Pick one variable only.
    */
   PickVariable:
   printf ("Which variable (-1(exit), 0 ... %d): ", nvar-1);
   scanf ("%d", &ivar);
   if (ivar < 0) exit (0);
   if (ivar >= nvar) goto PickVariable;
   }

/*
 * Data space for the contouring routine.
 * We assume that the first two variables are x,y.
 */
printf ("Copy data...\n");
for (jb = 0; jb < nblock; ++jb)
   {
   ssp_AllocateArray (&fn[jb], nnx[jb], nny[jb]);
   for (ix = 0; ix < nnx[jb]; ++ix)
      for (iy = 0; iy < nny[jb]; ++iy)
         {
         fn[jb].x[ix][iy] = Var[jb][0][ix][iy];
         fn[jb].y[ix][iy] = Var[jb][1][ix][iy];
         if ( do_vectors )
            {
            u = Var[jb][ivar][ix][iy];
            v = Var[jb][ivar2][ix][iy];
            mag = sqrt( u * u + v * v );
            fn[jb].u[ix][iy] = u;
            fn[jb].v[ix][iy] = v;
            fn[jb].value[ix][iy] = mag;
            }
         else
            {
            fn[jb].value[ix][iy] = Var[jb][ivar][ix][iy];
            }
	 /*
	 printf ("[%d][%d] x=%e, y=%e, value=%e\n", ix, iy,
	    fn[jb].x[ix][iy], fn[jb].y[ix][iy], fn[jb].value[ix][iy]);
         */
         }
   fn[jb].ixmin = 0;
   fn[jb].ixmax = nnx[jb] - 1;
   fn[jb].iymin = 0;
   fn[jb].iymax = nny[jb] - 1;
   }


/*
 ********************************
 * Set the Plotting Parameters. *
 ********************************
 *
 * Default plotting parameters.
 * Plotting units are millimetres (on an A4 page).
 */
if ( land_scape == 1) 
   {
   xplotsize = 200.0;
   yplotsize = 120.0;
   xorig = 45.0;
   yorig = 20.0;
   }
else
   {
   xplotsize = 150.0;
   yplotsize = 100.0;
   xorig = 30.0;
   yorig = 20.0;
   } /* endif */
nlevel = 16;
trueshape = 1;

/*
 * Find extreme values and decide on scales.
 */
printf ("Find extremes...\n");
for (jb = 0; jb < nblock; ++jb)
   ssp_ArrayExtremes ( &(fn[jb]) );

/*
 * Set up nominal ranges.
 */
vstart = fn[0].vmin;
vend = fn[0].vmax;
xstart = fn[0].xmin;
xend = fn[0].xmax;
ystart = fn[0].ymin;
yend = fn[0].ymax;

if (nblock > 1)
   {
   for (jb = 1; jb < nblock; ++jb)
      {
      if (fn[jb].vmin < vstart) vstart = fn[jb].vmin;
      if (fn[jb].vmax > vend) vend = fn[jb].vmax;
      if (fn[jb].xmin < xstart) xstart = fn[jb].xmin;
      if (fn[jb].xmax > xend) xend = fn[jb].xmax;
      if (fn[jb].ymin < ystart) ystart = fn[jb].ymin;
      if (fn[jb].ymax > yend) yend = fn[jb].ymax;
      }
   }

/* Allow space for the mirror image in y = 0. */
if (mirror_in_yaxis == 1)
   {
   ystart = -yend;
   }

vrange = vend - vstart;
vtick = vrange / ((double) nlevel);
if ( do_vectors != 1 )
   {
   /* Bring the ends in by half an interval */
   vstart += 0.5 * vtick;
   vend   -= 0.5 * vtick;
   }

xrange = xend - xstart;
xtick = xrange / 5.0;
yrange = yend - ystart;
ytick = yrange / 5.0;

printf ("Nominal ranges...\n");
printf ("xstart = %e, xend = %e, xtick = %e\n", xstart, xend, xtick);
printf ("ystart = %e, yend = %e, ytick = %e\n", ystart, yend, ytick);
printf ("vstart = %e, vend = %e, vtick = %e\n", vstart, vend, vtick);

if (fabs(vtick) < 1.0e-10 || fabs(xtick) < 1.0e-10 ||
    fabs(ytick) < 1.0e-10)
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
   printf ("vstart vend vtick: ");
   scanf ("%lf %lf %lf", &vstart, &vend, &vtick);
   xrange = xend - xstart;
   yrange = yend - ystart;
   vrange = vend - vstart;
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

if ( do_vectors )
   {
   printf ("Scale factor for vectors (mm per vector unit) : ");
   scanf ("%lf", &vector_scale);
   }

/*
 * Offset for creating mirror images in the y-axis.
 */
y_temp = 2.0 * ystart * yscale;

/*
 * Scale the spatial variables.
 */
printf ("Scale the spatial data...\n");
for (jb = 0; jb < nblock; ++jb)
   {
   for (ix = 0; ix < nnx[jb]; ++ix)
      for (iy = 0; iy < nny[jb]; ++iy)
         {
         fn[jb].x[ix][iy] = (fn[jb].x[ix][iy] - xstart) * xscale;
         fn[jb].y[ix][iy] = (fn[jb].y[ix][iy] - ystart) * yscale;
         }
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
   if ((pfp = fopen(PlotFileName, "w")) == NULL)
      {
      printf ("Could not open PLOT file: %s\n", PlotFileName);
      exit (-1);
      }
   }

if ( land_scape == 1 )
   {
   ssp_BeginPlot (screen, pfile, PlotFileName, "ps", "landscape");
   }
else
   {
   ssp_BeginPlot (screen, pfile, PlotFileName, "ps", "portrait");
   } /* endif */
ssp_SetOrigin (xorig, yorig);

ssp_SetLineThickness (0.2);
if ( add_labels == 1 )
   {
   ssp_PlotText (-15.0, yplotsize+17.0, 2.5, 0.0, DataFileName, 16);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotText (xx+5.0, yy, 2.5, 0.0, VarName[ivar], 16);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotText (xx+5.0, yy, 2.5, 0.0, title, 32);

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

   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotText (xx+2.0, yy, 2.5, 0.0, "v1 v2 dv", 8);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, vstart, 2);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, vend, 2);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, vtick, 2);
   }

ssp_SetLineThickness (0.35);
xaxislength = xrange * xscale;
if (xplotsize > xaxislength) xaxislength = xplotsize;
ssp_PlotXAxis (0.0, 0.0, xaxislength, 0.0, xscale,
   xstart, xtick, VarName[0], -1);
yaxislength = yrange * yscale;
if (yplotsize > yaxislength) yaxislength = yplotsize;
ssp_PlotYAxis (0.0, 0.0, yaxislength, 0.0, yscale,
   ystart, ytick, VarName[1], -1);

ssp_SetClipLimits (0.0, 0.0, xplotsize, yplotsize);
ssp_SetClipOn ();

if ( (icontour && do_colours) || do_filled_colours )
   {
   /*
    * Put a colour map in here by writing the contour
    * values down the right-hand side of the page in the
    * appropriate colours.
    */
   ssp_SetClipOff ();
   for (i = 15; i >= 0; --i)
      {
      level = vstart + i * (vend - vstart) / 15.0;
      hue = ssp_MapToColour (level, vstart, vend);
      bright = 1.0;
      sat = 1.0;
      ssp_SetHSBColour (hue, sat, bright);
      xx = 210.0;
      yy = 20.0 + 5.0 * i;
      ssp_PlotENumber (xx, yy, 2.5, 0.0, level, 3);
      }
   ssp_SetNoColour ();
   ssp_SetClipOn ();
   }

if ( do_filled_grey )
   {
   /*
    * Put a grey-scale map in here by writing the contour
    * values down the right-hand side of the page in the
    * appropriate colours.
    */
   ssp_SetClipOff ();
   for (i = 15; i >= 0; --i)
      {
      level = vstart + i * (vend - vstart) / 15.0;
      hue = 0.0;
      sat = 0.0;  /* zero saturation should force a grey scale */
      bright = (level - vstart) / (vend - vstart);
      ssp_SetHSBColour (hue, sat, bright);
      xx = 210.0;
      yy = 20.0 + 5.0 * i;
      ssp_PlotENumber (xx, yy, 2.5, 0.0, level, 3);
      }
   ssp_SetNoColour ();
   ssp_SetClipOn ();
   }

if ( icontour )
   {
   /*
    * The contours themselves.
    */
   for (jb = 0; jb < nblock; ++jb)
      {
      ssp_SetLineThickness (0.25);
      for (level = vstart; level <= vend+0.001*vtick; level += vtick)
         {
         sprintf (txt, "Block[%d], Contour level = %f", jb, level);
         ssp_Message (txt);
         if (do_colours)
            {
            hue = ssp_MapToColour (level, vstart, vend);
            bright = 1.0;
            sat = 1.0;
            ssp_SetHSBColour (hue, sat, bright);
            }
         ssp_ContourArray (&fn[jb], level);
         if (do_colours) ssp_SetNoColour();
         }
      ssp_SetLineThickness (0.12);
      if (iedge == 1) ssp_PlotEdge ( &fn[jb] );
      }  /* for (jb... */
   }  /* if (icontour ... */

if ( do_filled_colours || do_filled_grey )
   {
   /*
    * Fill the cells with colour or grey-shade.
    */
   for (jb = 0; jb < nblock; ++jb)
      {
      ssp_ColourFillArray ( &fn[jb], vstart, vend, do_filled_grey );
      ssp_SetNoColour();
      ssp_SetLineThickness (0.12);
      if (iedge == 1) ssp_PlotEdge ( &fn[jb] );
      }  /* for (jb... */
   }  /* if (do_filled_colours... */

if ( do_vectors )
   {
   /*
    * Draw an array of arrows to represent the vector quantities.
    */
   for (jb = 0; jb < nblock; ++jb)
      {
      ssp_SetLineThickness (0.18);
      count = 0;
      climit = 100;
      for (ix = 0; ix < nnx[jb]; ++ix)
         {
         for (iy = 0; iy < nny[jb]; ++iy)
            {
            ++count;
            if (count >= climit)
               {
               sprintf (txt, "Block[%d], ix= %d, iy=%d", jb, ix, iy);
               ssp_Message (txt);
               count = 0;
               }

            xs = fn[jb].x[ix][iy];
            ys = fn[jb].y[ix][iy];
            u = fn[jb].u[ix][iy];
            v = fn[jb].v[ix][iy];
            mag = fn[jb].value[ix][iy];

            x1 = xs;
            y1 = ys;
            x2 = xs + u * vector_scale;
            y2 = ys + v * vector_scale;
            len = mag * vector_scale;
            headlen = 0.2 * len;
            if (headlen > 3.0) headlen = 3.0;
            width = 0.5 * headlen;

            if (do_colours)
               {
               /* Set the colour on the vector magnitude. */
               hue = ssp_MapToColour (mag, vstart, vend);
               bright = 1.0;
               sat = 1.0;
               ssp_SetHSBColour (hue, sat, bright);
               }
            ssp_Arrow( x1, y1, x2, y2, headlen, width, 1 );
            } /* for (ix... */
         } /* for (iy... */

      if (do_colours) ssp_SetNoColour();
      ssp_SetLineThickness (0.12);
      if (iedge == 1) ssp_PlotEdge ( &fn[jb] );
      }  /* for (jb... */
   }  /* if ( do_vectors ... */

if (imesh == 1)
   {
   ssp_SetLineThickness (0.10);
   for (jb = 0; jb < nblock; ++jb) ssp_PlotMesh ( &fn[jb] );
   }

/*
 * For the mirror image,
 * invert the y-coordinate and redo the contours.
 */
if (mirror_in_yaxis == 1)
   {
   /*
    * Negate the y-coordinates.
    */
   for (jb = 0; jb < nblock; ++jb)
      {
      for (ix = 0; ix < nnx[jb]; ++ix)
         for (iy = 0; iy < nny[jb]; ++iy)
            {
            fn[jb].y[ix][iy] = -fn[jb].y[ix][iy] - y_temp;
            if ( do_vectors ) fn[jb].v[ix][iy] = -fn[jb].v[ix][iy];
            }
      }  /* for (jb = 0.... */

   if (icontour == 1)
      {
      /*
       * Redo the contours.
       */
      for (jb = 0; jb < nblock; ++jb)
         {
         ssp_SetLineThickness (0.25);
         for (level = vstart; level <= vend+0.001*vtick; level += vtick)
            {
            sprintf (txt, "Block[%d], Contour level = %f", jb, level);
            ssp_Message (txt);
            if (do_colours)
               {
               hue = ssp_MapToColour (level, vstart, vend);
               bright = 1.0;
               sat = 1.0;
               ssp_SetHSBColour (hue, sat, bright);
               }
            ssp_ContourArray (&fn[jb], level);
            if (do_colours) ssp_SetNoColour();
            }
         ssp_SetLineThickness (0.12);
         if (iedge == 1) ssp_PlotEdge ( &fn[jb] );
         }
      }  /* if (icontour == 1 ... */

   if ( do_filled_colours || do_filled_grey )
      {
      /*
       * Fill the cells with colour or grey-shade.
       */
      for (jb = 0; jb < nblock; ++jb)
         {
         ssp_ColourFillArray ( &fn[jb], vstart, vend, do_filled_grey );
         ssp_SetNoColour();
         ssp_SetLineThickness (0.12);
         if (iedge == 1) ssp_PlotEdge ( &fn[jb] );
         }  /* for (jb... */
      }  /* if (do_filled_colours... */

   if ( do_vectors )
      {
      /*
       * Draw an array of arrows to represent the vector quantities.
       */
      for (jb = 0; jb < nblock; ++jb)
         {
         ssp_SetLineThickness (0.18);
         count = 0;
         climit = 100;
         for (ix = 0; ix < nnx[jb]; ++ix)
            {
            for (iy = 0; iy < nny[jb]; ++iy)
               {
               ++count;
               if (count >= climit)
                  {
                  sprintf (txt, "Block[%d], ix= %d, iy=%d", jb, ix, iy);
                  ssp_Message (txt);
                  count = 0;
                  }

               xs = fn[jb].x[ix][iy];
               ys = fn[jb].y[ix][iy];
               u = fn[jb].u[ix][iy];
               v = fn[jb].v[ix][iy];
               mag = fn[jb].value[ix][iy];

               x1 = xs;
               y1 = ys;
               x2 = xs + u * vector_scale;
               y2 = ys + v * vector_scale;
               len = mag * vector_scale;
               headlen = 0.2 * len;
               if (headlen > 3.0) headlen = 3.0;
               width = 0.5 * headlen;

               if (do_colours)
                  {
                  /* Set the colour on the vector magnitude. */
                  hue = ssp_MapToColour (mag, vstart, vend);
                  bright = 1.0;
                  sat = 1.0;
                  ssp_SetHSBColour (hue, sat, bright);
                  }
               ssp_Arrow( x1, y1, x2, y2, headlen, width, 1 );
               } /* for (ix... */
            } /* for (iy... */

         if (do_colours) ssp_SetNoColour();
         ssp_SetLineThickness (0.12);
         if (iedge == 1) ssp_PlotEdge ( &fn[jb] );
         }  /* for (jb... */
      }  /* if ( do_vectors ... */

   if (imesh == 1)
      {
      ssp_SetLineThickness (0.10);
      for (jb = 0; jb < nblock; ++jb) ssp_PlotMesh ( &fn[jb] );
      }
   }  /* if (mirror...  */


/*
 * Read and act on any trailing geometry records.
 * This will be done on the first pass so, any subsequent
 * plots will miss out.
 * I guess that this can be improved by storing the segment
 * data points and replaying for each plot.
 */
ssp_SetLineThickness (0.35);
finished = 0;
while (finished != 1)
   {
   /* Attempt to read 4 coordinates. */
   if ( fscanf(dfp, "%lf %lf %lf %lf", &x1, &y1, &x2, &y2) == EOF )
      finished = 1;
   /* printf ("x1 = %f, y1 = %f, x2 = %f, y2 = %f\n", x1, y1, x2, y2); */
   if (finished != 1)
      {
      if (mirror_in_yaxis == 1)
         {
         /*
          * If the segment is not along the y-axis,
          * scale the line segment and plot it and its image.
          */
         if (fabs(y1) > 1.0e-6 || fabs(y2) > 1.0e-6)
            {
            x1 = (x1 - xstart) * xscale;
            y1 = (y1 - ystart) * yscale;
            x2 = (x2 - xstart) * xscale;
            y2 = (y2 - ystart) * yscale;
            ssp_Move (x1, y1);
            ssp_Plot (x2, y2);
            ssp_Move (x1, -y1 - y_temp);
            ssp_Plot (x2, -y2 - y_temp);
            }
         }
      else
         {
         /*
          * Scale the line segment and plot it.
          */
         x1 = (x1 - xstart) * xscale;
         y1 = (y1 - ystart) * yscale;
         x2 = (x2 - xstart) * xscale;
         y2 = (y2 - ystart) * yscale;
         ssp_Move (x1, y1);
         ssp_Plot (x2, y2);
         }

      }  /* if (finished != 1... */
   }     /* while (!finished... */

ssp_SetClipOff ();
ssp_EndPlot ();

printf ("\nEnd of plot; do another (1/0)? ");
scanf ("%d", &go_again);
if ( go_again == 1) goto Make_a_Plot;

return 0;
}  /* end of gcp.c */


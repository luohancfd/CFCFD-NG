/** \file contour.c
 * \ingroup plot
 * \brief Generic Contour Plotter for 2-Dimensional Data.
 *
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
 * \author PA Jacobs
 *
 * \version 1.0 :    Oct-91 : First cut
 * \version 2.0 : 30-Sep-92 : geometry records added
 * \version 3.0 : 12-Jan-94 : Added vector plotting as originally coded by
 *                Keith Weinman.
 * \version 4.0 : 17-Nov-97 : UNIX style command line and GIF file generation
 *                Renamed to "contour.c"
 * \version 4.01: 20-Dec-97 : Finer lines for drawing smaller postscript plots.
 *                Fix too many tic-marks for true-shape plots.
 * \version 4.02: 04-Feb-97 : Default column is set to 2 (after columns 0 and 1).
 *                This should be a sensible default in the cases 
 *                where there is only one dependent variable.
 * \version 4.03: 24-Feb-00 : ixskip and iyskip added to ssp_PlotMesh()
 * \version 4.04: 29-May-05 : Added option for increasing number of contour
 *                lines without specifying value range and increments.
 */

/*-------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include "../../util/source/compiler.h"
#if STDLIBH
#  include <stdlib.h>
#  include <string.h>
#endif

#if CONIOH
#  include <conio.h>
#endif

#include "../../geometry/source/geom.h"
#include "ssp.h"
#include "../../util/source/logfile.h"
#include "../../util/source/useful.h"

/*-------------------------------------------------------*/

int main ( int argc, char *argv[] )
{
char   DataFileName[32];
FILE   *dfp;
#define  BLOCKS  100
int    ivar, ivar1, ivar2, nvar;
int    ix, iy, nnx[BLOCKS], nny[BLOCKS];
int    jb, nblock;
#define  DIMVAR  50
char   VarName[DIMVAR][32];
double ***Var[BLOCKS];

#define  NCHAR  512
char   title[NCHAR];
char   *token, txt[NCHAR], msg_text[NCHAR], PlotFileName[32];
char   orientation[12], output_type[12];
int    do_contours, do_mesh, do_edge;
int    ixskip, iyskip;
int    draw_axes, draw_table;
int    do_vectors, with_colours, do_filled;
int    true_shape, mirror_flag;
int    add_labels;
int    nlevel, i;
double level, u, v;
double hue, bright, sat;
double pw_axis, pw_labels, pw_contours;
double pw_geometry, pw_mesh, pw_vectors;

int    command_line_error;
int    xrange_given;
double xmin_given, xmax_given, dx_given;
int    yrange_given;
double ymin_given, ymax_given, dy_given;
int    vrange_given;
double vmin_given, vmax_given, dv_given;

struct ssp_data_array fn[BLOCKS];
double xstart, xend, xtick;
double ystart, yend, ytick;
double vstart, vend, vtick;
double xscale, yscale, xrange, yrange, vrange;
double xplotsize, yplotsize, pixels_per_mm, xorig, yorig;
char pixels_per_mm_str[32];
double xaxislength, yaxislength;
double xx, yy, sc, ang;
double y_temp;
double vector_scale, mag, xs, ys, len, headlen, width;
int    count, climit;
int    finished;
double x1, y1, x2, y2;

open_log_file("mb_cont.log");
log_message ("------------------------\n", 1);
log_message ("Generic Contour Plotting\n", 1);
log_message ("------------------------\n", 1);
log_message ("$Id: contour.c,v 1.9 2003/08/28 06:13:13 peterj Exp $\n", 1);
log_message ("\n", 1);

/*
 * Set defaults...
 */
strcpy (DataFileName, "default.gen");
strcpy (PlotFileName, "default.ps");
strcpy (output_type,  "ps");
strcpy (orientation,  "portrait");
add_labels     = 1;
do_filled      = 0;
do_contours    = 1;
do_mesh        = 0;
ixskip         = 1;
iyskip         = 1;
do_edge        = 0;
with_colours   = 0;
draw_table     = 1;
draw_axes      = 1;
mirror_flag    = 0;
xplotsize      = 120.0;  /* millimetres */
yplotsize      = 120.0;  /* millimetres */
pixels_per_mm  = 16.0;   /* Stefan Hess' value to get finer gif contour plots. */ 
xorig          = 30.0;   /* mm          */
yorig          = 20.0;   /* mm          */
pw_axis        = 0.20;   /* pen width axis, mm */
pw_labels      = 0.20;
pw_contours    = 0.12;
pw_mesh        = 0.12;
pw_geometry    = 0.20;
pw_vectors     = 0.12;
ivar1          = 2;   /* guess at index of contour variable      */
                      /* This shoule be safe as x and y are      */
                      /* generally columns 0 and 1 respectively. */
nlevel         = 16;
true_shape     = 1;
vrange_given   = 0;
xrange_given   = 0;
yrange_given   = 0;
/*
 * One day, we need to return and fix the following items.
 */
do_vectors     = 0;
vector_scale   = 1.0;  /* mm per vector unit */
ivar2          = 6;


/*
 * Decode command line arguments.
 */
command_line_error = 0;
if (argc < 2)
   {
   /* no arguments supplied */
   command_line_error = 1; 
   goto usage;
   }

i = 1;
while (i < argc)
   {
   /* process the next command-line argument */

   if (strcmp(argv[i],"-fi") == 0)
      {
      /* Set the input data filename */
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      strcpy(DataFileName, argv[i]);
      i++;
      sprintf (msg_text, "Setting DataFileName = %s\n", DataFileName);
      log_message (msg_text, 0);
      }
   else if (strcmp(argv[i],"-fo") == 0)
      {
      /* Set the filename for output. */
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      strcpy(PlotFileName, argv[i]);
      i++;
      sprintf (msg_text, "Setting PlotFileName = %s\n", PlotFileName);
      log_message (msg_text, 0);
      }
   else if (strcmp(argv[i],"-gif") == 0)
      {
      strcpy(output_type, "gif");
      i++;
      log_message ("Selecting GIF output.\n", 0);
      }
   else if (strcmp(argv[i],"-pixelspermm") == 0)
      {
      /* Set the increment for plotting the mesh. */
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &pixels_per_mm);
      i++;
      sprintf (msg_text, "Setting pixel_sper_mm = %f.2\n", pixels_per_mm);
      log_message (msg_text, 0);
      }
   else if (strcmp(argv[i],"-ps") == 0)
      {
      strcpy(output_type, "ps");
      i++;
      log_message ("Selecting POSTSCRIPT output.\n", 0);
      }
   else if (strcmp(argv[i],"-hpgl") == 0)
      {
      strcpy(output_type, "hpgl");
      i++;
      log_message ("Selecting HPGL output.\n", 0);
      }
   else if (strcmp(argv[i],"-nolabel") == 0)
      {
      add_labels = 0;
      i++;
      log_message ("Setting add_labels = 0.\n", 0);
      }
   else if (strcmp(argv[i],"-noaxes") == 0)
      {
      draw_axes = 0;
      i++;
      log_message ("Suppress drawing of axes.\n", 0);
      }
   else if (strcmp(argv[i],"-notable") == 0)
      {
      draw_table = 0;
      i++;
      log_message ("Suppress drawing of table of colour values.\n", 0);
      }
   else if (strcmp(argv[i],"-mirror") == 0)
      {
      mirror_flag = 1;
      i++;
      log_message ("Setting mirror image in y-axis.\n", 0);
      }
   else if (strcmp(argv[i],"-notrueshape") == 0)
      {
      true_shape = 0;
      i++;
      log_message ("Setting true_shape = 0.\n", 0);
      }
   else if (strcmp(argv[i],"-colour") == 0)
      {
      with_colours = 1;
      i++;
      log_message ("Setting with_colours = 1.\n", 0);
      }
   else if (strcmp(argv[i],"-fill") == 0)
      {
      do_filled   = 1;
      do_contours = 0;
      i++;
      log_message ("Setting filled rather than contours.\n", 0);
      }
   else if (strcmp(argv[i],"-nocontours") == 0)
      {
      do_contours = 0;
      i++;
      log_message ("Suppress contours.\n", 0);
      }
   else if (strcmp(argv[i],"-edge") == 0)
      {
      do_edge = 1;
      i++;
      log_message ("Setting do_edge = 1.\n", 0);
      }
   else if (strcmp(argv[i],"-mesh") == 0)
      {
      do_mesh = 1;
      i++;
      log_message ("Setting do_mesh = 1.\n", 0);
      }
   else if (strcmp(argv[i],"-var") == 0)
      {
      /* Set variable index. */
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%d", &ivar1);
      i++;
      sprintf (msg_text, "Setting var1 index = %d\n", ivar1);
      log_message (msg_text, 0);
      }
   else if (strcmp(argv[i],"-ixskip") == 0)
      {
      /* Set the increment for plotting the mesh. */
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%d", &ixskip);
      i++;
      sprintf (msg_text, "Setting ixskip = %d\n", ixskip);
      log_message (msg_text, 0);
      }
   else if (strcmp(argv[i],"-iyskip") == 0)
      {
      /* Set the increment for plotting the mesh. */
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%d", &iyskip);
      i++;
      sprintf (msg_text, "Setting iyskip = %d\n", iyskip);
      log_message (msg_text, 0);
      }
   else if (strcmp(argv[i],"-levels") == 0)
      {
      /* Set vrange_given. */
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &vmin_given);
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &vmax_given);
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &dv_given);
      i++;
      vrange_given = 1;
      sprintf (msg_text, "Setting vrange = %g, %g, %g\n", 
         vmin_given, vmax_given, dv_given);
      log_message (msg_text, 0);
      }
   else if (strcmp(argv[i],"-xrange") == 0)
      {
      /* Set xrange_given. */
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &xmin_given);
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &xmax_given);
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &dx_given);
      i++;
      xrange_given = 1;
      sprintf (msg_text, "Setting xrange = %g, %g, %g\n", 
         xmin_given, xmax_given, dx_given);
      log_message (msg_text, 0);
      }
   else if (strcmp(argv[i],"-yrange") == 0)
      {
      /* Set yrange_given. */
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &ymin_given);
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &ymax_given);
      i++;
      if (i >= argc) { command_line_error = 1; goto usage; }
      sscanf (argv[i], "%lf", &dy_given);
      i++;
      yrange_given = 1;
      sprintf (msg_text, "Setting yrange = %g, %g, %g\n", 
         ymin_given, ymax_given, dy_given);
      log_message (msg_text, 0);
      }
   else if (strcmp(argv[i],"-numconts") == 0)
      {
      /* Set number of contours, independent of value range */
      i++;
      vrange_given = 0;
      sscanf (argv[i], "%d", &nlevel);
      i++;
      sprintf (msg_text, "Setting number of contours = %d\n", nlevel);
      log_message (msg_text, 0);
      }
   else 
      {
      sprintf (msg_text, "Unknown option: %s\n", argv[i]);
      log_message (msg_text, 1);
      command_line_error = 1;
      goto usage;
      }
   } /* while... */

/*
 * If the command line arguments are incorrect,
 * write some advice to the console then finish.
 */
usage:
if (command_line_error == 1)
   {
   log_message ("Usage: (defaults are shown in parentheses)\n", 1);

   log_message ("-fi <data_file>              (default.gen)\n", 1);
   log_message ("-fo <output_file>            (default.ps)\n", 1);
   log_message ("-gif|-ps|-hpgl               (-ps)\n", 1);
   log_message ("-pixelspermm <n>             (16.0)\n", 1);

   log_message ("-var <index>                 (2)\n", 1);
   log_message ("-levels <vmin> <vmax> <dv>   (search data)\n", 1);
   log_message ("-xrange <xmin> <xmax> <dx>   (search data)\n", 1);
   log_message ("-yrange <ymin> <ymax> <dy>   (search data)\n", 1);

   log_message ("-nolabel                     (label)\n", 1);
   log_message ("-noaxes                      (draw axes)\n", 1);
   log_message ("-notable                     (draw table)\n", 1);
   log_message ("-mirror                      (no mirror image)\n", 1);
   log_message ("-notrueshape                 (trueshape)\n", 1);
   log_message ("-colour                      (no colour)\n", 1);
   log_message ("-fill                        (contour)\n", 1);
   log_message ("-nocontours                  (contour)\n", 1);
   log_message ("-edge                        (no edge)\n", 1);
   log_message ("-mesh                        (no mesh)\n", 1);
   log_message ("-ixskip <i>                  (1)\n", 1);
   log_message ("-iyskip <i>                  (1)\n", 1);

   log_message ("\nNotes...\n", 1);
   log_message ("The -var option selects the variable to contour.\n", 1);
   exit (-1);
   }


/*
 * ***********************
 * * Read the data file. *
 * ***********************
 */

if ((dfp = fopen(DataFileName, "r")) == NULL) {
    printf ("Could not open data file: %s\n", DataFileName);
    exit (-1);
 }

/* The first contains a title string. */
if ( fgets (title, NCHAR, dfp) == NULL ) {
    printf("Problem reading file.\n");
    exit(BAD_INPUT_ERROR);
}
printf ("\nTitle: %s\n", title);

/* The second line contains the number of variables. */
if ( fgets (txt, NCHAR, dfp) == NULL ) {
    printf("Problem reading file.\n");
    exit(BAD_INPUT_ERROR);
}
sscanf (txt, "%d", &nvar);
printf ("Number of variables: %d\n", nvar);
if (nvar > DIMVAR)
   {
   printf ("DIMVAR is too small; rebuild code.\n");
   exit (-1);
   }

/* The variable names follow, one per line. */
for (ivar = 0; ivar < nvar; ++ivar)
   {
   if ( fgets (txt, NCHAR, dfp) == NULL ) {
      printf("Problem reading file.\n");
      exit(BAD_INPUT_ERROR);
   }
   sscanf (txt, "%s", VarName[ivar]);
   printf ("Variable[%d]: %s\n", ivar, VarName[ivar]);
   }

/*
 * After the variable names, read the number of blocks.
 */
if ( fgets (txt, NCHAR, dfp) == NULL ) {
    printf("Problem reading file.\n");
    exit(BAD_INPUT_ERROR);
}
sscanf (txt, "%d", &nblock);
printf ("Number of blocks: %d\n", nblock);
if (nblock > BLOCKS)
   {
   printf ("BLOCKS is too small; rebuild code.\n");
   exit (-1);
   }

/* 
 * For each block, allocate memory and read data. 
 */
for (jb = 0; jb < nblock; ++jb)
   {
   /* Read the number of points in the y- and x-index directions. */
   if ( fgets (txt, NCHAR, dfp) == NULL ) {
      printf("Problem reading file.\n");
      exit(BAD_INPUT_ERROR);
   }
   sscanf (txt, "%d %d", &nny[jb], &nnx[jb] );
   printf ("nny = %d,  nnx = %d\n", nny[jb], nnx[jb]);

   printf ("Allocating memory for block[%d]...\n", jb);
   if ((Var[jb] = (double ***) malloc (nvar * sizeof(double **))) == NULL)
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
   for (ix = 0; ix < nnx[jb]; ++ix) {
      for (iy = 0; iy < nny[jb]; ++iy) {
         /* 
          * Each line represents on data point on the mesh.
	  * Data values are assumed to be separated by white-space. 
	  */
         if ( fgets (txt, NCHAR, dfp) == NULL ) {
            printf("Problem reading file.\n");
            exit(BAD_INPUT_ERROR);
         }
	 token = strtok( txt, " " );
	 sscanf(token, "%lf", &Var[jb][0][ix][iy]);
         for (ivar = 1; ivar < nvar; ++ivar) {
	    token = strtok( NULL, " " );
	    sscanf(token, "%lf", &Var[jb][ivar][ix][iy]);
	 }
      } /* end of iy loop */
   } /* end ix loop */

   }   /* end of jb loop  */


/*
 * Check the indices of the selected variables.
 */
if (ivar1 >= 0 || ivar1 < nvar)
   {
   printf ("Selected variable 1 = %s\n", VarName[ivar1]);
   }
else
   {
   printf ("Variable 1 index invalid: %d\n", ivar1);
   exit(-1);
   }
if ( do_vectors )
   {
   if (ivar2 >= 0 || ivar2 < nvar)
      {
      printf ("Selected variable 2 = %s\n", VarName[ivar2]);
      }
   else
      {
      printf ("Variable 2 index invalid: %d\n", ivar2);
      exit(-1);
      }
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
            u = Var[jb][ivar1][ix][iy];
            v = Var[jb][ivar2][ix][iy];
            mag = sqrt( u * u + v * v );
            fn[jb].u[ix][iy] = u;
            fn[jb].v[ix][iy] = v;
            fn[jb].value[ix][iy] = mag;
            }
         else
            {
            fn[jb].value[ix][iy] = Var[jb][ivar1][ix][iy];
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
 */

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
if (mirror_flag == 1)
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

sprintf (msg_text, "Nominal ranges...\n");
log_message (msg_text, 0);
sprintf (msg_text, "xstart = %e, xend = %e, xtick = %e\n", xstart, xend, xtick);
log_message (msg_text, 0);
sprintf (msg_text, "ystart = %e, yend = %e, ytick = %e\n", ystart, yend, ytick);
log_message (msg_text, 0);
sprintf (msg_text, "vstart = %e, vend = %e, vtick = %e\n", vstart, vend, vtick);
log_message (msg_text, 0);

if (fabs(vtick) < 1.0e-10 || fabs(xtick) < 1.0e-10 ||
    fabs(ytick) < 1.0e-10)
   {
   log_message ("Data looks strange; will not do contours!\n", 1);
   do_contours = 0;
   }

/*
 * Alter the nominal ranges if desired.
 */
if ( vrange_given )
   {
   vstart = vmin_given;
   vend   = vmax_given;
   vtick  = dv_given;
   vrange = vend - vstart;
   }
if ( xrange_given )
   {
   xstart = xmin_given;
   xend   = xmax_given;
   xtick  = dx_given;
   xrange = xend - xstart;
   }
if ( yrange_given )
   {
   ystart = ymin_given;
   yend   = ymax_given;
   ytick  = dy_given;
   yrange = yend - ystart;
   }


/*
 * Plot scales.
 */
xscale = xplotsize / xrange;
yscale = yplotsize / yrange;
if (true_shape == 1)
   {
   /* Choose the smallest scale */
   if (yscale < xscale) 
      xscale = yscale;
   else
      yscale = xscale;
   /* Choose the largest tic marks */
   if (ytick < xtick) 
      ytick = xtick;
   else
      xtick = ytick;
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

ssp_Init();
ssp_SetParameter( "orientation", orientation );
ssp_SetParameter( "device", output_type );
sprintf(pixels_per_mm_str, "%f.2", pixels_per_mm);
ssp_SetParameter( "pixels_per_mm", pixels_per_mm_str );
ssp_BeginPlot( "", PlotFileName );
ssp_SetOrigin( xorig, yorig );

ssp_SetLineThickness (pw_labels);
if ( add_labels )
   {
   ssp_PlotText (-15.0, yplotsize+22.0, 2.5, 0.0, DataFileName, 16);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotText (xx+5.0, yy, 2.5, 0.0, VarName[ivar1], 16);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotText (xx+5.0, yy, 3.0, 0.0, title, strlen(title));

   ssp_PlotText (-15.0, yplotsize+17.0, 2.5, 0.0, "x1 x2 dx", 8);
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

   ssp_PlotText (-15.0, yplotsize+12.0, 2.5, 0.0, "v1 v2 dv", 8);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, vstart, 2);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, vend, 2);
   ssp_Where (&xx, &yy, &sc, &ang);
   ssp_PlotENumber (xx+2.0, yy, 2.5, 0.0, vtick, 2);
   }

ssp_SetLineThickness (pw_axis);
xaxislength = xrange * xscale;
if (xplotsize > xaxislength) xaxislength = xplotsize;
if (draw_axes == 1) {
   ssp_PlotXAxis (0.0, 0.0, xaxislength, 0.0, xscale,
      xstart, xtick, VarName[0], -1);
} /* end if */
yaxislength = yrange * yscale;
if (yplotsize > yaxislength) yaxislength = yplotsize;
if (draw_axes == 1) {
   ssp_PlotYAxis (0.0, 0.0, yaxislength, 0.0, yscale,
      ystart, ytick, VarName[1], -1);
} /* end if */

ssp_SetClipLimits (0.0, 0.0, xplotsize, yplotsize);
ssp_SetClipOn ();

if ( draw_table == 1 && 
     ((do_contours && with_colours) || (do_filled && with_colours)) )
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
      xx = xplotsize+2.0;
      yy = 20.0 + 5.0 * i;
      ssp_PlotENumber (xx, yy, 2.5, 0.0, level, 3);
      }
   ssp_SetNoColour ();
   ssp_SetClipOn ();
   }

if ( do_filled && !with_colours )
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
      xx = xplotsize + 2.0;
      yy = 20.0 + 5.0 * i;
      ssp_PlotENumber (xx, yy, 2.5, 0.0, level, 3);
      }
   ssp_SetNoColour ();
   ssp_SetClipOn ();
   }

if ( do_contours )
   {
   /*
    * The contours themselves.
    */
   for (jb = 0; jb < nblock; ++jb)
      {
      ssp_SetLineThickness (pw_contours);
      for (level = vstart; level <= vend+0.001*vtick; level += vtick)
         {
         sprintf (txt, "Block[%d], Contour level = %f", jb, level);
         log_message (txt, 0);
         if ( with_colours )
            {
            hue = ssp_MapToColour (level, vstart, vend);
            bright = 1.0;
            sat = 1.0;
            ssp_SetHSBColour (hue, sat, bright);
            }
         ssp_ContourArray (&fn[jb], level);
         if ( with_colours ) ssp_SetNoColour();
         }
      ssp_SetLineThickness (pw_mesh);
      if ( do_edge ) ssp_PlotEdge ( &fn[jb] );
      }  /* for (jb... */
   }  /* if (icontour ... */

if ( do_filled )
   {
   /*
    * Fill the cells with colour or grey-shade.
    */
   for (jb = 0; jb < nblock; ++jb)
      {
      ssp_ColourFillArray ( &fn[jb], vstart, vend, !with_colours );
      ssp_SetNoColour();
      ssp_SetLineThickness (pw_mesh);
      if ( do_edge ) ssp_PlotEdge ( &fn[jb] );
      }  /* for (jb... */
   }  /* if (do_filled_colours... */

if ( do_vectors )
   {
   /*
    * Draw an array of arrows to represent the vector quantities.
    */
   for (jb = 0; jb < nblock; ++jb)
      {
      ssp_SetLineThickness (pw_vectors);
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
               log_message (txt, 0);
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

            if ( with_colours )
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

      if ( with_colours ) ssp_SetNoColour();
      ssp_SetLineThickness (pw_mesh);
      if ( do_edge ) ssp_PlotEdge ( &fn[jb] );
      }  /* for (jb... */
   }  /* if ( do_vectors ... */

if ( do_mesh )
   {
   ssp_SetLineThickness (pw_mesh);
   for (jb = 0; jb < nblock; ++jb) {
      ssp_PlotMesh ( &fn[jb], ixskip, iyskip );
   } /* end for */
   }

/*
 * For the mirror image,
 * invert the y-coordinate and redo the contours.
 */
if ( mirror_flag )
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

   if ( do_contours == 1)
      {
      /*
       * Redo the contours.
       */
      for (jb = 0; jb < nblock; ++jb)
         {
         ssp_SetLineThickness (pw_contours);
         for (level = vstart; level <= vend+0.001*vtick; level += vtick)
            {
            sprintf (txt, "Block[%d], Contour level = %f", jb, level);
            log_message (txt, 0);
            if ( with_colours )
               {
               hue = ssp_MapToColour (level, vstart, vend);
               bright = 1.0;
               sat = 1.0;
               ssp_SetHSBColour (hue, sat, bright);
               }
            ssp_ContourArray (&fn[jb], level);
            if ( with_colours ) ssp_SetNoColour();
            }
         ssp_SetLineThickness (pw_mesh);
         if ( do_edge ) ssp_PlotEdge ( &fn[jb] );
         }
      }  /* if ( do_contours ... */

   if ( do_filled )
      {
      /*
       * Fill the cells with colour or grey-shade.
       */
      for (jb = 0; jb < nblock; ++jb)
         {
         ssp_ColourFillArray ( &fn[jb], vstart, vend, !with_colours );
         ssp_SetNoColour();
         ssp_SetLineThickness (pw_mesh);
         if ( do_edge ) ssp_PlotEdge ( &fn[jb] );
         }  /* for (jb... */
      }  /* if ( do_filled... */

   if ( do_vectors )
      {
      /*
       * Draw an array of arrows to represent the vector quantities.
       */
      for (jb = 0; jb < nblock; ++jb)
         {
         ssp_SetLineThickness (pw_vectors);
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
                  log_message (txt, 0);
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

               if ( with_colours )
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

         if ( with_colours ) ssp_SetNoColour();
         ssp_SetLineThickness (pw_mesh);
         if ( do_edge ) ssp_PlotEdge ( &fn[jb] );
         }  /* for (jb... */
      }  /* if ( do_vectors ... */

   if ( do_mesh )
      {
      ssp_SetLineThickness (pw_mesh);
      for (jb = 0; jb < nblock; ++jb) {
         ssp_PlotMesh ( &fn[jb], ixskip, iyskip );
      } /* end for */
      }
   }  /* if (mirror_flag...  */


/*
 * Read and act on any trailing geometry records.
 * This will be done on the first pass so, any subsequent
 * plots will miss out.
 * I guess that this can be improved by storing the segment
 * data points and replaying for each plot.
 */
ssp_SetLineThickness (pw_geometry);
finished = 0;
while (finished != 1)
   {
   /* Attempt to read 4 coordinates. */
   if ( fscanf(dfp, "%lf %lf %lf %lf", &x1, &y1, &x2, &y2) == EOF )
      finished = 1;
   /* printf ("x1 = %f, y1 = %f, x2 = %f, y2 = %f\n", x1, y1, x2, y2); */
   if (finished != 1)
      {
      if (mirror_flag == 1)
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
ssp_Final();

fclose(dfp);

for (jb = 0; jb < nblock; ++jb) {
   ssp_FreeArray (&fn[jb]);
}

/* 
 * For each block, free memory. 
 */
for (jb = 0; jb < nblock; ++jb) {
    printf ("Freeing memory for block[%d]...\n", jb);
    for (ivar = 0; ivar < nvar; ++ivar) {
	for (ix = 0; ix < nnx[jb]; ++ix) {
	    free(Var[jb][ivar][ix]);
	}
	free(Var[jb][ivar]);
      }
   free(Var[jb]);
   }   /* end of jb loop  */

close_log_file ();

return 0;
}  /* end of main() */


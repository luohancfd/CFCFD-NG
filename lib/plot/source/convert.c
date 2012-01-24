/** \file convert.c
 * \ingroup plot
 * \brief  Convert a plot-move file from plotxy to
 *	   actual plotter commands using ssp routines -- no longer used.
 * 
 *         The input file contains lines of the form
 *         ipen, x, y  (option 0) or x, y, ipen  (option 1)
 *
 * \author PA Jacobs
 *
 * \version 1.0  :
 * \version 1.1  :   -Nov-92 : add an option to swap the order of the input data
 * \version 1.2  : 11-Jan-94 : Make the addition of the date/time optional.
 * \version 1.3  : 31-Aug-94 : Fix the line clipping to work for portrait mode
 *                 plotting.  Take advantage of the more compact
 *                 postscript.
 */

#include <stdio.h>
#include <time.h>

#include "compiler.h"
#include "geom.h"
#include "ssp.h"

#if (STDLIBH)
#  include <stdlib.h>
#endif


int main ()
{
int    ipen, count, climit;
int    pen_option, date_option, portrait_option;
double x, y, ang, fact;
char   in_name[32], out_name[32];
char   string[132];
FILE   *in_file;
time_t bintime;

printf ("\n");
printf ("Convert a generic plot-move file to plotter commands.\n");
printf ("Rev. 1.3, 31-Aug-94\n");
printf ("\n");

/*
 * My default order is to put the pen status first followed by the
 * coordinates.  This allows later extension to plotting character
 * strings and the like.
 * The old Method-of-Characteristics program (MOC) would generate
 * a plot-move file with the pen status last.
 */
printf ("Order 0=(ipen,x,y) 1=(x,y,ipen)   : ");
scanf ("%d", &pen_option);
printf ("Add date to RHS/Top of page (1/0) : ");
scanf ("%d", &date_option);
printf ("Orientation 0=landscape 1=portrait: ");
scanf ("%d", &portrait_option);

/*
 * Open the plot-move file.
 */
printf ("input (plot-move) file: ");
scanf ("%s", in_name);
if ( (in_file = fopen (in_name, "r")) == NULL )
   {
   printf ("Could not open file : %s\n", in_name);
   exit (-1);
   }

/*
 * Initialize the drawing
 */
printf ("output file           : ");
scanf ("%s", out_name);

if ( portrait_option == 1 )
   {
   ssp_BeginPlot (0, 1, out_name, "ps", "portrait");
   ssp_SetClipLimits (0.0, 0.0, 210.0, 300.0);
   }
else
   {
   ssp_BeginPlot (0, 1, out_name, "ps", "landscape");
   ssp_SetClipLimits (0.0, 0.0, 300.0, 210.0);
   }

ssp_SetLineThickness (0.18);
ssp_SetClipOn ();
ssp_SetOrigin (0.0,0.0);

if ( date_option )
   {
   time (&bintime);
   if ( portrait_option == 1 )
      {
      /* Add date to the top of the page. */
      ssp_PlotText (5.0, 270.0, 3.0, 0.0, out_name, 32);
      ssp_Where (&x, &y, &fact, &ang);
      ssp_PlotText (x+5.0, y, 3.0, 0.0, ctime(&bintime), 26);
      }
   else
      {
      /* Add date and time to RHS of page. */
      ssp_PlotText (270.0, 5.0, 3.0, 90.0, out_name, 32);
      ssp_Where (&x, &y, &fact, &ang);
      ssp_PlotText (x, y+5.0, 3.0, 90.0, ctime(&bintime), 26);
      }
   }

/*
 * Process the plotting moves.
 */
count = 0;
climit = 1000;
printf ("\n");
while ( !feof(in_file) )
   {
   /* Read a line until we hit the end of the file. */
   if ( NULL == fgets (string, 132, in_file) ) break;

   if (pen_option == 0)
      {
      /*
       * Try to read the default format.
       * If we fail, go on to the next line in the plot-move file.
       */
      if (0 == sscanf (string, "%d %lf %lf", &ipen, &x, &y)) continue;
      }
   else
      {
      if (0 == fscanf (in_file, "%lf %lf %d", &x, &y, &ipen)) continue;
      }
   if (ipen == 2) ssp_Plot (x, y);
   if (ipen == 3) ssp_Move (x, y);

   ++count;
   if (count >= climit)
      {
      /* Write out a star every climit moves. */
      printf ("*");
      fflush (stdout);
      count = 0;
      }
   }
printf ("\n");

ssp_EndPlot ();


return (0);
}  /* end of convert.c */


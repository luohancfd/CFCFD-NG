/* showplot.c  -- Display a plot-move file from plotxy on the screen.
 *                The input file contains lines of the form
 *                ipen, x, y  (option 0) or
 *                x, y, ipen  (option 1)
 * -------------------------------------------------------------
 *
 * Written by ... P. A. Jacobs
 * ----------     38 Ellington St
 *                Tarragindi, Qld 4121.
 *
 * Version...
 * -------
 * 1.0  : 27-Sep-94 : Adapted from convert.c
 *
 */

#include <stdio.h>
#include <time.h>

#include "compiler.h"
#include "geom.h"
#include "ssp.h"

#if (STDLIBH)
#  include <stdlib.h>
#endif


main ()
{
int    ipen;
int    pen_option, portrait_option;
double x, y, ang, fact;
char   in_name[32];
char   string[132];
FILE   *in_file;

printf ("\n");
printf ("Display a generic plot-move file.\n");
printf ("Rev. 1.0, 27-Sep-94\n");
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
if ( portrait_option == 1 )
   {
   ssp_BeginPlot (1, 0, "junk", "ps", "portrait");
   ssp_SetClipLimits (0.0, 0.0, 210.0, 300.0);
   }
else
   {
   ssp_BeginPlot (1, 0, "junk", "ps", "landscape");
   ssp_SetClipLimits (0.0, 0.0, 300.0, 210.0);
   }

ssp_SetLineThickness (0.18);
ssp_SetClipOn ();
ssp_SetOrigin (0.0,0.0);
if ( portrait_option == 1 )
   {
   /* Scale the plot so that it all fits. */
   ssp_Factor ( 210.0/300.0 );
   }


/*
 * Process the plotting moves.
 */
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
   }

ssp_EndPlot ();


return (0);
}  /* end of convert.c */


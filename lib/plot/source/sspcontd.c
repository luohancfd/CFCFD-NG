/* sspcontd.c  -- exercise contouring routines */

#include "compiler.h"
#include <math.h>
#include <stdio.h>
#include "geom.h"
#include "ssp.h"

main ()
{
int    i,j;
double x, y, dx, dy, level, d_level;
struct ssp_data_point t1, t2, t3, t4;
struct ssp_data_array fn;
int    ix, iy, ixdim, iydim;
int    option;
char   txtstring[72];

printf ("Ssp demonstration program...\n");
printf ("0 = contour single elements, 1 = contour array\n");
scanf ("%d", &option);

/* Initialize the drawing */
ssp_BeginPlot (1, 1, "contour.ps", "ps", "landscape");

ssp_PlotText (30.0, 152.0, 8.0, 0.0,
	      "Contour Plotting Demonstration.", 31);

ssp_SetOrigin (30.0,30.0);              /* set origin on page */
ssp_Factor (30.0);                      /* suitable scale for data */

/* Axes ... arrows with solid heads, tic marks etc */
ssp_Message ("Draw Axes");

ssp_SetLineThickness (0.5);
ssp_PlotXAxis (0.0, 0.0, 6.0, 0.0, 1.0, 0.0, 1.0, "X-Axis", -1);
ssp_PlotYAxis (0.0, 0.0, 4.0, 0.0, 1.0, 0.0, 1.0, "Y-Axis", -1);

ssp_SetClipLimits (0.0, 0.0, 6.0, 4.0);
ssp_SetClipOn ();

if (option == 0)
   {
   /*
    * sample data
    */
   ssp_Message ("Sample data for single quadrilateral.");

   t1.x = 2.5; t1.y = 1.0; t1.value = 0.5;
   t2.x = 2.0; t2.y = 3.0; t2.value = 1.0;
   t3.x = 1.0; t3.y = 0.5; t3.value = 3.0;

   ssp_SetLineThickness (0.35);

   ssp_Move (t1.x, t1.y);
   ssp_Plot (t2.x, t2.y);
   ssp_Plot (t3.x, t3.y);
   ssp_Plot (t1.x, t1.y);

   ssp_PlotSymbol (t1.x, t1.y, 0.12, CIRCLE_SYM);
   ssp_PlotSymbol (t2.x, t2.y, 0.12, DIAMOND_SYM);
   ssp_PlotSymbol (t3.x, t3.y, 0.12, DEL_SYM);

   ssp_ContourTriangle (&t1, &t2, &t3, 1.2);
   ssp_ContourTriangle (&t1, &t2, &t3, 1.5);
   ssp_ContourTriangle (&t1, &t2, &t3, 1.8);

   /*
    * Try out quadrilateral contouring.
    */
   t1.x = 5.0; t1.y = 1.0; t1.value = 0.5;
   t2.x = 4.5; t2.y = 3.0; t2.value = 1.0;
   t3.x = 3.3; t3.y = 3.5; t3.value = 3.0;
   t4.x = 3.0; t4.y = 0.5; t4.value = 0.45;

   ssp_Move (t1.x, t1.y);
   ssp_Plot (t2.x, t2.y);
   ssp_Plot (t3.x, t3.y);
   ssp_Plot (t4.x, t4.y);
   ssp_Plot (t1.x, t1.y);

   ssp_ContourQuad (&t1, &t2, &t3, &t4, 1.2);
   ssp_ContourQuad (&t1, &t2, &t3, &t4, 1.5);
   ssp_ContourQuad (&t1, &t2, &t3, &t4, 2.5);
   ssp_ContourQuad (&t1, &t2, &t3, &t4, 0.75);
   }
else if (option == 1)
   {
   ssp_Message ("Set up data for an array.");

   ixdim = 20;
   iydim = 15;
   dx = 0.2;
   dy = 0.15;
   ssp_AllocateArray (&fn, ixdim, iydim);
   for (ix = 0; ix < ixdim; ++ix)
      {
      for (iy = 0; iy < iydim; ++iy)
	 {
	 x = 0.5 + dx * ix;
	 y = dy * iy + 0.3 * sin(x);
	 fn.x[ix][iy] = x;
	 fn.y[ix][iy] = y;
	 x = x - 2.0;
	 y = y - 1.0;
	 fn.value[ix][iy] = sin(x*x + y*y) / (x*x + y*y + 0.1);
	 }
      }

   fn.ixmin = 0;
   fn.ixmax = ixdim - 1;
   fn.iymin = 0;
   fn.iymax = iydim - 1;

   ssp_ArrayExtremes (&fn);
   d_level = (fn.vmax - fn.vmin) / 16.0;

   for (level = fn.vmin + 0.5 * d_level;
	level <= fn.vmax - 0.5 * d_level;
	level += d_level)
      {
      sprintf (txtstring, "Contour level = %f", level);
      ssp_Message (txtstring);
      ssp_SetLineThickness (0.35);
      ssp_ContourArray (&fn, level);
      }

   ssp_SetLineThickness (0.10);
   ssp_PlotMesh (&fn);
   }

ssp_EndPlot ();    /* finish up screen demo */

return (0);
}  /* end of sspd */


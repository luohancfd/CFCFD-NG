/** \file ssp.c
 * \ingroup plot
 * \brief Device independent routines for Simple Scientific Plotting.
 *
 * A set of Calcomp-like plotting routines for C.
 * This file contains the device independent routines; the
 * kernel routines may be found in sspkern.c.
 *
 * All plotting is done in a section of the screen with
 * proportions similar to an A4 page.  Plotting dimensions
 * are in mm with the x-axis aligned with the long edge
 * of the page.
 *
 * Routines supported in this file...
 * \verbatim
 *    ssp_PlotNumber (x, y, size, directn, value, dec_places);
 *    ssp_PlotENumber (x, y, size, directn, value, dec_places);
 *    ssp_PlotSymbol (x, y, size, option);
 *    ssp_Arrow (x1, y1, x2, y2, headlen, width, option);
 *    ssp_PlotXAxis (x0, y0, length, angle, scale,
 *	 	     first_value, delta_value, title, option);
 *    ssp_PlotYAxis (x0, y0, length, angle, scale,
 *		     first_value, delta_value, title, option);
 *
 * Notes ...
 * -----
 * (1) Graphics window definition:
 *     Map an a4 page to this window in mm units.
 *     The x-axis is along the long edge of the page.
 * (2) Definitions of the plotting units:
 *     Page Coordinates : these are the units in which marks are made
 *                        on the screen or page.
 *     Plot Coordinates : these are the units in which the user specifies
 *                        the plotting moves.
 * \endverbatim
 *
 * \author PA Jacobs
 *
 * \version 1.0, September 1988
 * \version 2.0, October   1991
 * \version 2.1, 10-Nov-93 : allow wider characters in the axis functions
 *               Assume that width = 0.8 * height rather than 2/3.
 * \version 2.11,20-Dec-97 : PlotXAxis no longer sets its own pen thickness
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

/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_PlotNumber (double x, double y,
		    double size, double directn,
		    double value,
		    int dec_places)

#else

int ssp_PlotNumber (x, y, size, directn, value, dec_places)
double x, y, size, directn, value;
int    dec_places;

#endif

/*
 * Purpose...
 * -------
 * Write out a decimal number onto the plotting page.
 * This is accomplished by first converting it to a
 * string representation and then calling ssp_PlotText().
 *
 * Input...
 * -----
 * x, y    : location of the lower left corner of the number
 * size    : height of the characters
 * directn : angle of test direction wrt the x-axis
 * value   : actual value of the number to be written
 * dec_places : number of decimal places to include
 *
 */

{  /* begin ssp_PlotNumber() */
char string[32];
int  count, i;

/*
 * First, write out the number with the correct number
 * of decimal places.
 */
if (dec_places < 0) dec_places = 0;
if (dec_places > 10) dec_places = 10;

switch (dec_places) {
   case 0 : {
      sprintf (string, "%10.0f", value);
      break;
      }
   case 1 : {
      sprintf (string, "%11.1f", value);
      break;
      }
   case 2 : {
      sprintf (string, "%12.2f", value);
      break;
      }
   case 3 : {
      sprintf (string, "%13.3f", value);
      break;
      }
   case 4 : {
      sprintf (string, "%14.4f", value);
      break;
      }
   case 5 : {
      sprintf (string, "%15.5f", value);
      break;
      }
   case 6 : {
      sprintf (string, "%16.6f", value);
      break;
      }
   case 7 : {
      sprintf (string, "%17.7f", value);
      break;
      }
   case 8 : {
      sprintf (string, "%18.8f", value);
      break;
      }
   case 9 : {
      sprintf (string, "%19.9f", value);
      break;
      }
   case 10 : {
      sprintf (string, "%20.10f", value);
      break;
      }
   }  /* end of switch... */

/*
 * Count the number of leading blanks and
 * shuffle characters to the left.
 */
count = 0;
while (string[count] == ' ' && count < 32) ++count;
if (count > 0)
   {
   for (i = 0; i + count < 32; ++i)
      string[i] = string[i + count];
   }

/*
 * Count the number of characters and then send to
 * the text plotting routine.
 */
count = 0;
while (string[count] != 0 && count < 32) ++count;
ssp_PlotText (x, y, size, directn, string, count);

return (0);
}  /* end of ssp_PlotNumber() */

/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_PlotENumber (double x, double y,
		     double size, double directn,
		     double value,
		     int dec_places)

#else

int ssp_PlotENumber (x, y, size, directn, value, dec_places)
double x, y, size, directn, value;
int    dec_places;

#endif

/*
 * Purpose...
 * -------
 * Write out a decimal number onto the plotting page
 * using EXPONENTIAL format.
 * This is accomplished by first converting it to a
 * string representation and then calling ssp_PlotText().
 *
 * Input...
 * -----
 * x, y    : location of the lower left corner of the number
 * size    : height of the characters
 * directn : angle of test direction wrt the x-axis
 * value   : actual value of the number to be written
 * dec_places : number of decimal places to include
 *
 */

{  /* begin ssp_PlotENumber() */
char string[32];
int  count, i;

/*
 * First, write out the number with the correct number
 * of decimal places.
 */
if (dec_places < 0) dec_places = 0;
if (dec_places > 10) dec_places = 10;

switch (dec_places) {
   case 0 : {
      sprintf (string, "%10.0e", value);
      break;
      }
   case 1 : {
      sprintf (string, "%11.1e", value);
      break;
      }
   case 2 : {
      sprintf (string, "%12.2e", value);
      break;
      }
   case 3 : {
      sprintf (string, "%13.3e", value);
      break;
      }
   case 4 : {
      sprintf (string, "%14.4e", value);
      break;
      }
   case 5 : {
      sprintf (string, "%15.5e", value);
      break;
      }
   case 6 : {
      sprintf (string, "%16.6e", value);
      break;
      }
   case 7 : {
      sprintf (string, "%17.7e", value);
      break;
      }
   case 8 : {
      sprintf (string, "%18.8e", value);
      break;
      }
   case 9 : {
      sprintf (string, "%19.9e", value);
      break;
      }
   case 10 : {
      sprintf (string, "%20.10e", value);
      break;
      }
   }  /* end of switch... */

/*
 * Count the number of leading blanks and
 * shuffle characters to the left.
 */
count = 0;
while (string[count] == ' ' && count < 32) ++count;
if (count > 0)
   {
   for (i = 0; i + count < 32; ++i)
      string[i] = string[i + count];
   }

/*
 * Count the number of characters and then send to
 * the text plotting routine.
 */
count = 0;
while (string[count] != 0 && count < 32) ++count;
ssp_PlotText (x, y, size, directn, string, count);

return (0);
}  /* end of ssp_PlotENumber() */

/*-------------------------------------------------------------------*/

#if (PROTO)

int ssp_PlotSymbol (double x, double y, double size, int option)

#else

int ssp_PlotSymbol (x, y, size, option)
double x, y, size;
int    option;

#endif

/*
 * Purpose...
 * -------
 * Draw a centred symbol at (x, y).
 *
 * Input...
 * -----
 * x, y    : coordinates of the centre of the symbol
 * size    : diameter of the symbol (plot units)
 * option  : index to the particular symbol
 *           (see the header file "ssp.h" for a list)
 */

{
double x1, y1, x2, y2;
double dtheta, theta, radius;
int    i;

x1 = x - size / 2.0;      /* Bottom left point */
y1 = y - size / 2.0;
x2 = x + size / 2.0;      /* Top right point   */
y2 = y + size / 2.0;

switch (option)
  {
  case CROSS_SYM : {  /* Diagonal cross */
      ssp_Move (x1, y1);
      ssp_Plot (x2, y2);
      ssp_Move (x1, y2);
      ssp_Plot (x2, y1);
      break;
      }
  case PLUS_SYM : {  /* Straight cross */
      ssp_Move (x, y1);
      ssp_Plot (x, y2);
      ssp_Move (x1, y);
      ssp_Plot (x2, y);
      break;
      }
  case DIAMOND_SYM : {  /* Diamond */
      ssp_Move (x, y);
      ssp_Plot (x, y2);
      ssp_Plot (x2, y);
      ssp_Plot (x, y1);
      ssp_Plot (x1, y);
      ssp_Plot (x, y2);
      break;
      }
  case TRIANGLE_SYM : {  /* Triangle */
      ssp_Move (x, y);
      ssp_Plot (x, y2);
      ssp_Plot (x2, y1);
      ssp_Plot (x1, y1);
      ssp_Plot (x, y2);
      break;
      }
  case DEL_SYM : {  /* Inverted triangle */
      ssp_Move (x, y);
      ssp_Plot (x, y1);
      ssp_Plot (x2, y2);
      ssp_Plot (x1, y2);
      ssp_Plot (x, y1);
      break;
      }
  case SQUARE_SYM : {  /* Square */
      ssp_Move (x, y);
      ssp_Move (x, y2);
      ssp_Plot (x1, y2);
      ssp_Plot (x1, y1);
      ssp_Plot (x2, y1);
      ssp_Plot (x2, y2);
      ssp_Plot (x, y2);
      break;
      }
  case CIRCLE_SYM : {  /* Rough circle */
     radius = size / 2.0;
     dtheta = ssp_PI / 10.0;
     ssp_Move (x, y);
     ssp_Plot (x, y+radius);
     for (i = 0; i <= 20; ++i)
	{
	theta = ssp_PI * 0.5 + dtheta * i;
	ssp_Plot (x + radius * cos(theta), y + radius * sin(theta));
	}
     break;
     }
  }    /* case option */
ssp_Move (x, y);   /* leave at centre of symbol */

return (0);
}  /* end of ssp_PlotSymbol() */
/*-------------------------------------------------------------------------*/

#if (PROTO)

int ssp_Arrow (double x1, double y1, double x2, double y2,
	       double HeadLen, double HeadWidth, int option)

#else

int ssp_Arrow (x1, y1, x2, y2, HeadLen, HeadWidth, option)
double x1, y1, x2, y2;
double HeadLen, HeadWidth;
int    option;

#endif

/*
 * Purpose...
 * -------
 * Draw an arrow with its tail at x1, y1 and head at x2, y2.
 *
 * Input...
 * -----
 * x1, y1   : coordinates of the tail (plot units)
 * x2, y2   : coordinates of the head (plot units)
 * HeadLen  : length of the arrowhead (plot units)
 * HeadWidth : width of the arrowhead (plot units)
 * option   : type of arrowhead
 *            0 = open head
 *            1 = closed head
 *            2 = filled-in head
 *
 * Revisions...
 * ---------
 * 12-Jan-94 : Trap zero-length arrows.
 *
 */

{
double x3, y3, x4, y4, x5, y5, xx, yy, delx, dely;
double xunit, yunit, xnorm, ynorm, Len;
int    nstep, ii;

/*
 * Draw shaft
 */
ssp_Move (x1, y1);
ssp_Plot (x2, y2);
Len = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));   /* line length */

/*
 * Check for zero-length arrow.
 */
if (Len < 1.0e-10) return (0);

/*
 * unit vector components
 */
xunit = (x2 - x1) / Len;
yunit = (y2 - y1) / Len;

/*
 * Prepare to add the arrowhead
 */
x3 = x2 - HeadLen * xunit;
y3 = y2 - HeadLen * yunit;
/*
 * unit vector normal to shaft
 */
xnorm = -yunit;
ynorm = xunit;
/*
 * positions for the wide part of the arrowhead
 */
x4 = x3 + HeadWidth * xnorm / 2.0;
y4 = y3 + HeadWidth * ynorm / 2.0;
x5 = x3 - HeadWidth * xnorm / 2.0;
y5 = y3 - HeadWidth * ynorm / 2.0;

switch (option)
  {
  case 0 : {  /* open head */
      ssp_Move (x4, y4);
      ssp_Plot (x2, y2);
      ssp_Plot (x5, y5);
      break;
      }
  case 1 : {  /* closed head */
      ssp_Move (x4, y4);
      ssp_Plot (x2, y2);
      ssp_Plot (x5, y5);
      ssp_Plot (x4, y4);
      break;
      }
  case 2 : {  /* Filled in Head */
      nstep = 20;
      delx = (x4 - x5) / nstep;
      dely = (y4 - y5) / nstep;
      ssp_Move (x4, y4);       /* draw a closed head */
      ssp_Plot (x2, y2);
      ssp_Plot (x5, y5);
      ssp_Plot (x4, y4);
      for (ii = 1; ii <= nstep; ++ii)      /* fill in the head */
        {
        xx = x5 + delx * ii;
        yy = y5 + dely * ii;
	ssp_Move (xx, yy);
	ssp_Plot (x2, y2);
        }
      break;
      }  /* filled in head */
  }  /* case option */

/*
 * Move to the tip of the arrowhead.
 */
ssp_Move (x2, y2);

return (0);
}   /* end of ssp_Arrow() */

/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_PlotXAxis (double x0, double y0,
		   double length, double angle,
		   double scale,
		   double first_value, double delta_value,
		   char *title,
		   int option)

#else

int ssp_PlotXAxis (x0, y0, length, angle, scale,
		   first_value, delta_value, title, option)
double x0, y0;
double length, angle, scale;
double first_value, delta_value;
char   *title;
int    option;

#endif

/*
 * Purpose...
 * -------
 * Draw an x-axis and title it.
 *
 * Input...
 * -----
 * x0, y0  : starting point for the axis in current plotting units
 * length  : length of the axis in current plotting units
 * angle   : angle of the axis with the frame of reference (degree)
 * scale   : scale factor required to multiply data units to get
 *           plotting units
 * first_value : data value at the start of the axis
 * delta_value : jump in data value between tick marks
 * title   : pointer to the title string
 * option  : +1 : put titles on the positive side of the axis
 *           -1 : put titles on the negative side of the axis
 *
 * Revisions...
 * ---------
 * 20-Dec-97 : Use the pen thickness set outside the routing.
 *
 */

{  /* begin ssp_PlotXAxis() */
double number_height, title_height, tick_length;
double current_factor, current_angle;
double title_length, title_offset;
double number_length, number_offset;
int    count, nexp;
double sine, cosine;
double x1, y1, x2, y2, dist, tick_value;
double max_label, fnexp;

/*
 * Default parameters in plotting units (mm).
 */
number_height = 3.5;
title_height = 5.0;
if (option == 1)
   {
   /*
    * Tic marks on the positive side of the axis
    * Labels on the negative side of the axis.
    */
   tick_length = 2.0;
   number_offset = -6.0;
   title_offset = -13.0;
   }
else
   {
   /*
    * Tic marks on the negative side of the axis
    * Labels on the negative side of the axis.
    */
   tick_length = -2.0;
   number_offset = -8.0;
   title_offset = -15.0;
   }

/*
 * Present value for the global scale factor.
 */
ssp_Where (&x1, &y1, &current_factor, &current_angle);

/*
 * Convert to plotting units.
 */
number_height /= current_factor;
title_height /= current_factor;
tick_length /= current_factor;
number_offset /= current_factor;
title_offset /= current_factor;

/*
 * Angle of the X-axis.
 */
cosine = cos(ssp_PI * angle / 180.0);
sine   = sin(ssp_PI * angle / 180.0);

/*
 * Length of title.
 */
count = 0;
while (title[count] != 0 && count < 32) ++count;
title_length = 0.8 * title_height * count;

/*
 * Assume that numbers drawn on the axis have 3 characters.
 */
number_length = 0.8 * number_height * 3;

/*
 * Work out the value of the largest label
 * The factor 1.001 should avoid rounding problems when putting
 * down the tic marks.
 */
tick_value = first_value;
max_label = fabs(tick_value);
dist = 0;
while (dist <= length * 1.001)
   {
   dist += delta_value * scale;
   tick_value += delta_value;
   }
if (fabs(tick_value) > max_label) max_label = fabs(tick_value);

/*
 * Now compute an exponent so that the maximum label is less
 * than or equal to 100.0 but greater than 1
 */
nexp = 0;
while (max_label > 100.0)
   {
   --nexp;
   max_label /= 10.0;
   }
while (max_label < 1.0)
   {
   ++nexp;
   max_label *= 10.0;
   }
fnexp = (double) nexp;

/*
 * Draw the axis itself.
 */
/* ssp_SetLineThickness (0.5);*/
x1 = x0 + cosine * length;
y1 = y0 + sine * length;
ssp_Move (x0, y0);
ssp_Plot (x1, y1);

/*
 * Draw the tick marks and label them with scaled numbers.
 */
dist = 0;
tick_value = first_value;
while (dist <= length * 1.001)
   {
   x1 = x0 + cosine * dist;
   y1 = y0 + sine * dist;
   x2 = x1 - sine * tick_length;
   y2 = y1 + cosine * tick_length;
   /* ssp_SetLineThickness (0.5); */
   ssp_Move (x1, y1);
   ssp_Plot (x2, y2);

   /* ssp_SetLineThickness (0.2); */
   x2 = x1 - sine * number_offset - cosine * 0.5 * number_length;
   y2 = y1 + cosine * number_offset - sine * 0.5 * number_length;
   ssp_PlotNumber (x2, y2, number_height, angle,
		   tick_value * pow(10.0,fnexp), 2);

   dist += delta_value * scale;
   tick_value += delta_value;
   }

/*
 * Add the title at the midpoint.
 */
x1 = x0 + cosine * 0.5 * length - sine * title_offset
     - cosine * 0.5 * title_length;
y1 = y0 + sine * 0.5 * length + cosine * title_offset
     - sine * 0.5 * title_length;
/* ssp_SetLineThickness (0.2); */
ssp_PlotText (x1, y1, title_height, angle, title, count);
/*
 * Add the exponent if it is not zero.
 */
ssp_Where (&x1, &y1, &current_factor, &current_angle);
if (nexp != 0)
   {
   ssp_Where (&x1, &y1, &current_factor, &current_angle);
   ssp_PlotText (x1, y1, title_height, angle, " * 10", 5);
   ssp_Where (&x1, &y1, &current_factor, &current_angle);
   ssp_PlotNumber (x1, y1 + title_height * 0.5,
		   title_height * 0.7, angle, fnexp, 0);
   }

return (0);
}  /* end of ssp_PlotXAxis() */

/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_PlotYAxis (double x0, double y0,
		   double length, double angle,
		   double scale,
		   double first_value, double delta_value,
		   char *title,
		   int option)

#else

int ssp_PlotYAxis (x0, y0, length, angle, scale,
		   first_value, delta_value, title, option)
double x0, y0;
double length, angle, scale;
double first_value, delta_value;
char   *title;
int    option;

#endif

/*
 * Purpose...
 * -------
 * Draw a y-axis and title it.
 *
 * Input...
 * -----
 * x0, y0  : starting point for the axis in current plotting units
 * length  : length of the axis in current plotting units
 * angle   : angle of the axis with the frame of reference (degree)
 * scale   : scale factor required to multiply data units to get
 *           plotting units
 * first_value : data value at the start of the axis
 * delta_value : jump in data value between tick marks
 * title   : pointer to the title string
 * option  : +1 : put titles on the positive side of the axis
 *           -1 : put titles on the negative side of the axis
 *
 * Revisions...
 * ---------
 * 20-Dec-97 : use alreay set line thickness
 *
 */

{  /* begin ssp_PlotYAxis() */
double number_height, title_height, tick_length;
double current_factor, current_angle;
double title_length, title_offset;
double number_length, number_offset;
int    count, nexp;
double sine, cosine;
double x1, y1, x2, y2, dist, tick_value;
double max_label, fnexp;

/*
 * Default parameters in plotting units (mm).
 */
number_height = 3.5;
title_height = 5.0;
if (option == 1)
   {
   /*
    * Tic marks on the positive side of the axis
    * Labels on the negative side of the axis.
    */
   tick_length = -2.0;
   number_offset = 6.0;
   title_offset = 13.0;
   }
else
   {
   /*
    * Tic marks on the negative side of the axis
    * Labels on the negative side of the axis.
    */
   tick_length = 2.0;
   number_offset = 8.0;
   title_offset = 15.0;
   }

/*
 * Present value for the global scale factor.
 */
ssp_Where (&x1, &y1, &current_factor, &current_angle);

/*
 * Convert to plotting units.
 */
number_height /= current_factor;
title_height /= current_factor;
tick_length /= current_factor;
number_offset /= current_factor;
title_offset /= current_factor;

/*
 * Angle of the Y-axis.
 */
cosine = cos(ssp_PI * (angle + 90.0) / 180.0);
sine   = sin(ssp_PI * (angle + 90.0) / 180.0);

/*
 * Length of title.
 */
count = 0;
while (title[count] != 0 && count < 32) ++count;
title_length = 0.8 * title_height * count;

/*
 * Assume that numbers drawn on the axis have 3 characters.
 */
number_length = 0.8 * number_height * 3;

/*
 * Work out the value of the largest label
 */
tick_value = first_value;
max_label = fabs(tick_value);
dist = 0;
while (dist <= length * 1.001)
   {
   dist += delta_value * scale;
   tick_value += delta_value;
   }
if (fabs(tick_value) > max_label) max_label = fabs(tick_value);

/*
 * Now compute an exponent so that the maximum label is less
 * than or equal to 100.0 but greater than 1
 */
nexp = 0;
while (max_label > 100.0)
   {
   --nexp;
   max_label /= 10.0;
   }
while (max_label < 1.0)
   {
   ++nexp;
   max_label *= 10.0;
   }
fnexp = (double) nexp;

/*
 * Draw the axis itself.
 */
/* ssp_SetLineThickness (0.5); */
x1 = x0 + cosine * length;
y1 = y0 + sine * length;
ssp_Move (x0, y0);
ssp_Plot (x1, y1);

/*
 * Draw the tick marks and label them with scaled numbers.
 * The factor 1.001 should avoid rounding problems when putting
 * down the tic marks.
 */
dist = 0;
tick_value = first_value;
while (dist <= length * 1.001)
   {
   x1 = x0 + cosine * dist;
   y1 = y0 + sine * dist;
   x2 = x1 - sine * tick_length;
   y2 = y1 + cosine * tick_length;
   /* ssp_SetLineThickness (0.5); */
   ssp_Move (x1, y1);
   ssp_Plot (x2, y2);

   x2 = x1 - sine * number_offset - cosine * 0.5 * number_length;
   y2 = y1 + cosine * number_offset - sine * 0.5 * number_length;
   /* ssp_SetLineThickness (0.2); */
   ssp_PlotNumber (x2, y2, number_height, (angle+90.0),
		   tick_value * pow(10.0,fnexp), 2);

   dist += delta_value * scale;
   tick_value += delta_value;
   }

/*
 * Add the title at the midpoint.
 */
x1 = x0 + cosine * 0.5 * length - sine * title_offset
     - cosine * 0.5 * title_length;
y1 = y0 + sine * 0.5 * length + cosine * title_offset
     - sine * 0.5 * title_length;
/* ssp_SetLineThickness (0.2); */
ssp_PlotText (x1, y1, title_height, (angle+90.0), title, count);
/*
 * Add the exponent if it is not zero.
 */
ssp_Where (&x1, &y1, &current_factor, &current_angle);
if (nexp != 0)
   {
   ssp_Where (&x1, &y1, &current_factor, &current_angle);
   ssp_PlotText (x1, y1, title_height, (angle+90.0), " * 10", 5);
   ssp_Where (&x1, &y1, &current_factor, &current_angle);
   ssp_PlotNumber (x1 - (0.5 * title_height), y1,
		   title_height * 0.7, (angle+90.0), fnexp, 0);
   }

return (0);
}  /* end of ssp_PlotYAxis() */

/*-----------------------------------------------------------------*/

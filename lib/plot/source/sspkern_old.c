/* sspkern.c
 *
 * Simple Scientific Plotting Routines: Kernel Routines
 *
 */

/*
 * Purpose ...
 * -------
 * A set of Calcomp-like plotting routines for C.
 * This file contains all of the environment dependent routines.
 * All arguments specifying lengths in the following 
 * calls are assumed to be in millimetres.
 *
 * Routines supported ...
 *
 * Basic plotting routines
 *    ssp_BeginPlot (screen_flag, record_flag, file, device, orientation)
 *    ssp_EndPlot ()
 *    ssp_StartNewPlot ()
 *
 *    ssp_SetClipOn ()
 *    ssp_SetClipOff ()
 *    ssp_SetClipLimits (x_min, y_min, x_max, y_max)
 *
 *    ssp_SetLineThickness (thickness)
 *    ssp_SetFillGrey (lightness)
 *    ssp_SetHSBColour (hue, saturation, brightness)
 *    ssp_SetNoColour ()
 *    ssp_MapToColour ()
 *
 *    ssp_Plot (x, y)
 *    ssp_Move (x, y)
 *    ssp_Polygon (n, x[], y[], grey_fill, colour_fill, border)
 *
 *    ssp_Factor (zoom_factor)
 *    ssp_Rotate (angle)
 *    ssp_Where (*x, *y, *scale_factor, *angle)
 *
 *    ssp_PlotText (x, y, size, directn, *text, n)
 *
 * Swapping text and graphics screens
 *    ssp_TextOn ()
 *    ssp_TextOff ()
 *
 * Clearing the part of the screen reserved for text
 *    ssp_ClearLine (line)
 *    ssp_ClearText ()
 *
 * Auxiliary routines
 *    ssp_ReportError (*error_string)
 *    *ssp_Compact (x)
 *
 * Version ...
 * -------
 * 1.0, September 1988
 * 2.0, October   1991
 * 2.1, November  1991  : Starbase included
 * 2.2, December  1991  : raw hpgl included
 * 3.0, June      1992  : move text to top of IBM-PC screen
 * 3.1, 12 July   1992  : screen buffer allocated early and retained
 *                        ssp_WaitPrompt() and ssp_StartNewPlot()
 *                        added
 * 3.2, 26 July   1992  : Line clipping added
 *                        Move instruction removed from the postscript
 *                        output.
 * 3.3, 1  Sep    1992  : Tektronix screen control added
 * 3.4, 30 Sep    1992  : Precheck contouring data before contouring
 *                        each quadrilateral or triangle (Keith Weinman).
 *                        Change point size to the (correct) value of
 *                        72.27 points per inch.  Previously we have
 *                        been using the value of 72 big-points per inch.
 * 3.5, 07 Nov    1992  : Full Function-Prototypes
 * 4.0, 14 Nov    1992  : TopSpeed C++ compiler
 * 4.1, 09 April  1993  : Polygon drawing and filling
 * 4.11, 25-May-93 : geom.h -- a new header file
 * 4.12, 31-May-93 : tried to get more solid lines on the TopSpeed
 *                   screen when drawing filled polygons. The problem
 *                   seems to be in the filling of adjacent polygons
 *                   writing over the already drawn borders of other
 *                   polygons.  This does not seem to be a problem
 *                   on the postscript printers, so just leave it.
 * 4.13, 06-Jun-93 : plot polygons with optional borders
 * 4.2,  07-Jun-93 : conforming postscript
 *                   portrait/landscape option in ssp_BeginPlot()
 * 4.3,  20-Oct-93 : HSB colour function (and the NoColour fuunction)
 * 4.31, 24-Oct-93 : Added ssp_MapToColour() function.
 *                   Increase the size of the bounding box.
 * 5.0,  05-Nov-93 : Removed the pen functions and the Starbase
 *                   capability.
 * 6.0,  10-Nov-93 : Replaced the text functions with glyph-based
 *                   functions.
 * 6.1,  12-Nov-93 : Terminate writing of strings when a Null char
 *                   is encountered.
 * 6.2,  21-Dec-93 : Turbo-C routines working under TopSpeed
 *                   They should also work under Turbo-C (of course).
 * 6.3,  31-Dec-93 : Added colour_fill to polygon drawing.
 * 6.4,  11-Jan-94 : Removed some of the "dead wood" such as META...
 *                   Fixed the postscript colour polygons.
 * 7.0,  02-Apr-94 : Postsript plotting now done in points rather
 *                   than millimetres - smaller, faster ps files.
 * 7.1,  13-Jun-94 : Don't stroke postscript paths until necessary.
 *                   This should make the PS file shorter but adds
 *                   a bit of logic to check for active_ps_path each
 *                   time some postscript is output.
 * 7.2,  29-Aug-94 : Added more postscript macros to try to reduce
 *                   the size of a colour postscript file.
 *                   Also used a format of %4.2f for the floating
 *                   -point colour arguments.
 * 7.3   25-Apr-95 : Replaced NULL with 0 when referring to the Null
 *                   character.  The BCOS2 compiler seems to like to
 *                   reserve NULL for pointer values.
 * 8.0   17-Nov-97 : GIF output added via the gd 1.2 library by
 *                   Thomas Boutell of the Quest Protein Database Center.
 *                   Change the plotting area to a 200mm square.
 *                   This will suit the portrait mode of the postscript
 *                   printers and will suit the GIF generation.
 * 8.1   08-Aug-98 : Add a better colour selection for GIF images.
 *                   Included utility function to convert HSB to RGB
 *                   Plotting area changed to a 200mm by 170mm box.
 *
 *
 * Written by ... P. A. Jacobs
 * ----------     38 Ellington St,
 *                Tarragindi, Qld 4121.
 *
 * Notes ...
 * -----
 * (1) Graphics window definition:
 *     Map a 200mm by 170mm box to the display window.
 * (2) Definitions of the plotting units:
 *     Page Coordinates : these are the units (mm) in which marks
 *                        are made on the screen or page.
 *     Plot Coordinates : these are the units in which the user
 *                        specifies the plotting moves.
 *     Device Coordinates : The device/display plotting units.
 *                          These may be mm, pixels, points, or
 *                          1/40 mm (for HP plotters)
 *
 *
 *------------------------------------------------------------------
 */

/*
 * ------------------
 * Global Definitions
 * ------------------
 */

#include "../../util/source/compiler.h"
#include "geom.h"
#include "ssp.h"
#include "gd.h"

static char rcsid[]="$Id: sspkern_old.c,v 1.4 1998/10/01 05:38:27 peterj Exp $";

#if (GRAPH_ENV == TURBOC_G)
  /*
   * The Turbo-C graphics can be used from Borland's Turbo-C environment
   * (of course) or from TopSpeed C using the Clip-window library or
   * from X-Windows using Neville Richter's emulator.
   *
   * To use the TopSpeed compiler running the Turbo-C emulation,
   * the _CLIP_WIN_ macro needs to be defined before including conio.h.
   * This should not upset the Turbo-C compiler.
   */
#  if (TOPSPEED_C == 1)
#      define _CLIP_WIN_
#  endif
   /*
    * For DOS systems only; Unix does not have dos.h or conio.h.
    */
#  if (TOPSPEED_C == 1 || TURBO_C == 1)
#     include <dos.h>
#     include <conio.h>
#  endif
#  include <graphics.h>
#endif

#if (GRAPH_ENV == TOPSPEED_G)
#  include <dos.h>
#  include <conio.h>
#  include <graph.h>
#endif

#include <stdio.h>
#include <math.h>
#include <time.h>

#if STDLIBH == 1
#  include <stdlib.h>
#  include <string.h>
#  include <stdarg.h>
#endif

/*
 * -------------------------
 * Screen related quantities
 * -------------------------
 */

double XScrnDim = 0.90;       /* Fraction of physical screen covered */
double YScrnDim = 0.61;       /* by the graphics window. This leaves */
			      /* lines 1 and 2 for text if the plot  */
			      /* is put in the bottom left of the    */
			      /* screen.                             */
int    FirstLine = 1;
int    LastLine  = 2;

double XPlotDim = 200.0;      /* Real plot dimensions                */
double YPlotDim = 170.0;
int    XCursor = 20;          /* Size of cursor in pixels            */
int    YCursor = 10;

int    XPlotPixels;           /* Plot window size in pixels          */
int    YPlotPixels;
int    XMaxPixels;            /* Screen size in pixels               */
int    YMaxPixels;


/*
 * ---------------------------
 * Generic plotting quantities
 * ---------------------------
 */

double plot_scale;           /* Page Units = Plot Units * plot_scale  */
double global_angle;         /* Angle of the frame-of-ref with the    */
			     /* horizontal axis                       */

double x_old, y_old;         /* Previous location (page units)        */
double x_origin, y_origin;   /* Plotting origin (page units)          */

int    do_clipping;            /* clipping flag (=1 to clip)          */
double x_min_clip, y_min_clip; /* lower-left clip window (page units) */
double x_max_clip, y_max_clip; /* upper-right clip       (page units) */

int    display_on_screen;    /* =1 : display on screen                */
int    record_in_file;       /* =1 : record commands in file          */
int    outdevice;            /* output device                         */
int    portrait;             /* =0 : landscape mode                   */
                             /* =1 : portrait mode                    */
FILE   *plotfile;            /* text file for plotting commands       */
int    GraphicsOn;           /* true if graphics is displayed         */

double dash_length_mm;       /* length of dash in page units, mm      */
double space_length_mm;      /* length od space in page units, mm     */
double combined_length_mm;   /* dash+space                            */

int    in_colour;            /* ==1 the line work is in colour        */
                             /* ==0 the line work is in black/white   */
double red, blue, green;     /* Colour values in the range 0.0 to 1.0 */
double hue, brightness, saturation;                 /* line colours   */
int    vga_colour;           /* VGA colour value (for lines)          */

double fill_hue, fill_brightness, fill_saturation;  /* fill colours   */
double fill_grey_level;      /* Grey level for shading polygons etc   */
                             /*             0.0=black, 1.0 = white    */
int    fill_vga_colour;      /* VGA colour value for filling polygons */

/*
 * The following definitions are basically those in the Turbo-C
 * graphics environment.
 */
#define  VGA_BLACK         0
#define  VGA_BLUE          1
#define  VGA_GREEN         2
#define  VGA_CYAN          3
#define  VGA_RED           4
#define  VGA_MAGENTA       5
#define  VGA_BROWN         6
#define  VGA_LIGHTGRAY     7
#define  VGA_DARKGRAY      8
#define  VGA_LIGHTBLUE     9
#define  VGA_LIGHTGREEN    10
#define  VGA_LIGHTCYAN     11
#define  VGA_LIGHTRED      12
#define  VGA_LIGHTMAGENTA  13
#define  VGA_YELLOW        14
#define  VGA_WHITE         15

/*
 * ------------------
 * Glyph related data
 * ------------------
 */

char **glyph_stroke;
int glyph_code[96];
int nstroke[96];

/*
 * --------------------------
 * Turbo-C related quantities
 * --------------------------
 */

int    GraphDriver;
int    GraphMode;
int    ErrorCode;
unsigned buffer_size;
void   *buffer;           /* buffer for saving screens */
/*
 * If we wish to not have bgi drivers in the program area, use ...
 * char   PathToDriver[] = "c:/tc/bgi";
 * else the following line is the most basic arrangement
 * char   PathToDriver[] = "";
 */
char   PathToDriver[] = "";

/*
 * ---------------------------
 * TopSpeed related quantities
 * ---------------------------
 */
#if (GRAPH_ENV == TOPSPEED_G)
   struct videoconfig videoc;
#endif
int  pixels_per_col, pixels_per_row;
int  col_offset, row_offset;
int  col_text, row_text;
unsigned char fill_mask_00[8], fill_mask_05[8], fill_mask_10[8];

/*
 * ----------------------------
 * Microsoft related quantities
 * ----------------------------
 */

/*
 * ---------------------------
 * Xwindows related quantities
 * ---------------------------
 */

/*
 * -----------------------------
 * Postscript related quantities
 * -----------------------------
 */

#define  Lstr  10
#define  points_per_mm  2.8453
int      active_ps_path;

/*
 * ----------------------
 * GIF related quantities
 * ----------------------
 * Also uses some of the pixel quantities listed above.
 */
char       gif_file_name[132];
gdImagePtr gif_image;
struct     hsv_colour hsv_values;
struct     rgb_colour rgb_values;
int        background_colour, gif_colour, fill_gif_colour;
#define    MAX_GIF_HUE 24
int        gif_hue[MAX_GIF_HUE];
#define    MAX_GIF_GREY 8
int        gif_greyshade[MAX_GIF_GREY];
float      delta_hue, delta_grey;
int        ihue, igrey;
int        gif_black, gif_blue, gif_green, gif_cyan;
int        gif_red, gif_magenta, gif_brown, gif_lightgray;
int        gif_darkgray, gif_lightblue, gif_lightgreen, gif_lightcyan;
int        gif_lightred, gif_lightmagenta, gif_yellow, gif_white;

/*
 * -----------------------
 * HPGL related quantities
 * -----------------------
 */
#define  hp2mm  40.0
#define  ETX    3
#define  ESC    27

/*
 * ---------
 * Tektronix
 * ---------
 */
#define  BEL   0x0007
#define  BS    0x0008
#define  LF    0x000A
#define  FF    0x000C
#define  CR    0x000D
/* ESC 1B above */
#define  GS    0x001D
#define  US    0x001F
#define  SPACE 0x0020

char hi_y, lo_y, hi_x, lo_x;

/* ----------------------------------------------------------
 * For information on Tektronix 4010 plotting see...
 * Mahlon Kelly
 * Mainframe graphics on a microcomputer.
 * Byte October 1983, Vol 8 (10), pp 439--442
 *
 * The order of characters is hi_y, lo_y, hi_x, lo_x.
 * A line or move is made after the arrival of lo_x.
 * An (ix,iy) coordinate is reconstituted using
 *   ix = (hi_x - 32) * 32 + lo_x - 64
 *   iy = (hi_y - 32) * 32 + lo_y - 96
 * where
 *   32 <= hi_y <= 63
 *   96 <= lo_y <= 127
 *   32 <= hi_x <= 63
 *   64 <= lo_x <= 95
 *------------------------------------------------------------
 */


/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_BeginPlot (int screen_flag, int record_flag,
		   char *out_file, char *device, char *vmode)

#else

int ssp_BeginPlot (screen_flag, record_flag, out_file, device, vmode)
int screen_flag, record_flag;
char *out_file;
char *device;
char *vmode;

#endif

/*
 * Purpose...
 * -------
 * Initialize the plot.
 *
 * Input...
 * -----
 * screen_flag : =1 : display on screen
 *               =0 : DO NOT display on screen
 * record_flag : =1 : record commands in a file
 *               =0 : DO NOT record commands
 * out_file : a character string giving the name of the recording file
 * device   : a character string indicating the recording format
 *            (only the first letter is significant)
 *            "meta"   : metafile to be later translated
 *            "hpgl"   : (Raw) Hewlett-Packard Graphics Language
 *            "ps"     : postscript
 *            "gif"    : GIF for web documents
 * vmode    : Viewing mode for postscript output
 *            (only the first letter is significant)
 *            "portrait"
 *            "landscape"
 */

{
char letter;
time_t bintime;

/*
 * Set the device flag.
 */

display_on_screen = (screen_flag == 1);
record_in_file = (record_flag == 1);

/*
 * Determine the recording device.
 */
outdevice = PSCRIPT;             /* Postscript is default */
letter = device[0];
if (letter == 'h' || letter == 'H')
   outdevice = HPGL;
else if (letter == 'p' || letter == 'P')
   outdevice = PSCRIPT;
else if (letter == 'g' || letter == 'G')
   outdevice = GIF;

letter = vmode[0];
if (letter == 'p' || letter == 'P')
   portrait = 1;
else
   portrait = 0;

/*
 * Initialize the Generic Plotting Variables.
 */
plot_scale = 1.0;    /* plot in mm units */
global_angle = 0.0;

x_origin = 0.0;      /* origin in bottom left corner */
y_origin = 0.0;

x_min_clip = 0.0;   /* lower-left of the clip window  */
y_min_clip = 0.0;
x_max_clip = 200.0; /* upper-right of the clip window */
y_max_clip = 170.0;

do_clipping = 0;    /* By default, do NOT clip */

x_old = x_origin;    /* start at origin */
y_old = y_origin;

/* Default Grey level and colours */
in_colour = 0;   /* default black/white lines */
red   = 0.5;
blue  = 0.5;
green = 0.5;
hue        = 0.0;
saturation = 0.0;
brightness = 1.0;
vga_colour = VGA_WHITE;

fill_grey_level = 0.5;
fill_hue        = 0.0;
fill_saturation = 0.0;
fill_brightness = 1.0;
fill_vga_colour = VGA_WHITE;

/* Allocate memory and assign the glyph encoded strings. */
ssp_AllocGlyph ();
ssp_InitGlyph ();

/* Default dashed line in millimetres */
ssp_SetDash (4.0, 2.0);


if (record_in_file)
  {
  /*
   * Open my own file to record the plotting commands.
   */

  if (outdevice == HPGL)
     {
     /* HP plotters or Laserjet printers
      * First, use Printer Control Language (PCL) to reset the
      * laserjet and to go to HPGL mode.
      */
     if ((plotfile = fopen(out_file, "w")) == NULL)
        {
        printf ("ssp_BeginPlot -- cannot open plot file %s\n", out_file);
        exit (1);
        }
     fprintf (plotfile, "%cE%c&l1O%c%%0B\n", ESC, ESC, ESC);
     /* HPGL: reset the plotter and select pen 1 */
     fprintf (plotfile, "IN;PU;PA0,0;SP1;\n");
     fprintf (plotfile, "PW%f;PT%f;\n", 0.35, 0.35);
     }
  else if (outdevice == PSCRIPT)
     {
     active_ps_path = 0;

     if ((plotfile = fopen(out_file, "w")) == NULL)
        {
        printf ("ssp_BeginPlot -- cannot open plot file %s\n", out_file);
        exit (1);
        }
     fprintf (plotfile, "%%!PS-Adobe-1.0\n");
     fprintf (plotfile, "%%%%Title: %s\n", out_file);
     fprintf (plotfile, "%%%%Creator: Simple Scientific-Plotting Package\n");
     time (&bintime);
     fprintf (plotfile, "%%%%CreationDate: %s", ctime(&bintime) );
     fprintf (plotfile, "%%%%Pages: 1\n");
     fprintf (plotfile, "%%%%DocumentFonts: Courier\n");
     /* We shall assume that the figure fits a 200mm by 170mm box. */
     fprintf (plotfile, "%%%%BoundingBox: 0 0 567 482\n");
     fprintf (plotfile, "%%%%EndComments\n");
     fprintf (plotfile, "%%\n");
     fprintf (plotfile, "save\n");
     fprintf (plotfile, "%% Convert from millimetres to points\n");
     fprintf (plotfile, "/mm {72.27 mul 25.4 div} def\n");
     fprintf (plotfile, "%%\n");
     fprintf (plotfile, "%% Shorthand for frequent commands\n");
     fprintf (plotfile, "/l {lineto} def\n");
     fprintf (plotfile, "/m {moveto} def\n");
     fprintf (plotfile, "/s {stroke} def\n");
     fprintf (plotfile, "/n {newpath} def\n");
     fprintf (plotfile, "/sg {setgray} def\n");
     fprintf (plotfile, "/sc {sethsbcolor} def\n");
     fprintf (plotfile, "/f {fill} def\n");
     fprintf (plotfile, "/gr {grestore} def\n");
     fprintf (plotfile, "/gs {gsave} def\n");
     fprintf (plotfile, "/cp {closepath} def\n");
     fprintf (plotfile, "%%\n");
     fprintf (plotfile, "%% Default settings\n");
     fprintf (plotfile, "1 setlinecap %% semicircular ends on lines\n");
     fprintf (plotfile, "%9.2f mm setlinewidth\n", 0.35);
     fprintf (plotfile, "%%%%EndProlog\n");
     fprintf (plotfile, "%%%%Page: 1 1\n");

     if (portrait == 1)
        {
        fprintf (plotfile, "%%Portrait mode.\n");
        fprintf (plotfile, "10 10 translate\n");
        fprintf (plotfile, "0 rotate\n");
        }
     else
        {
        fprintf (plotfile, "%%Landscape mode.\n");
        fprintf (plotfile, "7.32 72.27 mul 0 translate\n");
        fprintf (plotfile, "90 rotate\n");
        }

     fprintf (plotfile, "%%\n");
     fprintf (plotfile, "%%----------------------------------------\n");
     }  /* End of Postscript prologue. */
  else if (outdevice == GIF)
     {
     /* 
      * Set up to build the image in memory.
      * Don't write the file itself until the end of plotting.
      */
     strcpy(gif_file_name, out_file);
     XPlotPixels = 500;  /* 2.5 times 200 */
     YPlotPixels = 425;  /* 2.5 times 170 */
     gif_image = gdImageCreate(XPlotPixels, YPlotPixels);
     /* White background */
     gif_white = gdImageColorAllocate(gif_image, 255, 255, 255);
     background_colour = gif_white;
     /* Black foreground */
     gif_black = gdImageColorAllocate(gif_image, 0, 0, 0);
     gif_colour = gif_black;
     /* Colour wheel for contouring; bright colours */
     delta_hue = 1.0 / MAX_GIF_HUE;
     for (ihue = 0; ihue < MAX_GIF_HUE; ihue++) {
         hsv_values.h = delta_hue * ihue;
         hsv_values.s = 1.0;
         hsv_values.v = 1.0;
         hsv2rgb( &hsv_values, &rgb_values );
         gif_hue[ihue] = gdImageColorAllocate( gif_image, 
            (int) 255 * rgb_values.r, 
            (int) 255 * rgb_values.g, 
            (int) 255 * rgb_values.b );
     } /* end for */
     /* levels of grey for contouring */
     delta_grey = 1.0 / MAX_GIF_GREY;
     for (igrey = 0; igrey < MAX_GIF_GREY; igrey++) {
         hsv_values.h = 0.0;
         hsv_values.s = 0.0;
         hsv_values.v = delta_grey * igrey;
         hsv2rgb( &hsv_values, &rgb_values );
         gif_greyshade[igrey] = gdImageColorAllocate( gif_image, 
            (int) 255 * rgb_values.r, 
            (int) 255 * rgb_values.g, 
            (int) 255 * rgb_values.b );
     } /* end for */
     /* Other colours to match 16 VGA colours. */
     gif_blue         = gdImageColorAllocate(gif_image,   0,   0, 255);
     gif_green        = gdImageColorAllocate(gif_image,   0, 255,   0); 
     gif_cyan         = gdImageColorAllocate(gif_image,   0, 255, 255);
     gif_red          = gdImageColorAllocate(gif_image, 255,   0,   0);
     gif_magenta      = gdImageColorAllocate(gif_image, 255,   0, 255);
     gif_brown        = gdImageColorAllocate(gif_image, 255,  93,  63);
     gif_lightgray    = gdImageColorAllocate(gif_image, 200, 200, 200);
     gif_darkgray     = gdImageColorAllocate(gif_image, 100, 100, 100);
     gif_lightblue    = gdImageColorAllocate(gif_image, 128, 128, 255);
     gif_lightgreen   = gdImageColorAllocate(gif_image, 128, 255, 128);
     gif_lightcyan    = gdImageColorAllocate(gif_image, 128, 255, 255);
     gif_lightred     = gdImageColorAllocate(gif_image, 255, 128, 128);
     gif_lightmagenta = gdImageColorAllocate(gif_image, 255, 128, 255);
     gif_yellow       = gdImageColorAllocate(gif_image, 255, 255,   0);
     }
  } /* end of record_in_file */ 

if (display_on_screen)
   {
   GraphicsOn = 1;

#  if (GRAPH_ENV == TURBOC_G)
   /*
    * Initialize the Turbo-C Graphics environment.
    */
#  if (OLIVETTI)
      GraphDriver = ATT400;
      GraphMode = ATT400HI;
#  else
      GraphDriver = DETECT;
#  endif
   initgraph (&GraphDriver, &GraphMode, PathToDriver);
   ErrorCode = graphresult ();
   if (ErrorCode != grOk)
      {
      printf (" Graphics Error in ssp_BeginPlot : %s\n",
      grapherrormsg (ErrorCode));
      exit (1);
      }
   XMaxPixels = getmaxx ();
   YMaxPixels = getmaxy ();
   cleardevice ();
   setcolor (getmaxcolor());

   /*
    * Place the graphics window in the lower part of the screen.
    */
   XPlotPixels = (int) (XMaxPixels * XScrnDim);
   YPlotPixels = (int) (YMaxPixels * YScrnDim);
   setviewport (0, YMaxPixels - YPlotPixels,
		XPlotPixels, YMaxPixels, 1);
   rectangle (0, 0, XPlotPixels, YPlotPixels);
   ErrorCode = graphresult ();
   if (ErrorCode != grOk)
      {
      printf (" Graphics Error in ssp_BeginPlot : %s\n",
	 grapherrormsg (ErrorCode));
      exit (1);
      }
   buffer_size = imagesize (0, 0, XPlotPixels, YPlotPixels);
   ErrorCode = graphresult ();
   if (ErrorCode != grOk)
      {
      printf ("Graphics Error in ssp_TextOn() : %s\n",
	      grapherrormsg (ErrorCode));
      exit (1);
      }

   buffer = NULL;
   if (buffer_size != -1)
      buffer = malloc (buffer_size);
   else
      buffer = NULL;

   /*
    * Lines FirstLine through LastLine are available for text.
    */
   gotoxy (1,FirstLine);
#  endif

#  if (GRAPH_ENV == TOPSPEED_G)
   /*
    * Initialize the TopSpeed Graphics environment.
    */
   _setvideomode (_VRES16COLOR);  /* VGA 16 colour 640x480 pixels */
   _getvideoconfig (&videoc);
   XMaxPixels = videoc.numxpixels - 1;
   YMaxPixels = videoc.numypixels - 1;

   /*
    * Place the graphics window in the lower part of the screen.
    */
   XPlotPixels = (int) (XMaxPixels * XScrnDim);
   YPlotPixels = (int) (YMaxPixels * YScrnDim);
   _setviewport (0, YMaxPixels - YPlotPixels,
		XPlotPixels, YMaxPixels);
   _rectangle (_GBORDER, 0, 0, XPlotPixels, YPlotPixels);

   /*
    * The grey fill mask for shading polygons etc.
    */
   fill_mask_05[0] = 0xAA; fill_mask_05[1] = 0x55;
   fill_mask_05[2] = 0xAA; fill_mask_05[3] = 0x55;
   fill_mask_05[4] = 0xAA; fill_mask_05[5] = 0x55;
   fill_mask_05[6] = 0xAA; fill_mask_05[7] = 0x55;

   /*
    * The black fill mask.
    */
   fill_mask_00[0] = 0x00; fill_mask_00[1] = 0x00;
   fill_mask_00[2] = 0x00; fill_mask_00[3] = 0x00;
   fill_mask_00[4] = 0x00; fill_mask_00[5] = 0x00;
   fill_mask_00[6] = 0x00; fill_mask_05[7] = 0x00;

   /*
    * The white fill mask.
    */
   fill_mask_10[0] = 0xFF; fill_mask_10[1] = 0xFF;
   fill_mask_10[2] = 0xFF; fill_mask_10[3] = 0xFF;
   fill_mask_10[4] = 0xFF; fill_mask_10[5] = 0xFF;
   fill_mask_10[6] = 0xFF; fill_mask_10[7] = 0xFF;

   /*
    * Default to the grey mask.
    */
   _setfillmask (fill_mask_05);

   /*
    * Characters are put down in text positions.
    */
   pixels_per_col = videoc.numxpixels / videoc.numtextcols;
   pixels_per_row = videoc.numypixels / videoc.numtextrows;
   col_offset = 1;
   row_offset = (YMaxPixels - YPlotPixels) / pixels_per_row + 1;

   /*
    * For now, don't worry about saving screen.
    */
   buffer = NULL;

   /*
    * Lines FirstLine through LastLine are available for text.
    */
   _settextposition (FirstLine, 1);
#  endif

#  if (GRAPH_ENV == TEKTERM)
   printf ("%c%c", ESC, FF);
   hi_y = 0 + 32;
   lo_y = 0 + 96;
   hi_x = 0 + 32;
   lo_x = 0 + 64;
   printf ("%c%c%c%c%c", GS, hi_y, lo_y, hi_x, lo_x);
#  endif

   }  /* end of screen initialization */

return (0);
}  /* end of ssp_BeginPlot() */



/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_StartNewPlot (void)

#else

int ssp_StartNewPlot ()

#endif

/*
 * Purpose...
 * -------
 * Clear the screen or change to a new page to start a new plot.
 *
 * Input...
 * -----
 * None
 */

{  /* begin ssp_StartNewPlot() */
/*
 * Initialize the Generic Plotting Variables.
 */
plot_scale = 1.0;    /* plot in mm units */
global_angle = 0.0;

x_origin = 0.0;      /* origin in bottom left corner */
y_origin = 0.0;

x_min_clip = 0.0;   /* lower-left of the clip window  */
y_min_clip = 0.0;
x_max_clip = 200.0; /* upper-right of the clip window */
y_max_clip = 170.0;

do_clipping = 0;    /* By default, do NOT clip */

x_old = x_origin;    /* start at origin */
y_old = y_origin;

/* Default dashed line in millimetres */
ssp_SetDash (4.0, 2.0);

/* Default Grey level and colours */
in_colour = 0;   /* default black/white lines */
red   = 0.5;
blue  = 0.5;
green = 0.5;
hue        = 0.0;
saturation = 0.0;
brightness = 1.0;
vga_colour = VGA_WHITE;

fill_grey_level = 0.5;
fill_hue        = 0.0;
fill_saturation = 0.0;
fill_brightness = 1.0;
fill_vga_colour = VGA_WHITE;

if (record_in_file)
  {
  if (outdevice == HPGL)
     {
     /* HP plotters or Laserjet printers */
     fprintf (plotfile, "IN;PU;PA0,0;SP1;\n");
     fprintf (plotfile, "PW%f;PT%f;\n", 0.35, 0.35);
     }
  else if (outdevice == PSCRIPT)
     {
     if (active_ps_path == 1)
        {
        /* Stroke the active path to finish it. */
        fprintf (plotfile, "s\n");
        active_ps_path = 0;
        }
     fprintf (plotfile, "showpage\n");
     fprintf (plotfile, "%%----------------------------------------\n");

     if (portrait == 1)
        {
        fprintf (plotfile, "%%Portrait mode.\n");
        fprintf (plotfile, "10 10 translate\n");
        fprintf (plotfile, "0 rotate\n");
        }
     else
        {
        fprintf (plotfile, "%%Landscape mode.\n");
        fprintf (plotfile, "7.32 72.27 mul 0 translate\n");
        fprintf (plotfile, "90 rotate\n");
        }

     fprintf (plotfile, "%%\n");
     fprintf (plotfile, "%%----------------------------------------\n");
     } /* end of postscript restart */
  else if ( outdevice == GIF )
     {
     /*
      * Destroy the old image and start over.
      */
     gdImageDestroy(gif_image);
     XPlotPixels = 500;  /* 2.5 times 200 */
     YPlotPixels = 425;  /* 2.5 times 170 */
     gif_image = gdImageCreate(XPlotPixels, YPlotPixels);
     /* White background */
     gif_white = gdImageColorAllocate(gif_image, 255, 255, 255);
     background_colour = gif_white;
     /* Black foreground */
     gif_black = gdImageColorAllocate(gif_image, 0, 0, 0);
     gif_colour = gif_black;
     /* Colour wheel for contouring; bright colours */
     delta_hue = 1.0 / MAX_GIF_HUE;
     for (ihue = 0; ihue < MAX_GIF_HUE; ihue++) {
         hsv_values.h = delta_hue * ihue;
         hsv_values.s = 1.0;
         hsv_values.v = 1.0;
         hsv2rgb( &hsv_values, &rgb_values );
         gif_hue[ihue] = gdImageColorAllocate( gif_image, 
            (int) 255 * rgb_values.r, 
            (int) 255 * rgb_values.g, 
            (int) 255 * rgb_values.b );
     } /* end for */
     /* levels of grey for contouring */
     delta_grey = 1.0 / MAX_GIF_GREY;
     for (igrey = 0; igrey < MAX_GIF_GREY; igrey++) {
         hsv_values.h = 0.0;
         hsv_values.s = 0.0;
         hsv_values.v = delta_grey * igrey;
         hsv2rgb( &hsv_values, &rgb_values );
         gif_greyshade[igrey] = gdImageColorAllocate( gif_image, 
            (int) 255 * rgb_values.r, 
            (int) 255 * rgb_values.g, 
            (int) 255 * rgb_values.b );
     } /* end for */
     /* Other colours to match 16 VGA colours. */
     gif_blue         = gdImageColorAllocate(gif_image,   0,   0, 255);
     gif_green        = gdImageColorAllocate(gif_image,   0, 255,   0); 
     gif_cyan         = gdImageColorAllocate(gif_image,   0, 255, 255);
     gif_red          = gdImageColorAllocate(gif_image, 255,   0,   0);
     gif_magenta      = gdImageColorAllocate(gif_image, 255,   0, 255);
     gif_brown        = gdImageColorAllocate(gif_image, 255,  93,  62);
     gif_lightgray    = gdImageColorAllocate(gif_image, 200, 200, 200);
     gif_darkgray     = gdImageColorAllocate(gif_image, 100, 100, 100);
     gif_lightblue    = gdImageColorAllocate(gif_image, 128, 128, 255);
     gif_lightgreen   = gdImageColorAllocate(gif_image, 128, 255, 128);
     gif_lightcyan    = gdImageColorAllocate(gif_image, 128, 255, 255);
     gif_lightred     = gdImageColorAllocate(gif_image, 255, 128, 128);
     gif_lightmagenta = gdImageColorAllocate(gif_image, 255, 128, 255);
     gif_yellow       = gdImageColorAllocate(gif_image, 255, 255,   0);
     }
  else
     ;  /* all other devices, do nothing */
  }

if (display_on_screen)
   {
   GraphicsOn = 1;

#  if (GRAPH_ENV == TURBOC_G)
   /*
    * Clear the screen and restart.
    */
   setgraphmode ( getgraphmode() );
   ErrorCode = graphresult ();
   if (ErrorCode != grOk)
      {
      printf (" Graphics Error in ssp_StartNewPlot : %s\n",
      grapherrormsg (ErrorCode));
      exit (1);
      }
   XMaxPixels = getmaxx ();
   YMaxPixels = getmaxy ();
   cleardevice ();
   setcolor (getmaxcolor());

   /*
    * Place the graphics window in the lower part of the screen.
    */
   XPlotPixels = (int) (XMaxPixels * XScrnDim);
   YPlotPixels = (int) (YMaxPixels * YScrnDim);
   setviewport (0, YMaxPixels - YPlotPixels,
		XPlotPixels, YMaxPixels, 1);
   rectangle (0, 0, XPlotPixels, YPlotPixels);
   ErrorCode = graphresult ();
   if (ErrorCode != grOk)
      {
      printf (" Graphics Error in ssp_StartNewPlot : %s\n",
	 grapherrormsg (ErrorCode));
      exit (1);
      }
   /*
    * Lines FirstLine through LastLine are available for text.
    */
   gotoxy (1,FirstLine);
#  endif

#  if (GRAPH_ENV == TOPSPEED_G)
   _clearscreen (_GCLEARSCREEN);
   /*
    * Place the graphics window in the lower part of the screen.
    */
   XPlotPixels = (int) (XMaxPixels * XScrnDim);
   YPlotPixels = (int) (YMaxPixels * YScrnDim);
   _setviewport (0, YMaxPixels - YPlotPixels,
		XPlotPixels, YMaxPixels);
   _rectangle (_GBORDER, 0, 0, XPlotPixels, YPlotPixels);

   /*
    * Lines FirstLine through LastLine are available for text.
    */
   _settextposition (FirstLine, 1);
#  endif

#  if (GRAPH_ENV == TEKTERM)
   printf ("%c%c", ESC, FF);
   hi_y = 0 + 32;
   lo_y = 0 + 96;
   hi_x = 0 + 32;
   lo_x = 0 + 64;
   printf ("%c%c%c%c%c", GS, hi_y, lo_y, hi_x, lo_x);
#  endif

   }  /* end of screen initialization */

return (0);
}  /* end of ssp_StartNewPlot() */

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_EndPlot (void)

#else

int ssp_EndPlot ()

#endif

/*
 * Purpose...
 * -------
 * Clean up the graphics and wipe the screen clean.
 * Also free memory for the glyph data.
 */

{
ssp_FreeGlyph ();

/*
 * Close the recording file if it is open.
 */
if (record_in_file)
   {
   if (outdevice == HPGL)
      fprintf (plotfile, "SP0;%c%%0A%cE\n", ESC, ESC);

   else if (outdevice == PSCRIPT)
      {
      if (active_ps_path == 1)
         {
         /* Stroke the active path to finish it. */
         fprintf (plotfile, "s\n");
         active_ps_path = 0;
         }
      fprintf (plotfile, "showpage\n");
      fprintf (plotfile, "%%%%Trailer\n");
      fprintf (plotfile, "restore\n");
      fclose (plotfile);
      }
   else if (outdevice == GIF)
      {
      /*
       * Now is the time to write the image to a GIF file.
       * The file needs to be opened as binary.
       */
      if ((plotfile = fopen(gif_file_name, "wb")) == NULL)
         {
         printf ("ssp_EndPlot -- cannot open GIF file %s\n", gif_file_name);
         exit (1);
         }
      gdImageGif(gif_image, plotfile);
      fclose(plotfile);
      gdImageDestroy(gif_image);
      }
   }  /* end of record_in_file */

/*
 * Now cleanup the graphics screen.
 */
if (display_on_screen)
   {
#  if (GRAPH_ENV == TURBOC_G)
   ssp_WaitMessage ("End of Graphics");
   closegraph ();
   GraphicsOn = 0;
   if (buffer != NULL) free (buffer);
#  endif

#  if (GRAPH_ENV == TOPSPEED_G)
   ssp_WaitMessage ("End of Graphics");
   _setvideomode (_DEFAULTMODE);
   GraphicsOn = 0;
   if (buffer != NULL) free (buffer);
#  endif

#  if (GRAPH_ENV == TEKTERM)
   printf ("%c%c", ESC, FF);
   printf ("%c", US);
#  endif

   }

return (0);
} /* end of ssp_EndPlot() */

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetClipOn (void)

#else

int ssp_SetClipOn ()

#endif

/* Purpose...
 * -------
 * Set the clipping flag ON.
 */

{
do_clipping = 1;
return (0);
}


/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetClipOff (void)

#else

int ssp_SetClipOff ()

#endif

/* Purpose...
 * -------
 * Set the clipping flag OFF.
 */

{
do_clipping = 0;
return (0);
}


/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetClipLimits (double x_min, double y_min,
		       double x_max, double y_max)

#else

int ssp_SetClipLimits (x_min, y_min, x_max, y_max)
double x_min, y_min, x_max, y_max;

#endif

/* Purpose...
 * -------
 * Set the current clipping limits.
 *
 * Input...
 * -----
 * x_min, y_min : lower-left corner of the clip window in current
 *                plot coordinates.  Note that these are converted
 *                to page coordinates (mm) for actual clipping.
 * x_max, y_max : upper-right corner of the clip page.
 *
 * Output...
 * ------
 * ssp_SetClipLimits() returns 0 for a normal return; -1 if there
 * are problems with the values of the clipping limits.
 *
 */

{
/* Check for a badly defined window. */
if (x_max_clip <= x_min_clip) return (-1);
if (y_max_clip <= y_min_clip) return (-1);

/* Convert to page coordinates. */
x_min_clip = x_min * plot_scale + x_origin;
x_max_clip = x_max * plot_scale + x_origin;

y_min_clip = y_min * plot_scale + y_origin;
y_max_clip = y_max * plot_scale + y_origin;

return (0);
}


/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_CheckVisible (struct point_2D *P1,
		      struct point_2D *P2)

#else

int ssp_CheckVisible (P1, P2)
struct point_2D  *P1, *P2;

#endif

/* Purpose...
 * -------
 * Check the visibility of a line segment.
 *
 * Input...
 * -----
 * *P1  : pointer to the first end point
 * *P2  : pointer to the second end point
 * The window clip limits are global variables.
 *
 * Output...
 * ------
 * ssp_CheckVisible() returns an integer value with one of
 * the following values: ALL_VISIBLE, NOT_VISIBLE, PART_VISIBLE.
 * On return with a PART_VISIBLE result, P1 is always outside
 * the clip window.
 *
 * Reference...
 * ---------
 * Chapter 3 of Procedural Elements for Computer Graphics by
 * David F. Rogers.
 */

{  /* begin ssp_CheckVisible() */
unsigned char ctemp;
double   temp;

/*
 * Calculate the end-point codes.
 */
P1->code = 0;
if ( P1->x < (x_min_clip - CLIP_TOL) ) P1->code |= 1;  /* set bit 0 */
if ( P1->x > (x_max_clip + CLIP_TOL) ) P1->code |= 2;  /* set bit 1 */
if ( P1->y < (y_min_clip - CLIP_TOL) ) P1->code |= 4;  /* set bit 2 */
if ( P1->y > (y_max_clip + CLIP_TOL) ) P1->code |= 8;  /* set bit 3 */

P2->code = 0;
if ( P2->x < (x_min_clip - CLIP_TOL) ) P2->code |= 1;
if ( P2->x > (x_max_clip + CLIP_TOL) ) P2->code |= 2;
if ( P2->y < (y_min_clip - CLIP_TOL) ) P2->code |= 4;
if ( P2->y > (y_max_clip + CLIP_TOL) ) P2->code |= 8;

/*
 * Now decide the visibility of the line.
 */
if ( (P1->code | P2->code) == 0 )
   {
   /* The line is completely visible. */
   return (ALL_VISIBLE);
   }
else if ( (P1->code & P2->code) != 0 )
  {
   /* The line is trivially invisible. */
   return (NOT_VISIBLE);
   }

/*
 * If, we get this far, the line may be partially visible.
 * Make sure that P1 is outside the window.
 */
if (P1->code == 0)
   {
   temp = P1->x; P1->x = P2->x; P2->x = temp;
   temp = P1->y; P1->y = P2->y; P2->y = temp;
   ctemp = P1->code; P1->code = P2->code; P2->code = ctemp;
   }
return (PART_VISIBLE);

}  /* end of ssp_CheckVisible() */


/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_ClipLine (struct point_2D *P1,
		  struct point_2D *P2)

#else

int ssp_ClipLine (P1, P2)
struct point_2D  *P1, *P2;

#endif

/* Purpose...
 * -------
 * Adjust the end-point coordinates of a line so that
 * the displayed line remains within the clip window.
 * The Sutherland-Cohen two-dimensional clipping algorithm
 * is used.
 *
 * Input...
 * -----
 * *P1  : pointer to the first end point
 * *P2  : pointer to the second end point
 * The window clip limits are global variables.
 *
 * Output...
 * ------
 * ssp_ClipLine() returns an integer value YES or NO indicating
 * whether any of the line is visible.
 * The coordinates of P1 and P2 may be altered.
 *
 * Reference...
 * ---------
 * Chapter 3 of Procedural Elements for Computer Graphics by
 * David F. Rogers.
 */

{ /* begin ssp_ClipLine() */
int    i, vertical, visible;
double slope, dx;
double y_intersect, x_intersect;


/*
 * Determine slope.
 */
dx = P2->x - P1->x;
if ( fabs(dx) < 1.0e-10 )
   {
   vertical = 1;
   }
else
   {
   vertical = 0;
   slope = (P2->y - P1->y) / dx;
   }

/*
 * Search each window edge for intersections.
 * Go around twice to make sure that we cut enough bits off
 * both end points.
 */

for (i = 1; i <= 2; ++i)
{
/*
 * First, check the left window edge.
 */
visible = ssp_CheckVisible (P1, P2);
if (visible == ALL_VISIBLE) return (ALL_VISIBLE);
if (visible == NOT_VISIBLE) return (NOT_VISIBLE);

/* Partially visible line; cut off a piece. */
if (vertical != 1 && P1->x < x_min_clip)
   {
   /*
    * Find the intersection with the left-edge of the window.
    * P1 is assumed to lie outside the window.
    */
   y_intersect = slope * (x_min_clip - P1->x) + P1->y;
   P1->x = x_min_clip;
   P1->y = y_intersect;
   }

visible = ssp_CheckVisible (P1, P2);
if (visible == ALL_VISIBLE) return (ALL_VISIBLE);
if (visible == NOT_VISIBLE) return (NOT_VISIBLE);

/* Partially visible line; cut off a piece. */
if (vertical != 1 && P1->x > x_max_clip)
   {
   /*
    * Find the intersection with the right-edge of the window.
    * P1 is assumed to lie outside the window.
    */
   y_intersect = slope * (x_max_clip - P1->x) + P1->y;
   P1->x = x_max_clip;
   P1->y = y_intersect;
   }

visible = ssp_CheckVisible (P1, P2);
if (visible == ALL_VISIBLE) return (ALL_VISIBLE);
if (visible == NOT_VISIBLE) return (NOT_VISIBLE);

/*
 * A horizontal line should not have made it past this point.
 */

/* Partially visible line; cut off a piece. */
if (P1->y < y_min_clip)
   {
   if (vertical == 1)
      {
      P1->y = y_min_clip;
      }
   else
      {
      /*
       * Find the intersection with the left-edge of the window.
       * P1 is assumed to lie outside the window.
       */
      x_intersect = (1.0 / slope) * (y_min_clip - P1->y) + P1->x;
      P1->y = y_min_clip;
      P1->x = x_intersect;
      }
   }


visible = ssp_CheckVisible (P1, P2);
if (visible == ALL_VISIBLE) return (ALL_VISIBLE);
if (visible == NOT_VISIBLE) return (NOT_VISIBLE);

/* Partially visible line; cut off a piece. */
if (P1->y > y_max_clip)
   {
   if (vertical == 1)
      {
      P1->y = y_max_clip;
      }
   else
      {
      /*
       * Find the intersection with the left-edge of the window.
       * P1 is assumed to lie outside the window.
       */
      x_intersect = (1.0 / slope) * (y_max_clip - P1->y) + P1->x;
      P1->y = y_max_clip;
      P1->x = x_intersect;
      }
   }

}  /* for (i = 1; i <= 2...  */


/*
 * At this point, we should have cut off all of the extra bits
 * and have a completely visible segment remaining.
 */
return (ALL_VISIBLE);

}  /* end of ssp_ClipLine() */


/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_ClearLine (int linenum)

#else

int ssp_ClearLine (linenum)
int linenum;

#endif

/*
 * Purpose...
 * -------
 * Clear a line in the text region of the screen.
 *
 * Input...
 * -----
 * linenum : number of the line to be cleared
 *           (should be in the range firstline <= linenum <= 25)
 */

{
int i;

if (linenum < FirstLine) return (1);
if (linenum > LastLine) return (1);

#if ( GRAPH_ENV == TURBOC_G && (TURBO_C || TOPSPEED_C) )
gotoxy (1,linenum);
for (i = 1; i <= 79; ++i) printf (" ");
#endif

#  if (GRAPH_ENV == TOPSPEED_G)
_settextposition (linenum, 1);
for (i = 1; i <= 79; ++i) printf (" ");
#  endif

#if (GRAPH_ENV == TEKTERM)
#  endif

return (0);
}

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_ClearText (void)

#else

int ssp_ClearText ()

#endif

/*
 * Purpose...
 * -------
 * Clear all of the lines in the text area.
 */

{
int i;
for (i = FirstLine; i <= LastLine; ++i) ssp_ClearLine (i);
return (0);
}

/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_Message (char *string)

#else

int ssp_Message (string)
char *string;

#endif

/*
 * Purpose...
 * -------
 * Write a message in the lines reserved for text.
 *
 * Input...
 * -----
 * string  : pointer to the message string
 *
 */

{
#if (GRAPH_ENV == TURBOC_G)
ssp_ClearLine (1);
ssp_ClearLine (2);
gotoxy (1,1);
#endif

#  if (GRAPH_ENV == TOPSPEED_G)
ssp_ClearLine (1);
ssp_ClearLine (2);
_settextposition (1, 1);
#  endif

#if (GRAPH_ENV == TEKTERM)
;
#endif

printf ("%s\n", string);

return (0);
}


/*----------------------------------------------------------------*/

#if (PROTO)

int ssp_WaitMessage (char *string)

#else

int ssp_WaitMessage (string)
char *string;

#endif

/*
 * Purpose...
 * -------
 * Write a message in the lines reserved for text and wait
 * for a keypress.
 *
 * Input...
 * -----
 * string  : pointer to the message string
 *
 */

{
#if (GRAPH_ENV == TURBOC_G)
ssp_ClearLine (1);
ssp_ClearLine (2);
gotoxy (1,1);
#endif

#if (GRAPH_ENV == TOPSPEED_G)
ssp_ClearLine (1);
ssp_ClearLine (2);
_settextposition (1, 1);
#endif

#if (GRAPH_ENV == TEKTERM)
#endif

printf ("%s\n", string);  /* The message itself... */
ssp_WaitPrompt ();

return (0);
}


/*----------------------------------------------------------------*/

#if (PROTO)

int ssp_WaitPrompt (void)

#else

int ssp_WaitPrompt ()

#endif

/*
 * Purpose...
 * -------
 * Write a prompt and wait for a keypress.
 *
 * Input...
 * -----
 * None
 *
 */

{
printf ("Press <c> to continue...");

#if (GRAPH_ENV == TURBOC_G || GRAPH_ENV == TOPSPEED_G)
   /* This does not need a carriage return */
   do {} while (getch() != 'c');
#else
   /* This requires a carriage return to put it into the buffer */
   do {} while (getchar() != 'c');
#endif

return (0);
}

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_Compact (double x, char *string)

#else

int ssp_Compact (x, string)
double x;
char *string;

#endif

/*
 * Purpose...
 * -------
 * Produce a compact string representation of a number.
 * This routine is mainly used for writing postscript coordinates
 * in points.
 *
 * Input...
 * -----
 * x    : value to be converted
 *
 * Output...
 * ------
 * string : pointer to the compacted string
 *
 * Revisions...
 * ---------
 * 02-Apr-94  : write only one decimal place for points.
 *
 */

{
int i, j;

for (i = 0; i < Lstr; ++i) string[i] = 0;
sprintf (string, "%9.1f", x);

/*
 * Strip off the leading spaces and fill trailing characters
 * with NULLs
 */
for (i = Lstr-1; i >= 1; --i)
   {
   if (string[0] == ' ')
      {
      /* shift characters left by 1 place */
      for (j = 0; j < Lstr - 1; ++j) string[j] = string[j+1];
      string[Lstr - 1] = 0;
      }
   else
      break;
   }

return (0);
}

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetDash (double dash, double space)

#else

int ssp_SetDash (dash, space)
double dash, space;

#endif

/*
 * Purpose...
 * -------
 * Set the dashed-line parameters.
 *
 * Input...
 * -----
 * dash  : length of dashes in current page units, mm
 * space : length of the spaces between dashes in mm
 *
 */

{
dash_length_mm = dash;
space_length_mm = space;
combined_length_mm = dash_length_mm + space_length_mm;

return (0);
}

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetFillGrey (double lightness)

#else

int ssp_SetFillGrey (lightness)
double lightness;

#endif

/*
 * Purpose...
 * -------
 * Set the grey level for subsequent shading
 *
 * Input...
 * -----
 * lightness : 0.0=black, 1.0=white.
 *
 */

{
if (lightness <= 0.0)
   fill_grey_level = 0.0;
else if (lightness >= 1.0)
   fill_grey_level = 1.0;
else
   fill_grey_level = lightness;

return (0);
}


/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetHSBColour (double h, double s, double b)

#else

int ssp_SetHSBColour (h, s, b)
double h, s, b;

#endif

/*
 * Purpose...
 * -------
 * Set the parameters for the HSB colour model as used
 * in postscript.  The VGA colour is then selected to be the
 * nearest match.
 *
 * Input...
 * -----
 * h   : 0.0 <= h <= 1.0, hue
 *       0.0 = pure red, 1/3 = pure green, 2/3 = pure blue
 * s   : 0.0 <= s <= 1.0, saturation
 *       0.0 = grey (or no colour), 1.0 = maximum colour concentration
 * b   : 0.0 <= b <= 1.0, brightness
 *       0.0 = black, 0.5 = colours, 1.0 = white
 *
 */

{
/*
 * Bring values to within range; 0.0 <= h < 1.0
 */
if (h < 0.0) {
   hue = h + 1.0 + ((int) h);
} else if (h >= 1.0) {
   hue = h - ((int) h);
} else {
   hue = h;
} /* end if */

if (s <= 0.0) {
   saturation = 0.0;
} else if (s >= 1.0) {
   saturation = 1.0;
} else {
   saturation = s;
} /* end if */

if (b <= 0.0) {
   brightness = 0.0;
} else if (b >= 1.0) {
   brightness = 1.0;
} else {
   brightness = b;
} /* end if */

/*
 * Select VGA colour values.
 */
if (saturation < 0.1) {
   /* Set a grey */
   if (brightness < 0.3) {
      vga_colour = VGA_DARKGRAY;
   } else if (brightness < 0.7) {
      vga_colour = VGA_LIGHTGRAY;
   } else {
      vga_colour = VGA_WHITE;
   } /* end if */
} else {
   /* Set a colour */
   if (hue <= 30.0/360.0 || hue > 330.0/360.0) {
      vga_colour = VGA_RED;
   } else if (hue <= 90.0/360.0) {
      vga_colour = VGA_YELLOW;
   } else if (hue <= 150.0/360.0) {
      vga_colour = VGA_GREEN;
   } else if (hue <= 210.0/360.0) {
      vga_colour = VGA_CYAN;
   } else if (hue <= 270.0/360.0) {
      vga_colour = VGA_BLUE;
   } else {
      vga_colour = VGA_MAGENTA;
   } /* end if */
} /* end if */

#  if (GRAPH_ENV == TURBOC_G)
   setcolor ( vga_colour );
#  endif

#  if (GRAPH_ENV == TOPSPEED_G)
   _setcolor ( vga_colour );
#  endif

/*
 * Select GIF colour values.
 */
if (saturation < 0.1) {
   /* Set a grey */
   igrey = (int) (brightness * MAX_GIF_GREY);      
   gif_colour = gif_greyshade[igrey];
} else {
   /* Set a colour */
   ihue = (int) (hue * MAX_GIF_HUE);      
   gif_colour = gif_hue[ihue];
} /* end if */

if (record_in_file)
   {
   if (outdevice == HPGL)
      { ; }
   else if (outdevice == PSCRIPT)
      {
      if (active_ps_path == 1)
         {
         /* Stroke the active path to finish it. */
         fprintf (plotfile, "s\n");
         active_ps_path = 0;
         }
      fprintf (plotfile, "%4.2f %4.2f %4.2f sc\n",
               hue, saturation, brightness);
      } /* postscript */
   }

in_colour = 1;  /* remember that we are now doing lines in colour */
return (0);
}


/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetHSBFillColour (double h, double s, double b)

#else

int ssp_SetHSBFillColour (h, s, b)
double h, s, b;

#endif

/*
 * Purpose...
 * -------
 * Set the parameters for the HSB colour model for polygon filling.
 *
 * Input...
 * -----
 * h   : 0.0 <= h <= 1.0, hue
 *       0.0 = pure red, 1/3 = pure green, 2/3 = pure blue
 * s   : 0.0 <= s <= 1.0, saturation
 *       0.0 = grey (or no colour), 1.0 = maximum colour concentration
 * b   : 0.0 <= b <= 1.0, brightness
 *       0.0 = black, 0.5 = colours, 1.0 = white
 *
 */

{
/*
 * Bring values to within range; 0.0 <= h < 1.0
 */
if (h < 0.0) {
   fill_hue = h + 1.0 + ((int) h);
} else if (h >= 1.0) {
   fill_hue = h - ((int) h);
} else {
   fill_hue = h;
} /* end if */

if (s <= 0.0) {
   fill_saturation = 0.0;
} else if (s >= 1.0) {
   fill_saturation = 1.0;
} else {
   fill_saturation = s;
} /* end if */

if (b <= 0.0) {
   fill_brightness = 0.0;
} else if (b >= 1.0) {
   fill_brightness = 1.0;
} else {
   fill_brightness = b;
} /* end if */

/*
 * Select VGA colour value.
 */
if (fill_saturation < 0.1) {
   /* Set a grey */
   if (fill_brightness < 0.3) {
      fill_vga_colour = VGA_DARKGRAY;
   } else if (fill_brightness < 0.7) {
      fill_vga_colour = VGA_LIGHTGRAY;
   } else {
      fill_vga_colour = VGA_WHITE;
   }
} else {
   /* Set a colour */
   if (fill_hue <= 30.0/360.0 || fill_hue > 330.0/360.0) {
      fill_vga_colour = VGA_RED;
   } else if (fill_hue <= 90.0/360.0) {
      fill_vga_colour = VGA_YELLOW;
   } else if (fill_hue <= 150.0/360.0) {
      fill_vga_colour = VGA_GREEN;
   } else if (fill_hue <= 210.0/360.0) {
      fill_vga_colour = VGA_CYAN;
   } else if (fill_hue <= 270.0/360.0) {
      fill_vga_colour = VGA_BLUE;
   } else {
      fill_vga_colour = VGA_MAGENTA;
   } /* end if */
} /* end if */

/*
 * Select GIF colour values.
 */
if (fill_saturation < 0.1) {
   /* Set a grey */
   igrey = (int) (brightness * MAX_GIF_GREY);      
   fill_gif_colour = gif_greyshade[igrey];
} else {
   /* Set a colour */
   ihue = (int) (fill_hue * MAX_GIF_HUE);      
   fill_gif_colour = gif_hue[ihue];
} /* end if */

return (0);
}


/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetNoColour (void)

#else

int ssp_SetNoColour ()

#endif

/*
 * Purpose...
 * -------
 * Set the output colour back to black for postscript and
 * white for VGA displays.
 *
 */

{

hue        = 0.0;
saturation = 0.0;
brightness = 1.0;
vga_colour = VGA_WHITE;

fill_hue        = 0.0;
fill_saturation = 0.0;
fill_brightness = 1.0;
fill_vga_colour = VGA_WHITE;

#  if (GRAPH_ENV == TURBOC_G)
   setcolor ( vga_colour );
#  endif

#  if (GRAPH_ENV == TOPSPEED_G)
   _setcolor ( vga_colour );
#  endif

if (record_in_file)
   {
   if (outdevice == HPGL)
      { ; }
   else if (outdevice == PSCRIPT)
      {
      if (active_ps_path == 1)
         {
         /* Stroke the active path to finish it. */
         fprintf (plotfile, "s\n");
         active_ps_path = 0;
         }
      fprintf (plotfile, "0.0 sg\n");
      } /* postscript */
   else if (outdevice == GIF)
      {
      gif_colour = gif_black;
      }
   }

in_colour = 0;  /* remember that we are now doing lines in B/W */
return (0);
}


/*------------------------------------------------------------------*/

#if (PROTO)

double ssp_MapToColour (double level, double minimum, double maximum)

#else

double ssp_MapToColour (level, minimum, maximum)
double level, minimum, maximum;

#endif

/*
 * Purpose...
 * -------
 * Map a floating point value (between minimum and maximum values)
 * to a colour hue using the HSV colour model.
 * The minimum value corresponds to magenta and the maximum value
 * corresponds to red (moving from HUE_MAX back around to 0.0).
 *
 * Input...
 * -----
 * level   : double value  minimum <= level <= maximum
 * minimum : double value
 * maximum : double value
 *
 * Output...
 * ------
 * ssp_MapToColour() returns a double value specifying the colour
 * hue where 0.0 <= hue <= 1.0.
 * This range of values is suitable for postscript.
 *
 */

{
double hue_local;
/* 0.7 should correspond to blue-magenta */
#define  HUE_MAX  0.7

if ( fabs(maximum - minimum) < 1.0e-12 ) return (0.0);
if (level < minimum) return (HUE_MAX);;
if (level > maximum) return (0.0);

hue_local = HUE_MAX * (1.0 - (level - minimum) / (maximum - minimum));
return (hue_local);
}

/*-----------------------------------------------------------------*/
 
int hsv2rgb( struct hsv_colour *hsv, struct rgb_colour *rgb ) {
   /*
    * Purpose:
    * Convert Hue Saturation Value to Red Green Blue
    * All values are in the range [0.0 .. 1.0]
    *
    * Reference:
    * D. F. Rogers
    * Procedural Elements for Computer Graphics
    * McGraw Hill 1985
    */
   float S, H, V, F, M, N, K;
   int   I;
   
   S = hsv->s;  /* Saturation */
   H = hsv->h;  /* Hue */
   V = hsv->v;  /* value or brightness */
   
   if ( S == 0.0 ) {
      /* 
       * Achromatic case, set level of grey 
       */
      rgb->r = V;
      rgb->g = V;
      rgb->b = V;
   } else {
      /* 
       * Determine levels of primary colours. 
       */
      if (H >= 1.0) {
         H = 0.0;
      } else {
         H = H * 6;
      } /* end if */
      I = (int) H;   /* should be in the range 0..5 */
      F = H - I;     /* fractional part */

      M = V * (1 - S);
      N = V * (1 - S * F);
      K = V * (1 - S * (1 - F));

      if (I == 0) { rgb->r = V; rgb->g = K; rgb->b = M; }
      if (I == 1) { rgb->r = N; rgb->g = V; rgb->b = M; }
      if (I == 2) { rgb->r = M; rgb->g = V; rgb->b = K; }
      if (I == 3) { rgb->r = M; rgb->g = N; rgb->b = V; }
      if (I == 4) { rgb->r = K; rgb->g = M; rgb->b = V; }
      if (I == 5) { rgb->r = V; rgb->g = M; rgb->b = N; }
   } /* end if */

   return 0;
} /* end function hsv2rgb */

/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_SetLineThickness (double thickness)

#else

int ssp_SetLineThickness (thickness)
double thickness;

#endif

/*
 * Purpose...
 * -------
 * Set the line thickness (in mm).
 *
 * Input...
 * -----
 * thickness : the desired thickness in mm
 *
 */

{
if (thickness < 0.1) thickness = 0.1;
if (thickness > 10.0) thickness = 10.0;

if (display_on_screen)
   {
   ;
   }

if (record_in_file)
   {
   if (outdevice == HPGL)
      {
      /* For the laserjet, only pen 1 is functional */
      fprintf (plotfile, "SP%d;PT%f;\n", 1, thickness);
      fprintf (plotfile, "PW%f;\n", thickness);
      }
   else if (outdevice == PSCRIPT)
      {
      if (active_ps_path == 1)
         {
         /* Stroke the active path to finish it. */
         fprintf (plotfile, "s\n");
         active_ps_path = 0;
         }
      fprintf (plotfile, "%9.2f mm setlinewidth\n", thickness);
      }
   }

return (0);
}

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_Where (double *x, double *y,
	       double *scale_factor,
	       double *angle)

#else

int ssp_Where (x, y, scale_factor, angle)
double *x, *y, *scale_factor, *angle;

#endif

/*
 * Purpose...
 * -------
 * Tell the user where the current position is located.
 *
 * Output...
 * ------
 * x, y     : the present position (plot units)
 * scale_factor : scale_factor from plot units to page units
 * angle    : angle of plot frame-of-reference wrt the x-axis of
 *            the page
 *
 */

{
*x = (x_old - x_origin) / plot_scale;
*y = (y_old - y_origin) / plot_scale;

*scale_factor = plot_scale;
*angle = 180.0 / ssp_PI * global_angle;

return (0);
}


/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_Plot (double x, double y)

#else

int ssp_Plot (x, y)
double x, y;

#endif

/*
 * Purpose...
 * -------
 * Plot a line on the screen.  The line is plotted from
 * the previous point to the new point (x, y).
 *
 * Input...
 * -----
 * x, y  : coordinates of the new point (in plot units)
 *
 */

{
double x1, y1, x2, y2;
int    xp1, yp1, xp2, yp2;
char   xstr[Lstr], ystr[Lstr];
struct point_2D P1, P2;
int    visible;

if ( display_on_screen && (GraphicsOn != 1) )
   ssp_WaitMessage ("We are not currently if graphics mode...");

x1 = x_old;
y1 = y_old;    /* the old coordinates are already scaled */

/*
 * Convert the new coordinates to page units
 */
x2 = x * plot_scale + x_origin;
y2 = y * plot_scale + y_origin;

/*
 * Move the current point
 */
x_old = x2;
y_old = y2;

if (do_clipping == 1)
   {
   P1.x = x1; P1.y = y1;
   P2.x = x2; P2.y = y2;
   visible = ssp_ClipLine (&P1, &P2);
   if (visible == NOT_VISIBLE) return (0);
   x1 = P1.x; y1 = P1.y;
   x2 = P2.x; y2 = P2.y;
   }

if (display_on_screen)
   {
   xp1 = (int) (x1 / XPlotDim * XPlotPixels);
   yp1 = (int) ((1.0 - y1 / YPlotDim) * YPlotPixels);
   xp2 = (int) (x2 / XPlotDim * XPlotPixels);
   yp2 = (int) ((1.0 - y2 / YPlotDim) * YPlotPixels);

#  if (GRAPH_ENV == TURBOC_G)
   line (xp1, yp1, xp2, yp2);
   ErrorCode = graphresult ();
   if (ErrorCode != grOk)
      {
      printf (" Graphics Error in Plot : %s\n",
      grapherrormsg (ErrorCode));
      exit (1);
      }
#  endif

#  if (GRAPH_ENV == TOPSPEED_G)
   _moveto (xp1, yp1);
   _lineto (xp2, yp2);
#  endif

#  if (GRAPH_ENV == TEKTERM)
   xp1 = x1 / XPlotDim * 1024;
   yp1 = y1 / YPlotDim * 780;
   xp2 = x2 / XPlotDim * 1024;
   yp2 = y2 / YPlotDim * 780;

   if (xp1 > 1023) xp1 = 1023;
   if (xp1 < 0) xp1 = 0;
   if (yp1 > 779) yp1 = 779;
   if (yp1 < 0) yp1 = 0;

   hi_y = yp1 / 32;
   lo_y = yp1 - (hi_y * 32);
   hi_x = xp1 / 32;
   lo_x = xp1 - (hi_x * 32);
   hi_y += 32;
   lo_y += 96;
   hi_x += 32;
   lo_x += 64;
   printf ("%c%c%c%c%c", GS, hi_y, lo_y, hi_x, lo_x);

   if (xp2 > 1023) xp2 = 1023;
   if (xp2 < 0) xp2 = 0;
   if (yp2 > 779) yp2 = 779;
   if (yp2 < 0) yp2 = 0;

   hi_y = yp2 / 32;
   lo_y = yp2 - (hi_y * 32);
   hi_x = xp2 / 32;
   lo_x = xp2 - (hi_x * 32);

   hi_y += 32;
   lo_y += 96;
   hi_x += 32;
   lo_x += 64;
   printf ("%c%c%c%c", hi_y, lo_y, hi_x, lo_x);

#  endif

   }


if (record_in_file)
   {
   if (outdevice == HPGL)
      {
      fprintf (plotfile, "PD;PA%d,%d;\n",
	       (int)(x2 * hp2mm), (int)(y2 * hp2mm) );
      }
   else if (outdevice == PSCRIPT)
      {
      if (active_ps_path == 0)
         {
         /* Start a new active path. */
         active_ps_path = 1;
         /* First, move to the old point. */
         ssp_Compact (x1 * points_per_mm, xstr);
         ssp_Compact (y1 * points_per_mm, ystr);
         fprintf (plotfile, "n %s %s m ", xstr, ystr);
         }
      /* Now, line_to the new point. */
      ssp_Compact (x2 * points_per_mm, xstr);
      ssp_Compact (y2 * points_per_mm, ystr);
      fprintf (plotfile, "%s %s l\n", xstr, ystr);
      } /* end of postscript */
   else if (outdevice == GIF)
      {
      xp1 = (int) (x1 / XPlotDim * XPlotPixels);
      yp1 = (int) ((1.0 - y1 / YPlotDim) * YPlotPixels);
      xp2 = (int) (x2 / XPlotDim * XPlotPixels);
      yp2 = (int) ((1.0 - y2 / YPlotDim) * YPlotPixels);
      gdImageLine (gif_image, xp1, yp1, xp2, yp2, gif_colour);   
      } /* end of GIF */
   }

return (0);
}   /* end of ssp_Plot() */

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_Move (double x, double y)

#else

int ssp_Move (x, y)
double x, y;

#endif

/*
 * Purpose...
 * -------
 * Move to a new point on the screen.
 *
 * Input...
 * -----
 * x, y  : coordinates of the new point (in plot units)
 *
 */

{
double x1, y1, x2, y2;
int    xp1, yp1, xp2, yp2;
char   xstr[Lstr], ystr[Lstr];
struct point_2D P1, P2;
int    visible;

/*
 * Pick up the old point.
 */
x1 = x_old;
y1 = y_old;

/*
 * Convert the new coordinates to page units
 */
x2 = x * plot_scale + x_origin;
y2 = y * plot_scale + y_origin;

/*
 * Move the current point
 */
x_old = x2;
y_old = y2;


if (do_clipping == 1)
   {
   P1.x = x1; P1.y = y1;
   P2.x = x2; P2.y = y2;
   visible = ssp_ClipLine (&P1, &P2);
   if (visible == NOT_VISIBLE) return (0);
   x1 = P1.x; y1 = P1.y;
   x2 = P2.x; y2 = P2.y;
   }

if (record_in_file)
   {
   if (outdevice == HPGL)
      {
      fprintf (plotfile, "PU;PA%d,%d;\n",
	       (int)(x2 * hp2mm), (int)(y2 * hp2mm) );
      }
   else if (outdevice == PSCRIPT)
      {
      if (active_ps_path == 1)
         {
         /* Stroke the active path to finish it. */
         fprintf (plotfile, "s\n");
         active_ps_path = 0;
         }
      /*
       * Otherwise do nothing as the current pen position is remembered by
       * the ssp package.
       */
      }
   }

if (display_on_screen)
   {
   ;
   }

return (0);
}   /* end of ssp_Move() */

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_Polygon (int n, double x[], double y[],
                 int grey_fill, int colour_fill, int border_flag)

#else

int ssp_Polygon (n, x, y, grey_fill, colour_fill, border_flag)
int    n;
double x[], y[];
int    grey_fill, colour_fill, border_flag;

#endif

/*
 * Purpose...
 * -------
 * Draw a polygon and optionally fill it.
 *
 * Input...
 * -----
 * n         : number of vertices (maximum 10)
 * x[], y[]  : coordinates of the vertices (in plot units)
 * grey_fill : = 0, don't fill the polygon with grey
 *             = 1, fill the interior of the polygon with grey
 * colour_fill : = 0, don't fill the polygon with colour
 *               = 1, fill the interior of the polygon with colour
 * border_flag : = 0, don't draw the border
 *               = 1, draw the border
 *
 */

{
int    i, flag;
double xpage[10], ypage[10];
short  xp[10], yp[10];
int    xyp[20];
gdPoint p[10];
char   xstr[Lstr], ystr[Lstr];
int    vga_grey_level, gif_grey_level;

if ( display_on_screen && (GraphicsOn != 1) )
   ssp_WaitMessage ("We are not currently if graphics mode...");

/*
 * Check the number of vertices.
 */
if (n <= 2) return (0);
if (n > 10) n = 10;

/*
 * Convert the new coordinates to page units
 */
for (i = 0; i < n; ++i)
   {
   xpage[i] = x[i] * plot_scale + x_origin;
   ypage[i] = y[i] * plot_scale + y_origin;
   }

/*
 * Move the current point
 */
x_old = xpage[0];
y_old = ypage[0];

if (do_clipping == 1)
   {
   /* Need to think carefully about polygon clipping */
   /* For the moment, let the hardware look after it */
   }

if (display_on_screen)
   {
   for (i = 0; i < n; ++i)
      {
      xp[i] = (int) (xpage[i] / XPlotDim * XPlotPixels);
      yp[i] = (int) ((1.0 - ypage[i] / YPlotDim) * YPlotPixels);
      xyp[2*i] = xp[i];
      xyp[2*i+1] = yp[i];
      }

#  if (GRAPH_ENV == TURBOC_G)
   /* The grey-filled polygon. */
   if (grey_fill == 1)
      {
      if (fill_grey_level < 0.10)
         vga_grey_level = VGA_BLACK;
      else if (fill_grey_level < 0.50)
         vga_grey_level =  VGA_DARKGRAY;
      else if (fill_grey_level < 0.90)
         vga_grey_level = VGA_LIGHTGRAY;
      else
         vga_grey_level = VGA_WHITE;

      setfillstyle ( SOLID_FILL, vga_grey_level );
      setcolor ( vga_grey_level );  /* matching boundary */
      fillpoly (n, xyp);
      setcolor ( vga_colour );
      }

   /* Colour-filled region */
   if (colour_fill == 1)
      {
      setfillstyle ( SOLID_FILL, fill_vga_colour );
      setcolor ( fill_vga_colour );  /* matching boundary */
      fillpoly (n, xyp);
      setcolor ( vga_colour );
      }

   /* The border. */
   if (border_flag == 1)
      {
      setcolor ( vga_colour );
      for (i = 1; i <n; ++i)
         line (xp[i-1], yp[i-1], xp[i], yp[i]);
      line (xp[n-1], yp[n-1], xp[0], yp[0]);
      ErrorCode = graphresult ();
      if (ErrorCode != grOk)
         {
         printf (" Graphics Error in ssp_Polygon : %s\n",
         grapherrormsg (ErrorCode));
         exit (1);
         }
      }
#  endif

#  if (GRAPH_ENV == TOPSPEED_G)
   /* The grey-filled polygon. */
   if (grey_fill == 1)
      {
      /* Set the grey colour */
      if (fill_grey_level < 0.10)
         vga_grey_level = VGA_BLACK;
      else if (fill_grey_level < 0.50)
         vga_grey_level =  VGA_DARKGRAY;
      else if (fill_grey_level < 0.90)
         vga_grey_level = VGA_LIGHTGRAY;
      else
         vga_grey_level = VGA_WHITE;

      /* Now, fill the polygon. */
      _setfillmask ( fill_mask_10 );
      _setcolor ( vga_grey_level );
      flag = _GFILLINTERIOR;
      _polygon ( flag, n, xp, yp );
      _setcolor ( vga_colour );
      }

   /* The colour-filled polygon. */
   if (colour_fill == 1)
      {
      /* Set the intensity to solid fill */
      _setfillmask (fill_mask_10);
      /* Now, fill the polygon. */
      _setcolor ( fill_vga_colour );
      flag = _GFILLINTERIOR;
      _polygon (flag, n, xp, yp);
      _setcolor ( vga_colour );
      }

   /* The border. */
   if (border_flag == 1)
      {
      _setcolor ( vga_colour );
      _moveto (xp[0], yp[0]);
      for (i = 1; i < n; ++i)
         _lineto (xp[i], yp[i]);
      _lineto (xp[0], yp[0]);
      }
#  endif

#  if (GRAPH_ENV == TEKTERM)
#  endif

   } /* if display_on-screen */


if (record_in_file)
   {
   if (outdevice == HPGL)
      {
      fprintf (plotfile, "PU;PA%d,%d;\n",
	       (int)(xpage[0] * hp2mm), (int)(ypage[0] * hp2mm) );
      for (i = 1; i < n; ++i)
         {
         fprintf (plotfile, "PD;PA%d,%d;\n",
	       (int)(xpage[i] * hp2mm), (int)(ypage[i] * hp2mm) );
	 }
      fprintf (plotfile, "PD;PA%d,%d;\n",
	       (int)(xpage[0] * hp2mm), (int)(ypage[0] * hp2mm) );
      }  /* if (outdevice == HPGL... */

   else if (outdevice == PSCRIPT)
      {
      if (active_ps_path == 1)
         {
         /* Stroke the active path to finish it. */
         fprintf (plotfile, "s\n");
         active_ps_path = 0;
         }

      /* Set up the path defining the border. */
      ssp_Compact (xpage[0] * points_per_mm, xstr);
      ssp_Compact (ypage[0] * points_per_mm, ystr);
      fprintf (plotfile, "n %s %s m\n", xstr, ystr);
      for (i = 1; i < n; ++i)
         {
         ssp_Compact (xpage[i] * points_per_mm, xstr);
         ssp_Compact (ypage[i] * points_per_mm, ystr);
         fprintf (plotfile, "%s %s l\n", xstr, ystr);
         }
      fprintf (plotfile, "cp\n");
      fprintf (plotfile, "gs\n");

      if (grey_fill == 1)
         {
         fprintf (plotfile, "%4.2f sg\n", fill_grey_level);
         fprintf (plotfile, "f\n");
         }

      if (colour_fill == 1)
         {
         fprintf (plotfile, "%4.2f %4.2f %4.2f sc\n",
                  fill_hue, fill_saturation, fill_brightness);
         fprintf (plotfile, "f\n");
         }

      if (border_flag == 1)
         {
         /* Restore the stroke path from the stack. */
         if (in_colour == 1)
            fprintf (plotfile, "%4.2f %4.2f %4.2f sc\n",
                     hue, saturation, brightness);
         else
            fprintf (plotfile, "%4.2f sg\n", 0.0);
         fprintf (plotfile, "gr\n");
         fprintf (plotfile, "s\n");
         }
      else
         {
         /* pop the stroke path off the stack and discard */
         fprintf (plotfile, "gr\n");
         }
      }  /* if ... postscript... */

   else if (outdevice == GIF)
      {
      for (i = 0; i < n; ++i)
         {
         p[i].x = (int) (xpage[i] / XPlotDim * XPlotPixels);
         p[i].y = (int) ((1.0 - ypage[i] / YPlotDim) * YPlotPixels);
         }

      /* The grey-filled polygon. */

      if (grey_fill == 1)
         {
         /* Set the grey colour */
         gif_grey_level = (int) (fill_grey_level * MAX_GIF_GREY);

         /* Now, fill the polygon. */
         gdImageFilledPolygon(gif_image, p, n, 
            gif_greyshade[gif_grey_level]);
         }

      /* The colour-filled polygon. */
      if (colour_fill == 1 )
         {
         gdImageFilledPolygon(gif_image, p, n, fill_gif_colour);
         }

      /* The border. */
      if (border_flag == 1)
         {
         for (i = 1; i < n; ++i)
            gdImageLine (gif_image, p[i-1].x, p[i-1].y, 
                         p[i].x, p[i].y, gif_colour);
         gdImageLine (gif_image, p[n-1].x, p[n-1].y, 
                      p[0].x, p[0].y, gif_colour);
         }
      } /* if ... GIF */

   }  /* if (record_in_file... */

return (0);
}   /* end of ssp_Polygon() */

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_SetOrigin (double x, double y)

#else

int ssp_SetOrigin (x, y)
double x, y;

#endif

/*
 * Purpose ...
 * -------
 * Move to a new point on the screen and reset the origin to it.
 */

{
double x2, y2;

/*
 * Convert the new coordinates to page units
 */
x2 = x * plot_scale + x_origin;
y2 = y * plot_scale + y_origin;

/*
 * Move the current point
 */
x_old = x2;
y_old = y2;

/*
 * Reset the origin
 */
x_origin = x2;
y_origin = y2;


if (record_in_file)
   {
   if (outdevice == HPGL)
      {
      fprintf (plotfile, "PU;PA%d,%d;\n",
	       (int)(x2 * hp2mm), (int)(y2 * hp2mm) );
      }
   else if (outdevice == PSCRIPT)
      {
      if (active_ps_path == 1)
         {
         /* Stroke the active path to finish it. */
         fprintf (plotfile, "s\n");
         active_ps_path = 0;
         }
      }
   }

return (0);
}   /* end of ssp_SetOrigin() */

/*----------------------------------------------------------------*/

#if (PROTO)

int ssp_Rotate (double angle)

#else

int ssp_Rotate (angle)
double angle;

#endif

/*
 * Purpose...
 * -------
 * Rotate the present frame of reference (with respect to the
 * x-axis of the plotted page).
 *
 * Input...
 * -----
 * angle : rotation angle in degrees
 *
 */

{
global_angle += (ssp_PI / 180.0 * angle);
return (0);
}

/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_Factor (double zoom_factor)

#else

int ssp_Factor (zoom_factor)
double zoom_factor;

#endif

/*
 * Purpose...
 * -------
 * Change the global plot scale.
 *
 * Input...
 * -----
 * zoom_factor : multiplier for the global scale
 */

{
if (zoom_factor <= 0.0) return (1);   /* nonsense value */

plot_scale *= zoom_factor;
return (0);
}
/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_TextOn (void)

#else

int ssp_TextOn ()

#endif

/*
 * Purpose...
 * -------
 * Put the current graphics screen away in memory buffer.
 * Go from graphics mode to text mode.
 */

{
int xp1, yp1;

if (!display_on_screen) return (0);

if (GraphicsOn != 1)
   {
   ssp_WaitMessage ("We are not currently if graphics mode...");
   return (0);
   }

#if (GRAPH_ENV == TURBOC_G)
if (buffer != NULL)
   getimage (0, 0, XPlotPixels, YPlotPixels, buffer);
restorecrtmode();
if (buffer == NULL) printf ("Beware -- graphics screen not saved.");
#endif

#if (GRAPH_ENV == TOPSPEED_G)
_setvideomode (_DEFAULTMODE);
if (buffer == NULL) printf ("Beware -- graphics screen not saved.");
#endif

#if (GRAPH_ENV == TEKTERM)
printf ("%c%c", ESC, FF);

xp1 = 0;
yp1 = 780 - 1;

hi_y = yp1 / 32;
lo_y = yp1 - (hi_y * 32);
hi_x = xp1 / 32;
lo_x = xp1 - (hi_x * 32);

hi_y += 32;
lo_y += 96;
hi_x += 32;
lo_x += 64;

printf ("%c%c%c%c%c", GS, hi_y, lo_y, hi_x, lo_x);
printf ("%c", US);
#  endif


GraphicsOn = 0;
return (0);
}
/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_TextOff(void)

#else

int ssp_TextOff()

#endif

/*
 * Purpose...
 * -------
 * Go from text mode back to graphics mode.
 * Restore the current graphics screen from memory buffer.
 */

{
int xp1, yp1;

if (!display_on_screen) return (0);

#if (GRAPH_ENV == TURBOC_G)
   setgraphmode ( getgraphmode() );
   ErrorCode = graphresult ();
   if (ErrorCode != grOk)
      {
      printf ("Graphics Error in ssp_TextOff() : %s\n",
   	   grapherrormsg (ErrorCode));
      exit (1);
      }
   cleardevice ();
   setcolor (getmaxcolor());
   setviewport (0, YMaxPixels - YPlotPixels,
	        XPlotPixels, YMaxPixels, 1);
   rectangle (0, 0, XPlotPixels, YPlotPixels);
   if (buffer != NULL)
      {
      putimage (0, 0, buffer, COPY_PUT);
      }
   if (ErrorCode != grOk)
      {
      printf (" Graphics Error in ssp_TextOff() : %s\n",
	      grapherrormsg (ErrorCode));
      exit (1);
      }
#endif

#if (GRAPH_ENV == TOPSPEED_G)
   _setvideomode (_VRES16COLOR);
   XPlotPixels = (int) (XMaxPixels * XScrnDim);
   YPlotPixels = (int) (YMaxPixels * YScrnDim);
   _setviewport (0, YMaxPixels - YPlotPixels,
		XPlotPixels, YMaxPixels);
   _rectangle (_GBORDER, 0, 0, XPlotPixels, YPlotPixels);

   /*
    * Lines FirstLine through LastLine are available for text.
    */
   _settextposition (FirstLine, 1);
#endif

#if (GRAPH_ENV == TEKTERM)
   printf ("%c%c", ESC, FF);

   xp1 = 0;
   yp1 = 779;

   hi_y = yp1 / 32;
   lo_y = yp1 - (hi_y * 32);
   hi_x = xp1 / 32;
   lo_x = xp1 - (hi_x * 32);

   hi_y += 32;
   lo_y += 96;
   hi_x += 32;
   lo_x += 64;
   printf ("%c%c%c%c%c", GS, hi_y, lo_y, hi_x, lo_x);
#endif

GraphicsOn = 1;
return (0);
}
/*----------------------------------------------------------------*/

#if (PROTO)

int ssp_PlotText (double x, double y,
		  double size, double directn,
		  char *txt,
		  int n)

#else

int ssp_PlotText (x, y, size, directn, txt, n)
double x, y, size;
double directn;
char   *txt;
int    n;

#endif

/*
 * Purpose...
 * -------
 * Put some text onto the plot.
 *
 * Input...
 * -----
 * x, y   : starting coordinates of the lower left corner of the text
 * size   : height of the capitalized characters
 * directn : angle with the x-axis (degrees)
 * txt     : pointer to the text string
 * n       : number of characters to plot
 */

{
int    i;
double angle, scale;

/*
 * Write out the first n characters only.
 */
for (i = 0; i < n; ++i)
   {
   if (txt[i] == 0) break;
   ssp_DrawGlyph (txt[i], x, y, size, directn);
   ssp_Where (&x, &y, &scale, &angle);
   }

/*
 * On return, the current position should be at the end of the text.
 */
return (0);
}  /* end of ssp_PlotText() */

/*------------------------------------------------------------------*/

#if (PROTO)

int ssp_DashLine (double *x, double *y, int n)

#else

int ssp_DashLine (x, y, n)
double *x, *y;
int    n;

#endif

/*
 * Purpose...
 * -------
 * Draw a dashed line through the array of points.
 *
 * Input...
 * -----
 * x    : pointer to the array of x-coordinates, 0...n-1
 * y    : pointer to the array of y-coordinates, 0...n-1
 * n    : number of points in the arrays
 *
 */

{  /*  begin ssp_DashLine()  */
double dash_length, space_length, combined_length;
double segment_left;
double remainder;
double remaining_dash, remaining_space;
double x1, x2, y1, y2, dx, dy, xx, yy;
double current_scale, current_angle;
double sine, cosine;
int    i;

/*
 * Get the dash pattern in current plotting units.
 */
ssp_Where (&x1, &y1, &current_scale, &current_angle);
dash_length     = dash_length_mm / current_scale;
space_length    = space_length_mm / current_scale;
combined_length = combined_length_mm / current_scale;

/*
 * Start with a fresh dash.
 */
remainder = 0.0;

for (i = 0; i < n-1; ++i)
   {
   /*
    * Draw a dashed line between two points.
    */
   x1 = x[i];   y1 = y[i];
   x2 = x[i+1]; y2 = y[i+1];
   dx = x2 - x1;
   dy = y2 - y1;
   segment_left = sqrt(dx * dx + dy * dy);
   sine = dy / segment_left;
   cosine = dx / segment_left;
   xx = x1;
   yy = y1;
   ssp_Move (xx, yy);

   while (segment_left > 0.0)
      {
      if (remainder <= 0.0)
	 {
	 /*
	  * Start a new dash-space combination.
	  */
	 if (segment_left >= combined_length)
	    {
	    /*
	     * Draw a dash and space
	     */
	    xx += cosine * dash_length;
	    yy += sine * dash_length;
	    ssp_Plot (xx, yy);
	    segment_left -= dash_length;
	    xx += cosine * space_length;
	    yy += sine * space_length;
	    ssp_Move (xx, yy);
	    segment_left -= space_length;
	    remainder = 0.0;
	    }
	 else if (segment_left >= dash_length)
	    {
	    /*
	     * Draw a full dash and part of a space
	     */
	    xx += cosine * dash_length;
	    yy += sine * dash_length;
	    ssp_Plot (xx, yy);
	    segment_left -= dash_length;
	    xx += cosine * segment_left;
	    yy += sine * segment_left;
	    ssp_Move (xx, yy);
	    remainder = space_length - segment_left;
	    segment_left = 0.0;
	    }
	 else
	    {
	    /*
	     * Draw part of the dash only
	     */
	    xx += cosine * segment_left;
	    yy += sine * segment_left;
	    ssp_Plot (xx, yy);
	    remainder = combined_length - segment_left;
	    segment_left = 0.0;
	    }
	 }  /* end of new combination */

      else  /* remainder > 0 */
	 {
	 /*
	  * Finish off an old dash-space combination.
	  */
	 if (segment_left > remainder)
	    {
	    /*
	     * Plot all of the remaining part of the dash-space.
	     */
	    remaining_dash = remainder - space_length;
	    if (remaining_dash > 0.0)
	       {
	       xx += cosine * remaining_dash;
	       yy += sine * remaining_dash;
	       ssp_Plot (xx, yy);
	       segment_left -= remaining_dash;
	       remainder -= remaining_dash;
	       }
	    remaining_space = remainder;
	    xx += cosine * remaining_space;
	    yy += sine * remaining_space;
	    ssp_Move (xx, yy);
	    segment_left -= remaining_space;
	    remainder = 0.0;
	    }
	 else if (segment_left < (remainder - space_length) )
	    {
	    /*
	     * Plot part of the dash only, there is no room for
	     * anything else.
	     */
	    xx += cosine * segment_left;
	    yy += sine * segment_left;
	    ssp_Plot (xx, yy);
	    remainder -= segment_left;
	    segment_left = 0.0;
	    }
	 else
	    {
	    /*
	     * Plot the rest of the dash and some of the space.
	     */
	    remaining_dash = remainder - space_length;
	    if (remaining_dash > 0.0)
	       {
	       xx += cosine * remaining_dash;
	       yy += sine * remaining_dash;
	       ssp_Plot (xx, yy);
	       segment_left -= remaining_dash;
	       remainder -= remaining_dash;
	       }
	    xx += cosine * segment_left;
	    yy += sine * segment_left;
	    ssp_Move (xx, yy);
	    remainder -= segment_left;
	    segment_left = 0.0;
	    }
	 }

      }  /* end of while() */
   }     /* end of for(i...) */

return (0);
}  /*  end of ssp_DashLine()  */

/*-----------------------------------------------------------------*/

#if PROTO

int ssp_AllocGlyph (void)

#else

int ssp_AllocGlyph ()

#endif

/*
 * Purpose...
 * -------
 * Allocating memory for the encoded glyph strings.
 *
 */

{  /* begin ssp_AllocGlyph() */
int i;

glyph_stroke = (char **) malloc (96 * sizeof(char *));
if (glyph_stroke == NULL)
   {
   printf ("could not allocate memory for glyph_stroke\n");
   exit (-1);
   }
for (i = 0; i < 96; ++i)
   {
   glyph_stroke[i] = (char *) malloc (150 * sizeof(char));
   if (glyph_stroke[i] == NULL)
      {
      printf ("Could not allocate memory for glyph_stroke[%d]\n", i);
      }
   }

return 0;
}  /* end of ssp_AllocGlyph() */

/*---------------------------------------------------------------*/

#if PROTO

int ssp_FreeGlyph (void)

#else

int ssp_FreeGlyph ()

#endif

/*
 * Purpose...
 * -------
 * Free memory for the encoded glyph strings.
 *
 */

{  /* begin ssp_FreeGlyph() */
int i;

for (i = 0; i < 96; ++i)
   {
   if (glyph_stroke[i] != NULL) free (glyph_stroke[i]);
   }
if (glyph_stroke != NULL) free (glyph_stroke);

return 0;
}

/*---------------------------------------------------------------*/

#if PROTO

int ssp_InitGlyph (void)

#else

int ssp_InitGlyph ()

#endif

/*
 * Purpose...
 * -------
 * Initialize the stroke data for the characters (or glyphs).
 * Memory must have been allocated already.
 *
 * This data was extracted from the Hershey font files originally
 * created by A. V. Hershey at the U.S. Government National Bureau
 * of Standards and subsequently encoded by James Hurt, Cognition, Inc.
 *
 */

{  /* begin ssp_InitGlyph() */

glyph_code[0]  = 2199; nstroke[0]  =  1; strcpy (glyph_stroke[0], "JZ");
glyph_code[1]  = 2214; nstroke[1]  = 15; strcpy (glyph_stroke[1], "MWRFQHRTSHRF RRHRN RRYQZR[SZRY");
glyph_code[2]  = 2213; nstroke[2]  = 14; strcpy (glyph_stroke[2], "MWRMQNROSNRM RR[QZRYSZS\\R^Q_");
glyph_code[3]  = 2275; nstroke[3]  = 12; strcpy (glyph_stroke[3], "H]SFLb RYFRb RLQZQ RKWYW");
glyph_code[4]  = 2274; nstroke[4]  = 42; strcpy (glyph_stroke[4], "H\\PBP_ RTBT_ RXIWJXKYJYIWGTFPFMGKIKKLMMNOOUQWRYT RKKMMONUPWQXRYTYXWZT[P[MZKXKWLVMWLX");
glyph_code[5]  = 2271; nstroke[5]  = 32; strcpy (glyph_stroke[5], "F^[FI[ RNFPHPJOLMMKMIKIIJGLFNFPGSHVHYG[F RWTUUTWTYV[X[ZZ[X[VYTWT");
glyph_code[6]  = 2272; nstroke[6]  = 49; strcpy (glyph_stroke[6], "F_[NZO[P\\O\\N[MZMYNXPVUTXRZP[M[JZIXIUJSPORMSKSIRGPFNGMIMKNNPQUXWZZ[[[\\Z\\Y RM[KZJXJUKSMQ RMKNMVXXZZ[");
glyph_code[7]  = 2251; nstroke[7]  =  8; strcpy (glyph_stroke[7], "MWRHQGRFSGSIRKQL");
glyph_code[8]  = 2221; nstroke[8]  = 20; strcpy (glyph_stroke[8], "KYVBTDRGPKOPOTPYR]T`Vb RTDRHQKPPPTQYR\\T`");
glyph_code[9]  = 2222; nstroke[9]  = 20; strcpy (glyph_stroke[9], "KYNBPDRGTKUPUTTYR]P`Nb RPDRHSKTPTTSYR\\P`");
glyph_code[10] = 2219; nstroke[10] =  9; strcpy (glyph_stroke[10], "JZRFRR RMIWO RWIMO");
glyph_code[11] = 2232; nstroke[11] =  6; strcpy (glyph_stroke[11], "E_RIR[ RIR[R");
glyph_code[12] = 2211; nstroke[12] =  8; strcpy (glyph_stroke[12], "MWR[QZRYSZS\\R^Q_");
glyph_code[13] = 2231; nstroke[13] =  3; strcpy (glyph_stroke[13], "E_IR[R");
glyph_code[14] = 2210; nstroke[14] =  6; strcpy (glyph_stroke[14], "MWRYQZR[SZRY");
glyph_code[15] = 2220; nstroke[15] =  3; strcpy (glyph_stroke[15], "G][BIb");
glyph_code[16] = 2200; nstroke[16] = 40; strcpy (glyph_stroke[16], "H\\QFNGLJKOKRLWNZQ[S[VZXWYRYOXJVGSFQF RQFOGNHMJLOLRMWNYOZQ[ RS[UZVYWWXRXOWJVHUGSF");
glyph_code[17] = 2201; nstroke[17] = 11; strcpy (glyph_stroke[17], "H\\NJPISFS[ RRGR[ RN[W[");
glyph_code[18] = 2202; nstroke[18] = 45; strcpy (glyph_stroke[18], "H\\LJMKLLKKKJLHMGPFTFWGXHYJYLXNUPPRNSLUKXK[ RTFVGWHXJXLWNTPPR RKYLXNXSZVZXYYX RNXS[W[XZYXYV");
glyph_code[19] = 2203; nstroke[19] = 47; strcpy (glyph_stroke[19], "H\\LJMKLLKKKJLHMGPFTFWGXIXLWNTOQO RTFVGWIWLVNTO RTOVPXRYTYWXYWZT[P[MZLYKWKVLUMVLW RWQXTXWWYVZT[");
glyph_code[20] = 2204; nstroke[20] = 13; strcpy (glyph_stroke[20], "H\\THT[ RUFU[ RUFJUZU RQ[X[");
glyph_code[21] = 2205; nstroke[21] = 39; strcpy (glyph_stroke[21], "H\\MFKP RKPMNPMSMVNXPYSYUXXVZS[P[MZLYKWKVLUMVLW RSMUNWPXSXUWXUZS[ RMFWF RMGRGWF");
glyph_code[22] = 2206; nstroke[22] = 48; strcpy (glyph_stroke[22], "H\\WIVJWKXJXIWGUFRFOGMILKKOKULXNZQ[S[VZXXYUYTXQVOSNRNOOMQLT RRFPGNIMKLOLUMXOZQ[ RS[UZWXXUXTWQUOSN");
glyph_code[23] = 2207; nstroke[23] = 31; strcpy (glyph_stroke[23], "H\\KFKL RKJLHNFPFUIWIXHYF RLHNGPGUI RYFYIXLTQSSRVR[ RXLSQRSQVQ[");
glyph_code[24] = 2208; nstroke[24] = 63; strcpy (glyph_stroke[24], "H\\PFMGLILLMNPOTOWNXLXIWGTFPF RPFNGMIMLNNPO RTOVNWLWIVGTF RPOMPLQKSKWLYMZP[T[WZXYYWYSXQWPTO RPONPMQLSLWMYNZP[ RT[VZWYXWXSWQVPTO");
glyph_code[25] = 2209; nstroke[25] = 48; strcpy (glyph_stroke[25], "H\\XMWPURRSQSNRLPKMKLLINGQFSFVGXIYLYRXVWXUZR[O[MZLXLWMVNWMX RQSORMPLMLLMIOGQF RSFUGWIXLXRWVVXTZR[");
glyph_code[26] = 2212; nstroke[26] = 12; strcpy (glyph_stroke[26], "MWRMQNROSNRM RRYQZR[SZRY");
glyph_code[27] = 2213; nstroke[27] = 14; strcpy (glyph_stroke[27], "MWRMQNROSNRM RR[QZRYSZS\\R^Q_");
glyph_code[28] = 2241; nstroke[28] =  4; strcpy (glyph_stroke[28], "F^ZIJRZ[");
glyph_code[29] = 2238; nstroke[29] =  6; strcpy (glyph_stroke[29], "E_IO[O RIU[U");
glyph_code[30] = 2242; nstroke[30] =  4; strcpy (glyph_stroke[30], "F^JIZRJ[");
glyph_code[31] = 2215; nstroke[31] = 32; strcpy (glyph_stroke[31], "I[MJNKMLLKLJMHNGPFSFVGWHXJXLWNVORQRT RSFUGVHWJWLVNTP RRYQZR[SZRY");
glyph_code[32] = 2273; nstroke[32] = 56; strcpy (glyph_stroke[32], "E`WNVLTKQKOLNMMPMSNUPVSVUUVS RQKOMNPNSOUPV RWKVSVUXVZV\\T]Q]O\\L[JYHWGTFQFNGLHJJILHOHRIUJWLYNZQ[T[WZYYZX RXKWSWUXV");
glyph_code[33] = 2001; nstroke[33] = 18; strcpy (glyph_stroke[33], "H\\RFK[ RRFY[ RRIX[ RMUVU RI[O[ RU[[[");
glyph_code[34] = 2002; nstroke[34] = 45; strcpy (glyph_stroke[34], "G]LFL[ RMFM[ RIFUFXGYHZJZLYNXOUP RUFWGXHYJYLXNWOUP RMPUPXQYRZTZWYYXZU[I[ RUPWQXRYTYWXYWZU[");
glyph_code[35] = 2003; nstroke[35] = 32; strcpy (glyph_stroke[35], "G\\XIYLYFXIVGSFQFNGLIKKJNJSKVLXNZQ[S[VZXXYV RQFOGMILKKNKSLVMXOZQ[");
glyph_code[36] = 2004; nstroke[36] = 30; strcpy (glyph_stroke[36], "G]LFL[ RMFM[ RIFSFVGXIYKZNZSYVXXVZS[I[ RSFUGWIXKYNYSXVWXUZS[");
glyph_code[37] = 2005; nstroke[37] = 22; strcpy (glyph_stroke[37], "G\\LFL[ RMFM[ RSLST RIFYFYLXF RMPSP RI[Y[YUX[");
glyph_code[38] = 2006; nstroke[38] = 20; strcpy (glyph_stroke[38], "G[LFL[ RMFM[ RSLST RIFYFYLXF RMPSP RI[P[");
glyph_code[39] = 2007; nstroke[39] = 40; strcpy (glyph_stroke[39], "G^XIYLYFXIVGSFQFNGLIKKJNJSKVLXNZQ[S[VZXX RQFOGMILKKNKSLVMXOZQ[ RXSX[ RYSY[ RUS\\S");
glyph_code[40] = 2008; nstroke[40] = 27; strcpy (glyph_stroke[40], "F^KFK[ RLFL[ RXFX[ RYFY[ RHFOF RUF\\F RLPXP RH[O[ RU[\\[");
glyph_code[41] = 2009; nstroke[41] = 12; strcpy (glyph_stroke[41], "MXRFR[ RSFS[ ROFVF RO[V[");
glyph_code[42] = 2010; nstroke[42] = 20; strcpy (glyph_stroke[42], "KZUFUWTZR[P[NZMXMVNUOVNW RTFTWSZR[ RQFXF");
glyph_code[43] = 2011; nstroke[43] = 27; strcpy (glyph_stroke[43], "F\\KFK[ RLFL[ RYFLS RQOY[ RPOX[ RHFOF RUF[F RH[O[ RU[[[");
glyph_code[44] = 2012; nstroke[44] = 14; strcpy (glyph_stroke[44], "I[NFN[ ROFO[ RKFRF RK[Z[ZUY[");
glyph_code[45] = 2013; nstroke[45] = 30; strcpy (glyph_stroke[45], "F_KFK[ RLFRX RKFR[ RYFR[ RYFY[ RZFZ[ RHFLF RYF]F RH[N[ RV[][");
glyph_code[46] = 2014; nstroke[46] = 21; strcpy (glyph_stroke[46], "G^LFL[ RMFYY RMHY[ RYFY[ RIFMF RVF\\F RI[O[");
glyph_code[47] = 2015; nstroke[47] = 44; strcpy (glyph_stroke[47], "G]QFNGLIKKJOJRKVLXNZQ[S[VZXXYVZRZOYKXIVGSFQF RQFOGMILKKOKRLVMXOZQ[ RS[UZWXXVYRYOXKWIUGSF");
glyph_code[48] = 2016; nstroke[48] = 29; strcpy (glyph_stroke[48], "G]LFL[ RMFM[ RIFUFXGYHZJZMYOXPUQMQ RUFWGXHYJYMXOWPUQ RI[P[");
glyph_code[49] = 2017; nstroke[49] = 64; strcpy (glyph_stroke[49], "G]QFNGLIKKJOJRKVLXNZQ[S[VZXXYVZRZOYKXIVGSFQF RQFOGMILKKOKRLVMXOZQ[ RS[UZWXXVYRYOXKWIUGSF RNYNXOVQURUTVUXV_W`Y`Z^Z] RUXV\\W^X_Y_Z^");
glyph_code[50] = 2018; nstroke[50] = 45; strcpy (glyph_stroke[50], "G]LFL[ RMFM[ RIFUFXGYHZJZLYNXOUPMP RUFWGXHYJYLXNWOUP RI[P[ RRPTQURXYYZZZ[Y RTQUSWZX[Z[[Y[X");
glyph_code[51] = 2019; nstroke[51] = 34; strcpy (glyph_stroke[51], "H\\XIYFYLXIVGSFPFMGKIKKLMMNOOUQWRYT RKKMMONUPWQXRYTYXWZT[Q[NZLXKUK[LX");
glyph_code[52] = 2020; nstroke[52] = 16; strcpy (glyph_stroke[52], "I\\RFR[ RSFS[ RLFKLKFZFZLYF RO[V[");
glyph_code[53] = 2021; nstroke[53] = 23; strcpy (glyph_stroke[53], "F^KFKULXNZQ[S[VZXXYUYF RLFLUMXOZQ[ RHFOF RVF\\F");
glyph_code[54] = 2022; nstroke[54] = 15; strcpy (glyph_stroke[54], "H\\KFR[ RLFRX RYFR[ RIFOF RUF[F");
glyph_code[55] = 2023; nstroke[55] = 24; strcpy (glyph_stroke[55], "F^JFN[ RKFNV RRFN[ RRFV[ RSFVV RZFV[ RGFNF RWF]F");
glyph_code[56] = 2024; nstroke[56] = 21; strcpy (glyph_stroke[56], "H\\KFX[ RLFY[ RYFK[ RIFOF RUF[F RI[O[ RU[[[");
glyph_code[57] = 2025; nstroke[57] = 20; strcpy (glyph_stroke[57], "H]KFRQR[ RLFSQS[ RZFSQ RIFOF RVF\\F RO[V[");
glyph_code[58] = 2026; nstroke[58] = 16; strcpy (glyph_stroke[58], "H\\XFK[ RYFL[ RLFKLKFYF RK[Y[YUX[");
glyph_code[59] = 2223; nstroke[59] = 12; strcpy (glyph_stroke[59], "KYOBOb RPBPb ROBVB RObVb");
glyph_code[60] =  804; nstroke[60] =  3; strcpy (glyph_stroke[60], "KYKFY^");
glyph_code[61] = 2224; nstroke[61] = 12; strcpy (glyph_stroke[61], "KYTBTb RUBUb RNBUB RNbUb");
glyph_code[62] = 2262; nstroke[62] = 11; strcpy (glyph_stroke[62], "JZPLRITL RMORJWO RRJR[");
glyph_code[63] =  999; nstroke[63] =  3; strcpy (glyph_stroke[63], "JZJ]Z]");
glyph_code[64] = 2252; nstroke[64] =  8; strcpy (glyph_stroke[64], "MWSFRGQIQKRLSKRJ");
glyph_code[65] = 2101; nstroke[65] = 39; strcpy (glyph_stroke[65], "I]NONPMPMONNPMTMVNWOXQXXYZZ[ RWOWXXZZ[[[ RWQVRPSMTLVLXMZP[S[UZWX RPSNTMVMXNZP[");
glyph_code[66] = 2102; nstroke[66] = 33; strcpy (glyph_stroke[66], "G\\LFL[ RMFM[ RMPONQMSMVNXPYSYUXXVZS[Q[OZMX RSMUNWPXSXUWXUZS[ RIFMF");
glyph_code[67] = 2103; nstroke[67] = 28; strcpy (glyph_stroke[67], "H[WPVQWRXQXPVNTMQMNNLPKSKULXNZQ[S[VZXX RQMONMPLSLUMXOZQ[");
glyph_code[68] = 2104; nstroke[68] = 36; strcpy (glyph_stroke[68], "H]WFW[ RXFX[ RWPUNSMQMNNLPKSKULXNZQ[S[UZWX RQMONMPLSLUMXOZQ[ RTFXF RW[[[");
glyph_code[69] = 2105; nstroke[69] = 31; strcpy (glyph_stroke[69], "H[LSXSXQWOVNTMQMNNLPKSKULXNZQ[S[VZXX RWSWPVN RQMONMPLSLUMXOZQ[");
glyph_code[70] = 2106; nstroke[70] = 22; strcpy (glyph_stroke[70], "KXUGTHUIVHVGUFSFQGPIP[ RSFRGQIQ[ RMMUM RM[T[");
glyph_code[71] = 2107; nstroke[71] = 60; strcpy (glyph_stroke[71], "I\\QMONNOMQMSNUOVQWSWUVVUWSWQVOUNSMQM RONNPNTOV RUVVTVPUN RVOWNYMYNWN RNUMVLXLYM[P\\U\\X]Y^ RLYMZP[U[X\\Y^Y_XaUbObLaK_K^L\\O[");
glyph_code[72] = 2108; nstroke[72] = 28; strcpy (glyph_stroke[72], "G]LFL[ RMFM[ RMPONRMTMWNXPX[ RTMVNWPW[ RIFMF RI[P[ RT[[[");
glyph_code[73] = 2109; nstroke[73] = 18; strcpy (glyph_stroke[73], "MXRFQGRHSGRF RRMR[ RSMS[ ROMSM RO[V[");
glyph_code[74] = 2110; nstroke[74] = 25; strcpy (glyph_stroke[74], "MXSFRGSHTGSF RTMT_SaQbObNaN`O_P`Oa RSMS_RaQb RPMTM");
glyph_code[75] = 2111; nstroke[75] = 27; strcpy (glyph_stroke[75], "G\\LFL[ RMFM[ RWMMW RRSX[ RQSW[ RIFMF RTMZM RI[P[ RT[Z[");
glyph_code[76] = 2112; nstroke[76] = 12; strcpy (glyph_stroke[76], "MXRFR[ RSFS[ ROFSF RO[V[");
glyph_code[77] = 2113; nstroke[77] = 44; strcpy (glyph_stroke[77], "BcGMG[ RHMH[ RHPJNMMOMRNSPS[ ROMQNRPR[ RSPUNXMZM]N^P^[ RZM\\N]P][ RDMHM RD[K[ RO[V[ RZ[a[");
glyph_code[78] = 2114; nstroke[78] = 28; strcpy (glyph_stroke[78], "G]LML[ RMMM[ RMPONRMTMWNXPX[ RTMVNWPW[ RIMMM RI[P[ RT[[[");
glyph_code[79] = 2115; nstroke[79] = 36; strcpy (glyph_stroke[79], "H\\QMNNLPKSKULXNZQ[S[VZXXYUYSXPVNSMQM RQMONMPLSLUMXOZQ[ RS[UZWXXUXSWPUNSM");
glyph_code[80] = 2116; nstroke[80] = 36; strcpy (glyph_stroke[80], "G\\LMLb RMMMb RMPONQMSMVNXPYSYUXXVZS[Q[OZMX RSMUNWPXSXUWXUZS[ RIMMM RIbPb");
glyph_code[81] = 2117; nstroke[81] = 33; strcpy (glyph_stroke[81], "H\\WMWb RXMXb RWPUNSMQMNNLPKSKULXNZQ[S[UZWX RQMONMPLSLUMXOZQ[ RTb[b");
glyph_code[82] = 2118; nstroke[82] = 23; strcpy (glyph_stroke[82], "IZNMN[ ROMO[ ROSPPRNTMWMXNXOWPVOWN RKMOM RK[R[");
glyph_code[83] = 2119; nstroke[83] = 32; strcpy (glyph_stroke[83], "J[WOXMXQWOVNTMPMNNMOMQNRPSUUWVXW RMPNQPRUTWUXVXYWZU[Q[OZNYMWM[NY");
glyph_code[84] = 2120; nstroke[84] = 16; strcpy (glyph_stroke[84], "KZPFPWQZS[U[WZXX RQFQWRZS[ RMMUM");
glyph_code[85] = 2121; nstroke[85] = 28; strcpy (glyph_stroke[85], "G]LMLXMZP[R[UZWX RMMMXNZP[ RWMW[ RXMX[ RIMMM RTMXM RW[[[");
glyph_code[86] = 2122; nstroke[86] = 15; strcpy (glyph_stroke[86], "I[LMR[ RMMRY RXMR[ RJMPM RTMZM");
glyph_code[87] = 2123; nstroke[87] = 24; strcpy (glyph_stroke[87], "F^JMN[ RKMNX RRMN[ RRMV[ RSMVX RZMV[ RGMNM RWM]M");
glyph_code[88] = 2124; nstroke[88] = 21; strcpy (glyph_stroke[88], "H\\LMW[ RMMX[ RXML[ RJMPM RTMZM RJ[P[ RT[Z[");
glyph_code[89] = 2125; nstroke[89] = 22; strcpy (glyph_stroke[89], "H[LMR[ RMMRY RXMR[P_NaLbKbJaK`La RJMPM RTMZM");
glyph_code[90] = 2126; nstroke[90] = 16; strcpy (glyph_stroke[90], "I[WML[ RXMM[ RMMLQLMXM RL[X[XWW[");
glyph_code[91] = 2225; nstroke[91] = 40; strcpy (glyph_stroke[91], "KYTBRCQDPFPHQJRKSMSOQQ RRCQEQGRISJTLTNSPORSTTVTXSZR[Q]Q_Ra RQSSUSWRYQZP\\P^Q`RaTb");
glyph_code[92] = 2229; nstroke[92] =  3; strcpy (glyph_stroke[92], "NVRBRb");
glyph_code[93] = 2226; nstroke[93] = 40; strcpy (glyph_stroke[93], "KYPBRCSDTFTHSJRKQMQOSQ RRCSESGRIQJPLPNQPURQTPVPXQZR[S]S_Ra RSSQUQWRYSZT\\T^S`RaPb");
glyph_code[94] = 2246; nstroke[94] = 24; strcpy (glyph_stroke[94], "F^IUISJPLONOPPTSVTXTZS[Q RISJQLPNPPQTTVUXUZT[Q[O");
glyph_code[95] = 2218; nstroke[95] = 14; strcpy (glyph_stroke[95], "KYQFOGNINKOMQNSNUMVKVIUGSFQF");

return 0;
}  /* end of ssp_InitGlyph() */

/*---------------------------------------------------------------*/

#if PROTO

int ssp_DrawGlyph (char ascii_code,
                   double x_LL, double y_LL,
                   double h, double angle)

#else

int ssp_DrawGlyph (ascii_code, x_LL, y_LL, h, angle)
char ascii_code;
double x_LL, y_LL, h, angle;

#endif

/*
 * Purpose...
 * -------
 * Draw a character (glyph) with its lower-left corner at (x_LL, y_LL)
 * and height h
 *
 * Input...
 * -----
 * ascii_code  : ASCII code of the character
 * x_LL, y_LL  : lower left corner of the box containing the glyph
 * h           : height of the glyph in current plotting units
 * angle       : direction of the text (in degrees from the x-axis)
 *
 */

{  /* begin ssp_DrawGlyph() */
char   line[2][256];
int    iglyph, ix, iy, ileft, iright, BaseLine, CapLine;
double x, y, dx, dy;
int    ns, is;
double x_Ctr, y_Ctr, xscale, yscale, cosine, sine;

BaseLine = -9;                        /* for the normal size font */
CapLine = 12;

if (ascii_code > 127) ascii_code = ' ';
if (ascii_code < ' ') ascii_code = ' ';
iglyph = ascii_code - 32;   /* index for the glyph storage arrays */

/* Set the size of the characters. */
yscale = h / (CapLine - BaseLine);
xscale = yscale;

/* Set the direction. */
sine   = sin (3.1415927 * angle / 180.0);
cosine = cos (3.1415927 * angle / 180.0);

/*
 * Unpack the stroke data.
 */
ns = nstroke[iglyph];
for (is = 0; is < ns; ++is)
   {
   line[0][is] = glyph_stroke[iglyph][is * 2];
   line[1][is] = glyph_stroke[iglyph][is * 2 + 1];
   }

/*
 * The left and right limits.
 *
 * Note: this data can be used for proportional spacing
 */
ileft  = (int)line[0][0] - (int)'R';
iright = (int)line[1][0] - (int)'R';

/*
 * The Hershey moves are relative to the "centre" of the character.
 */
dx = -xscale * ileft;
dy = -yscale * BaseLine;
x_Ctr = x_LL + dx * cosine - dy * sine;
y_Ctr = y_LL + dx * sine + dy * cosine;

#if 0
   /*
    * Draw vertical sides to the bounding box.
    */
   dx = xscale * ileft; dy = yscale * BaseLine;
   x = x_Ctr + dx * cosine - dy * sine;
   y = y_Ctr + dx * sine + dy * cosine;
   ssp_Move (x, y);
   dx = xscale * ileft; dy = yscale * CapLine;
   x = x_Ctr + dx * cosine - dy * sine;
   y = y_Ctr + dx * sine + dy * cosine;
   ssp_Plot (x, y);
   dx = xscale * iright; dy = yscale * BaseLine;
   x = x_Ctr + dx * cosine - dy * sine;
   y = y_Ctr + dx * sine + dy * cosine;
   ssp_Move (x, y);
   dx = xscale * iright; dy = yscale * CapLine;
   x = x_Ctr + dx * cosine - dy * sine;
   y = y_Ctr + dx * sine + dy * cosine;
   ssp_Plot (x, y);
#endif


/*
 * The first data point is always a move.
 *
 * Note that Hershey Font data is in TV coordinate system
 * with y increasing down the screen and the origin is in the
 * centre of the character's bounding box.
 *
 */
ix = (int)line[0][1] - (int) 'R';
iy = (int)line[1][1] - (int) 'R';
dx = xscale * ix;
dy = -yscale * iy;
x = x_Ctr + dx * cosine - dy * sine;
y = y_Ctr + dx * sine + dy * cosine;
ssp_Move (x, y);

/*
 * Now, Stroke the glyph.
 */
for (is = 2; is < ns; is++)
   {
   if (line[0][is] == ' ') continue;

   /* .. process vector number ipnt */
   ix = (int)line[0][is] - (int) 'R';
   iy = (int)line[1][is] - (int) 'R';
   dx = xscale * ix;
   dy = -yscale * iy;
   x = x_Ctr + dx * cosine - dy * sine;
   y = y_Ctr + dx * sine + dy * cosine;

   if (line[0][is-1] == ' ') ssp_Move(x, y); else ssp_Plot(x, y);
   } /* for (is = ... */

/*
 * Move to the lower right corner of the box
 */
dx = xscale * (iright - ileft);
dy = 0.0;
x = x_LL + dx * cosine - dy * sine;
y = y_LL + dx * sine + dy * cosine;
ssp_Move (x, y);

return 0;
}  /* end of ssp_DrawGlyph() */

/*---------------------------------------------------------------*/


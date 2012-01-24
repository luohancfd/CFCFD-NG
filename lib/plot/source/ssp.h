/** \file ssp.h
 * \ingroup plot
 * \brief Header for Simple Scientific Plotting routines (for Tk).
 *
 * \author PA Jacobs
 * \version 1.0    : 26-Apr-99 : Adapted from ssp.h
 *
 */

/*
 * The header files "compiler.h", "geom.h"  and "tcl.h" must be included
 * before this file.
 */

#ifndef COMPILER_H
#  include "../../util/sourcecompiler.h"
#endif

#ifndef GEOM_H
#  include "geom.h"
#endif

#ifndef SSP_H
#define SSP_H  1

/*
 * Hard copy devices supported.
 */
#define  HPGL      1
#define  PSCRIPT   2
#define  GIF       3

/*
 * Miscellaneous global definitions...
 */
#define  ssp_PI  3.14159265359

/*
 * Centred symbols available.
 */
#define  CIRCLE_SYM    0
#define  SQUARE_SYM    1
#define  DIAMOND_SYM   2
#define  TRIANGLE_SYM  3
#define  DEL_SYM       4
#define  PLUS_SYM      5
#define  CROSS_SYM     6

/*
 * Clipping definitions.
 */
#define   ALL_VISIBLE   1
#define   NOT_VISIBLE   0
#define   PART_VISIBLE  2

#define   CLIP_TOL      1.0e-10


/*
 * Data structures for points.
 */

struct poly_3D { int A; int B; int C; int D; };
/*
 * The A, B, C, D contain the indices to the vertices of the
 * polygon patch.  Vertices are ordered counter-clockwise.
 * If C == -1, we have a line.
 * If C >= 0 and D == -1, we have a triangle.
 */

struct ssp_data_point { double x; double y; double value; };

struct ssp_data_vector { double *x; double *value;
			 double imin; double imax;
			 double xmin; double xmax;
			 double vmin; double vmax;
			 double xscale; double vscale; };

struct ssp_data_array { double **x; double **y; double **value;
                        double **u; double **v;
			int ixmin; int ixmax;
			int iymin; int iymax;
			double xmin; double xmax;
			double ymin; double ymax;
			double vmin; double vmax;
			double xscale; double yscale; double vscale; };
/*
 * Array indexing is value[ix][iy], ixmin <= ix <= ixmax,
 *                                  iymin <= iy <= iymax.
 */


/*
 * -------------------
 * Function Prototypes...
 * -------------------
 */

#ifdef WITHTCLTK
    int ssp_Init(Tcl_Interp *interp_ptr);
#else
    int ssp_Init( void );
#endif
int ssp_InitInternal( void );
int ssp_SetParameter( char *name, char *value_string );
int ssp_BeginPlot( char *canvas, char *out_file );
int ssp_Compact(double x, char *xstr);
int ssp_EndPlot(void);
int ssp_Final(void);
int ssp_UpdateDisplay( void );

int ssp_SetClipOn (void);
int ssp_SetClipOff (void);
int ssp_SetClipLimits (double x_min, double y_min,
		       double x_max, double y_max);
int ssp_ClipLine (struct point_2D *P1, struct point_2D *P2);
int ssp_InsideClipRegion( double x, double y );
int ssp_CheckVisible (struct point_2D *P1, struct point_2D *P2);

int ssp_Message (char *string);

int ssp_Plot (double x, double y);
int ssp_Move (double x, double y);
int ssp_Polygon (int n, double x[], double y[],
                 int grey_fill, int colour_fill, int border_flag);
int ssp_SetOrigin (double x, double y);
int ssp_Factor (double new_scale);
int ssp_Where (double *x, double *y,
	       double *scale_factor, double *angle);

int ssp_SetDash (double dash, double space);
int ssp_DashLine (double *x, double *y, int n);

int ssp_SetLineThickness (double thickness);
int ssp_SetFillGrey (double lightness);
int ssp_SetHSBColour (double h, double s, double b);
int ssp_SetHSBFillColour (double h, double s, double b);
int ssp_SetNoColour (void);
double ssp_MapToColour (double level, double minimum, double maximum);

int ssp_DrawGlyph (char ascii_code,
                    double x_LL, double y_LL,
                    double h, double angle);
int ssp_AllocGlyph (void);
int ssp_InitGlyph (void);
int ssp_FreeGlyph (void);

int ssp_PlotText (double x, double y,
		  double size, double directn,
		  char *txt, int n);
int ssp_PlotSymbol (double x, double y,
		    double size, int option);
int ssp_PlotNumber (double x, double y,
		    double size, double directn,
		    double value, int dec_places);
int ssp_PlotENumber (double x, double y,
		    double size, double directn,
		    double value, int dec_places);

int ssp_Arrow (double x1, double y1,
	       double x2, double y2,
	       double HeadLen, double HeadWidth,
	       int option);

int ssp_PlotXAxis (double x0, double y0,
		   double length, double angle, double scale,
		   double first_value, double delta_value,
		   char *title, int option);

int ssp_PlotYAxis (double x0, double y0,
		   double length, double angle, double scale,
		   double first_value, double delta_value,
		   char *title, int option);

int ssp_ContourTriangle (struct ssp_data_point *A,
			 struct ssp_data_point *B,
			 struct ssp_data_point *C,
			 double level);

int ssp_ContourTriangle3D (struct point_3D *A,
		           struct point_3D *B,
		           struct point_3D *C,
		           double A_value,
		           double B_value,
		           double C_value,
		           double level);

int ssp_ContourQuad (struct ssp_data_point *A,
		     struct ssp_data_point *B,
		     struct ssp_data_point *C,
		     struct ssp_data_point *D,
		     double level);

int ssp_ContourQuad3D (struct point_3D *A,
		       struct point_3D *B,
		       struct point_3D *C,
		       struct point_3D *D,
		       double A_value,
		       double B_value,
		       double C_value,
		       double D_value,
		       double level);

int ssp_ColourFillTriangle (struct ssp_data_point *A,
		            struct ssp_data_point *B,
  		            struct ssp_data_point *C,
		            double level_min,
		            double level_max,
		            int    grey_fill );

int ssp_ColourFillQuad (struct ssp_data_point *A,
		        struct ssp_data_point *B,
		        struct ssp_data_point *C,
  		        struct ssp_data_point *D,
		        double level_min,
		        double level_max,
		        int    grey_fill );

int ssp_ColourFillArray (struct ssp_data_array *F,
		         double level_min, double level_max,
		         int grey_fill );
		
int ssp_ColourFillTriangle3D (struct point_3D *A,
		              struct point_3D *B,
  		              struct point_3D *C,
		              double A_value,
		              double B_value,
		              double C_value,
		              double level_min,
		              double level_max,
		              int    grey_fill );

int ssp_ColourFillQuad3D (struct point_3D *A,
		          struct point_3D *B,
		          struct point_3D *C,
		          struct point_3D *D,
		          double A_value,
		          double B_value,
		          double C_value,
		          double D_value,
		          double level_min,
		          double level_max,
		          int    grey_fill );

int ssp_AllocateVector (struct ssp_data_vector *V, int idim);
int ssp_FreeVector (struct ssp_data_vector *V);

int ssp_AllocateArray (struct ssp_data_array *F, int ixdim, int iydim);
int ssp_FreeArray(struct ssp_data_array *F);

int ssp_VectorExtremes (struct ssp_data_vector *V);

int ssp_ArrayExtremes (struct ssp_data_array *F);

int ssp_AutoPlot (struct ssp_data_vector *V,
		  double x0, double y0,
		  double x_length, double y_length);

int ssp_ContourArray (struct ssp_data_array *F, double level);

int ssp_PlotMesh (struct ssp_data_array *F, int ixskip, int iyskip);
int ssp_PlotEdge (struct ssp_data_array *F);


int ssp_SetOrigin3D (double x, double y, double z);

int ssp_SetTrimetricProjection (double phi, double theta);

int ssp_SetIsometricProjection (void);

int ssp_SetOrthographicProjection (int iplane);

int ssp_Plot3D (double x, double y, double z);

int ssp_Move3D (double x, double y, double z);

int ssp_ApplyProjectionMatrix (struct point_3D *pt1,
                               struct point_3D *pt2);

int ssp_Polygon3D (int n, struct point_3D vtx[],
                   int grey_fill, int colour_fill, int border_flag);

int ssp_axes3D (double xaxis, double yaxis, double zaxis,
                double xtick, double ytick, double ztick,
                double tick);

#endif

/*-----------------------------------------------------------*/


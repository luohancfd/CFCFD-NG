/* sspd.c  -- exercise Simple Scientific Plotting 
 * #ifdef WITHTCLTK
 *    This is a loadable module for Tcl/Tk
 * #else
 *    This is a stand-alone program with no screen graphics.
 * #endif
 *
 * P. J. 26-Apr-1999, 28-Sep-1999
 */

#include <stdio.h>
#include <math.h>

#include "../../util/source/compiler.h"
#include "geom.h"
#ifdef WITHTCLTK
    #include <tcl.h>
#endif
#include "ssp_colours.h"
#include "ssp.h"


void run_demo( char *canvas_name ) {
    int    j;
    double xs[5], ys[5];
    double theta;
   
    xs[0] = 1000.0; ys[0] = 0.1e-6;   /* sample data */
    xs[1] = 2000.0; ys[1] = 0.4e-6;
    xs[2] = 3000.0; ys[2] = 0.9e-6;
    xs[3] = 4000.0; ys[3] = 1.6e-6;
    xs[4] = 5000.0; ys[4] = 2.5e-6;

    /*
     * Initialize the drawing
     */
    ssp_SetParameter( "orientation", "landscape" );
    ssp_SetParameter( "device", "ps" );
    ssp_BeginPlot( canvas_name, "sspd.ps" );

    ssp_Message("Ssp plotting demo.");

    ssp_PlotText(35.0, 150.0, 8.0, 0.0,
	          "Ssp Plotting Demonstration.", 27);

    ssp_SetOrigin(30.0,20.0);              /* set origin on page */
    ssp_Factor(30.0);                      /* suitable scale */
    
    ssp_SetLineThickness(0.5);
    ssp_PlotXAxis(0.0, 0.0, 6.0, 0.0, 1.0e-3, 0.0, 1.0e3, "X-Axis", -1);
    ssp_PlotYAxis(0.0, 0.0, 4.0, 0.0, 1.0e6, 0.0, 1.0e-6, "Y-Axis", -1);
    ssp_UpdateDisplay();
    
    ssp_Message("Here is the sample data ...");
    for (j = 0; j <= 4; ++j)
        printf("x[%d] = %f, y[%d] = %f\n", j, xs[j], j, ys[j]);
    
    ssp_SetLineThickness(0.18);
    ssp_Move (0.0, 0.0);
    for (j = 0; j < 5; ++j) {
        ssp_Plot(xs[j]/1000.0, ys[j]*1.0e6);
        ssp_PlotSymbol(xs[j]/1000.0, ys[j]*1.0e6, 0.1, CIRCLE_SYM);
    }
    
    /*
     * Try out a dashed line.
     */
    xs[0] = 0.0;  ys[0] = 2.0;
    xs[1] = 2.4;  ys[1] = 3.0;
    xs[2] = 3.0;  ys[2] = 1.0;
    xs[3] = 5.1;  ys[3] = 4.0;
    
    ssp_DashLine (xs, ys, 4);
    
    /*
     * Try out the line clipping.
     */
    ssp_SetOrigin(6.0, 1.0);
    ssp_SetClipLimits(-0.25, -0.25, 0.75, 0.75);
    ssp_SetClipOn();
    for (theta = 0.0; theta <= 3.14; theta += 0.4) {
        ssp_SetHSBColour(theta/3.14, 1.0, 0.5);
        ssp_Move( sin(theta), cos(theta) );
        ssp_Plot( -sin(theta), -cos(theta) );
    }
    ssp_SetClipOff();

    /*
     * Fill in a grey-shaded polygon with a coloured border.
     */
    xs[0] = 0.0; ys[0] = 2.0;
    xs[1] = 1.0; ys[1] = 2.3;
    xs[2] = 0.5; ys[2] = 3.0;
    ssp_SetHSBColour(HSB_GREEN, 1.0, 0.5);
    ssp_SetFillGrey(0.5);
    ssp_Polygon(3, xs, ys, 1, 0, 1);

    /*
     * Fill in a coloured polygon.
     */
    xs[0] = 0.75; ys[0] = 2.0;
    xs[1] = 1.75; ys[1] = 2.3;
    xs[2] = 1.25; ys[2] = 3.0;
    ssp_SetHSBColour(HSB_MAGENTA, 1.0, 0.5);
    ssp_SetHSBFillColour(HSB_YELLOW, 1.0, 0.5);
    ssp_Polygon(3, xs, ys, 0, 1, 1);

    ssp_Message( "End of plotting.");
    ssp_EndPlot();

} /* end function run_demo */


#ifdef WITHTCLTK

    /* 
     * Declarations for the application-specific command procedures.
     */
    int ssp_tk_demo( ClientData clientdata,
                  Tcl_Interp *interp,
                  int argc, char *argv[] );

    /* 
     * The initialization procedure that is called when the 
     * package is loaded.  It seems that we need to have the 
     * first letter capitalised.
     */
    int Ssp_tk_Init( Tcl_Interp *interp ) {
        /* 
         * Register the new command 
         */
        Tcl_CreateCommand( interp, "ssp_tk_demo", ssp_tk_demo, 
            (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
        /*
         * Declare that we implement the gtest package
         * so that we may load the package via 
         * "package require gtest"
         */
        Tcl_PkgProvide( interp, "ssp_tk", "0.1");
        return TCL_OK;
    } /* end of function Gtest_Init */


    int ssp_tk_demo( ClientData clientdata,
                  Tcl_Interp *interp,
                  int argc, char *argv[] ) {
        char result_string[132];
        char canvas_name[132];

        if ( argc != 2 ) {
            Tcl_SetResult(interp, 
               "Usage: ssp_tk_demo canvas_name", TCL_STATIC);
            return TCL_ERROR;
        } /* end if */

        strcpy( canvas_name, argv[1] );
        ssp_Init( interp );
        run_demo( canvas_name );
        ssp_Final();

        sprintf( result_string, "Finished plotting demonstration" );
        Tcl_SetResult(interp, result_string, TCL_STATIC);
        return TCL_OK;

    } /* end function ssp_tk_demo */

#else

    int main () {
        printf( "Run demo without Tcl/Tk graphics.\n");
        ssp_Init();
        run_demo( "" );
        ssp_Final();
        printf( "Finished demo.\n");
        return 0;
    } /* end function main */

#endif


/* ssp3dd.c  -- exercise 3D Simple Scientific Plotting routines */

#include <stdio.h>
#include <math.h>

#include "../../util/source/compiler.h"
#include "../../geometry/source/geom.h"
#include "ssp.h"


main ()
{
struct point_3D A, B, C, D, E, F, G, H, I, J;

/*
 * Set up the vertices of a cube with the corner removed.
 * This is example 3-10 from Rogers and Adams (1990).
 */
A.x = 2.0; A.y = 1.0; A.z = 2.0;
B.x = 3.0; B.y = 1.0; B.z = 2.0;
C.x = 3.0; C.y = 1.5; C.z = 2.0;
D.x = 2.5; D.y = 2.0; D.z = 2.0;
E.x = 2.0; E.y = 2.0; E.z = 2.0;
F.x = 2.0; F.y = 1.0; F.z = 1.0;
G.x = 3.0; G.y = 1.0; G.z = 1.0;
H.x = 3.0; H.y = 2.0; H.z = 1.0;
I.x = 2.0; I.y = 2.0; I.z = 1.0;
J.x = 3.0; J.y = 2.0; J.z = 1.5;

/*
 * Initialize the drawing
 */
ssp_BeginPlot (1, 1, "ssp3dd.ps", "ps", "landscape");

ssp_Message ("Ssp plotting demo.");

ssp_PlotText (35.0, 150.0, 12.0, 0.0,
	      "3D Plotting Demonstration.", 26);

ssp_SetOrigin (100.0,80.0);              /* set origin on page */
ssp_Factor (20.0);                      /* suitable scale */

ssp_SetLineThickness (0.5);

/*
 * Set the 3D origin and type of projection.
 */
ssp_SetOrigin3D (0.0, 0.0, 0.0);
ssp_SetIsometricProjection ();

/*
 * Plot Axes
 */
ssp_Move3D (0.0, 0.0, 0.0);
ssp_Plot3D (2.0, 0.0, 0.0);

ssp_Move3D (0.0, 0.0, 0.0);
ssp_Plot3D (0.0, 2.0, 0.0);

ssp_Move3D (0.0, 0.0, 0.0);
ssp_Plot3D (0.0, 0.0, 2.0);

/*
 * Now, draw the visible faces of the cube.
 */
ssp_Move3D (A.x, A.y, A.z);
ssp_Plot3D (B.x, B.y, B.z);
ssp_Plot3D (C.x, C.y, C.z);
ssp_Plot3D (D.x, D.y, D.z);
ssp_Plot3D (E.x, E.y, E.z);
ssp_Plot3D (A.x, A.y, A.z);

ssp_Plot3D (F.x, F.y, F.z);
ssp_Plot3D (I.x, I.y, I.z);
ssp_Plot3D (E.x, E.y, E.z);

ssp_Move3D (I.x, I.y, I.z);
ssp_Plot3D (H.x, H.y, H.z);
ssp_Plot3D (J.x, J.y, J.z);
ssp_Plot3D (D.x, D.y, D.z);

ssp_Move3D (J.x, J.y, J.z);
ssp_Plot3D (C.x, C.y, C.z);


ssp_EndPlot ();    /* finish up screen demo */


return (0);
}  /* end of ssp3dd */


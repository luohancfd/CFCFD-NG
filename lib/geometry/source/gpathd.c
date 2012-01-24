/** \file gpathd.c
 * \ingroup geom
 * \brief Sample driver for the gpath functions for geometry specification.
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "geom.h"
#include "bezier.h"
#include "gpath.h"

/*--------------------------------------------------------------*/

int main ()
{
    struct point_3D pos, loca[10], locb[30];
    double t;
    int n, jseg;
    double s, r;
    struct GPathPolyLine A, AA, c1, c2, c3, c4, arc, combined;

    printf ("# Test drive gpath functions...\n");
    gpath_init( &A );
    printf ("#--------------------------------------\n");
    printf ("# A Bezier polyline...\n");
    loca[0].x = 0.0; loca[0].y = 0.0; loca[0].z = 0.0;
    loca[1].x = 0.2; loca[1].y = 1.0; loca[1].z = 0.0;
    loca[2].x = 0.8; loca[2].y = 1.0; loca[2].z = 0.0;
    loca[3].x = 1.0; loca[3].y = 0.0; loca[3].z = 0.0;
    gpath_add_element(&A, GPATH_BEZIER, 4, loca);

    loca[0].x = 1.0; loca[0].y =  0.0; loca[0].z = 0.0;
    loca[1].x = 1.2; loca[1].y = -1.0; loca[1].z = 0.0;
    loca[2].x = 1.8; loca[2].y = -1.0; loca[2].z = 0.0;
    loca[3].x = 2.0; loca[3].y =  0.0; loca[3].z = 0.0;
    gpath_add_element(&A, GPATH_BEZIER, 4, loca);

    printf("# Length of polyline = %g\n", A.length );
    printf("#  t  x   y   z   sin(x)\n");
    for (t = 0; t <= 1.0; t += 0.1) {
	gpath_polyline_eval( &A, t, &pos );
	printf ("%f %f %f %f\n", t, pos.x, pos.y, pos.z);
    }
    printf ("#--------------------------------------\n");


    printf ("# A Bezier spline with 4 segments.\n");
    n = 4;
    loca[0].x = 0.0; loca[0].y = sin(0.0); loca[0].z = 0.0;
    loca[1].x = 0.5; loca[1].y = sin(0.5); loca[1].z = 0.0;
    loca[2].x = 1.0; loca[2].y = sin(1.0); loca[2].z = 0.0;
    loca[3].x = 1.5; loca[3].y = sin(1.5); loca[3].z = 0.0;
    loca[4].x = 2.0; loca[4].y = sin(2.0); loca[4].z = 0.0;
    bezier_3_spline( locb, n, loca );
    gpath_init( &AA );
    for (jseg = 0; jseg <= n; ++jseg ) {
	gpath_add_element( &AA, GPATH_BEZIER, 4, &(locb[jseg*4]) );
    }
    printf("#  t  x   y   z   sin(x)\n");
    for (t = 0; t <= 1.001; t += 0.05) {
	gpath_polyline_eval( &AA, t, &pos );
	printf ("%f %f %f %f %f\n", t, pos.x, pos.y, pos.z, sin(pos.x));
    } /* end for */
    printf ("#--------------------------------------\n\n");


    printf ("# Bezier Coons patch...\n");
    gpath_init( &c1 );
    loca[0].x = 0.0; loca[0].y = 0.0; loca[0].z = 0.0;
    loca[1].x = 0.2; loca[1].y = 0.2; loca[1].z = 0.0;
    loca[2].x = 0.8; loca[2].y = 0.2; loca[2].z = 0.0;
    loca[3].x = 1.0; loca[3].y = 0.0; loca[3].z = 0.0;
    gpath_add_element(&c1, GPATH_BEZIER, 4, loca);
    loca[0].x = 1.0; loca[0].y =  0.0; loca[0].z = 0.0;
    loca[1].x = 1.2; loca[1].y = -0.2; loca[1].z = 0.0;
    loca[2].x = 1.8; loca[2].y = -0.2; loca[2].z = 0.0;
    loca[3].x = 2.0; loca[3].y =  0.0; loca[3].z = 0.0;
    gpath_add_element(&c1, GPATH_BEZIER, 4, loca);

    gpath_init( &c2 );
    loca[0].x = 0.0; loca[0].y = 2.0; loca[0].z = 0.0;
    loca[1].x = 0.2; loca[1].y = 2.2; loca[1].z = 0.0;
    loca[2].x = 0.8; loca[2].y = 2.2; loca[2].z = 0.0;
    loca[3].x = 1.0; loca[3].y = 2.0; loca[3].z = 0.0;
    gpath_add_element(&c2, GPATH_BEZIER, 4, loca);
    loca[0].x = 1.0; loca[0].y = 2.0; loca[0].z = 0.0;
    loca[1].x = 1.2; loca[1].y = 1.8; loca[1].z = 0.0;
    loca[2].x = 1.8; loca[2].y = 1.8; loca[2].z = 0.0;
    loca[3].x = 2.0; loca[3].y = 2.0; loca[3].z = 0.0;
    gpath_add_element(&c2, GPATH_BEZIER, 4, loca);

    gpath_init( &c3 );
    loca[0].x = 0.0; loca[0].y = 0.0; loca[0].z = 0.0;
    loca[1].x = 0.2; loca[1].y = 0.2; loca[1].z = 0.0;
    loca[2].x = 0.2; loca[2].y = 0.8; loca[2].z = 0.0;
    loca[3].x = 0.0; loca[3].y = 1.0; loca[3].z = 0.0;
    gpath_add_element(&c3, GPATH_BEZIER, 4, loca);
    loca[0].x =  0.0; loca[0].y = 1.0; loca[0].z = 0.0;
    loca[1].x = -0.2; loca[1].y = 1.2; loca[1].z = 0.0;
    loca[2].x = -0.2; loca[2].y = 1.8; loca[2].z = 0.0;
    loca[3].x =  0.0; loca[3].y = 2.0; loca[3].z = 0.0;
    gpath_add_element(&c3, GPATH_BEZIER, 4, loca);

    gpath_init( &c4 );
    loca[0].x = 2.0; loca[0].y = 0.0; loca[0].z = 0.0;
    loca[1].x = 2.2; loca[1].y = 0.2; loca[1].z = 0.0;
    loca[2].x = 2.2; loca[2].y = 0.8; loca[2].z = 0.0;
    loca[3].x = 2.0; loca[3].y = 1.0; loca[3].z = 0.0;
    gpath_add_element(&c4, GPATH_BEZIER, 4, loca);
    loca[0].x = 2.0; loca[0].y = 1.0; loca[0].z = 0.0;
    loca[1].x = 1.8; loca[1].y = 1.2; loca[1].z = 0.0;
    loca[2].x = 1.8; loca[2].y = 1.8; loca[2].z = 0.0;
    loca[3].x = 2.0; loca[3].y = 2.0; loca[3].z = 0.0;
    gpath_add_element(&c4, GPATH_BEZIER, 4, loca);

    printf("# Coons-patch surface, ready for GNU-Plot\n");
    printf("# s  r  x  y  z\n");
    for (s = 0.0; s <= 1.0; s += 0.333333) {
	for (r = 0.0; r <= 1.0; r += 0.333333) {
	    coons_patch(&c1, &c2, &c3, &c4, r, s, &pos);
	    printf ("%12.6f %12.6f %12.6f %12.6f %12.6f\n",
		    s, r, pos.x, pos.y, pos.z);
	}
	printf ("\n");  /* blank line to separate rows */
    }

    printf ("# A circular arc...\n");
    loca[0].x = 0.0; loca[0].y = 1.0; loca[0].z = 0.0;
    loca[1].x = 0.0; loca[1].y = 0.0; loca[1].z = 1.0;
    loca[2].x = 0.0; loca[2].y = 0.0; loca[2].z = 0.0;
    printf("# Length of arc = %g\n", arc_length(loca) );
    printf("#  t  x   y   z\n");
    for (t = 0; t <= 1.0; t += 0.1) {
	arc_eval( loca, t, &pos );
	printf ("%f %f %f %f\n", t, pos.x, pos.y, pos.z);
    }
    printf ("#--------------------------------------\n");

    printf ("# A circular arc plus line...\n");
    gpath_init( &arc );
    loca[0].x = 0.0; loca[0].y = 1.0; loca[0].z = 0.0;
    loca[1].x = 0.0; loca[1].y = 0.0; loca[1].z = 1.0;
    loca[2].x = 0.0; loca[2].y = 0.0; loca[2].z = 0.0;
    gpath_add_element(&arc, GPATH_ARC, 3, loca);
    loca[0].x = 0.0; loca[0].y = 0.0; loca[0].z = 1.0;
    loca[1].x = 0.0; loca[1].y = -1.5708; loca[1].z = 1.0;
    /* The pi/2 line length was selected to match the arc length. 
     * That way we keep the same arc sampling as the single-arc test. */
    gpath_add_element(&arc, GPATH_LINE, 2, loca);
    printf("#  t  x   y   z\n");
    for (t = 0; t <= 1.0001; t += 0.05) {
	gpath_polyline_eval( &arc, t, &pos );
	printf ("%f %f %f %f\n", t, pos.x, pos.y, pos.z);
    } /* end for */
    printf ("#--------------------------------------\n");

    printf ("# A combined path...\n");
    gpath_init( &combined );
    gpath_append_polyline( &combined, &c1, 1 );
    gpath_append_polyline( &combined, &c4, 1 );
    gpath_append_polyline( &combined, &c2, -1 );
    gpath_append_polyline( &combined, &c3, -1 );
    printf("#  t  x   y   z\n");
    for (t = 0; t <= 1.0001; t += 0.05) {
	gpath_polyline_eval( &combined, t, &pos );
	printf ("%f %f %f %f\n", t, pos.x, pos.y, pos.z);
    } /* end for */
    printf ("#--------------------------------------\n");

    printf ("# Done.\n");
    return 0;
}  /* end main() */


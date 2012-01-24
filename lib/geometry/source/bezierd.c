/** \file bezierd.c
 * \ingroup geom
 * \brief Sample driver for the Bezier curves of third order
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
    struct point_3D B[4], pos, loca[10];
    double t;

    printf ("Test drive Bezier curves...\n");

    printf ("A single third-degree Bezier curve...\n");
    B[0].x = 0.0; B[0].y = 0.0; B[0].z = 0.0;
    B[1].x = 0.2; B[1].y = 1.0; B[1].z = 0.0;
    B[2].x = 0.8; B[2].y = 1.0; B[2].z = 0.0;
    B[3].x = 1.0; B[3].y = 0.0; B[3].z = 0.0;
    for (t = 0; t <= 1.001; t += 0.1) {
	bezier_3_eval(B, t, &pos);
	printf ("%f %f %f\n", t, pos.x, pos.y);
    }
    printf ("length of Bezier curve = %g\n", bezier_length(3, B));

    printf ("A single n-degree Bezier curve...\n");
    loca[0].x = 0.0; loca[0].y = 0.0; loca[0].z = 0.0;
    loca[1].x = 0.2; loca[1].y = 1.0; loca[1].z = 0.0;
    loca[2].x = 0.8; loca[2].y = 1.0; loca[2].z = 0.0;
    loca[3].x = 1.0; loca[3].y = 0.0; loca[3].z = 0.0;
    loca[4].x = 1.2; loca[4].y = 1.0; loca[4].z = 0.0;
    loca[5].x = 1.8; loca[5].y = 1.0; loca[5].z = 0.0;
    loca[6].x = 2.0; loca[6].y = 0.0; loca[6].z = 0.0;

    for (t = 0; t <= 1.001; t += 0.1) {
	bezier_eval(6, loca, t, &pos);
	printf ("%f %f %f %f\n", t, pos.x, pos.y, pos.z);
    }

    printf ("length of Bezier curve = %g\n", bezier_length(6, loca));
    return 0;
}  /* end main() */

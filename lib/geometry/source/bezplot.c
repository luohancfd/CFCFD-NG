/** \file bezplot.c
 * \ingroup nm
 *
 * \brief  Read a Bezier curve description and generate a file of
 * control point coordinates and a file of curve coordinates.
 *
 * \author PA Jacobs
 *
 * \date 21-Feb-95
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "../../util/source/compiler.h"
#include "geom.h"
#include "bezier.h"
#include "gpath.h"

int main()
{
    int     n, i;
    double  t;
    struct  point_3D loca[MAX_BEZ_DEGREE], pnt;
    FILE    *bez_file, *curve_file, *control_file;
#   define NCHAR 132
    char    textline[NCHAR];

    bez_file = fopen("bezier.dat", "r");
    if (bez_file == NULL) {
	printf ("Could not open Bezier data file.\n");
	exit(-1);
    }

    /*
     * Read the Bezier description file.
     */
    fgets (textline, NCHAR, bez_file);
    sscanf (textline, "%d", &n);
    if (n > MAX_BEZ_DEGREE) n = MAX_BEZ_DEGREE;

    for (i = 0; i <= n; ++i) {
	fgets (textline, NCHAR, bez_file);
	sscanf (textline, "%lf %lf %lf", &(loca[i].x),
		&(loca[i].y), &(loca[i].z) );
    }

    fclose (bez_file);

    /*
     * Now compute an array of points along the curve and write
     * them to a file for plotting.
     */
    curve_file = fopen("curve.dat", "w");
    if (curve_file == NULL) {
	printf ("Could not write curve file.\n");
	exit(-1);
    }

    for (t = 0.0; t <= 1.0001; t += 0.005) {
	bezier_eval(n, loca, t, &pnt);
	fprintf(curve_file, "%f %f\n", pnt.x, pnt.y);
    }

    fclose (curve_file);

    /*
     * Write out the control-point coordinates so that they may be
     * plotted separately.
     */
    control_file = fopen("control.dat", "w");
    if (control_file == NULL) {
	printf ("Could not write control point file.\n");
	exit(-1);
    }

    for (i = 0; i <= n; ++i)
	fprintf (control_file, "%f %f\n", loca[i].x, loca[i].y);

    fclose (control_file);

    return 0;
}  /* end of main() */


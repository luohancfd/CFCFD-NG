// l_tube.cxx
// Refactored from l1d code 25-Sep-2012

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "l_tube.hh"

TubeModel::TubeModel(std::string config_file_name, int echo_input)
{
    ConfigParser dict = ConfigParser(config_file_name);
    /*
     * Tube parameters: 
     * n : number of steps for internal discretization (say 10000)
     * nseg : number of tube segments in the user description
     * xb[0],    Diamb[0],    linear[0]
     * ...
     * xb[nseg], Diamb[nseg], linear[nseg]
     */
    dict.parse_int("global_data", "tube_n", n, 10000);
    dict.parse_int("global_data", "tube_nseg", nseg, 1);
    std::vector<double> vdbl_default;
    vdbl_default.resize(nseg+1);
    for ( size_t i = 0; i < vdbl_default.size(); ++i ) vdbl_default[i] = 0.0;
    dict.parse_vector_of_doubles("global_data", "tube_xb", xb, vdbl_default);
    dict.parse_vector_of_doubles("global_data", "tube_d", Diamb, vdbl_default);
    std::vector<int> vint_default;
    vint_default.resize(nseg+1);
    for ( size_t i = 0; i < vint_default.size(); ++i ) vint_default[i] = 0;
    dict.parse_vector_of_ints("global_data", "tube_linear", linear, vint_default);
    if (echo_input == 1) {
	cout << "    tube_n = " << n << endl;
	cout << "    tube_nseg = " << nseg << endl;
	cout << "    tube_xb =";
	for ( int i = 0; i < nseg+1; ++i ) cout << " " << xb[i];
	cout << endl;
	cout << "    tube_d =";
	for ( int i = 0; i < nseg+1; ++i ) cout << " " << Diamb[i];
	cout << endl;
	cout << "    tube_linear =";
	for ( int i = 0; i < nseg+1; ++i ) cout << " " << linear[i];
	cout << endl;
    }
    // Now, read the loss-factor patches.
    dict.parse_int("global_data", "KL_n", nKL, 0);
    vdbl_default.resize(nKL);
    for ( size_t i = 0; i < vdbl_default.size(); ++i ) vdbl_default[i] = 0.0;
    dict.parse_vector_of_doubles("global_data", "KL_xL", xbeginK, vdbl_default);
    dict.parse_vector_of_doubles("global_data", "KL_xR", xendK, vdbl_default);
    dict.parse_vector_of_doubles("global_data", "KL_K", K, vdbl_default);
    if (echo_input == 1) {
	cout << "    KL_nseg = " << nKL << endl;
	cout << "    KL_xL =";
	for ( int i = 0; i < nKL; ++i ) cout << " " << xbeginK[i];
	cout << endl;
	cout << "    KL_xR =";
	for ( int i = 0; i < nKL; ++i ) cout << " " << xendK[i];
	cout << endl;
	cout << "    KL_K =";
	for ( int i = 0; i < nKL; ++i ) cout << " " << K[i];
	cout << endl;
    }
    // Now, read the wall-temperature patches.
    dict.parse_double("global_data", "T_nominal", Tnominal, 300.0);
    dict.parse_int("global_data", "Tpatch_n", nT, 0);
    vdbl_default.resize(nT);
    for ( size_t i = 0; i < vdbl_default.size(); ++i ) vdbl_default[i] = 0.0;
    dict.parse_vector_of_doubles("global_data", "Tpatch_xL", xbeginT, vdbl_default);
    dict.parse_vector_of_doubles("global_data", "Tpatch_xR", xendT, vdbl_default);
    for ( size_t i = 0; i < vdbl_default.size(); ++i ) vdbl_default[i] = 300.0;
    dict.parse_vector_of_doubles("global_data", "Tpatch_T", Tlocal, vdbl_default);
    if (echo_input == 1) {
	cout << "    T_nominal = " << Tnominal << endl;
	cout << "    Tpatch_n = " << nT << endl;
	cout << "    Tpatch_xL =";
	for ( int i = 0; i < nT; ++i ) cout << " " << xbeginT[i];
	cout << endl;
	cout << "    Tpatch_xR =";
	for ( int i = 0; i < nT; ++i ) cout << " " << xendT[i];
	cout << endl;
	cout << "    Tpatch_T =";
	for ( int i = 0; i < nT; ++i ) cout << " " << Tlocal[i];
	cout << endl;
    }
    printf("Set tube area specification\n");
    double myPI = 4.0*atan(1.0);
    diam.resize(n);
    area.resize(n);
    T_Wall.resize(n);
    K_over_L.resize(n);
    x1 = xb[0];
    dx = (xb[nseg] - xb[0]) / n;
    for ( int ix = 0; ix < n; ++ix ) {
        double real_x = x1 + dx * ix;
	// Locate the appropriate tube segment and interpolate area.
	// Start the search from the left and stop when the xb[iseg]
	// exceeds the given x-location.
	int iseg = 0;
        int found_segment = 0;
        for ( iseg = 0; iseg <= nseg; ++iseg ) {
            if (real_x < xb[iseg]) {
                found_segment = 1;
                /* on leaving this loop, iseg indicates the segment */
                break;
            } 
        } // end for iseg...

        if ( iseg == 0 ) {
            /* We are upstream of the tube. */
            diam[ix] = Diamb[iseg];
        } else if (found_segment == 0) {
            /* We are downstream of the tube. */
            diam[ix] = Diamb[nseg];
        } else {
            /* We are between xb[iseg-1] and xb[iseg]. */
            if (linear[iseg - 1] == 1) {
                /* Linear interpolation */
                double xx = (real_x - xb[iseg-1]) / (xb[iseg] - xb[iseg-1]);
                diam[ix] = Diamb[iseg-1] * (1.0 - xx) + Diamb[iseg] * xx;
            } else {
                /* Cubic interpolation to give dArea/dx == 0 at ends */
                double xx = (real_x - xb[iseg - 1]) / (xb[iseg] - xb[iseg-1]);
                double alpha2 = 3.0 * (Diamb[iseg] - Diamb[iseg-1]);
                double alpha3 = -2.0 / 3.0 * alpha2;
                diam[ix] = Diamb[iseg - 1] + (alpha2 + alpha3 * xx) * xx * xx;
            }   /* end if */
        }   /* end if */
        area[ix] = myPI * 0.25 * diam[ix] * diam[ix];

        /*
         * Pipe-fitting loss coefficients:
         * Most of the tube is "smooth", so assume zero, then
         * search the loss patches to see if we are within one.
         */
        double K_over_L_value = 0.0;
        for ( int iseg = 0; iseg < nKL; ++iseg ) {
            if (real_x >= xbeginK[iseg] && real_x <= xendK[iseg]) {
                K_over_L_value = K[iseg] / (xendK[iseg] - xbeginK[iseg]);
            }
        } /* end for iseg */
        K_over_L[ix] = K_over_L_value;

        /*
         * Local variations of wall temperature:
         * Assume a nominal temperature then
         * search the loss patches to see if we are within one.
         */
        double T_Wall_value = Tnominal;
        for ( int iseg = 0; iseg < nT; ++iseg ) {
            if (real_x >= xbeginT[iseg] && real_x <= xendT[iseg]) {
                T_Wall_value = Tlocal[iseg];
            }
        } /* end for iseg */
        T_Wall[ix] = T_Wall_value;
    } /* end for ix... */
}


TubeModel::~TubeModel()
{
    // FIX-ME should clean up the vectors,
    // however, it probably doesn't matter because we only ever make one of these.
}


int TubeModel::read_area(std::string file_name)
// Read the Area(x) specification from a file.
{
    FILE* gf = fopen(file_name.c_str(), "r");
    if ( !gf ) {
        printf("\nCould not open %s; BAILING OUT\n", file_name.c_str());
        return FAILURE;
    } 
#   define NCHAR 320
    char line[NCHAR];
    if (fgets(line, NCHAR, gf) == NULL) {
        printf("Empty area specification file.\n");
        return FAILURE;
    }
    sscanf(line, "%d", &(n));
    if (fgets(line, NCHAR, gf) == NULL) {
        printf("Empty area specification file.\n");
        return FAILURE;
    }
    sscanf(line, "%lf %lf", &(x1), &(dx));
    diam.resize(n);
    area.resize(n);
    T_Wall.resize(n);
    K_over_L.resize(n);
    for ( int ix = 0; ix < n; ++ix ) {
        if (fgets(line, NCHAR, gf) == NULL) {
            printf("Premature end of area file, tube segment[%d]\n", ix);
            return FAILURE;
        }
        sscanf(line, "%lf %lf %lf %lf", &(diam[ix]), &(area[ix]), &(T_Wall[ix]), &(K_over_L[ix]));
    } // end for ix...
#   undef NCHAR
    if ( gf ) fclose(gf);
    return SUCCESS;
} // end read_area()


int TubeModel::write_area(std::string file_name)
// Write the area(x) specification out to a file.
{
    printf("Write area file.\n");
    FILE* gf = fopen(file_name.c_str(), "w");
    if ( !gf ) {
        printf("\nCould not open %s\n", file_name.c_str());
        return FAILURE;
    } 
    fprintf(gf, "%d\n", n);
    fprintf(gf, "%e %e\n", x1, dx);
    for ( int ix = 0; ix < n; ++ix ) {
        fprintf(gf, "%e %e %e %e\n", diam[ix], area[ix], T_Wall[ix], K_over_L[ix]);
    } 
    if ( gf ) fclose(gf);
    return SUCCESS;
} // end write_area()


int TubeModel::write_dump_file(std::string file_name)
{
    printf("Writing dump file for diameter, area, L_over_L, T_wall...\n");
    FILE* dumpf = fopen(file_name.c_str(), "w");
    if ( !dumpf ) {
        printf("\nCould not open %s\n", file_name.c_str());
        return FAILURE;
    }
    fprintf(dumpf, "# x  diam  area  K_over_L  T_Wall\n");
    for ( int ix = 0; ix < n; ++ix ) {
        double xx = dx * ix + x1;
        fprintf(dumpf, "%e %e %e %e %e\n", xx, diam[ix], area[ix], K_over_L[ix], T_Wall[ix]);
    }
    if ( dumpf ) fclose(dumpf);
    return SUCCESS;
} // end write_dump_file()

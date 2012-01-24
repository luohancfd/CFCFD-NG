/** \file roberts.h
 * ingroup nm
 * \brief Header file for the coordinate-stretching functions.
 *
 *-------------------------------------------------------------------
 */

#ifndef COMPILER_H
#  include "../../util/source/compiler.h"
#endif


/* Full function prototypes... */

int distribute_points (double x1, double x2, int n, double x[],
		       double beta_end1, double beta_end2);

int distribute_points_1 (double x1, double x2, int n, double x[],
			 int end1, int end2, double beta);

int distribute_points_2 (double x1, double x2, int n, double x[],
			 double beta_end1, double beta_end2);

int distribute_points_3 (double x1, double x2, int n, double x[],
                         double beta, double yc);

double roberts (double eta, double alpha, double beta);

double roberts_rev (double eta, double beta);

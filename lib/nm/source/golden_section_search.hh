// Author: Rowan J. Gollan
// Version: 06-Jun-2008
//            Initial version.

#ifndef GOLDEN_SECTION_SEARCH_HH
#define GOLDEN_SECTION_SEARCH_HH

double golden_section_search(double (*f)(double),
			     double a, double c,
			     int &result_flag,
			     double tolerance=1.0e-11,
			     int max_iterations=100);

#endif

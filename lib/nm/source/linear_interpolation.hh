#ifndef LINEAR_INTERPOLATION_HH
#define LINEAR_INTERPOLATION_HH

#include <vector>

int 
linear_eval(double xval,
	    double &yval,
	    const std::vector<double> x,
	    const std::vector<double> y);

// python friendly
double 
linear_eval_py(double xval,
	       std::vector<double> x,
	       std::vector<double> y);
#endif

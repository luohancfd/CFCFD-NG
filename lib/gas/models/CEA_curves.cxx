// Author: Rowan J. Gollan
// Date: 30-July-2008

#include <cmath>

#include "CEA_curves.hh"

bool check_T_transport_range(double T, CEA_transport_params &p)
{
    if ( T >= p.T_low && T <= p.T_high )
	return true;
    else
	return false;
}
double eval_CEA_transport_curve(double T, CEA_transport_params &p)
{
    double log_val = p.A*log(T) + p.B/T + p.C/(T*T) + p.D;
    double val = exp(log_val);
    return val;
}

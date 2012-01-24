// Author: Rowan J. Gollan
// Date: 30-July-2007

#ifndef CEA_CURVES_HH
#define CEA_CURVES_HH

struct CEA_transport_params {
    double T_low;
    double T_high;
    double A;
    double B;
    double C;
    double D;
};

bool check_T_transport_range(double T, CEA_transport_params &p);
double eval_CEA_transport_curve(double T, CEA_transport_params &p);

#endif

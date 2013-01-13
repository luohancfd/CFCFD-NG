/** \file line_shapes.cxx
 *  \ingroup radiation
 *
 *  \brief Functions for describing line profiles
 *
 *  \author Daniel F. Potter
 *  \version 19-Oct-12: initial implementation
 *            
 **/

#include <cstdlib>
#include <math.h>

using namespace std;

double eval_Gaussian_profile( double delta_nu, double gamma )
{
    // delta_nu: distance from center of profile
    //    gamma: half-width at half-maximum
    double A = (1.0/gamma)*sqrt(log(2.0)/M_PI);
    double B = -log(2.0)*(delta_nu/gamma)*(delta_nu/gamma);
    
    return A*exp(B);
}

double eval_Lorentzian_profile( double delta_nu, double gamma )
{
    /* Line shape as a function of delta_nu for a Lorentz profile */
    double bp_nu = gamma / M_PI / ( ( delta_nu * delta_nu ) + gamma * gamma );

    return bp_nu;
}

double calculate_Voigt_width( double gamma_L, double gamma_G )
{
    double d = ( gamma_L - gamma_G ) / ( gamma_L + gamma_G );
    double beta = 0.023665 * exp ( 0.6 * d) + 0.0418 * exp ( -1.9 * d);
    double alpha = 0.18121;
    double R_d = 1.0 - alpha * ( 1.0 - d * d ) - beta * sin(M_PI * d);

    return R_d * ( gamma_L + gamma_G);
}

double eval_Voigt_profile( double delta_nu, double gamma_V, double gamma_L, double gamma_G )
{
    // Ref: Whiting (1968) JQRST Vol. 8 pp 1379-1384
    double R_l = delta_nu / ( 2.0 * gamma_V );
    double R_d = gamma_L / gamma_V;

    // Approximate expression
    double tmpA = ( 1.0 - R_d ) * exp( -2.772 * R_l * R_l ) + R_d / ( 1.0 + 4.0 * R_l * R_l );
    double tmpC = 2.0 * gamma_V * (1.065 + 0.447 * R_d + 0.058 * R_d * R_d);

    double b_nu = tmpA / tmpC;

    return b_nu;
}


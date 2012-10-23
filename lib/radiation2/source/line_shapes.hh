/** \file line_shapes.hh
 *  \ingroup radiation2
 *
 *  \brief Functions for describing line profiles
 *
 *  \author Daniel F. Potter
 *  \version 19-Oct-12: initial implementation
 *
 **/

#ifndef LINE_SHAPES_HH
#define LINE_SHAPES_HH

double eval_Gaussian_profile( double delta_nu, double gamma );

double eval_Lorentzian_profile( double delta_nu, double gamma );

double calculate_Voigt_width( double gamma_L, double gamma_G );

double eval_Voigt_profile( double delta_nu, double gamma_V, double gamma_L, double gamma_G );

#endif

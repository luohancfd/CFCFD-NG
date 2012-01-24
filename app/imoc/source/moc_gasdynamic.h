/** \file moc_gasdynamic.h
 * \ingroup imoc
 * \brief Header file for basic gas dynamic relations. 
 */

double T0_T(double M, double g);
double P0_P(double M, double g);

double NuFromM(double M, double g);
double MFromNu( double Nu, double g );
double MFromNu_approximate( double Nu, double g );
double MachAngle( double M );


/** \file radiation_constants.hh
 *  \ingroup radiation2
 *
 *  \author Daniel F Potter
 *  \version 02-July-09
 *
 **/

#ifndef RAD_CONST_HH
#define RAD_CONST_HH

const double RC_sigma_SI   = 5.670400e-8;   /* Stefan-Boltzmann constant W/(m^2 K^4)     */
const double RC_c          = 2.99792458e10; /* cm/s, speed of light                      */
const double RC_c_SI       = 2.99792458e8;  /* m/s, speed of light (in S.I.)             */
const double RC_h          = 6.626076e-27;  /* erg.s, Planck's constant                  */
const double RC_h_SI       = 6.626076e-34;  /* J.s, Planck's constant (in S.I.)          */
const double RC_k          = 1.380658e-16;  /* erg/deg, Boltzmann's constant             */
const double RC_k_SI       = 1.380658e-23;  /* J/K , Boltzmann's constant (in S.I.)      */
const double RC_e          = 4.803e-10;     /* esu, charge on electron                   */
const double RC_e_SI       = 1.60218e-19;   /* Coulombs, charge on electron (in S.I.)    */
const double RC_eps0_SI    = 8.85419e-12;   /* F/m, Permittivity of a vacuum (in S.I.)   */
const double RC_m          = 9.109390e-28;  /* g, mass of an electron                    */
const double RC_m_SI       = 9.109390e-31;  /* kg, mass of an electron (in S.I.)         */
const double RC_Na         = 6.022137e23;   /* particles/mol, Avogadro's number          */
const double RC_a0         = 0.52918e-08;   /* Bohr radius (cm)                          */
const double RC_a0_SI      = 0.52918e-10;   /* Bohr radius (m)                           */
const double RC_alpha      = 7.292e-3;      /* Fine structure constant (ND)              */
const double RC_H_ionise   = 1.0967877e5;   /* Hydrogen atom ionisation potential (cm-1) */
const double RC_H_ionise_J = 2.1787113e-18; /* Hydrogen atom ionisation potential (J)    */
const double RC_R_u        = 8.314472;      /* J/(mol.K)                                 */
const double RC_E0         = 1.3875117e-11; /* erg, energy potential of first Bohr orbit */

#endif


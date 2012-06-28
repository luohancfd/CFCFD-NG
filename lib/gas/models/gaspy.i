/** \file gaspy.i
 *  \ingroup libgas 
 *  \brief SWIG interface header file for C++ gas model library.
 *  \author RJG
 *  
 *  NOTE: This differs from libgas.i in that this is for a 
 *        standalone gas Python module.  It does not hook into
 *        the larger framework of the CFD code.   
 */
%define DOCSTRING
"Python interface to the gas model classes, now implemented in C++"
%enddef
%module(docstring=DOCSTRING) gaspy

%include "std_string.i"
%include "std_vector.i"

// The following magic allows us to pass a list of numbers
// to be collected as a C++ vector.
#ifndef STD_VECTOR_TEMPLATES_ALREADY_DEFINED
%template(vectori) std::vector<int>;
%template(vectord) std::vector<double>;
#define STD_VECTOR_TEMPLATES_ALREADY_DEFINED
#endif

%include "libgas.i"

%{
#include "../../util/source/useful.h"
#include "composite-gas-model.hh"
#include "look-up-table.hh"
#include "LUT-plus-composite-gas-model.hh"
#include "../kinetics/reaction-update.hh"
#include "../kinetics/chemical-kinetic-ODE-update.hh"
#include "../kinetics/chemical-kinetic-system.hh"
#include "../kinetics/reaction-rate-coeff.hh"
#include "../kinetics/generalised-Arrhenius.hh"
#include "../kinetics/pressure-dependent-rate.hh"
#include "../kinetics/reaction.hh"
#include "thermal-behaviour-model.hh"
#include "noneq-thermal-behaviour.hh"
#include "chemical-species.hh"
#include "chemical-equilibrium-system.hh"
#include "chemical-species-library.hh"
#include "species-energy-modes.hh"
#include "../kinetics/energy-exchange-mechanism.hh"
#include "CI-functor.hh"
#include "GuptaYos-mixing-rule.hh"
%}

%include "../../nm/source/ode_system.hh"
%include "../../nm/source/zero_system.hh"
%include "../kinetics/reaction-update.hh"
%include "../kinetics/chemical-kinetic-ODE-update.hh"
%include "../kinetics/chemical-kinetic-system.hh"
%include "../kinetics/reaction-rate-coeff.hh"
%include "../kinetics/pressure-dependent-rate.hh"
%include "../kinetics/reaction.hh"
%include "thermal-behaviour-model.hh"
%include "noneq-thermal-behaviour.hh"
%include "chemical-species.hh"
%include "chemical-equilibrium-system.hh"
%include "chemical-species-library.hh"
%include "species-energy-modes.hh"
%include "../kinetics/energy-exchange-mechanism.hh"
%include "../../nm/source/functor.hh"
%include "CI-functor.hh"
%include "transport-coefficients-model.hh"
%include "GuptaYos-mixing-rule.hh"
%include "collision-integral.hh"
%include "binary-interaction.hh"
%include "diatom-electronic-level.hh"

%pythoncode %{
%}

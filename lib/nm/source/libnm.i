/** \file libnm.i
 *  \ingroup nm
 *  \brief SWIG interface file for C/C++ numerical methods
 *
 *  \author Rowan J Gollan
 *  \date 26-Apr-2006
 *
 **/    

%define NM_DOCSTRING
"Python interface to C/C++ numerical methods."
%enddef
%module(docstring=NM_DOCSTRING) libnm

%include "std_string.i"
%include "std_vector.i"
/** %include "cpointer.i"
    %pointer_class(double, doublep)
**/

// The following magic allows us to pass a list of numbers
// to be collected as a C++ vector.
#ifndef STD_VECTOR_TEMPLATES_ALREADY_DEFINED
%template(vectori) std::vector<int>;
%template(vectord) std::vector<double>;
#define STD_VECTOR_TEMPLATES_ALREADY_DEFINED
#endif

%include "typemaps.i"
%apply double *INOUT { double *h };

//%template(vectord) std::vector<double>;

%{
#include <sstream>
#include <string>
#include <vector>
#include "../../util/source/useful.h"
#include "no_fuss_linear_algebra.hh"
#include "ode_system.hh"
#include "ode_solver.hh"
#include "zero_system.hh"
#include "zero_finders.hh"
#include "ode_step.hh"
#include "exponential_integrals.hh"
#include "linear_interpolation.hh"
%}

%ignore Valmatrix::operator=;
%ignore str;
%ignore set;
%include "no_fuss_linear_algebra.hh"
%ignore eval;
%include "zero_system.hh"
%include "zero_finders.hh"
%include "zero_system.hh"
%ignore str;
%include "ode_solver.hh"
%include "exponential_integrals.hh"
%include "linear_interpolation.hh"

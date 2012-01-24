/** \file libfluid.i
 *  \ingroup libgas 
 *  \brief SWIG interface header file for C++ ideal fluid model library.
 *  \author PAJ
 */
%define DOCSTRING
"Python interface to the ideal fluid classes, implemented in C++"
%enddef
%module(docstring=DOCSTRING) libfluid

%include "std_string.i"
%include "std_vector.i"

// The following magic allows us to pass a list of numbers
// to be collected as a C++ vector.
#ifndef STD_VECTOR_TEMPLATES_ALREADY_DEFINED
%template(vectori) std::vector<int>;
%template(vectord) std::vector<double>;
#define STD_VECTOR_TEMPLATES_ALREADY_DEFINED
#endif

%{
#include <string>
#include <vector>
#include <cstdio>
#include "fluid_thermo.hh"
%}

%include "fluid_thermo.hh"




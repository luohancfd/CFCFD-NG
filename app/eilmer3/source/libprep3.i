/// \file libprep3.i
/// \brief SWIG interface header for the combined libgas2 and libgeom2 module.
/// \author PJ
/// \version 23-May-2007  Adapted from Rowan's libtpib.i.
///
/// To be used in e3prep.py so that only one binary
/// library is used to contain all of the C++ functions.
/// That way we end up with just one copy of the numerical method functions
/// and the like.  Although the old approach of loading libgas2 and libgeom2
/// was fine on Linux and Cygwin, the Mac OSX system seemed to be sensitive
/// to the order in which the separate modules were loaded.

%define PREP_DOCSTRING
"Python interface to the combined libgas and libgeom2 modules."
%enddef
%module(docstring=PREP_DOCSTRING) libprep3

%include "std_string.i"
%include "std_vector.i"
%include "../../../lib/geometry2/source/libgeom2.i"
%include "../../../lib/gas/models/libgas.i"
%include "../../../lib/radiation/source/librad.i"
%include "radiation_transport.hh"

// The following magic allows us to pass a list of numbers
// to be collected as a C++ vector in the CFlowCondition constructor.
#ifndef STD_VECTOR_TEMPLATES_ALREADY_DEFINED
%template(vectori) std::vector<int>;
%template(vectord) std::vector<double>;
#define STD_VECTOR_TEMPLATES_ALREADY_DEFINED
#endif

%{
#include <string>
#include <vector>
#include "c-flow-condition.hh"
#include "kernel.hh"
#include "radiation_transport.hh"
%}

%rename(cflowcondition_print) operator<<( std::ostream &os, const CFlowCondition &cfc );
%include "c-flow-condition.hh"
%extend CFlowCondition {
    char *__str__() {
        static char tmp[512];
	strncpy(tmp, self->str().c_str(), 511);
	tmp[511] = '\0';
	return tmp;
    }
};

// To get access to the gas model pointer and its associated services.
%include "kernel.hh"

